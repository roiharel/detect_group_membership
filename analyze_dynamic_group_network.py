"""Dynamic, weighted, directional between-group association network.

EDGE WEIGHT - corrected for shared collaring effort
---------------------------------------------------
Raw co-occurrence counts scale with how many animals each group has collared, so
well-collared dyads look more social. The correction is exact and cheap. For an
hour h, if group A has n_A(h) collars observed anywhere and n_A(h,c) of them in
spatial cluster c:

    opportunity(A,B,h) = n_A(h) * n_B(h)                  # co-observed animal pairs
    contact(A,B,h)     = sum_c n_A(h,c) * n_B(h,c)        # pairs sharing a cluster

    association(A,B)   = sum_h contact / sum_h opportunity

This is the fraction of realisable cross-group animal-pairings that actually
co-occurred, so a dyad with 2 collars each and one with 20 each are comparable.
Collar coverage across groups ranges 3%-42%, so this matters.

DIRECTION
---------
association() is symmetric. Direction comes from normalising within the sending
group: of all the intergroup association A has, what share goes to B?

    share(A->B) = association(A,B) / sum_C association(A,C)

share(A->B) != share(B->A) whenever the two groups have different numbers of
partners or different total outward association. A small group whose only
neighbour is a large one sends nearly all its association there, while the large
one spreads its own across many - a real asymmetry that a symmetric edge hides.

Nodes carry short codes; node area is demographic group size where available
(from EAS), never collar count.

Read-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
MEMBERSHIP = Path(
    r"C:\Users\rharel\Documents\New project\outputs"
    r"\canonical_robust_hourly_membership_shared_full_20260722\canonical_hourly_membership.parquet")
DEMOGRAPHICS = PROJECT / "outputs" / "authoritative_group_names_2026-09-03" / "group_demographics.csv"


def short_code(name: str, used: set[str]) -> str:
    """Three-letter node label, disambiguated on collision."""
    parts = [p for p in str(name).replace("_", " ").split() if p]
    if len(parts) > 1:                      # LapisSplinter-style handled below
        cand = "".join(p[0] for p in parts)[:3].title()
    else:
        s = parts[0]
        caps = [i for i, ch in enumerate(s) if ch.isupper()]
        cand = ("".join(s[i] for i in caps)[:3].title() if len(caps) > 1 else s[:3].title())
    base, i = cand, 2
    while cand in used:
        cand = f"{base[:2]}{i}"
        i += 1
    used.add(cand)
    return cand


def pair_products(counts: pd.DataFrame, keys: list[str], unit: str, n: str,
                  out: str) -> pd.DataFrame:
    left = counts.rename(columns={unit: "unit_a", n: "n_a"})
    right = counts.rename(columns={unit: "unit_b", n: "n_b"})
    m = left.merge(right, on=keys)
    m = m[m.unit_a < m.unit_b].copy()
    m[out] = m.n_a * m.n_b
    return m


def build_network(membership: Path, unit_col: str, freq: str,
                  min_opportunity: int) -> pd.DataFrame:
    m = pd.read_parquet(membership, columns=["window_start", "animal_id", "temp_group_id", unit_col])
    m = m[m[unit_col].notna()].copy()
    m["window_start"] = pd.to_datetime(m["window_start"])
    m[unit_col] = m[unit_col].astype(str)
    m["period"] = (m.window_start.dt.to_period("W").dt.start_time if freq == "weekly"
                   else m.window_start.values.astype("datetime64[M]"))
    print(f"rows {len(m):,}  hours {m.window_start.nunique():,}  units {m[unit_col].nunique()}")

    obs = (m.groupby(["period", "window_start", unit_col], observed=True)
             .animal_id.nunique().reset_index(name="n"))
    opp = pair_products(obs, ["period", "window_start"], unit_col, "n", "opportunity")
    opp = (opp.groupby(["period", "unit_a", "unit_b"], observed=True)
              .agg(opportunity=("opportunity", "sum"), hours=("window_start", "nunique"))
              .reset_index())

    cl = (m.groupby(["period", "window_start", "temp_group_id", unit_col], observed=True)
            .animal_id.nunique().reset_index(name="n"))
    con = pair_products(cl, ["period", "window_start", "temp_group_id"], unit_col, "n", "contact")
    con = (con.groupby(["period", "unit_a", "unit_b"], observed=True)
              .agg(contact=("contact", "sum"),
                   contact_hours=("window_start", "nunique")).reset_index())

    net = opp.merge(con, on=["period", "unit_a", "unit_b"], how="left").fillna(
        {"contact": 0, "contact_hours": 0})
    net = net[net.opportunity >= min_opportunity].copy()
    net["association"] = net.contact / net.opportunity
    return net


def add_direction(net: pd.DataFrame) -> pd.DataFrame:
    """Row-normalise within the sending group to get share(A->B)."""
    a = net.rename(columns={"unit_a": "source", "unit_b": "target"})
    b = net.rename(columns={"unit_b": "source", "unit_a": "target"})
    d = pd.concat([a, b], ignore_index=True)
    total = d.groupby(["period", "source"], observed=True).association.transform("sum")
    d["out_share"] = np.where(total > 0, d.association / total, np.nan)
    d["out_total"] = total
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership", type=Path, default=MEMBERSHIP)
    ap.add_argument("--unit-col", default="dynamic_social_unit")
    ap.add_argument("--freq", choices=["weekly", "monthly"], default="monthly")
    ap.add_argument("--min-opportunity", type=int, default=200,
                    help="minimum co-observed animal-pair-hours per dyad-period")
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "dynamic_group_network_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    net = build_network(a.membership, a.unit_col, a.freq, a.min_opportunity)
    directed = add_direction(net)
    net.to_csv(a.output_dir / f"group_network_{a.freq}_undirected.csv", index=False)
    directed.to_csv(a.output_dir / f"group_network_{a.freq}_directed.csv", index=False)

    units = sorted(set(net.unit_a) | set(net.unit_b))
    used: set[str] = set()
    codes = {u: short_code(u, used) for u in units}
    nodes = pd.DataFrame({"unit": units, "code": [codes[u] for u in units]})
    if DEMOGRAPHICS.exists():
        demo = pd.read_csv(DEMOGRAPHICS)
        demo["key"] = demo.group_id.str.replace(" ", "", regex=False).str.lower()
        nodes["key"] = nodes.unit.str.replace(" ", "", regex=False).str.lower()
        nodes = nodes.merge(demo[["key", "group_size", "no_collars"]], on="key", how="left")
    nodes.to_csv(a.output_dir / "node_codes.csv", index=False)
    print("\nnode codes:")
    print("  " + "   ".join(f"{codes[u]}={u}" for u in units))

    strongest = (net.groupby(["unit_a", "unit_b"], as_index=False)
                    .agg(opportunity=("opportunity", "sum"), contact=("contact", "sum"),
                         periods=("period", "nunique")))
    strongest["association"] = strongest.contact / strongest.opportunity
    strongest = strongest.sort_values("association", ascending=False)
    strongest.to_csv(a.output_dir / "dyad_association_overall.csv", index=False)
    print(f"\n=== strongest dyads (effort-corrected, pooled over {a.freq} periods) ===")
    print(f"{'dyad':<28}{'periods':>8}{'oppPairHrs':>12}{'assoc':>8}")
    for r in strongest.head(12).itertuples(index=False):
        print(f"{codes[r.unit_a]+' - '+codes[r.unit_b]:<28}{r.periods:>8}"
              f"{r.opportunity:>12,}{r.association:>8.3f}")

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "membership_source": str(a.membership), "unit_column": a.unit_col,
        "frequency": a.freq, "min_opportunity_pair_hours": a.min_opportunity,
        "edge_weight": "association = sum_h sum_c n_A(h,c)*n_B(h,c) / sum_h n_A(h)*n_B(h)",
        "effort_correction": "denominator is co-observed animal-pair-hours, so dyads with "
                             "different collar counts are comparable (coverage 3%-42%)",
        "direction": "out_share(A->B) = association(A,B) / sum_C association(A,C)",
        "node_size": "demographic group size from EAS, not collar count",
        "node_codes": codes,
        "periods": int(net.period.nunique()), "dyads": int(len(strongest)),
        "caveat": "Association is spatial cluster co-membership at the hourly scale "
                  "(adaptive 120-900 m edges). It is not fine-scale contact and not affiliation.",
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
