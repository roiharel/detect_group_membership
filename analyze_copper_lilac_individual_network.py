"""Dynamic individual-level association network for Copper and Lilac.

Nodes are individual animals, not groups. Copper sits on the left arc and Lilac
on the right, so cross-origin association reads as chords crossing the middle and
within-origin association stays on its own side. Integration over time is visible
as the middle filling in.

EDGE WEIGHT - corrected for shared observation effort
-----------------------------------------------------
    association(a,b) = contact bins within `radius` / bins where BOTH were observed

The denominator is co-observed 2-minute bins for that specific pair, so animals
with different tracking duration or collar reliability are comparable. This is the
same effort correction as the group-level network, applied per animal pair.

DIRECTION
---------
    out_share(a->b) = association(a,b) / sum_c association(a,c)

i.e. what fraction of everything animal a does socially it does with b. Asymmetric
because a and b have different totals: a peripheral animal with one partner sends
all its share there, while a well-connected animal spreads its own thin.

COHORT
------
Nine collars deployed 2025-07-31 are labelled Lilac but form a distinct block
(see analyze_lilac_fission_hypothesis.py). They are drawn with a heavy outline so
the two Lilac cohorts can be told apart.

Read-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
SRC = PROJECT / "outputs" / "copper_lilac_effort_corrected_integration"
PAIRS = SRC / "copper_lilac_pair_month_contact_rates.csv"
POSITIONS = SRC / "copper_lilac_fusion_2min_positions.parquet"
DEPLOYMENT = pd.Timestamp("2025-08-01")
NEW_LILAC = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
             "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"}


def short_name(animal_id: str, used: set[str]) -> str:
    """24AC11_5W6X -> AC11 ; disambiguated on collision."""
    stem = str(animal_id).split("_")[0]
    cand = stem[2:] if len(stem) > 2 else stem
    base, i = cand, 2
    while cand in used:
        cand = f"{base}.{i}"
        i += 1
    used.add(cand)
    return cand


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, default=5.0)
    ap.add_argument("--freq", choices=["monthly", "weekly"], default="monthly")
    ap.add_argument("--min-bins", type=int, default=60,
                    help="minimum co-observed 2-min bins per pair-period")
    ap.add_argument("--cohort-month", default=None,
                    help="restrict to animals tracked in this month, e.g. 2025-05. Fixes the "
                         "node set for the whole series, so a change in the network cannot be "
                         "caused by animals entering or leaving it. Chosen before the "
                         "2025-07-31 deployment, this excludes every new collar.")
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "copper_lilac_individual_network_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(PAIRS, parse_dates=["month"])
    pairs = pairs[pairs.radius_m == a.radius].copy()
    if pairs.empty:
        raise SystemExit(f"no rows at radius {a.radius}; available: "
                         f"{sorted(pd.read_csv(PAIRS).radius_m.unique())}")

    pos = pd.read_parquet(POSITIONS, columns=["bin_2min", "animal_id", "origin_group"])
    info = (pos.groupby(["animal_id", "origin_group"], as_index=False)
               .bin_2min.agg(["min", "max"]).reset_index(drop=True))
    info.columns = ["animal_id", "origin_group", "first_fix", "last_fix"]
    info["cohort"] = np.where(info.animal_id.isin(NEW_LILAC), "new_Aug2025",
                              np.where(info.first_fix < DEPLOYMENT, "original", "later"))

    if a.freq == "weekly":
        raise SystemExit("The saved pair table is monthly. Use --freq monthly, or extend "
                         "this script to rebuild pair rates from the 2-minute parquet.")
    pairs["period"] = pairs["month"]

    cohort_ids = None
    if a.cohort_month:
        anchor = pd.Timestamp(a.cohort_month)
        tracked = pos[(pos.bin_2min >= anchor) &
                      (pos.bin_2min < anchor + pd.offsets.MonthBegin(1))]
        cohort_ids = set(tracked.animal_id.unique())
        if not cohort_ids:
            raise SystemExit(f"no animals tracked in {a.cohort_month}")
        pairs = pairs[pairs.animal_a.isin(cohort_ids) & pairs.animal_b.isin(cohort_ids)]
        info = info[info.animal_id.isin(cohort_ids)]
        by_origin = info.groupby("origin_group").size().to_dict()
        print(f"cohort filter: animals tracked in {a.cohort_month} -> "
              f"{len(cohort_ids)} animals {by_origin}")

    d = pairs[pairs.opportunity_bins >= a.min_bins].copy()
    print(f"radius {a.radius:g} m | pair-periods kept {len(d):,} of {len(pairs):,} "
          f"| animals {len(set(d.animal_a) | set(d.animal_b))} | periods {d.period.nunique()}")

    d.rename(columns={"opportunity_bins": "opportunity", "contact_bins": "contact"}, inplace=True)
    und = (d.groupby(["period", "animal_a", "animal_b", "origin_a", "origin_b", "pair_type"],
                     as_index=False).agg(opportunity=("opportunity", "sum"),
                                         contact=("contact", "sum")))
    und["association"] = und.contact / und.opportunity

    left = und.rename(columns={"animal_a": "source", "animal_b": "target"})
    right = und.rename(columns={"animal_b": "source", "animal_a": "target"})
    dire = pd.concat([left, right], ignore_index=True)
    tot = dire.groupby(["period", "source"], observed=True).association.transform("sum")
    dire["out_share"] = np.where(tot > 0, dire.association / tot, np.nan)

    und.to_csv(a.output_dir / f"individual_network_{a.freq}_undirected.csv", index=False)
    dire.to_csv(a.output_dir / f"individual_network_{a.freq}_directed.csv", index=False)

    animals = sorted(set(und.animal_a) | set(und.animal_b))
    info = info[info.animal_id.isin(animals)].sort_values(["origin_group", "cohort", "animal_id"])
    used: set[str] = set()
    info["code"] = [short_name(x, used) for x in info.animal_id]
    info.to_csv(a.output_dir / "node_codes.csv", index=False)

    for g, sub in info.groupby("origin_group"):
        print(f"\n{g} ({len(sub)}):")
        for c, s2 in sub.groupby("cohort"):
            print(f"  {c:12s} " + " ".join(s2.code))

    overall = (und.groupby(["animal_a", "animal_b", "pair_type"], as_index=False)
                  .agg(opportunity=("opportunity", "sum"), contact=("contact", "sum"),
                       periods=("period", "nunique")))
    overall["association"] = overall.contact / overall.opportunity
    overall.to_csv(a.output_dir / "pair_association_overall.csv", index=False)

    byclass = (und.groupby(["period", "pair_type"], as_index=False)
                  .agg(opportunity=("opportunity", "sum"), contact=("contact", "sum")))
    byclass["association"] = byclass.contact / byclass.opportunity
    byclass.to_csv(a.output_dir / "association_by_pair_type_period.csv", index=False)
    print(f"\n=== pooled association by pair type ===")
    print(byclass.groupby("pair_type").apply(
        lambda x: pd.Series({"periods": x.period.nunique(),
                             "association": x.contact.sum() / x.opportunity.sum()}),
        include_groups=False).round(4).to_string())

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "source_pairs": str(PAIRS), "radius_m": a.radius, "frequency": a.freq,
        "min_co_observed_bins_per_pair_period": a.min_bins,
        "edge_weight": "contact bins within radius / co-observed 2-minute bins for that pair",
        "direction": "out_share(a->b) = association(a,b) / sum_c association(a,c)",
        "n_animals": len(animals), "n_periods": int(und.period.nunique()),
        "cohort_month": a.cohort_month,
        "cohort_animals": sorted(cohort_ids) if cohort_ids else None,
        "new_lilac_cohort": sorted(NEW_LILAC),
        "caveat": "Conditioned on canonical Copper-Lilac fusion hours, which saturate after "
                  "2025-08-01 (see analyze_fusion_hour_collar_dependence.py). Association is "
                  "GPS proximity, not affiliation.",
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
