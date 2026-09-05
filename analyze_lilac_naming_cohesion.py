"""Is 'Lilac' after 2025-08-01 one social unit, or two entities under one name?

Six collars first labelled Lilac appear on 2025-08-01. In the monthly effort-
corrected analysis, the within-Lilac contact rate halves at exactly that date
while within-Copper is unchanged - but only when the new collars are included.
Among consistently tracked animals the within-Lilac rate does NOT fall. That is
the signature of a naming problem: animals labelled Lilac that do not actually
associate with the Lilac core.

This script tests it directly. It builds the within-group pairwise contact matrix
after the deployment and asks whether it separates into blocks aligned with the
collar cohort (original vs newly deployed). Copper is the control: three collars
were added there on the same day, so if block structure is an artifact of
deployment timing rather than group identity, Copper should show it too.

Read-only. Writes to a new dated output directory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

PROJECT = Path(__file__).resolve().parent
SRC = PROJECT / "outputs" / "copper_lilac_effort_corrected_integration"
PAIRS = SRC / "copper_lilac_pair_month_contact_rates.csv"
POSITIONS = SRC / "copper_lilac_fusion_2min_positions.parquet"
DEPLOYMENT = pd.Timestamp("2025-08-01")


def cohorts(positions_path: Path) -> pd.DataFrame:
    pos = pd.read_parquet(positions_path, columns=["bin_2min", "animal_id", "origin_group"])
    first = pos.groupby(["animal_id", "origin_group"], as_index=False).bin_2min.agg(["min", "max"])
    first.columns = ["animal_id", "origin_group", "first_fix", "last_fix"]
    first["cohort"] = np.where(first.first_fix < DEPLOYMENT, "original", "new_Aug2025")
    return first


def contact_matrix(pairs: pd.DataFrame, animals: list[str], radius: float) -> pd.DataFrame:
    d = pairs[(pairs.radius_m == radius) & pairs.animal_a.isin(animals) & pairs.animal_b.isin(animals)]
    g = (d.groupby(["animal_a", "animal_b"], as_index=False)
           .agg(opportunity=("opportunity_bins", "sum"), contact=("contact_bins", "sum")))
    g = g[g.opportunity >= 60]
    g["rate"] = g.contact / g.opportunity
    m = pd.DataFrame(np.nan, index=animals, columns=animals, dtype=float)
    for r in g.itertuples(index=False):
        m.loc[r.animal_a, r.animal_b] = r.rate
        m.loc[r.animal_b, r.animal_a] = r.rate
    np.fill_diagonal(m.values, np.nan)
    return m.dropna(how="all").dropna(how="all", axis=1)


def block_split(m: pd.DataFrame) -> pd.Series:
    """Two-block hierarchical split on 1 - (rate / max rate)."""
    filled = m.fillna(0.0)
    mx = filled.values.max()
    dist = 1.0 - (filled.values / mx if mx > 0 else filled.values)
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0
    z = linkage(squareform(dist, checks=False), method="average")
    return pd.Series(fcluster(z, 2, criterion="maxclust"), index=m.index, name="block")


def summarise(m: pd.DataFrame, labels: pd.Series, name: str) -> dict:
    a = labels.reindex(m.index)
    vals, out = m.values, {}
    same, diff = [], []
    for i, ai in enumerate(m.index):
        for j, aj in enumerate(m.index):
            if j <= i or np.isnan(vals[i, j]):
                continue
            (same if a[ai] == a[aj] else diff).append(vals[i, j])
    out[f"{name}_within_label_mean"] = float(np.mean(same)) if same else np.nan
    out[f"{name}_across_label_mean"] = float(np.mean(diff)) if diff else np.nan
    out[f"{name}_n_within"] = len(same)
    out[f"{name}_n_across"] = len(diff)
    if same and diff and np.mean(same) > 0:
        out[f"{name}_across_over_within"] = float(np.mean(diff) / np.mean(same))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--radius", type=float, default=5.0)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "lilac_naming_cohesion_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    info = cohorts(POSITIONS)
    pairs = pd.read_csv(PAIRS, parse_dates=["month"])
    post = pairs[pairs.month >= DEPLOYMENT]

    report, rows = {}, []
    for group in ["Lilac", "Copper"]:
        members = info[info.origin_group == group]
        ids = sorted(members.animal_id)
        m = contact_matrix(post, ids, a.radius)
        if m.empty or len(m) < 4:
            print(f"{group}: insufficient post-deployment support")
            continue
        coh = members.set_index("animal_id").cohort.reindex(m.index)
        blk = block_split(m)

        report.update(summarise(m, coh, f"{group}_by_cohort"))
        report.update(summarise(m, blk, f"{group}_by_datablock"))
        report[f"{group}_n_animals"] = int(len(m))
        # agreement between the data-driven split and the deployment cohort
        ct = pd.crosstab(blk, coh)
        report[f"{group}_block_vs_cohort_crosstab"] = ct.to_dict()

        print(f"\n=== {group} ({len(m)} animals with post-deployment support) @ {a.radius:.0f} m ===")
        print("data-driven 2-block split vs collar cohort:")
        print(ct.to_string())
        wc = report.get(f"{group}_by_cohort_within_label_mean", np.nan)
        ac = report.get(f"{group}_by_cohort_across_label_mean", np.nan)
        print(f"  contact within same cohort : {wc:.4f}  (n={report[f'{group}_by_cohort_n_within']} pairs)")
        print(f"  contact across cohorts     : {ac:.4f}  (n={report[f'{group}_by_cohort_n_across']} pairs)")
        if wc and np.isfinite(wc) and wc > 0:
            print(f"  ratio across/within        : {ac / wc:.2f}")

        mm = m.copy()
        mm.insert(0, "cohort", coh)
        mm.insert(1, "block", blk)
        mm.to_csv(a.output_dir / f"{group.lower()}_post_deployment_contact_matrix.csv")
        for animal in m.index:
            rows.append({"group": group, "animal_id": animal, "cohort": coh[animal],
                         "data_block": int(blk[animal]),
                         "mean_contact_to_original": float(np.nanmean(
                             m.loc[animal, coh[coh == "original"].index.intersection(m.index)])),
                         "mean_contact_to_new": float(np.nanmean(
                             m.loc[animal, coh[coh == "new_Aug2025"].index.intersection(m.index)]))})

    pd.DataFrame(rows).to_csv(a.output_dir / "animal_cohort_affinity.csv", index=False)
    (a.output_dir / "lilac_naming_cohesion.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
