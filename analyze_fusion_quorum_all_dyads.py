"""Does the fusion-hour quorum saturate for every group pair, or only Copper-Lilac?

The fusion filter (plot_canonical_5m_shared_history.py, build_hourly_pair_rows)
scores an hour as fusion for a dyad when >=2 collars of EACH unit sit in one
spatial cluster and >=4 collars are in it in total. That rule's sensitivity rises
with collar count, so a dyad can appear permanently fused simply because it is
well collared.

For Copper-Lilac the control arm reaches a fusion rate of exactly 1.000 in 45 of
45 post-deployment weeks - a detection ceiling, not a biological fact. Since the
same filter gates every dyad in the project (including the 12-dyad 2 m analysis),
this script measures the ceiling across all of them.

For each dyad it reports:
  opportunity hours - both units had >=2 collars observed somewhere that hour
  fusion hours      - the quorum was met inside one cluster
  fusion rate       - fusion / opportunity
  saturated months  - months at rate >= 0.999 (the diagnostic for a ceiling)

Read-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
DEFAULT_MEMBERSHIP = Path(
    r"C:\Users\rharel\Documents\New project\outputs"
    r"\canonical_robust_hourly_membership_shared_full_20260722\canonical_hourly_membership.parquet"
)
MIN_PER_GROUP = 2
MIN_TOTAL = 4


def unordered_pairs(counts: pd.DataFrame, keys: list[str], unit: str, n: str) -> pd.DataFrame:
    """Self-join a (keys, unit, n) table into unordered unit pairs."""
    left = counts.rename(columns={unit: "unit_a", n: "n_a"})
    right = counts.rename(columns={unit: "unit_b", n: "n_b"})
    m = left.merge(right, on=keys)
    return m[m.unit_a < m.unit_b]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership", type=Path, default=DEFAULT_MEMBERSHIP)
    ap.add_argument("--unit-col", default="dynamic_social_unit")
    ap.add_argument("--min-opportunity-hours", type=int, default=200,
                    help="only report dyads with at least this much opportunity")
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "fusion_quorum_all_dyads_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    m = pd.read_parquet(a.membership,
                        columns=["window_start", "animal_id", "temp_group_id", a.unit_col])
    m = m[m[a.unit_col].notna()].copy()
    m["window_start"] = pd.to_datetime(m["window_start"])
    m[a.unit_col] = m[a.unit_col].astype(str)
    print(f"rows {len(m):,}  hours {m.window_start.nunique():,}  units {m[a.unit_col].nunique()}")

    # opportunity: >=2 collars of a unit observed anywhere in the hour
    obs = (m.groupby(["window_start", a.unit_col], observed=True).animal_id.nunique()
             .reset_index(name="n"))
    obs = obs[obs.n >= MIN_PER_GROUP]
    opp = unordered_pairs(obs, ["window_start"], a.unit_col, "n")
    opp["dyad"] = opp.unit_a + " - " + opp.unit_b
    print(f"opportunity dyad-hours: {len(opp):,}")

    # detection: the same quorum inside a single cluster, >=4 collars total
    cl = (m.groupby(["window_start", "temp_group_id", a.unit_col], observed=True)
            .animal_id.nunique().reset_index(name="n"))
    cl = cl[cl.n >= MIN_PER_GROUP]
    det = unordered_pairs(cl, ["window_start", "temp_group_id"], a.unit_col, "n")
    det = det[(det.n_a + det.n_b) >= MIN_TOTAL]
    det["dyad"] = det.unit_a + " - " + det.unit_b
    det = det[["window_start", "dyad"]].drop_duplicates()
    print(f"detected  dyad-hours: {len(det):,}")

    opp = opp[["window_start", "dyad"]].drop_duplicates()
    opp["fusion"] = opp.set_index(["window_start", "dyad"]).index.isin(
        det.set_index(["window_start", "dyad"]).index)
    opp["month"] = opp.window_start.values.astype("datetime64[M]")

    monthly = (opp.groupby(["dyad", "month"], as_index=False)
                  .agg(opportunity_hours=("fusion", "size"), fusion_hours=("fusion", "sum")))
    monthly["fusion_rate"] = monthly.fusion_hours / monthly.opportunity_hours
    monthly.to_csv(a.output_dir / "monthly_fusion_rate_by_dyad.csv", index=False)

    dyads = (monthly.groupby("dyad", as_index=False)
                    .agg(months=("month", "nunique"),
                         opportunity_hours=("opportunity_hours", "sum"),
                         fusion_hours=("fusion_hours", "sum")))
    dyads["fusion_rate"] = dyads.fusion_hours / dyads.opportunity_hours
    sat = (monthly[monthly.opportunity_hours >= 20]
           .assign(saturated=lambda d: d.fusion_rate >= 0.999)
           .groupby("dyad", as_index=False)
           .agg(months_scored=("saturated", "size"), months_saturated=("saturated", "sum")))
    dyads = dyads.merge(sat, on="dyad", how="left")
    dyads["pct_months_saturated"] = 100 * dyads.months_saturated / dyads.months_scored
    dyads = dyads[dyads.opportunity_hours >= a.min_opportunity_hours].sort_values(
        "fusion_rate", ascending=False)
    dyads.to_csv(a.output_dir / "fusion_rate_by_dyad.csv", index=False)

    print(f"\n=== dyads with >= {a.min_opportunity_hours} opportunity hours "
          f"(quorum: >={MIN_PER_GROUP}/unit, >={MIN_TOTAL} total) ===")
    print(f"{'dyad':<34} {'oppHrs':>8} {'fusHrs':>8} {'rate':>6} {'satMonths':>10}")
    print("-" * 72)
    for r in dyads.itertuples(index=False):
        flag = "  <-- ceiling" if (r.pct_months_saturated or 0) >= 50 else ""
        sm = f"{int(r.months_saturated or 0)}/{int(r.months_scored or 0)}"
        print(f"{r.dyad:<34} {r.opportunity_hours:>8,} {r.fusion_hours:>8,} "
              f"{r.fusion_rate:>6.3f} {sm:>10}{flag}")

    n_ceiling = int((dyads.pct_months_saturated >= 50).sum())
    print(f"\n{n_ceiling} of {len(dyads)} dyads spend most scored months at the detection ceiling.")

    (a.output_dir / "fusion_quorum_all_dyads.json").write_text(json.dumps({
        "membership_source": str(a.membership),
        "unit_column": a.unit_col,
        "quorum": {"min_per_unit": MIN_PER_GROUP, "min_total": MIN_TOTAL},
        "normalisation": "fusion dyad-hours / dyad-hours where both units had >=2 collars observed",
        "dyads_reported": int(len(dyads)),
        "dyads_mostly_saturated": n_ceiling,
        "caveat": "Saturation means the quorum was met in every scored hour. It does not by "
                  "itself prove the groups were apart - it means the filter cannot distinguish "
                  "degrees of association for that dyad, so anything conditioned on its fusion "
                  "hours measures detection rather than association.",
    }, indent=2), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
