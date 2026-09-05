"""Do the two pipelines say the same thing about the same animal-hour?

Two independent pipelines now assign a social state to every animal-hour, from different
GPS vintages, with different clustering rules that nonetheless agree at ARI 1.00 in the
median hour:

  THIS PROJECT   adaptive kNN linkage (120-900 m) -> `social_context`, one of
                 with_origin, mixed_with_origin_present, other, isolated,
                 mixed_without_origin_unit, mixed_tie_or_unclear
  NEW PROJECT    DBSCAN eps 300 m -> `hourly_state`, one of FULL_ORIGIN_GROUP,
                 WITHIN_GROUP_SPLIT, ISOLATED, DISPERSAL_WITH_OTHER_GROUP,
                 BETWEEN_GROUP_MERGE, MIXED_CLUSTER

The vocabularies are not the same shape, and that is the point. The other pipeline makes
a distinction this one carries only implicitly: a cluster of five Lilac and one Jade is
DISPERSAL_WITH_OTHER_GROUP, not a merge -- exactly the single-animal confound that had to
be stripped out of the encounter network here (639 of 2,867 detected-encounter dyad-days,
22.3%). If the two pipelines agree on that boundary, the correction is independently
confirmed; if they do not, one of them is wrong and this says which rows to look at.

The join is on (animal_id, hour) over the window both cover. Everything is reported as a
confusion matrix rather than a single agreement number, because a single number would
hide exactly the disagreements worth reading.

Read-only. Outputs: outputs/general_structure_2026_09/phase5_state_crosswalk/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HOME = Path(r"C:\Users\rharel\Documents")
OTHER = (HOME / "New project/outputs/hourly_grouping_validation_shared_20260703"
                "/animal_hourly_grouping.parquet")
NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
OUT = Path("outputs/general_structure_2026_09/phase5_state_crosswalk")

# how this project's context maps onto the other pipeline's vocabulary, where a mapping
# exists at all. `mixed_with_origin_present` deliberately spans two of their states --
# it says the animal is in a mixed cluster WITH its own group present, which is their
# BETWEEN_GROUP_MERGE when both sides are represented and their
# DISPERSAL_WITH_OTHER_GROUP seen from the majority side when only one animal crosses.
EXPECTED = {
    "with_origin": {"FULL_ORIGIN_GROUP", "WITHIN_GROUP_SPLIT"},
    "mixed_with_origin_present": {"BETWEEN_GROUP_MERGE", "DISPERSAL_WITH_OTHER_GROUP",
                                  "MIXED_CLUSTER"},
    "isolated": {"ISOLATED"},
    "other": {"DISPERSAL_WITH_OTHER_GROUP", "MIXED_CLUSTER"},
    "mixed_without_origin_unit": {"DISPERSAL_WITH_OTHER_GROUP", "BETWEEN_GROUP_MERGE",
                                  "MIXED_CLUSTER"},
    "mixed_tie_or_unclear": {"MIXED_CLUSTER", "BETWEEN_GROUP_MERGE"},
}


def load_other() -> pd.DataFrame:
    cols = ["timestamp", "animal_id", "origin_group", "hourly_state", "current_group",
            "cluster_size_total", "same_origin_cluster_size", "cluster_n_origin_groups",
            "cluster_origin_group_counts", "is_daytime_window"]
    d = pq.read_table(OTHER, columns=cols).to_pandas()
    d["animal_id"] = d["animal_id"].astype(str)
    d["origin_group"] = d["origin_group"].astype(str)
    d["hourly_state"] = d["hourly_state"].astype(str)
    d["hour"] = pd.to_datetime(d["timestamp"], utc=True).dt.tz_localize(None)
    return d.drop(columns=["timestamp"])


def load_ours() -> pd.DataFrame:
    cols = ["window_start", "animal_id", "origin_group", "dynamic_social_unit",
            "social_context", "association_event_id", "temp_group_size",
            "association_cluster_n_animals", "is_observed", "is_carried_night"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
    d["animal_id"] = d["animal_id"].astype(str)
    d = d.rename(columns={"window_start": "hour"})
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    o = load_other()
    m = load_ours()
    print("other pipeline %s rows, ours %s rows" % (format(len(o), ","),
                                                    format(len(m), ",")))
    lo = max(o["hour"].min(), m["hour"].min())
    hi = min(o["hour"].max(), m["hour"].max())
    o = o[o["hour"].between(lo, hi)]
    m = m[m["hour"].between(lo, hi)]
    print("shared window %s to %s" % (lo, hi))

    j = m.merge(o, on=["animal_id", "hour"], how="inner", suffixes=("_ours", "_other"))
    print("joined animal-hours: %s" % format(len(j), ","))

    # only rows this project actually observed are comparable: a carried-night row here
    # is an assumption, and the other pipeline's observed table has no such row
    obs = j[j["is_observed"]]
    ct = pd.crosstab(obs["social_context"], obs["hourly_state"])
    ct.to_csv(args.output_dir / "confusion_context_by_state.csv")
    row_pct = (ct.T / ct.sum(axis=1)).T.round(3)
    row_pct.to_csv(args.output_dir / "confusion_row_share.csv")

    agree = 0
    for ctx, states in EXPECTED.items():
        if ctx in ct.index:
            agree += int(ct.loc[ctx, [s for s in states if s in ct.columns]].sum())
    concordance = agree / max(1, int(ct.to_numpy().sum()))

    # the specific question: does their DISPERSAL_WITH_OTHER_GROUP land where this
    # project's single-animal encounter rows land?
    disp = obs[obs["hourly_state"].eq("DISPERSAL_WITH_OTHER_GROUP")]
    solo_like = obs[obs["same_origin_cluster_size"].eq(1)
                    & obs["cluster_n_origin_groups"].ge(2)]
    merge = obs[obs["hourly_state"].eq("BETWEEN_GROUP_MERGE")]

    report = {
        "shared_window": [str(lo), str(hi)],
        "rows_other": int(len(o)), "rows_ours": int(len(m)),
        "joined_animal_hours": int(len(j)),
        "joined_observed": int(len(obs)),
        "join_rate_of_ours": round(float(len(j) / max(1, len(m))), 4),
        "animals_in_join": int(j["animal_id"].nunique()),
        "concordance_under_expected_mapping": round(float(concordance), 4),
        "state_totals_other": obs["hourly_state"].value_counts().to_dict(),
        "context_totals_ours": obs["social_context"].value_counts().to_dict(),
        "dispersal_with_other_group": {
            "n": int(len(disp)),
            "our_context": disp["social_context"].value_counts().to_dict(),
            "median_same_origin_in_cluster": float(
                disp["same_origin_cluster_size"].median()) if len(disp) else None,
            "share_alone_on_their_side": round(float(
                (disp["same_origin_cluster_size"] == 1).mean()), 3) if len(disp) else None,
        },
        "single_animal_side_rows": {
            "n": int(len(solo_like)),
            "their_state": solo_like["hourly_state"].value_counts().to_dict(),
            "our_context": solo_like["social_context"].value_counts().to_dict(),
        },
        "between_group_merge": {
            "n": int(len(merge)),
            "our_context": merge["social_context"].value_counts().to_dict(),
            "median_same_origin_in_cluster": float(
                merge["same_origin_cluster_size"].median()) if len(merge) else None,
        },
    }
    with open(args.output_dir / "state_crosswalk_report.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("=" * 92)
    print("STATE CROSSWALK")
    print("=" * 92)
    print("joined %s animal-hours (%s observed here), %d animals; concordance %.1f%%"
          % (format(len(j), ","), format(len(obs), ","), j["animal_id"].nunique(),
             100 * concordance))
    print("\nrow shares: this project's social_context (rows) by their hourly_state")
    print(row_pct.to_string())
    print("\ncounts:")
    print(ct.to_string())
    print("\nTheir DISPERSAL_WITH_OTHER_GROUP (%s rows): %s"
          % (format(len(disp), ","), report["dispersal_with_other_group"]["our_context"]))
    print("  alone on their own side in %.0f%% of them"
          % (100 * (report["dispersal_with_other_group"]["share_alone_on_their_side"]
                    or 0)))
    print("\nRows where the animal is the ONLY one of its group in a mixed cluster "
          "(%s): %s" % (format(len(solo_like), ","),
                        report["single_animal_side_rows"]["their_state"]))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
