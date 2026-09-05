from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


BASE = Path(r"C:\Users\rharel\Documents")
DEFAULT_STATUS = (
    BASE
    / "New project"
    / "outputs"
    / "rerun_20260703_full_hierarchical_sampling_filtered_eps300"
    / "proximity_status_full.parquet"
)
DEFAULT_CANONICAL = (
    BASE
    / "New project"
    / "outputs"
    / "canonical_robust_hourly_membership_local_2h_support"
    / "canonical_hourly_membership.parquet"
)
DEFAULT_OUT = (
    BASE
    / "group_mebership"
    / "outputs"
    / "dynamic_social_unit_merge_gamm"
    / "proximity_status_dynamic_social_unit.parquet"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a proximity status table where sustained non-origin association "
            "rows use dynamic_social_unit as the analysis group."
        )
    )
    parser.add_argument("--status-file", type=Path, default=DEFAULT_STATUS)
    parser.add_argument("--canonical", type=Path, default=DEFAULT_CANONICAL)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    status = pd.read_parquet(args.status_file)
    status["timestamp"] = pd.to_datetime(status["timestamp"], utc=True).dt.tz_localize(None)
    status["animal_id"] = status["animal_id"].astype(str)
    status["basal_group"] = status["basal_group"].astype(str)

    canonical = pd.read_parquet(
        args.canonical,
        columns=[
            "animal_id",
            "window_start",
            "origin_group",
            "dynamic_social_unit",
            "dynamic_assignment",
            "dynamic_target_group",
        ],
    )
    canonical = canonical.rename(columns={"window_start": "timestamp"})
    canonical["timestamp"] = pd.to_datetime(canonical["timestamp"], utc=True).dt.tz_localize(None)
    canonical["animal_id"] = canonical["animal_id"].astype(str)
    canonical = canonical.drop_duplicates(["animal_id", "timestamp"])

    merged = status.merge(
        canonical,
        on=["animal_id", "timestamp"],
        how="left",
        suffixes=("", "_canonical"),
    )
    reassigned = merged["dynamic_assignment"].eq("sustained_non_origin_association")
    has_unit = merged["dynamic_social_unit"].notna()
    merged["origin_basal_group"] = merged["basal_group"]
    merged["basal_group"] = merged["basal_group"].where(
        ~(reassigned & has_unit),
        merged["dynamic_social_unit"].astype(str),
    )
    merged["analysis_group_source"] = "origin_or_default"
    merged.loc[reassigned & has_unit, "analysis_group_source"] = "sustained_non_origin_dynamic_social_unit"

    status_columns = list(status.columns)
    extra_columns = [
        "origin_basal_group",
        "origin_group",
        "dynamic_social_unit",
        "dynamic_assignment",
        "dynamic_target_group",
        "analysis_group_source",
    ]
    out = merged[status_columns + [col for col in extra_columns if col in merged.columns]].copy()
    out.to_parquet(args.out, index=False)

    summary = {
        "source_status": str(args.status_file),
        "source_canonical": str(args.canonical),
        "output": str(args.out),
        "rows": int(len(out)),
        "animals": int(out["animal_id"].nunique()),
        "original_groups": int(out["origin_basal_group"].nunique()),
        "dynamic_analysis_groups": int(out["basal_group"].nunique()),
        "reassigned_rows": int((reassigned & has_unit).sum()),
        "reassigned_animals": int(out.loc[reassigned & has_unit, "animal_id"].nunique()),
        "dynamic_assignment_counts": {
            str(k): int(v)
            for k, v in out["dynamic_assignment"].fillna("unmatched_status_row").value_counts().items()
        },
    }
    metadata = args.out.with_name(args.out.stem + "_metadata.json")
    metadata.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
