from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE_COLUMNS = [
    "window_start",
    "animal_id",
    "tag_id",
    "sex",
    "age",
    "origin_group",
    "assigned_social_unit",
    "dynamic_social_unit",
    "social_context",
    "association_type",
    "association_event_id",
    "association_event_label",
    "membership_confidence",
    "is_observed",
    "is_carried_night",
    "is_local_2h_supported",
    "longitude",
    "latitude",
]

KEYS = ["night_date", "animal_id"]
NIGHT_HOURS = 15  # 16:00-23:00 plus 00:00-06:00 UTC, inclusive.


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collapse canonical hourly membership to one row per animal and biological night."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-output", type=Path)
    parser.add_argument("--sample-rows", type=int, default=5000)
    parser.add_argument("--metadata-output", type=Path)
    return parser.parse_args()


def dominant_value(
    frame: pd.DataFrame,
    value_column: str,
    output_column: str,
) -> pd.DataFrame:
    valid = frame.dropna(subset=[value_column]).copy()
    valid[value_column] = valid[value_column].astype(str)
    counts = (
        valid.groupby(KEYS + [value_column], observed=True)
        .size()
        .rename("dominant_hours")
        .reset_index()
    )
    counts = counts.sort_values(
        KEYS + ["dominant_hours", value_column],
        ascending=[True, True, False, True],
        kind="stable",
    )
    result = counts.drop_duplicates(KEYS).rename(
        columns={value_column: output_column, "dominant_hours": f"{output_column}_hours"}
    )
    return result[KEYS + [output_column, f"{output_column}_hours"]]


def classify_immediate_membership(row: pd.Series) -> str:
    assigned = str(row["nightly_assigned_social_unit"])
    origin = str(row["origin_group"])
    if assigned == "ISOLATED":
        return "isolated"
    if assigned.startswith("mixed:"):
        return "mixed_groups"
    if assigned == origin:
        return "with_origin_group"
    return "with_other_group"


def main() -> None:
    args = parse_args()
    hourly = pd.read_parquet(args.source, columns=SOURCE_COLUMNS)
    hourly["window_start"] = pd.to_datetime(hourly["window_start"])
    night_mask = hourly["window_start"].dt.hour.ge(16) | hourly["window_start"].dt.hour.le(6)
    night = hourly.loc[night_mask].copy()
    night["night_date"] = night["window_start"].dt.normalize()
    morning = night["window_start"].dt.hour.le(6)
    night.loc[morning, "night_date"] -= pd.Timedelta(days=1)
    night["is_exact_hour_supported"] = (
        night["is_observed"].fillna(False) & ~night["is_local_2h_supported"].fillna(False)
    )
    night = night.sort_values(KEYS + ["window_start"], kind="stable")

    if night.empty:
        raise ValueError("No nighttime hourly rows were found")
    duplicates = int(night.duplicated(["window_start", "animal_id"]).sum())
    if duplicates:
        raise ValueError(f"Found {duplicates} duplicate animal-hour rows")

    summary = (
        night.groupby(KEYS, observed=True)
        .agg(
            tag_id=("tag_id", "first"),
            sex=("sex", "first"),
            age=("age", "first"),
            origin_group=("origin_group", "first"),
            n_hourly_rows=("window_start", "size"),
            first_hour_utc=("window_start", "min"),
            last_hour_utc=("window_start", "max"),
            n_observed_hours=("is_observed", "sum"),
            n_exact_hour_supported_hours=("is_exact_hour_supported", "sum"),
            n_carried_hours=("is_carried_night", "sum"),
            n_local_2h_supported_hours=("is_local_2h_supported", "sum"),
            median_membership_confidence=("membership_confidence", "median"),
            minimum_membership_confidence=("membership_confidence", "min"),
            nightly_longitude=("longitude", "median"),
            nightly_latitude=("latitude", "median"),
            n_assigned_social_units=("assigned_social_unit", "nunique"),
            n_dynamic_social_units=("dynamic_social_unit", "nunique"),
            n_social_contexts=("social_context", "nunique"),
            n_association_types=("association_type", "nunique"),
        )
        .reset_index()
    )

    dominant_specs = [
        ("assigned_social_unit", "nightly_assigned_social_unit"),
        ("dynamic_social_unit", "nightly_dynamic_social_unit"),
        ("social_context", "nightly_social_context"),
        ("association_type", "nightly_association_type"),
        ("association_event_id", "nightly_association_event_id"),
        ("association_event_label", "nightly_association_event_label"),
    ]
    for source_column, output_column in dominant_specs:
        summary = summary.merge(
            dominant_value(night, source_column, output_column),
            on=KEYS,
            how="left",
            validate="one_to_one",
        )

    summary["night_start_utc"] = summary["night_date"] + pd.Timedelta(hours=16)
    summary["night_end_utc_exclusive"] = summary["night_date"] + pd.Timedelta(days=1, hours=7)
    summary["expected_night_hours"] = NIGHT_HOURS
    summary["night_coverage_fraction"] = np.minimum(summary["n_hourly_rows"] / NIGHT_HOURS, 1.0)
    summary["is_complete_night"] = summary["n_hourly_rows"].eq(NIGHT_HOURS)
    summary["dominant_assignment_fraction"] = (
        summary["nightly_assigned_social_unit_hours"] / summary["n_hourly_rows"]
    )
    summary["dominant_dynamic_fraction"] = (
        summary["nightly_dynamic_social_unit_hours"] / summary["n_hourly_rows"]
    )
    summary["dominant_association_fraction"] = (
        summary["nightly_association_type_hours"] / summary["n_hourly_rows"]
    )
    summary["nightly_membership_class"] = summary.apply(classify_immediate_membership, axis=1)
    summary["has_within_night_membership_conflict"] = (
        summary["n_assigned_social_units"].gt(1) | summary["n_dynamic_social_units"].gt(1)
    )
    summary["is_ambiguous_night"] = (
        summary["dominant_assignment_fraction"].le(0.5)
        | summary["night_coverage_fraction"].lt(0.5)
    )
    summary["evidence_type"] = np.select(
        [
            summary["n_exact_hour_supported_hours"].gt(0),
            summary["n_local_2h_supported_hours"].gt(0),
            summary["n_carried_hours"].gt(0),
        ],
        ["exact_hour_supported", "local_2h_supported", "carried_night_only"],
        default="no_supported_hour",
    )

    output_columns = [
        "night_date",
        "night_start_utc",
        "night_end_utc_exclusive",
        "animal_id",
        "tag_id",
        "sex",
        "age",
        "origin_group",
        "nightly_membership_class",
        "nightly_assigned_social_unit",
        "nightly_dynamic_social_unit",
        "nightly_social_context",
        "nightly_association_type",
        "nightly_association_event_id",
        "nightly_association_event_label",
        "n_hourly_rows",
        "expected_night_hours",
        "night_coverage_fraction",
        "is_complete_night",
        "first_hour_utc",
        "last_hour_utc",
        "n_observed_hours",
        "n_exact_hour_supported_hours",
        "n_carried_hours",
        "n_local_2h_supported_hours",
        "evidence_type",
        "dominant_assignment_fraction",
        "dominant_dynamic_fraction",
        "dominant_association_fraction",
        "median_membership_confidence",
        "minimum_membership_confidence",
        "has_within_night_membership_conflict",
        "is_ambiguous_night",
        "n_assigned_social_units",
        "n_dynamic_social_units",
        "n_social_contexts",
        "n_association_types",
        "nightly_longitude",
        "nightly_latitude",
    ]
    summary = summary[output_columns].sort_values(
        ["night_date", "origin_group", "animal_id"], kind="stable"
    )

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    summary.to_csv(partial, index=False, date_format="%Y-%m-%d %H:%M:%S")
    os.replace(partial, output)

    if args.sample_output:
        sample_output = args.sample_output.resolve()
        sample_output.parent.mkdir(parents=True, exist_ok=True)
        summary.head(args.sample_rows).to_csv(
            sample_output, index=False, date_format="%Y-%m-%d %H:%M:%S"
        )

    metadata = {
        "source": str(args.source.resolve()),
        "output": str(output),
        "sample_output": str(args.sample_output.resolve()) if args.sample_output else None,
        "row_definition": "one animal in one biological night",
        "night_definition_utc": "16:00 through 06:00 inclusive; night_date is the evening date",
        "classification_rule": (
            "Dominant hourly assigned_social_unit defines immediate nightly membership; "
            "dominant dynamic_social_unit is retained as the conservative longer-term membership."
        ),
        "rows": int(len(summary)),
        "animals": int(summary["animal_id"].nunique()),
        "origin_groups": int(summary["origin_group"].nunique(dropna=True)),
        "start_night": str(summary["night_date"].min()),
        "end_night": str(summary["night_date"].max()),
        "complete_nights": int(summary["is_complete_night"].sum()),
        "membership_class_counts": summary["nightly_membership_class"].value_counts().to_dict(),
        "association_type_counts": summary["nightly_association_type"].value_counts().to_dict(),
        "within_night_membership_conflicts": int(summary["has_within_night_membership_conflict"].sum()),
        "ambiguous_nights": int(summary["is_ambiguous_night"].sum()),
        "columns": output_columns,
        "csv_bytes": int(output.stat().st_size),
    }
    metadata_output = args.metadata_output or output.with_suffix(".metadata.json")
    metadata_output.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
