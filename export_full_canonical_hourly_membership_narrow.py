from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd


NARROW_COLUMNS = [
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
    "temp_group_size",
    "association_cluster_n_animals",
    "membership_confidence",
    "is_observed",
    "is_carried_night",
    "is_local_2h_supported",
    "longitude",
    "latitude",
]

SORT_COLUMNS = ["window_start", "origin_group", "animal_id"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export every row of the canonical hourly membership table using the narrow sample schema."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")

    membership = pd.read_parquet(args.source, columns=NARROW_COLUMNS)
    membership = membership.sort_values(SORT_COLUMNS, kind="stable").reset_index(drop=True)

    if membership.empty:
        raise ValueError("Canonical membership source contains no rows")
    if membership[["window_start", "animal_id"]].isna().any().any():
        raise ValueError("Missing window_start or animal_id values")
    duplicate_count = int(membership.duplicated(["window_start", "animal_id"]).sum())
    if duplicate_count:
        raise ValueError(f"Found {duplicate_count} duplicate animal-hour rows")

    membership.to_csv(partial, index=False)
    os.replace(partial, output)

    metadata = {
        "source": str(args.source.resolve()),
        "output": str(output),
        "row_definition": "one animal in one hourly window",
        "rows": int(len(membership)),
        "columns": NARROW_COLUMNS,
        "column_count": len(NARROW_COLUMNS),
        "animals": int(membership["animal_id"].nunique()),
        "origin_groups": int(membership["origin_group"].nunique(dropna=True)),
        "start": str(membership["window_start"].min()),
        "end": str(membership["window_start"].max()),
        "duplicate_animal_hours": duplicate_count,
        "sort_order": SORT_COLUMNS,
        "csv_bytes": int(output.stat().st_size),
    }
    metadata_output = args.metadata_output or output.with_suffix(".metadata.json")
    metadata_output.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
