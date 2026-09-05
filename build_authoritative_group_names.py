"""Snapshot authoritative group names from EAS and build a pipeline override table.

WHY
---
`origin_group` in the canonical membership pipeline is derived from the `group_id`
column of `network_v1_cleaned_gps_v1.parquet`, which is a STATIC per-animal field
frozen when that file was written (2026-06-18). The authoritative animal->group
mapping lives on EAS and is actively maintained: between the local copy
(2026-04-24) and the EAS copy (2026-08-17), 16 animals changed group, including
two that move in or out of Copper/Lilac. Nothing in the pipeline reads it.

This script snapshots the EAS metadata, diffs it against the GPS file, and writes
an override table that `build_canonical_robust_hourly_membership.py` can consume
via --group-overrides.

It also copies the demographic group-size table (group_size, no_collars,
collar_coverage_percent), which is the independent source needed to stop using
observed collar counts as a proxy for group size.

EAS must be reachable (mapped drive or UNC). Read-only with respect to EAS.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

PROJECT = Path(__file__).resolve().parent
EAS_MOVEBANK = Path(r"Z:\baboon\working\data\processed\2025\acc\movebank_metadata.csv")
EAS_DEMOGRAPHICS = Path(r"Z:\baboon\working\data\processed\2025\metadata\GS_collars_demographics.xlsx")
EAS_REFERENCE = Path(r"Z:\baboon\working\data\processed\2025\metadata\Baboons MBRP Mpala Kenya-reference-data.csv")
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")


def gps_groups(path: Path) -> pd.DataFrame:
    f = pq.ParquetFile(path)
    seen = []
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=["animal_id", "group_id"]).to_pandas()
        seen.append(t.drop_duplicates())
    d = pd.concat(seen).drop_duplicates()
    multi = d.animal_id.value_counts()
    multi = multi[multi > 1]
    if len(multi):
        print(f"  NOTE: {len(multi)} animals carry >1 group_id in GPS; keeping first by sort")
    return (d.sort_values(["animal_id", "group_id"])
             .drop_duplicates("animal_id")
             .rename(columns={"group_id": "group_gps"}))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--movebank", type=Path, default=EAS_MOVEBANK)
    ap.add_argument("--demographics", type=Path, default=EAS_DEMOGRAPHICS)
    ap.add_argument("--gps", type=Path, default=GPS)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "authoritative_group_names_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    if not a.movebank.exists():
        raise SystemExit(f"EAS not reachable: {a.movebank}\n"
                         "Connect to the institute network/VPN, then re-run.")

    eas = pd.read_csv(a.movebank, dtype=str, low_memory=False)
    eas = (eas[["animal-id", "animal-group-id"]]
           .rename(columns={"animal-id": "animal_id", "animal-group-id": "group_eas"})
           .dropna(subset=["animal_id"]).drop_duplicates("animal_id"))
    eas["source_mtime"] = datetime.fromtimestamp(
        a.movebank.stat().st_mtime, tz=timezone.utc).isoformat(timespec="seconds")
    print(f"EAS metadata     : {len(eas)} animals  (mtime {eas.source_mtime.iloc[0]})")

    gps = gps_groups(a.gps)
    print(f"GPS group_id     : {len(gps)} animals")

    m = gps.merge(eas[["animal_id", "group_eas"]], on="animal_id", how="outer", indicator=True)
    both = m[m._merge == "both"]
    override = both[both.group_gps != both.group_eas].copy()
    override["reason"] = "EAS metadata is authoritative and newer than the GPS snapshot"

    print(f"\n=== animals whose authoritative group differs from GPS: {len(override)} ===")
    if len(override):
        print(override[["animal_id", "group_gps", "group_eas"]]
              .sort_values(["group_gps", "animal_id"]).to_string(index=False))

    out = override[["animal_id", "group_eas", "group_gps", "reason"]].rename(
        columns={"group_eas": "group_id"})
    out.to_csv(a.output_dir / "group_overrides.csv", index=False)

    eas[["animal_id", "group_eas"]].rename(columns={"group_eas": "group_id"}).to_csv(
        a.output_dir / "authoritative_animal_groups.csv", index=False)

    only_gps = m[m._merge == "left_only"].animal_id.tolist()
    only_eas = m[m._merge == "right_only"].animal_id.tolist()
    print(f"\nin GPS but not EAS metadata : {len(only_gps)}")
    print(f"in EAS metadata but not GPS : {len(only_eas)}")

    demo = None
    if a.demographics.exists():
        demo = pd.read_excel(a.demographics)
        demo.to_csv(a.output_dir / "group_demographics.csv", index=False)
        print(f"\n=== demographic group sizes ({len(demo)} groups) ===")
        print(demo.to_string(index=False))

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "prepared_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "eas_movebank": {"path": str(a.movebank), "animals": int(len(eas)),
                         "mtime_utc": eas.source_mtime.iloc[0]},
        "gps_snapshot": {"path": str(a.gps), "animals": int(len(gps)),
                         "mtime_utc": datetime.fromtimestamp(
                             a.gps.stat().st_mtime, tz=timezone.utc).isoformat(timespec="seconds")},
        "overrides": int(len(override)),
        "override_detail": out.to_dict("records"),
        "in_gps_not_eas": len(only_gps),
        "in_eas_not_gps": len(only_eas),
        "demographics_groups": int(len(demo)) if demo is not None else None,
        "not_established": [
            "The nine animals collared as Lilac on 2025-07-31 are labelled Lilac in EVERY "
            "available source (EAS 2026-08-17, EAS reference 2025-12-05, local 2026-04-24). "
            "The server does not correct them. Their failure to associate with the Lilac core "
            "is therefore not a stale-label artifact and needs a field explanation.",
            "group_id is static per animal in the GPS snapshot, so it cannot represent a group "
            "change over time; these overrides are also a single static value per animal.",
        ],
    }, indent=2), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
