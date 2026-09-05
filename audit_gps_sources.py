"""Which GPS file is canonical? An inventory, because three answers are in circulation.

Three pipelines are running against what everyone calls "the cleaned GPS", and they do
not agree on how much data that is:

  * this project's frozen export      350 animals, 26 origin groups, 1,924,104 rows
  * the parallel hourly-grouping run  372 animals, 25,550,285 raw fixes, to 2026-07-02
  * the shared drive itself           \\10.126.19.90\EAS_shared\...\v1_cleaned\gps_v1.parquet

A per-animal discrepancy is already documented between copies: 25AA07_4S5T's last fix is
2026-01-29 in one and 2026-03-18 in another, which is enough to move it in and out of the
permanent-disperser set. Optimising anything downstream before this is settled is wasted
work, so this script states the facts rather than assuming them.

For every candidate it reports rows, animals, origin groups, date range, file mtime and a
content hash of the (animal, timestamp) key space, so two copies can be compared without
reading them both into memory at once. Read-only.

Outputs: outputs/general_structure_2026_09/phase5_source_audit/
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HOME = Path(r"C:\Users\rharel\Documents")
SHARED = Path(r"\\10.126.19.90\EAS_shared\baboon\working\data\processed\2025\gps"
              r"\v1_cleaned\gps_v1.parquet")
OUT = Path("outputs/general_structure_2026_09/phase5_source_audit")

RAW_CANDIDATES = [
    ("shared_drive", SHARED),
    ("new_project_network_v1", HOME / "New project/outputs/network_v1_cleaned_gps_v1.parquet"),
]
DERIVED_CANDIDATES = [
    ("frozen_narrow_export",
     Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")),
    ("hourly_grouping_observed",
     HOME / "New project/outputs/hourly_grouping_validation_shared_20260703"
            "/animal_hourly_grouping.parquet"),
    ("hourly_grouping_night_filled",
     HOME / "New project/outputs/hourly_grouping_validation_shared_20260703"
            "/animal_hourly_grouping_night_filled.parquet"),
]
PROBE_ANIMALS = ["25AA07_4S5T", "24AC17_7J8K", "24AA01_5O8B"]


def describe_raw(path: Path) -> dict:
    """Rows, animals, span and a key-space hash, without loading the whole file."""
    f = pq.ParquetFile(path)
    names = set(f.schema_arrow.names)
    acol = "animal_id" if "animal_id" in names else None
    tcol = "timestamp" if "timestamp" in names else None
    gcol = "group_id" if "group_id" in names else None
    cols = [c for c in (acol, tcol, gcol) if c]
    t = pq.read_table(path, columns=cols).to_pandas()
    t[acol] = t[acol].astype(str)
    out = {
        "rows": int(f.metadata.num_rows),
        "size_mb": round(path.stat().st_size / 1e6, 1),
        "modified": pd.Timestamp(path.stat().st_mtime, unit="s").isoformat(),
        "animals": int(t[acol].nunique()),
        "columns": len(f.schema_arrow.names),
    }
    if tcol:
        out["start"] = str(t[tcol].min())
        out["end"] = str(t[tcol].max())
    if gcol:
        out["origin_groups"] = int(t[gcol].astype(str).nunique())
    # a hash of the sorted (animal, hour) key space: two copies of the same underlying
    # record agree here even if columns or ordering differ
    key = (t[acol] + "|" + t[tcol].dt.floor("h").astype(str)).drop_duplicates()
    out["animal_hours"] = int(len(key))
    out["keyspace_sha1"] = hashlib.sha1(
        "\n".join(sorted(key)).encode("utf-8")).hexdigest()[:16]
    probes = {}
    for a in PROBE_ANIMALS:
        s = t.loc[t[acol].eq(a), tcol]
        probes[a] = {"n": int(len(s)),
                     "first": str(s.min()) if len(s) else None,
                     "last": str(s.max()) if len(s) else None}
    out["probe_animals"] = probes
    return out


def describe_derived(path: Path) -> dict:
    """The same census for an animal-hour table, whatever its column names."""
    if path.suffix == ".csv":
        d = pd.read_csv(path, usecols=["window_start", "animal_id", "origin_group"],
                        parse_dates=["window_start"])
        acol, tcol, gcol = "animal_id", "window_start", "origin_group"
    else:
        f = pq.ParquetFile(path)
        names = set(f.schema_arrow.names)
        tcol = "timestamp" if "timestamp" in names else "window_start"
        gcol = "origin_group" if "origin_group" in names else "group_id"
        d = pq.read_table(path, columns=["animal_id", tcol, gcol]).to_pandas()
        acol = "animal_id"
    d[acol] = d[acol].astype(str)
    d[tcol] = pd.to_datetime(d[tcol], utc=True, errors="coerce").dt.tz_localize(None)
    key = (d[acol] + "|" + d[tcol].dt.floor("h").astype(str)).drop_duplicates()
    return {
        "rows": int(len(d)),
        "size_mb": round(path.stat().st_size / 1e6, 1),
        "modified": pd.Timestamp(path.stat().st_mtime, unit="s").isoformat(),
        "animals": int(d[acol].nunique()),
        "origin_groups": int(d[gcol].astype(str).nunique()),
        "start": str(d[tcol].min()), "end": str(d[tcol].max()),
        "animal_hours": int(len(key)),
        "keyspace_sha1": hashlib.sha1(
            "\n".join(sorted(key)).encode("utf-8")).hexdigest()[:16],
        "probe_animals": {a: {
            "n": int(d[acol].eq(a).sum()),
            "first": str(d.loc[d[acol].eq(a), tcol].min())
            if d[acol].eq(a).any() else None,
            "last": str(d.loc[d[acol].eq(a), tcol].max())
            if d[acol].eq(a).any() else None} for a in PROBE_ANIMALS},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    report = {"raw": {}, "derived": {}, "unreachable": []}
    for name, path in RAW_CANDIDATES:
        if not path.exists():
            report["unreachable"].append({"name": name, "path": str(path)})
            print("  unreachable  %-24s %s" % (name, path))
            continue
        print("  reading raw  %-24s %s" % (name, path))
        report["raw"][name] = {"path": str(path), **describe_raw(path)}
    for name, path in DERIVED_CANDIDATES:
        if not path.exists():
            report["unreachable"].append({"name": name, "path": str(path)})
            print("  unreachable  %-24s %s" % (name, path))
            continue
        print("  reading      %-24s %s" % (name, path))
        report["derived"][name] = {"path": str(path), **describe_derived(path)}

    # do any two sources share a key space?
    keys = {k: v["keyspace_sha1"] for k, v in
            {**report["raw"], **report["derived"]}.items()}
    same = {}
    for a in keys:
        for b in keys:
            if a < b and keys[a] == keys[b]:
                same.setdefault(a, []).append(b)
    report["identical_keyspaces"] = same

    with open(args.output_dir / "gps_source_audit.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("=" * 96)
    print("GPS SOURCE AUDIT")
    print("=" * 96)
    hdr = ("source", "rows", "animals", "groups", "animal-hours", "start", "end")
    print("%-30s %12s %8s %7s %13s %11s %11s" % hdr)
    for kind in ("raw", "derived"):
        for name, v in report[kind].items():
            print("%-30s %12s %8s %7s %13s %11s %11s"
                  % (name, format(v["rows"], ","), v["animals"],
                     v.get("origin_groups", "-"), format(v["animal_hours"], ","),
                     str(v.get("start", ""))[:10], str(v.get("end", ""))[:10]))
    print("\nprobe animals (first / last fix):")
    for kind in ("raw", "derived"):
        for name, v in report[kind].items():
            for a, p in v["probe_animals"].items():
                if p["n"]:
                    print("  %-30s %-14s n=%-8s %s .. %s"
                          % (name, a, format(p["n"], ","), str(p["first"])[:16],
                             str(p["last"])[:16]))
    if report["unreachable"]:
        print("\nunreachable:")
        for u in report["unreachable"]:
            print("  %-28s %s" % (u["name"], u["path"]))
    print("\nidentical key spaces: %s" % (report["identical_keyspaces"] or "none"))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
