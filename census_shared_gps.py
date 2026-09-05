"""Census the canonical GPS file on the EAS share, and say how it differs from our copy.

The frozen export's own run metadata names its source:

    \\\\10.126.19.90\\EAS_shared\\baboon\\working\\data\\processed\\2025\\gps\\v1_cleaned\\gps_v1.parquet

which is the same path the parallel hourly-grouping pipeline uses. There was never a
second canonical source -- the file is updated in place, so every "different source" was a
different vintage of one file. This measures the vintage that is on the share now against
the local copies, so the difference is a number rather than an assumption.

Must be run from a process that can see the UNC share. Git Bash cannot resolve it here;
PowerShell can, so invoke this with `powershell python census_shared_gps.py`.

Read-only. Outputs: outputs/general_structure_2026_09/phase5_source_audit/
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

SHARE = ("//10.126.19.90/EAS_shared/baboon/working/data/processed/2025/gps"
         "/v1_cleaned/gps_v1.parquet")
LOCAL = Path(r"C:\Users\rharel\Documents\New project\outputs"
             r"\network_v1_cleaned_gps_v1.parquet")
OUT = Path("outputs/general_structure_2026_09/phase5_source_audit")
PROBES = ["25AA07_4S5T", "24AC17_7J8K", "24AA01_5O8B"]


def census(path, label):
    """Rows, animals, groups, span and probe endpoints for one GPS parquet."""
    f = pq.ParquetFile(path)
    names = f.schema_arrow.names
    cols = [c for c in ("animal_id", "timestamp", "group_id") if c in names]
    t = pq.read_table(path, columns=cols).to_pandas()
    t["animal_id"] = t["animal_id"].astype(str)
    out = {
        "label": label, "path": str(path),
        "rows": int(f.metadata.num_rows),
        "size_gb": round(os.path.getsize(path) / 1e9, 3),
        "modified": str(pd.Timestamp(os.path.getmtime(path), unit="s")),
        "columns": len(names),
        "animals": int(t["animal_id"].nunique()),
        "start": str(t["timestamp"].min()), "end": str(t["timestamp"].max()),
    }
    if "group_id" in t:
        out["groups"] = int(t["group_id"].astype(str).nunique())
    out["probes"] = {}
    for a in PROBES:
        s = t.loc[t["animal_id"].eq(a), "timestamp"]
        out["probes"][a] = {"n": int(len(s)),
                            "first": str(s.min()) if len(s) else None,
                            "last": str(s.max()) if len(s) else None}
    per_animal = t.groupby("animal_id")["timestamp"].max()
    out["animals_ending_before_2026_06"] = int((per_animal < "2026-06-01").sum())
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    ap.add_argument("--share", default=SHARE)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    report = {}
    for label, path in (("share_current", args.share), ("local_network_v1", LOCAL)):
        if not os.path.exists(path):
            print("  unreachable: %s" % path)
            report[label] = {"unreachable": str(path)}
            continue
        print("  reading %s ..." % label, flush=True)
        report[label] = census(path, label)

    a, b = report.get("share_current"), report.get("local_network_v1")
    if a and "rows" in a and b and "rows" in b:
        report["delta_share_minus_local"] = {
            "rows": a["rows"] - b["rows"],
            "animals": a["animals"] - b["animals"],
            "end_shift_days": round(
                (pd.Timestamp(a["end"]) - pd.Timestamp(b["end"])).total_seconds()
                / 86400, 1),
        }
    with open(args.output_dir / "shared_gps_census.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("=" * 88)
    print("CANONICAL GPS CENSUS")
    print("=" * 88)
    for k, v in report.items():
        if not isinstance(v, dict) or "start" not in v:
            continue
        print("%-18s rows %14s  animals %4s  groups %3s  %s .. %s  (mod %s)"
              % (k, format(v["rows"], ","), v["animals"], v.get("groups", "-"),
                 v["start"][:10], v["end"][:10], v["modified"][:16]))
    if "delta_share_minus_local" in report:
        print("\nshare minus local: %s" % report["delta_share_minus_local"])
    print("\nprobe animals:")
    for k, v in report.items():
        if not isinstance(v, dict) or "probes" not in v:
            continue
        for aid, p in v["probes"].items():
            print("  %-18s %-14s n=%-9s %s .. %s"
                  % (k, aid, format(p["n"], ",") if p["n"] else 0,
                     str(p["first"])[:16], str(p["last"])[:16]))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
