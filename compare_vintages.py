"""What changed between the frozen 2026-07-22 chain and the 2026-09-05 rebuild?

A rebuild that changes both the data and the coverage rule will move almost every number
in the paper, and the only way to report that honestly is to measure it rather than
re-quote the new figures as if nothing happened. This compares the two vintages on the
quantities the figures actually use, and separates three sources of change:

  NEW DATA      79 more days and 42 more collars on the canonical share
  NEW COVERAGE  --min-fixes 1 rather than 3, recovering hours that had only one or two
                fixes -- real observations that a three-fix median threshold discarded
  INFERENCE     interpolated hours, which should be a rounding error; if they are not,
                that is the finding

The comparison is restricted to the OVERLAP WINDOW and, where noted, to the animals
present in both, because otherwise "more animals and more days" explains everything and
nothing.

Read-only. Outputs: outputs/general_structure_2026_09/phase6_vintage_comparison/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

NEW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
OLD = Path("outputs/membership_export_narrow_frozen20260722/"
           "canonical_hourly_membership_narrow.csv")
GEN_NEW = Path("outputs/general_structure_2026_09")
GEN_OLD = Path("outputs/general_structure_2026_09_frozen20260722")
OUT = GEN_NEW / "phase6_vintage_comparison"

COLS = ["window_start", "animal_id", "origin_group", "dynamic_social_unit",
        "social_context", "association_event_id", "is_observed", "is_carried_night"]


def load(path: Path, label: str) -> pd.DataFrame:
    use = COLS + (["position_support_type"] if label == "new" else [])
    head = pd.read_csv(path, nrows=1)
    use = [c for c in use if c in head.columns]
    d = pd.read_csv(path, usecols=use, parse_dates=["window_start"])
    d["animal_id"] = d["animal_id"].astype(str)
    return d


def census(d: pd.DataFrame) -> dict:
    return {
        "rows": int(len(d)),
        "animals": int(d["animal_id"].nunique()),
        "origin_groups": int(d["origin_group"].astype(str).nunique()),
        "units": int(d["dynamic_social_unit"].astype(str).nunique()),
        "start": str(d["window_start"].min()), "end": str(d["window_start"].max()),
        "observed": int(d["is_observed"].sum()),
        "carried_night": int(d["is_carried_night"].sum()),
        "context_shares": {k: round(float(v), 4) for k, v in
                           d["social_context"].value_counts(normalize=True).items()},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not NEW.exists() or not OLD.exists():
        raise SystemExit("need both vintages; missing %s"
                         % (NEW if not NEW.exists() else OLD))
    new, old = load(NEW, "new"), load(OLD, "old")
    report = {"new": census(new), "old": census(old)}

    if "position_support_type" in new.columns:
        report["new"]["support_types"] = {
            k: int(v) for k, v in new["position_support_type"].value_counts().items()}
        interp = new["position_support_type"].eq("interpolated")
        report["new"]["interpolated_rows"] = int(interp.sum())
        report["new"]["interpolated_share"] = round(float(interp.mean()), 5)
        report["new"]["interpolated_animals"] = int(
            new.loc[interp, "animal_id"].nunique())

    # like for like: the overlap window, and the animals both vintages know about
    lo = max(new["window_start"].min(), old["window_start"].min())
    hi = min(new["window_start"].max(), old["window_start"].max())
    shared = sorted(set(new["animal_id"]) & set(old["animal_id"]))
    n_ov = new[new["window_start"].between(lo, hi) & new["animal_id"].isin(shared)]
    o_ov = old[old["window_start"].between(lo, hi) & old["animal_id"].isin(shared)]
    report["overlap"] = {
        "window": [str(lo), str(hi)],
        "shared_animals": len(shared),
        "animals_only_new": int(len(set(new["animal_id"]) - set(old["animal_id"]))),
        "animals_only_old": int(len(set(old["animal_id"]) - set(new["animal_id"]))),
        "rows_new": int(len(n_ov)), "rows_old": int(len(o_ov)),
        "row_ratio": round(float(len(n_ov) / max(1, len(o_ov))), 3),
    }

    # the same animal-hours in both: does the state agree where both have a row?
    key = ["animal_id", "window_start"]
    j = (n_ov[key + ["social_context", "dynamic_social_unit"]]
         .merge(o_ov[key + ["social_context", "dynamic_social_unit"]],
                on=key, suffixes=("_new", "_old")))
    if len(j):
        report["shared_animal_hours"] = {
            "n": int(len(j)),
            "context_agreement": round(
                float(j["social_context_new"].eq(j["social_context_old"]).mean()), 4),
            "unit_agreement": round(
                float(j["dynamic_social_unit_new"]
                      .eq(j["dynamic_social_unit_old"]).mean()), 4),
        }
        ct = pd.crosstab(j["social_context_old"], j["social_context_new"])
        ct.to_csv(args.output_dir / "context_old_by_new.csv")
        report["shared_animal_hours"]["biggest_moves"] = {
            "%s -> %s" % (a, b): int(ct.loc[a, b])
            for a in ct.index for b in ct.columns if a != b
            and ct.loc[a, b] >= max(50, 0.01 * ct.loc[a].sum())}

    # headline products, where both vintages have them
    prods = {
        "opportunity_dyad_day": "phase1_opportunity/opportunity_dyad_day.csv",
        "excursions_dominant": "phase4b_individual_axis/excursions_dominant_gap0.csv",
        "weekly_network_metrics": "phase4d_axis_b_frozen/weekly_network_metrics_frozen.csv",
        "transfer_edges": "phase4g_excursion_destinations/transfer_edges_dominant.csv",
    }
    report["products"] = {}
    for name, rel in prods.items():
        entry = {}
        for label, base in (("new", GEN_NEW), ("old", GEN_OLD)):
            p = base / rel
            entry[label] = int(len(pd.read_csv(p))) if p.exists() else None
        if entry.get("new") and entry.get("old"):
            entry["ratio"] = round(entry["new"] / entry["old"], 3)
        report["products"][name] = entry

    with open(args.output_dir / "vintage_comparison.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("=" * 84)
    print("VINTAGE COMPARISON  frozen 2026-07-22  vs  rebuild 2026-09-05")
    print("=" * 84)
    for k in ("rows", "animals", "origin_groups", "units", "observed", "carried_night"):
        print("  %-16s old %14s   new %14s"
              % (k, format(report["old"][k], ","), format(report["new"][k], ",")))
    print("  %-16s old %14s   new %14s"
          % ("span end", report["old"]["end"][:10], report["new"]["end"][:10]))
    if "interpolated_rows" in report["new"]:
        print("  interpolated     %s rows (%.3f%% of new) across %d animals"
              % (format(report["new"]["interpolated_rows"], ","),
                 100 * report["new"]["interpolated_share"],
                 report["new"]["interpolated_animals"]))
    ov = report["overlap"]
    print("\n  overlap window %s .. %s; %d shared animals (+%d new, -%d gone)"
          % (ov["window"][0][:10], ov["window"][1][:10], ov["shared_animals"],
             ov["animals_only_new"], ov["animals_only_old"]))
    print("  rows on shared animals in window: old %s -> new %s  (x%.2f)"
          % (format(ov["rows_old"], ","), format(ov["rows_new"], ","),
             ov["row_ratio"]))
    if "shared_animal_hours" in report:
        s = report["shared_animal_hours"]
        print("\n  animal-hours present in BOTH: %s" % format(s["n"], ","))
        print("    social_context agreement %.1f%%   dynamic unit agreement %.1f%%"
              % (100 * s["context_agreement"], 100 * s["unit_agreement"]))
        if s.get("biggest_moves"):
            print("    largest reclassifications:")
            for k, v in sorted(s["biggest_moves"].items(), key=lambda x: -x[1])[:8]:
                print("      %-52s %s" % (k, format(v, ",")))
    print("\n  product row counts:")
    for name, e in report["products"].items():
        print("    %-24s old %-10s new %-10s %s"
              % (name, e.get("old"), e.get("new"),
                 ("x%.2f" % e["ratio"]) if e.get("ratio") else ""))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
