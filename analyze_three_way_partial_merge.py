"""Partial merges with a visible boundary: pure A, pure B and a mixed A+B, all at once.

A between-group encounter is usually described as if the two groups either meet or do
not. There is a third configuration, and it is the one that carries the most information
about how permeable a boundary is: in the SAME hour, a cluster holding only group A, a
cluster holding only group B, and a third cluster holding animals from both. The groups
are simultaneously apart and mixed, and the mixed cluster is a boundary you can count
animals across.

This measures it on the frozen export, mirroring a definition arrived at independently in
the parallel hourly-grouping pipeline (207 hour-cases, 137 events, 12 pairs there):

  * both units have more than MIN_SIDE_OBSERVED animals observed that hour
  * there is a cluster with 2+ animals of A and none of B
  * there is a cluster with 2+ animals of B and none of A
  * there is a cluster with at least one of each

Consecutive qualifying hours for the same pair collapse into one event. Read-only.

Outputs: outputs/general_structure_2026_09/phase5_three_way/
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
OUT = Path("outputs/general_structure_2026_09/phase5_three_way")

MIN_SIDE_OBSERVED = 5      # animals of each unit observed that hour, as in the other run
MIN_PURE = 2               # animals in a pure cluster for it to count as the group's body


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cols = ["window_start", "animal_id", "dynamic_social_unit", "association_event_id",
            "is_observed"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
    d = d[d["is_observed"] & d["association_event_id"].notna()]
    print("observed animal-hours with a cluster: %s" % format(len(d), ","))

    counts = (d.groupby(["window_start", "association_event_id", "dynamic_social_unit"])
              .size().rename("n").reset_index())
    present = (counts.groupby(["window_start", "dynamic_social_unit"])["n"].sum()
               .rename("unit_n").reset_index())
    big = present[present["unit_n"] > MIN_SIDE_OBSERVED]
    units_by_hour = big.groupby("window_start")["dynamic_social_unit"].apply(list)

    # cluster composition as {hour: {cluster: {unit: n}}}
    comp: dict = {}
    for r in counts.itertuples(index=False):
        comp.setdefault(r.window_start, {}).setdefault(
            r.association_event_id, {})[r.dynamic_social_unit] = int(r.n)

    rows = []
    for hour, units in units_by_hour.items():
        if len(units) < 2:
            continue
        clusters = comp.get(hour, {})
        for a, b in combinations(sorted(units), 2):
            pure_a = pure_b = mixed = None
            for cid, cc in clusters.items():
                na, nb = cc.get(a, 0), cc.get(b, 0)
                if na >= MIN_PURE and nb == 0 and pure_a is None:
                    pure_a = (cid, na)
                elif nb >= MIN_PURE and na == 0 and pure_b is None:
                    pure_b = (cid, nb)
                elif na >= 1 and nb >= 1 and mixed is None:
                    mixed = (cid, na, nb)
            if pure_a and pure_b and mixed:
                rows.append({"hour": hour, "unit_a": a, "unit_b": b,
                             "pure_a_n": pure_a[1], "pure_b_n": pure_b[1],
                             "mixed_a_n": mixed[1], "mixed_b_n": mixed[2],
                             "mixed_cluster": mixed[0]})
    hrs = pd.DataFrame(rows)
    if hrs.empty:
        print("no qualifying hours")
        return
    hrs["pair"] = hrs["unit_a"] + " + " + hrs["unit_b"]
    hrs = hrs.sort_values(["pair", "hour"])
    hrs.to_csv(args.output_dir / "three_way_hours.csv", index=False)

    # collapse consecutive hours of the same pair into events
    ev = []
    for pair, g in hrs.groupby("pair"):
        g = g.sort_values("hour")
        t = g["hour"].astype("int64") // 3_600_000_000_000
        start = prev = t.iloc[0]
        acc = [g.iloc[0]]
        for i in range(1, len(g)):
            if t.iloc[i] == prev + 1:
                acc.append(g.iloc[i])
            else:
                ev.append({"pair": pair, "start": pd.Timestamp(start * 3600, unit="s"),
                           "end": pd.Timestamp(prev * 3600, unit="s"),
                           "hours": prev - start + 1,
                           "max_mixed_a": max(r["mixed_a_n"] for r in acc),
                           "max_mixed_b": max(r["mixed_b_n"] for r in acc)})
                start, acc = t.iloc[i], [g.iloc[i]]
            prev = t.iloc[i]
        ev.append({"pair": pair, "start": pd.Timestamp(start * 3600, unit="s"),
                   "end": pd.Timestamp(prev * 3600, unit="s"),
                   "hours": prev - start + 1,
                   "max_mixed_a": max(r["mixed_a_n"] for r in acc),
                   "max_mixed_b": max(r["mixed_b_n"] for r in acc)})
    events = pd.DataFrame(ev).sort_values("hours", ascending=False)
    events.to_csv(args.output_dir / "three_way_events.csv", index=False)

    mixed_side = pd.concat([hrs["mixed_a_n"], hrs["mixed_b_n"]])
    report = {
        "min_side_observed": MIN_SIDE_OBSERVED, "min_pure": MIN_PURE,
        "hour_cases": int(len(hrs)),
        "events": int(len(events)),
        "pairs": int(hrs["pair"].nunique()),
        "median_event_hours": float(events["hours"].median()),
        "max_event_hours": int(events["hours"].max()),
        "mixed_side_median": float(mixed_side.median()),
        "mixed_side_is_one_animal": round(float((mixed_side == 1).mean()), 3),
        "top_pairs": hrs["pair"].value_counts().head(10).to_dict(),
    }
    with open(args.output_dir / "three_way_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("=" * 84)
    print("THREE-WAY PARTIAL MERGE  (pure A, pure B and a mixed A+B in one hour)")
    print("=" * 84)
    for k, v in report.items():
        if k != "top_pairs":
            print("  %-28s %s" % (k, v))
    print("\n  top pairs:")
    for k, v in report["top_pairs"].items():
        print("      %-34s %d hours" % (k, v))
    print("\n  longest events:")
    print(events.head(8).to_string(index=False))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
