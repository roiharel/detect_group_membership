"""Rebuild the whole chain on the 2026-09-05 canonical GPS, all 392 collars.

Twelve steps stand between a GPS file and Figures 1 and 2, and running them by hand is
how a vintage mismatch gets introduced. This drives them in order, logs each one, and
stops on the first failure.

WHAT CHANGES relative to the frozen 2026-07-22 chain

  data     the canonical share file as of 2026-09-05: 30,143,804 fixes, 392 animals,
           27 groups, through 2026-09-05, against 350 animals through 2026-07-22
  method   --min-fixes 1 instead of 3, which recovers 409,196 observed animal-hours
           (+36%) that a three-fix threshold was discarding, including essentially the
           entire record of the 22 collars that fix only every two hours
  method   interpolation of fix-free hours between bracketing observations, with the
           ceiling scaled to each collar's own fix rate (3 x its natural spacing, capped
           at 6 h) so a two-hourly collar's cadence hole is filled while a dense collar's
           six-hour failure is not. Filled rows are flagged
           `position_support_type='interpolated'` and are passengers: excluded from the
           kNN scale, unable to bridge two clusters, and attached to the cluster of their
           nearest observed neighbour

Everything else is pinned to the frozen run's settings -- local_2h support, 7-day
dynamic-unit rule, sparse-origin guard at 3 animals, adaptive k=2 factor 1.65 clipped
120-900 m -- so the comparison isolates data and coverage from the rest of the method.

THE OLD VINTAGE IS ARCHIVED, NOT OVERWRITTEN. `outputs/general_structure_2026_09` and
`outputs/membership_export_narrow` are renamed with a `_frozen20260722` suffix before
anything writes to those paths, so the frozen figures stay reproducible from the archive.

Run with `--dry-run` to print the plan without executing it.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT = Path(r"C:\Users\rharel\Documents\group_mebership")
UPSTREAM = Path(r"C:\Users\rharel\Documents\New project")
GPS = PROJECT / "data" / "gps_v1_canonical_20260905.parquet"

STAMP = "20260905"
CANON = PROJECT / "outputs" / ("canonical_%s_v2" % STAMP)
MEMBERSHIP = CANON / "canonical_hourly_membership.parquet"
WITH_EVENTS = CANON / "canonical_hourly_membership_with_association_events.parquet"
NARROW_DIR = PROJECT / "outputs" / "membership_export_narrow"
NARROW = NARROW_DIR / "canonical_hourly_membership_narrow.csv"
GEN = PROJECT / "outputs" / "general_structure_2026_09"
LOGS = PROJECT / "outputs" / ("rebuild_logs_%s" % STAMP)

FINE_5M = CANON / "canonical_5m_shared_history"
FINE_5M_EXP = CANON / "canonical_5m_shared_history_shuffle_expectation"
FINE_2M = CANON / "canonical_2m_shared_history"
FINE_2M_EXP = CANON / "canonical_2m_shared_history_shuffle_expectation"

# pinned to the frozen run so only data and coverage differ
MEMBERSHIP_ARGS = [
    "--window-minutes", "60",
    "--min-fixes", "1",
    # ceiling scaled to each collar's own fix rate: 6 h for the sparse collars whose
    # natural spacing is hours, 3 h for dense ones whose holes are failures
    "--interpolate-max-gap-hours", "6",
    "--interpolate-cadence-multiple", "3",
    "--temporal-support-mode", "local_2h",
    "--support-radius-minutes", "60", "--max-support-lag-minutes", "60",
    "--method", "adaptive", "--adaptive-k", "2", "--adaptive-factor", "1.65",
    "--adaptive-min-edge-m", "120", "--adaptive-max-edge-m", "900",
    "--night-anchor-time", "15:30", "--night-end-hour", "6",
    "--origin-presence-min", "2",
    "--dynamic-min-days", "7", "--dynamic-min-outside-hours", "24",
    "--dynamic-max-gap-hours", "36",
    "--dynamic-origin-sparse-radius-m", "500",
    "--dynamic-origin-sparse-min-hours", "12",
    "--dynamic-origin-sparse-min-fraction", "0.10",
    "--dynamic-origin-sparse-min-animals", "3",
]


def step(name, cwd, cmd, produces=None, env_extra=None):
    return {"name": name, "cwd": cwd, "cmd": cmd, "produces": produces,
            "env": env_extra or {}}


def plan() -> list[dict]:
    py = sys.executable
    return [
        step("membership", UPSTREAM,
             [py, "-u", "build_canonical_robust_hourly_membership.py",
              "--gps-file", str(GPS), "--output-dir", str(CANON)] + MEMBERSHIP_ARGS,
             produces=MEMBERSHIP),
        step("association_events", UPSTREAM,
             [py, "-u", "add_association_events_to_membership.py",
              "--input", str(MEMBERSHIP), "--output", str(WITH_EVENTS)],
             produces=WITH_EVENTS),
        step("narrow_export", PROJECT,
             [py, "-u", "export_full_canonical_hourly_membership_narrow.py",
              "--source", str(WITH_EVENTS), "--output", str(NARROW),
              "--metadata-output",
              str(NARROW_DIR / "canonical_hourly_membership_narrow.metadata.json")],
             produces=NARROW),
        step("fine_5m", UPSTREAM,
             [py, "-u", "plot_canonical_5m_shared_history.py",
              "--canonical", str(WITH_EVENTS), "--gps-2min", str(GPS),
              "--output-dir", str(FINE_5M), "--edge-radius-m", "5"],
             produces=FINE_5M / "canonical_5m_2min_merge_metric_rows.csv"),
        step("fine_5m_expectation", UPSTREAM,
             [py, "-u", "plot_canonical_5m_shuffle_expectation_fast.py",
              "--input-csv", str(FINE_5M / "canonical_5m_2min_merge_metric_rows.csv"),
              "--output-dir", str(FINE_5M_EXP)],
             produces=FINE_5M_EXP / "canonical_5m_shuffle_expectation_2min_rows.csv"),
        step("fine_2m", UPSTREAM,
             [py, "-u", "plot_canonical_5m_shared_history.py",
              "--canonical", str(WITH_EVENTS), "--gps-2min", str(GPS),
              "--output-dir", str(FINE_2M), "--edge-radius-m", "2"],
             produces=FINE_2M / "canonical_5m_2min_merge_metric_rows.csv"),
        step("fine_2m_expectation", UPSTREAM,
             [py, "-u", "plot_canonical_5m_shuffle_expectation_fast.py",
              "--input-csv", str(FINE_2M / "canonical_5m_2min_merge_metric_rows.csv"),
              "--output-dir", str(FINE_2M_EXP)],
             produces=FINE_2M_EXP / "canonical_5m_shuffle_expectation_2min_rows.csv"),
        step("phase1_opportunity", PROJECT,
             [py, "-u", "phase1_opportunity_table.py"],
             produces=GEN / "phase1_opportunity" / "opportunity_dyad_day.csv",
             env_extra={"REBUILD_CANON": str(CANON)}),
        step("phase2_two_stage", PROJECT,
             [py, "-u", "phase2_two_stage_events.py"],
             produces=(GEN / "phase2_two_stage_events"
                       / "stage1_events_with_stage2_mixing.csv"),
             env_extra={"REBUILD_CANON": str(CANON)}),
        step("excursions", PROJECT,
             [py, "-u", "analyze_individual_axis_identifiability.py"],
             produces=(GEN / "phase4b_individual_axis"
                       / "excursions_dominant_gap0.csv")),
        step("weekly_network_metrics", PROJECT,
             [py, "-u", "rederive_axis_b_on_frozen_cohort.py"],
             produces=(GEN / "phase4d_axis_b_frozen"
                       / "weekly_network_metrics_frozen.csv")),
        step("excursion_destinations", PROJECT,
             [py, "-u", "derive_excursion_destinations.py"],
             produces=(GEN / "phase4g_excursion_destinations"
                       / "transfer_edges_dominant.csv")),
        step("level2_inputs", PROJECT,
             [py, "-u", "derive_level2_inputs.py"],
             produces=(GEN / "phase4h_level2_inputs"
                       / "encounter_bouts_hourly.csv")),
        step("figure1", PROJECT,
             [py, "-u", "build_event_examples_figure.py"],
             produces=PROJECT / "docs/framing_2026-09-04/event_examples_figure.html"),
        step("figure2", PROJECT,
             [py, "-u", "build_figure3.py"],
             produces=PROJECT / "docs/framing_2026-09-04/figure3.html"),
    ]


def archive_old_vintage(dry: bool) -> list[str]:
    """Move the frozen vintage aside so the rebuild cannot silently overwrite it."""
    moved = []
    for d in (GEN, NARROW_DIR):
        if not d.exists():
            continue
        dest = d.with_name(d.name + "_frozen20260722")
        if dest.exists():
            moved.append("%s already archived" % d.name)
            continue
        moved.append("%s -> %s" % (d.name, dest.name))
        if not dry:
            shutil.move(str(d), str(dest))
    return moved


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--from-step", default=None,
                    help="skip everything before this step name")
    ap.add_argument("--only", nargs="*", default=None)
    # a rebuild should rebuild: skipping is opt-in, for resuming a failed run
    ap.add_argument("--skip-existing", action="store_true", default=False)
    args = ap.parse_args()

    steps = plan()
    names = [s["name"] for s in steps]
    if args.from_step:
        if args.from_step not in names:
            raise SystemExit("unknown step %r; choose from %s"
                             % (args.from_step, ", ".join(names)))
        steps = steps[names.index(args.from_step):]
    if args.only:
        steps = [s for s in steps if s["name"] in set(args.only)]

    if not GPS.exists():
        raise SystemExit("canonical GPS copy missing: %s" % GPS)
    LOGS.mkdir(parents=True, exist_ok=True)

    print("=" * 84)
    print("REBUILD ALL  (canonical GPS %s, all collars, min-fixes 1, interpolation 2 h)"
          % STAMP)
    print("=" * 84)
    for line in archive_old_vintage(args.dry_run):
        print("  archive: %s" % line)
    for s in steps:
        print("  step %-24s -> %s" % (s["name"], s["produces"]))
    if args.dry_run:
        print("\ndry run: nothing executed")
        return

    results = []
    for s in steps:
        prod = s["produces"]
        if args.skip_existing and prod is not None and Path(prod).exists():
            print("\n--- %s: already present, skipping (%s)" % (s["name"], prod))
            results.append({"step": s["name"], "status": "skipped"})
            continue
        log = LOGS / ("%s.log" % s["name"])
        env = dict(os.environ)
        env.update(s["env"])
        print("\n--- %s  (log: %s)" % (s["name"], log.name), flush=True)
        t0 = time.time()
        with open(log, "w", encoding="utf-8", errors="replace") as fh:
            proc = subprocess.run(s["cmd"], cwd=str(s["cwd"]), env=env,
                                  stdout=fh, stderr=subprocess.STDOUT)
        dt = time.time() - t0
        ok = proc.returncode == 0 and (prod is None or Path(prod).exists())
        results.append({"step": s["name"], "status": "ok" if ok else "failed",
                        "returncode": proc.returncode, "seconds": round(dt, 1)})
        print("    %s in %.0f s" % ("ok" if ok else "FAILED", dt), flush=True)
        with open(LOGS / "rebuild_status.json", "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=2)
        if not ok:
            print("\nstopped at %s. Tail of its log:" % s["name"])
            tail = log.read_text(encoding="utf-8", errors="replace").splitlines()[-25:]
            print("\n".join("    " + t for t in tail))
            raise SystemExit(1)

    print("\n" + "=" * 84)
    for r in results:
        print("  %-24s %-8s %s s" % (r["step"], r["status"], r.get("seconds", "-")))
    print("\nlogs in %s" % LOGS)


if __name__ == "__main__":
    main()
