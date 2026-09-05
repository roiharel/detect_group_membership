"""Phase 0: cohort ledger, event-ID crosswalk and definition audit.

Inventories every membership, structural-event, fine-scale and model-row product
in circulation, so that each saved result can be mapped to exactly one cohort.
Reads only; never refits or modifies a source.

Run from the project root:
    python phase0_cohort_ledger.py
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

PROJECT = Path(r"C:\Users\rharel\Documents\group_mebership")
UPSTREAM = Path(r"C:\Users\rharel\Documents\New project")
CANON = UPSTREAM / "outputs" / "canonical_robust_hourly_membership_shared_full_20260722"
INVENTORY_SRC = UPSTREAM / "outputs" / "canonical_robust_hourly_membership_local_2h_support"
OUT_DIR = PROJECT / "outputs" / "general_structure_2026_09" / "phase0_cohort_ledger"

HASH_LIMIT_BYTES = 600_000_000


def sha256_of(path: Path) -> str | None:
    """Hash a file in chunks; skip very large files to keep the audit fast."""
    if path.stat().st_size > HASH_LIMIT_BYTES:
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_stamp(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "exists": True,
        "bytes": stat.st_size,
        "modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(timespec="seconds"),
        "sha256": sha256_of(path),
    }


def describe_parquet(path: Path, animal_col: str, time_col: str, group_cols: list[str]) -> dict[str, object]:
    handle = pq.ParquetFile(path)
    names = set(handle.schema_arrow.names)
    wanted = [c for c in [animal_col, time_col, *group_cols] if c in names]
    frame = pd.read_parquet(path, columns=wanted)
    out: dict[str, object] = {"rows": int(handle.metadata.num_rows), "columns": len(names)}
    if animal_col in frame:
        out["animals"] = int(frame[animal_col].nunique())
    if time_col in frame:
        stamps = pd.to_datetime(frame[time_col])
        out["start"] = str(stamps.min())
        out["end"] = str(stamps.max())
    for col in group_cols:
        if col in frame:
            out[f"n_{col}"] = int(frame[col].dropna().nunique())
    return out


def describe_csv(path: Path, usecols: list[str], animal_col: str | None, time_col: str | None,
                 group_cols: list[str]) -> dict[str, object]:
    header = pd.read_csv(path, nrows=0)
    available = [c for c in usecols if c in header.columns]
    frame = pd.read_csv(path, usecols=available) if available else pd.read_csv(path, nrows=0)
    out: dict[str, object] = {"rows": int(len(frame)), "columns": int(header.shape[1])}
    if animal_col and animal_col in frame:
        out["animals"] = int(frame[animal_col].nunique())
    if time_col and time_col in frame:
        stamps = pd.to_datetime(frame[time_col], errors="coerce")
        out["start"] = str(stamps.min())
        out["end"] = str(stamps.max())
    for col in group_cols:
        if col in frame:
            out[f"n_{col}"] = int(frame[col].dropna().nunique())
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ledger: list[dict[str, object]] = []

    def add(product_id: str, tier: str, role: str, path: Path, unit: str,
            describe: dict[str, object] | None = None, note: str = "") -> None:
        record: dict[str, object] = {
            "product_id": product_id,
            "tier": tier,
            "role": role,
            "row_unit": unit,
            "path": str(path),
        }
        if not path.exists():
            record["exists"] = False
            record["note"] = f"missing; {note}".strip("; ")
            ledger.append(record)
            return
        record.update(file_stamp(path))
        if describe:
            record.update(describe)
        record["note"] = note
        ledger.append(record)

    # -- Tier A: membership products -----------------------------------------
    canonical = CANON / "canonical_hourly_membership_with_association_events.parquet"
    add(
        "membership.canonical_hourly_assoc", "A-membership",
        "FROZEN SOURCE: hourly membership with association events; parent of every"
        " Phase 1-3 product built from here on",
        canonical, "animal x hour",
        describe_parquet(canonical, "animal_id", "window_start",
                         ["origin_group", "dynamic_social_unit", "assigned_social_unit"]),
        note="canonical builder 20260722; night carry on; dynamic rule 7 days / 24 outside hours",
    )

    plain = CANON / "canonical_hourly_membership.parquet"
    add("membership.canonical_hourly_plain", "A-membership",
        "same builder without association-event columns", plain, "animal x hour",
        describe_parquet(plain, "animal_id", "window_start",
                         ["origin_group", "dynamic_social_unit"]),
        note="superseded by the association-event version for all new work")

    inventory_src = INVENTORY_SRC / "canonical_hourly_membership.parquet"
    add("membership.local_2h_support", "A-membership",
        "LEGACY: source of the 13,615-event inventory in canonical_group_merge_scale_log_scatter",
        inventory_src, "animal x hour",
        describe_parquet(inventory_src, "animal_id", "window_start",
                         ["origin_group", "dynamic_social_unit"]),
        note="different builder and start date from the frozen source; must be re-derived")

    prox = PROJECT / "outputs" / "dynamic_social_unit_merge_gamm" / "proximity_status_dynamic_social_unit.parquet"
    add("membership.proximity_status", "A-membership",
        "LEGACY: input to every saved encounter/hurdle model", prox, "animal x interval",
        describe_parquet(prox, "animal_id", "timestamp", ["dynamic_social_unit"]),
        note="21 analysis groups; metadata records 96,349 unmatched status rows")

    narrow = PROJECT / "outputs" / "membership_export_narrow" / "canonical_hourly_membership_narrow.csv"
    add("membership.export_narrow", "A-membership",
        "shared narrow export of the frozen source", narrow, "animal x hour",
        {"rows": 1_924_104, "columns": 20, "animals": 350, "n_origin_group": 26,
         "start": "2024-03-01 04:00:00", "end": "2026-07-22 06:00:00"},
        note="counts taken from its own metadata json; not re-scanned (368 MB)")

    nightly = PROJECT / "outputs" / "membership_export_nightly" / "canonical_nightly_membership.csv"
    add("membership.export_nightly", "A-membership",
        "nightly roll-up of the frozen source", nightly, "animal x night",
        describe_csv(nightly, ["animal_id", "night_date", "origin_group",
                               "nightly_dynamic_social_unit", "nightly_membership_class"],
                     "animal_id", "night_date",
                     ["origin_group", "nightly_dynamic_social_unit", "nightly_membership_class"]),
        note="night = 16:00 through 06:00 UTC, evening-dated")

    # -- Tier B: structural event products -----------------------------------
    assoc_index = CANON / "association_cluster_event_index.csv"
    add("events.association_cluster_index", "B-structural",
        "STAGE 1 SOURCE: hourly association clusters with continuing event ids",
        assoc_index, "hour x cluster",
        describe_csv(assoc_index, ["window_start", "association_event_id", "association_type",
                                   "association_dynamic_units", "temp_group_id"],
                     None, "window_start",
                     ["association_event_id", "association_type", "temp_group_id"]),
        note="derived from the frozen source; association_event_id continues across hours")

    merge_dir = PROJECT / "outputs" / "canonical_group_merge_scale_log_scatter"
    for name, unit, role in [
        ("canonical_group_merge_scale_events.csv", "merge event", "LEGACY merge inventory (74 dyads)"),
        ("canonical_within_group_split_events.csv", "split event", "LEGACY fission inventory"),
        ("canonical_isolated_events.csv", "isolated event", "LEGACY isolation inventory"),
        ("canonical_single_animal_separation_events.csv", "separation event", "LEGACY separation inventory"),
        ("canonical_disperser_events.csv", "dispersal segment", "LEGACY dispersal inventory"),
        ("canonical_event_size_duration_all_events.csv", "event", "LEGACY pooled inventory (13,615)"),
    ]:
        path = merge_dir / name
        add(f"events.{name[10:-4]}", "B-structural", role, path, unit,
            describe_csv(path, ["event_id", "start_time", "end_time", "pair", "origin_group",
                                "animal_id", "event_type"],
                         "animal_id", "start_time", ["pair", "origin_group", "event_type"]),
            note="built from membership.local_2h_support, NOT the frozen source")

    # -- Tier C: fine-scale products -----------------------------------------
    fine = {
        "fine.shuffle_5m": (CANON / "canonical_5m_shared_history_shuffle_expectation"
                            / "canonical_5m_shuffle_expectation_2min_rows.csv",
                            "STAGE 2 SOURCE at 5 m: widest eligible bin set; retains zero-cross bins"),
        "fine.shuffle_2m": (CANON / "canonical_2m_shared_history_shuffle_expectation"
                            / "canonical_5m_shuffle_expectation_2min_rows.csv",
                            "STAGE 2 SOURCE at 2 m: radius is set by the DIRECTORY, not the filename"),
    }
    for product_id, (path, role) in fine.items():
        add(product_id, "C-finescale", role, path, "dyad x 2-min bin",
            describe_csv(path, ["bin_2min", "pair_key", "association_cluster_id",
                                "event_number_for_pair"],
                         None, "bin_2min", ["pair_key", "association_cluster_id"]),
            note="eligible bin = at least one within-radius edge in a shared association cluster")

    bigmerge = (PROJECT / "outputs" / "dynamic_social_unit_merge_gamm" / "group_merge_mixing_dynamics"
                / "bigmerge_dynamic_social_unit_2min_vs_hourly_5m_no_copper_lilac_2min_metric_rows.csv")
    add("fine.bigmerge_no_copper_lilac", "C-finescale",
        "DEFECTIVE INPUT: the only fine-scale file the hurdle model reads; excludes Copper-Lilac by construction",
        bigmerge, "dyad x 2-min bin",
        describe_csv(bigmerge, ["bin_2min", "pair_key", "merge_episode_id"],
                     None, "bin_2min", ["pair_key", "merge_episode_id"]),
        note="filename states the exclusion; candidate days for that dyad remain coded as zero meetings")

    # -- Tier D: model row products ------------------------------------------
    hurdle = PROJECT / "outputs" / "dynamic_social_unit_merge_gamm" / "daily_interaction_hurdle"
    for name, unit, role in [
        ("daily_interaction_hurdle_daily_event_rows.csv", "dyad x day",
         "meeting-probability daily rows; 56,128 candidates, 117 positives"),
        ("daily_interaction_hurdle_model_rows.csv", "dyad x week",
         "meeting-probability weekly rows"),
        ("event_duration_model_rows.csv", "positive episode",
         "duration and integration event rows"),
        ("daily_interaction_positive_episodes.csv", "positive episode",
         "episodes stitched from contact-positive bins only"),
    ]:
        path = hurdle / name
        add(f"model.{name[:-4]}", "D-model", role, path, unit,
            describe_csv(path, ["pair_key", "period_start", "week_start", "episode_start_day",
                                "group_a", "group_b"],
                         None, "period_start", ["pair_key", "group_a", "group_b"]),
            note="rebuild required after Phase 1 and Phase 2")

    ledger_frame = pd.DataFrame(ledger)
    front = ["product_id", "tier", "role", "row_unit", "rows", "animals", "start", "end"]
    ordered = [c for c in front if c in ledger_frame.columns]
    ordered += [c for c in ledger_frame.columns if c not in ordered]
    ledger_frame = ledger_frame[ordered]
    ledger_frame.to_csv(OUT_DIR / "cohort_ledger.csv", index=False)

    # -- event-ID crosswalk ---------------------------------------------------
    crosswalk = pd.DataFrame(
        [
            {
                "identifier": "temp_group_id",
                "source_product": "membership.canonical_hourly_assoc",
                "scope": "one spatial cluster in one hour",
                "continues_across_hours": False,
                "shared_by_all_members": True,
                "use_for": "the exact cluster composition at a given hour",
                "do_not_use_for": "anything that must persist across hours",
            },
            {
                "identifier": "association_cluster_id",
                "source_product": "events.association_cluster_index",
                "scope": "one interpreted association cluster in one hour",
                "continues_across_hours": False,
                "shared_by_all_members": True,
                "use_for": "joining fine-scale 2-min rows to their structural cluster",
                "do_not_use_for": "counting events",
            },
            {
                "identifier": "association_event_id",
                "source_product": "events.association_cluster_index",
                "scope": "one association event, continuing across consecutive hours",
                "continues_across_hours": True,
                "shared_by_all_members": True,
                "use_for": "STAGE 1 encounter identity; the unit Phase 2 counts",
                "do_not_use_for": "fine-scale exposure, which needs the 2-min bin set",
            },
            {
                "identifier": "event_id (merge scale)",
                "source_product": "events.canonical_group_merge_scale_events",
                "scope": "pair x merge_scale x contiguous run, gap <= 2.5 h",
                "continues_across_hours": True,
                "shared_by_all_members": False,
                "use_for": "the legacy 13,615-event inventory",
                "do_not_use_for": "joining to the frozen source; built on a different membership product",
            },
            {
                "identifier": "event_number_for_pair",
                "source_product": "fine.shuffle_5m / fine.shuffle_2m",
                "scope": "sequential encounter counter within a dyad",
                "continues_across_hours": True,
                "shared_by_all_members": True,
                "use_for": "clustering fine-scale bins into encounters, and bootstrap clustering",
                "do_not_use_for": "matching to association_event_id without an explicit join",
            },
            {
                "identifier": "merge_episode_id",
                "source_product": "fine.bigmerge_no_copper_lilac",
                "scope": "big-merge episode under the eps100m / minpair2 / minep6h rule",
                "continues_across_hours": True,
                "shared_by_all_members": True,
                "use_for": "reproducing the saved hurdle results only",
                "do_not_use_for": "any new analysis",
            },
            {
                "identifier": "episode_id (positive)",
                "source_product": "model.daily_interaction_positive_episodes",
                "scope": "run of contact-positive bins bridged up to 14 h",
                "continues_across_hours": True,
                "shared_by_all_members": True,
                "use_for": "reproducing the saved duration and integration results only",
                "do_not_use_for": "duration or exposure statements; it omits supported zero-contact bins",
            },
        ]
    )
    crosswalk.to_csv(OUT_DIR / "event_id_crosswalk.csv", index=False)

    summary = {
        "phase": "0 - cohort ledger and definitions",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "frozen_source": str(canonical),
        "products_inventoried": int(len(ledger_frame)),
        "products_missing": int((~ledger_frame.get("exists", pd.Series(True, index=ledger_frame.index)).fillna(False)).sum()),
        "tiers": ledger_frame["tier"].value_counts().to_dict(),
        "outputs": {
            "cohort_ledger": str(OUT_DIR / "cohort_ledger.csv"),
            "event_id_crosswalk": str(OUT_DIR / "event_id_crosswalk.csv"),
        },
    }
    (OUT_DIR / "phase0_metadata.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(ledger_frame[[c for c in ["product_id", "tier", "rows", "animals", "start", "end"]
                        if c in ledger_frame.columns]].to_string(index=False))
    print()
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
