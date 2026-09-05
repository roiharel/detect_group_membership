"""Phase 2: the two-stage encounter unit.

Stage 1  structural encounter -- two dynamic social units co-present in one
         spatial cluster, stitched across consecutive hours. Detected without
         reference to close-range contact.
Stage 2  mixing within the encounter -- every eligible fine-scale 2-minute bin
         inside the stage-1 window is retained, including bins with zero
         cross-group contact.

Three durations are reported separately and must not be interchanged:
    structural_span_hours   last encounter hour minus first, plus one hour
    supported_exposure_h    eligible fine-scale bins x 2 minutes
    active_contact_hours    contact-positive bins x 2 minutes

Run from the project root:
    python phase2_two_stage_events.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(r"C:\Users\rharel\Documents\group_mebership")
UPSTREAM = Path(r"C:\Users\rharel\Documents\New project")
# the canonical run folder is overridable, so a rebuild does not need a code edit
CANON = Path(os.environ.get(
    "REBUILD_CANON",
    str(UPSTREAM / "outputs" / "canonical_robust_hourly_membership_shared_full_20260722")))
PHASE1 = PROJECT / "outputs" / "general_structure_2026_09" / "phase1_opportunity" / "opportunity_dyad_day.csv"
MEMBERSHIP = CANON / "canonical_hourly_membership_with_association_events.parquet"
FINE = {
    "5m": CANON / "canonical_5m_shared_history_shuffle_expectation" / "canonical_5m_shuffle_expectation_2min_rows.csv",
    "2m": CANON / "canonical_2m_shared_history_shuffle_expectation" / "canonical_5m_shuffle_expectation_2min_rows.csv",
}
OUT_DIR = PROJECT / "outputs" / "general_structure_2026_09" / "phase2_two_stage_events"

PRIMARY_GAP_HOURS = 3.0
SENSITIVITY_GAP_HOURS = [2.0, 3.0, 6.0, 14.0]
BIN_MINUTES = 2.0
LARGE_FRACTION = 0.5
SMALL_FRACTION = 0.25

# A stage-1 run this long is not a discrete encounter; it is a sustained
# association. The threshold matches the 7-day persistence rule the membership
# builder already uses for dynamic reassignment, so the two agree by design.
SUSTAINED_SPAN_HOURS = 168.0
EXAMPLE_MAX_SPAN_HOURS = 72.0


def parse_counts(value: object) -> dict[str, int]:
    if not isinstance(value, str) or not value:
        return {}
    counts: dict[str, int] = {}
    for piece in value.split(";"):
        if ":" not in piece:
            continue
        group, raw = piece.rsplit(":", 1)
        try:
            counts[str(group)] = int(raw)
        except ValueError:
            continue
    return counts


def load_encounter_hours() -> pd.DataFrame:
    """Dyad-hours in which both units shared a spatial cluster, with composition."""
    columns = ["animal_id", "window_start", "dynamic_social_unit", "is_observed",
               "temp_group_id", "temp_group_dynamic_counts"]
    data = pd.read_parquet(MEMBERSHIP, columns=columns)
    data["window_start"] = pd.to_datetime(data["window_start"])
    data = data.dropna(subset=["dynamic_social_unit"]).copy()
    data["dynamic_social_unit"] = data["dynamic_social_unit"].astype(str)

    observed = (
        data.groupby(["window_start", "dynamic_social_unit"], observed=True)["animal_id"]
        .nunique()
        .rename("unit_observed_n")
        .reset_index()
    )
    observed_lookup = {
        stamp: dict(zip(frame["dynamic_social_unit"], frame["unit_observed_n"]))
        for stamp, frame in observed.groupby("window_start", observed=True)
    }

    clusters = (
        data.dropna(subset=["temp_group_id"])
        .drop_duplicates(["window_start", "temp_group_id"])
        [["window_start", "temp_group_id", "temp_group_dynamic_counts"]]
    )

    from itertools import combinations

    rows: list[dict[str, object]] = []
    for record in clusters.itertuples(index=False):
        counts = parse_counts(record.temp_group_dynamic_counts)
        if len(counts) < 2:
            continue
        present = observed_lookup.get(record.window_start, {})
        cluster_size = int(sum(counts.values()))
        for unit_a, unit_b in combinations(sorted(counts), 2):
            observed_a = present.get(unit_a, np.nan)
            observed_b = present.get(unit_b, np.nan)
            if not observed_a or not observed_b:
                continue
            rows.append(
                {
                    "window_start": record.window_start,
                    "unit_a": unit_a,
                    "unit_b": unit_b,
                    "pair_key": f"{unit_a} - {unit_b}",
                    "n_a_in_cluster": counts[unit_a],
                    "n_b_in_cluster": counts[unit_b],
                    "observed_a": int(observed_a),
                    "observed_b": int(observed_b),
                    "frac_a": counts[unit_a] / observed_a,
                    "frac_b": counts[unit_b] / observed_b,
                    "cluster_size": cluster_size,
                }
            )
    hours = pd.DataFrame(rows)
    # one row per dyad-hour: keep the cluster in which the pair is best represented
    hours["pair_in_cluster"] = hours["n_a_in_cluster"] + hours["n_b_in_cluster"]
    hours = (
        hours.sort_values("pair_in_cluster", ascending=False)
        .drop_duplicates(["window_start", "pair_key"])
        .sort_values(["pair_key", "window_start"])
        .reset_index(drop=True)
    )
    return hours


def stitch_events(hours: pd.DataFrame, gap_hours: float) -> pd.DataFrame:
    """Group consecutive encounter hours of a dyad into stage-1 events."""
    rows = hours.sort_values(["pair_key", "window_start"]).copy()
    gaps = rows.groupby("pair_key", observed=True)["window_start"].diff().dt.total_seconds().div(3600)
    starts_new = gaps.isna() | gaps.gt(gap_hours)
    rows["event_number_in_dyad"] = starts_new.groupby(rows["pair_key"], observed=True).cumsum()
    rows["stage1_event_id"] = (
        rows["pair_key"].str.replace(" - ", "__", regex=False)
        + "__E"
        + rows["event_number_in_dyad"].astype(int).astype(str).str.zfill(4)
    )
    events = (
        rows.groupby("stage1_event_id", observed=True)
        .agg(
            pair_key=("pair_key", "first"),
            unit_a=("unit_a", "first"),
            unit_b=("unit_b", "first"),
            event_number_in_dyad=("event_number_in_dyad", "first"),
            start_hour=("window_start", "min"),
            end_hour=("window_start", "max"),
            encounter_hours=("window_start", "nunique"),
            max_frac_a=("frac_a", "max"),
            max_frac_b=("frac_b", "max"),
            median_frac_a=("frac_a", "median"),
            median_frac_b=("frac_b", "median"),
            max_n_a_in_cluster=("n_a_in_cluster", "max"),
            max_n_b_in_cluster=("n_b_in_cluster", "max"),
            max_cluster_size=("cluster_size", "max"),
            mean_observed_a=("observed_a", "mean"),
            mean_observed_b=("observed_b", "mean"),
        )
        .reset_index()
    )
    events["structural_span_hours"] = (
        (events["end_hour"] - events["start_hour"]).dt.total_seconds().div(3600) + 1.0
    )
    events["observed_hour_fraction"] = events["encounter_hours"] / events["structural_span_hours"]
    events["merge_scale"] = np.where(
        (events["median_frac_a"] >= LARGE_FRACTION) & (events["median_frac_b"] >= LARGE_FRACTION),
        "large_merge",
        np.where(
            (events["median_frac_a"] < SMALL_FRACTION) | (events["median_frac_b"] < SMALL_FRACTION),
            "small_subset_merge",
            "medium_partial_merge",
        ),
    )
    events["encounter_class"] = np.where(
        events["structural_span_hours"] >= SUSTAINED_SPAN_HOURS,
        "sustained_association",
        "discrete_encounter",
    )
    events["gap_rule_hours"] = gap_hours
    return events, rows


def attach_stage2(events: pd.DataFrame, radius: str) -> pd.DataFrame:
    """Aggregate every eligible fine-scale bin inside each stage-1 window."""
    usecols = ["bin_2min", "pair_key", "cross_edges", "total_edges",
               "observed_cross_edge_fraction", "shuffle_expected_cross_edge_fraction",
               "edge_modularity_q", "composition_entropy_norm", "pair_balance"]
    fine = pd.read_csv(FINE[radius], usecols=usecols, parse_dates=["bin_2min"])
    fine["pair_key"] = fine["pair_key"].astype(str)
    fine = fine.sort_values(["pair_key", "bin_2min"])

    windows = events[["stage1_event_id", "pair_key", "start_hour", "end_hour"]].copy()
    windows["window_end_exclusive"] = windows["end_hour"] + pd.Timedelta(hours=1)

    assigned: list[pd.DataFrame] = []
    for pair, bins in fine.groupby("pair_key", observed=True):
        pair_windows = windows[windows["pair_key"].eq(pair)]
        if pair_windows.empty:
            continue
        matched = pd.merge_asof(
            bins.sort_values("bin_2min"),
            pair_windows.sort_values("start_hour")[["start_hour", "window_end_exclusive", "stage1_event_id"]],
            left_on="bin_2min",
            right_on="start_hour",
            direction="backward",
        )
        matched = matched[matched["bin_2min"] < matched["window_end_exclusive"]]
        assigned.append(matched)
    if not assigned:
        return pd.DataFrame()
    matched = pd.concat(assigned, ignore_index=True)
    matched["is_positive_bin"] = matched["cross_edges"].gt(0)

    summary = (
        matched.groupby("stage1_event_id", observed=True)
        .agg(
            eligible_bins=("bin_2min", "nunique"),
            positive_bins=("is_positive_bin", "sum"),
            cross_edges=("cross_edges", "sum"),
            total_edges=("total_edges", "sum"),
            mean_shuffle_expected=("shuffle_expected_cross_edge_fraction", "mean"),
            mean_bin_obs_minus_shuffle=(
                "observed_cross_edge_fraction",
                lambda s: np.nan,
            ),
            mean_modularity=("edge_modularity_q", "mean"),
            mean_entropy=("composition_entropy_norm", "mean"),
            mean_balance=("pair_balance", "mean"),
        )
        .reset_index()
    )
    # bin-weighted difference, computed explicitly so the weighting is auditable
    matched["bin_obs_minus_shuffle"] = (
        matched["observed_cross_edge_fraction"] - matched["shuffle_expected_cross_edge_fraction"]
    )
    diffs = (
        matched.groupby("stage1_event_id", observed=True)["bin_obs_minus_shuffle"]
        .mean()
        .rename("mean_bin_obs_minus_shuffle")
        .reset_index()
    )
    summary = summary.drop(columns=["mean_bin_obs_minus_shuffle"]).merge(diffs, on="stage1_event_id")

    summary["positive_bins"] = summary["positive_bins"].astype(int)
    summary["supported_exposure_hours"] = summary["eligible_bins"] * BIN_MINUTES / 60.0
    summary["active_contact_hours"] = summary["positive_bins"] * BIN_MINUTES / 60.0
    summary["edge_weighted_cross_fraction"] = np.where(
        summary["total_edges"].gt(0), summary["cross_edges"] / summary["total_edges"], np.nan
    )
    summary["positive_bin_fraction"] = summary["positive_bins"] / summary["eligible_bins"]
    return summary.add_prefix(f"{radius}_").rename(columns={f"{radius}_stage1_event_id": "stage1_event_id"})


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("loading encounter hours ...")
    hours = load_encounter_hours()
    print(f"  encounter dyad-hours={len(hours):,}  dyads={hours['pair_key'].nunique()}")

    sensitivity: list[dict[str, object]] = []
    primary_events = None
    primary_hours = None
    for gap in SENSITIVITY_GAP_HOURS:
        events, tagged = stitch_events(hours, gap)
        sensitivity.append(
            {
                "gap_rule_hours": gap,
                "events": int(len(events)),
                "dyads": int(events["pair_key"].nunique()),
                "median_structural_span_hours": float(events["structural_span_hours"].median()),
                "mean_structural_span_hours": float(events["structural_span_hours"].mean()),
                "max_structural_span_hours": float(events["structural_span_hours"].max()),
                "median_encounter_hours": float(events["encounter_hours"].median()),
                "large_merge": int((events["merge_scale"] == "large_merge").sum()),
                "medium_partial_merge": int((events["merge_scale"] == "medium_partial_merge").sum()),
                "small_subset_merge": int((events["merge_scale"] == "small_subset_merge").sum()),
                "discrete_encounters": int((events["encounter_class"] == "discrete_encounter").sum()),
                "sustained_associations": int((events["encounter_class"] == "sustained_association").sum()),
                "sustained_dyads": int(
                    events.loc[events["encounter_class"] == "sustained_association", "pair_key"].nunique()
                ),
            }
        )
        if gap == PRIMARY_GAP_HOURS:
            primary_events, primary_hours = events, tagged
    sensitivity_frame = pd.DataFrame(sensitivity)
    sensitivity_frame.to_csv(OUT_DIR / "stage1_gap_rule_sensitivity.csv", index=False)

    events = primary_events
    print(f"  stage-1 events at gap={PRIMARY_GAP_HOURS} h: {len(events):,} across {events['pair_key'].nunique()} dyads")

    print("attaching stage-2 mixing ...")
    for radius in ["5m", "2m"]:
        stage2 = attach_stage2(events, radius)
        events = events.merge(stage2, on="stage1_event_id", how="left")
        covered = int(events[f"{radius}_eligible_bins"].notna().sum())
        print(f"  {radius}: {covered:,} of {len(events):,} events have fine-scale bins")

    events["has_finescale_5m"] = events["5m_eligible_bins"].notna()
    events["has_finescale_2m"] = events["2m_eligible_bins"].notna()

    front = [
        "stage1_event_id", "pair_key", "unit_a", "unit_b", "encounter_class", "merge_scale",
        "start_hour", "end_hour", "structural_span_hours", "encounter_hours",
        "observed_hour_fraction",
        "5m_supported_exposure_hours", "5m_active_contact_hours",
        "5m_eligible_bins", "5m_positive_bins", "5m_positive_bin_fraction",
        "5m_cross_edges", "5m_total_edges", "5m_edge_weighted_cross_fraction",
        "5m_mean_shuffle_expected", "5m_mean_bin_obs_minus_shuffle",
        "2m_eligible_bins", "2m_positive_bins", "2m_edge_weighted_cross_fraction",
        "has_finescale_5m", "has_finescale_2m",
    ]
    ordered = [c for c in front if c in events.columns] + [c for c in events.columns if c not in front]
    events = events[ordered].sort_values(["pair_key", "start_hour"]).reset_index(drop=True)
    events.to_csv(OUT_DIR / "stage1_events_with_stage2_mixing.csv", index=False)
    primary_hours.to_csv(OUT_DIR / "stage1_encounter_hours.csv", index=False)

    # -- worked example: retained versus omitted bins -------------------------
    fine5 = pd.read_csv(FINE["5m"], usecols=["bin_2min", "pair_key", "cross_edges", "total_edges"],
                        parse_dates=["bin_2min"])
    with_bins = events[
        events["has_finescale_5m"]
        & events["5m_positive_bins"].gt(0)
        & events["encounter_class"].eq("discrete_encounter")
        & events["structural_span_hours"].le(EXAMPLE_MAX_SPAN_HOURS)
    ]
    example = (
        with_bins.sort_values("5m_eligible_bins", ascending=False)
        .head(1)
        .iloc[0]
    )
    window_bins = fine5[
        fine5["pair_key"].eq(example["pair_key"])
        & fine5["bin_2min"].ge(example["start_hour"])
        & fine5["bin_2min"].lt(example["end_hour"] + pd.Timedelta(hours=1))
    ].copy()
    window_bins["retained_by_phase2"] = True
    window_bins["retained_by_saved_pipeline"] = window_bins["cross_edges"].gt(0)
    window_bins["stage1_event_id"] = example["stage1_event_id"]
    window_bins.sort_values("bin_2min").to_csv(OUT_DIR / "worked_example_bins.csv", index=False)

    example_note = {
        "stage1_event_id": str(example["stage1_event_id"]),
        "pair_key": str(example["pair_key"]),
        "start_hour": str(example["start_hour"]),
        "end_hour": str(example["end_hour"]),
        "structural_span_hours": float(example["structural_span_hours"]),
        "eligible_bins_retained_by_phase2": int(len(window_bins)),
        "bins_retained_by_saved_pipeline": int(window_bins["retained_by_saved_pipeline"].sum()),
        "bins_dropped_by_saved_pipeline": int((~window_bins["retained_by_saved_pipeline"]).sum()),
        "share_of_exposure_dropped_by_saved_pipeline": round(
            float((~window_bins["retained_by_saved_pipeline"]).mean()), 4
        ),
    }

    # -- dyad-level coverage --------------------------------------------------
    dyad = (
        events.groupby("pair_key", observed=True)
        .agg(
            stage1_events=("stage1_event_id", "size"),
            discrete_encounters=("encounter_class", lambda s: int((s == "discrete_encounter").sum())),
            sustained_associations=("encounter_class", lambda s: int((s == "sustained_association").sum())),
            encounter_hours=("encounter_hours", "sum"),
            structural_span_hours=("structural_span_hours", "sum"),
            max_structural_span_hours=("structural_span_hours", "max"),
            large_merge_events=("merge_scale", lambda s: int((s == "large_merge").sum())),
            events_with_5m_bins=("has_finescale_5m", "sum"),
            eligible_bins_5m=("5m_eligible_bins", "sum"),
            positive_bins_5m=("5m_positive_bins", "sum"),
            cross_edges_5m=("5m_cross_edges", "sum"),
            total_edges_5m=("5m_total_edges", "sum"),
        )
        .reset_index()
        .sort_values("stage1_events", ascending=False)
    )
    dyad["edge_weighted_cross_fraction_5m"] = np.where(
        dyad["total_edges_5m"].gt(0), dyad["cross_edges_5m"] / dyad["total_edges_5m"], np.nan
    )
    dyad.to_csv(OUT_DIR / "stage1_dyad_coverage.csv", index=False)

    covered5 = events["has_finescale_5m"]
    metadata = {
        "phase": "2 - two-stage encounter unit",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "membership_source": str(MEMBERSHIP),
        "finescale_sources": {k: str(v) for k, v in FINE.items()},
        "stage1_definition": "two dynamic social units co-present in one spatial cluster, "
                             f"consecutive encounter hours stitched when separated by <= {PRIMARY_GAP_HOURS} h",
        "stage2_definition": "every eligible fine-scale 2-minute bin inside the stage-1 window, "
                             "including bins with zero cross-group contact",
        "duration_fields": {
            "structural_span_hours": "last encounter hour minus first, plus one hour",
            "supported_exposure_hours": "eligible fine-scale bins x 2 minutes",
            "active_contact_hours": "contact-positive bins x 2 minutes",
        },
        "encounter_class_rule": f"structural span >= {SUSTAINED_SPAN_HOURS:.0f} h is a sustained "
                                "association, not a discrete encounter; threshold matches the "
                                "membership builder's 7-day dynamic-reassignment rule",
        "primary_gap_rule_hours": PRIMARY_GAP_HOURS,
        "stage1_events": int(len(events)),
        "stage1_dyads": int(events["pair_key"].nunique()),
        "stage1_events_by_class": events["encounter_class"].value_counts().to_dict(),
        "sustained_association_dyads": sorted(
            events.loc[events["encounter_class"].eq("sustained_association"), "pair_key"].unique().tolist()
        ),
        "discrete_encounter_dyads": int(
            events.loc[events["encounter_class"].eq("discrete_encounter"), "pair_key"].nunique()
        ),
        "median_discrete_encounter_span_hours": float(
            events.loc[events["encounter_class"].eq("discrete_encounter"), "structural_span_hours"].median()
        ),
        "stage1_events_by_scale": events["merge_scale"].value_counts().to_dict(),
        "events_with_5m_finescale": int(covered5.sum()),
        "events_with_2m_finescale": int(events["has_finescale_2m"].sum()),
        "dyads_with_5m_finescale": int(events.loc[covered5, "pair_key"].nunique()),
        "total_eligible_bins_5m": int(events["5m_eligible_bins"].fillna(0).sum()),
        "total_positive_bins_5m": int(events["5m_positive_bins"].fillna(0).sum()),
        "share_of_5m_exposure_with_zero_cross_contact": round(
            float(1.0 - events["5m_positive_bins"].fillna(0).sum() / events["5m_eligible_bins"].fillna(0).sum()), 4
        ),
        "median_structural_span_hours": float(events["structural_span_hours"].median()),
        "median_supported_exposure_hours_where_measured": float(
            events.loc[covered5, "5m_supported_exposure_hours"].median()
        ),
        "median_active_contact_hours_where_measured": float(
            events.loc[covered5, "5m_active_contact_hours"].median()
        ),
        "worked_example": example_note,
        "gap_rule_sensitivity": sensitivity,
        "outputs": {
            "events": str(OUT_DIR / "stage1_events_with_stage2_mixing.csv"),
            "encounter_hours": str(OUT_DIR / "stage1_encounter_hours.csv"),
            "dyad_coverage": str(OUT_DIR / "stage1_dyad_coverage.csv"),
            "gap_sensitivity": str(OUT_DIR / "stage1_gap_rule_sensitivity.csv"),
            "worked_example_bins": str(OUT_DIR / "worked_example_bins.csv"),
        },
    }
    (OUT_DIR / "phase2_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print()
    print(sensitivity_frame.to_string(index=False))
    print()
    print(json.dumps({k: v for k, v in metadata.items()
                      if k not in {"outputs", "gap_rule_sensitivity"}}, indent=2))


if __name__ == "__main__":
    main()
