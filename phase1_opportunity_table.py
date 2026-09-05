"""Phase 1: dyad-day opportunity table with explicit observation states.

Replaces the zero-filled candidate-day construction used by the saved hurdle
model. Every dyad-day carries shared observation duration and one of four
explicit states, so a supported negative is distinguishable from an unobserved
window and from an upstream exclusion.

Unit: one unordered pair of dynamic social units on one calendar day (UTC).

States
    detected_encounter    both units observed and co-present in one spatial
                          cluster in at least one shared hour
    observed_no_encounter both units observed in at least one shared hour with
                          exact-hour support, never in the same cluster
    insufficient_support  no shared hour, or no shared hour with exact-hour
                          support for both units
    excluded              dyad removed by an upstream rule of the pipeline
                          being audited (currently: the fine-scale input that
                          drops Copper-Lilac by construction)

Run from the project root:
    python phase1_opportunity_table.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from itertools import combinations
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
MEMBERSHIP = CANON / "canonical_hourly_membership_with_association_events.parquet"
FINE_5M = (CANON / "canonical_5m_shared_history_shuffle_expectation"
           / "canonical_5m_shuffle_expectation_2min_rows.csv")
BIGMERGE = (PROJECT / "outputs" / "dynamic_social_unit_merge_gamm" / "group_merge_mixing_dynamics"
            / "bigmerge_dynamic_social_unit_2min_vs_hourly_5m_no_copper_lilac_2min_metric_rows.csv")
SAVED_CANDIDATES = (PROJECT / "outputs" / "dynamic_social_unit_merge_gamm" / "daily_interaction_hurdle"
                    / "daily_interaction_hurdle_daily_event_rows.csv")
OUT_DIR = PROJECT / "outputs" / "general_structure_2026_09" / "phase1_opportunity"

EARTH_RADIUS_M = 6_371_000.0
DAY_START_HOUR = 3
DAY_END_HOUR = 16
UPSTREAM_EXCLUDED_DYADS = {"Copper - Lilac"}


def pair_key(unit_a: str, unit_b: str) -> str:
    left, right = sorted((str(unit_a), str(unit_b)))
    return f"{left} - {right}"


def haversine_m(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = phi2 - phi1
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2.0) ** 2
    return 2.0 * EARTH_RADIUS_M * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def parse_counts(value: object) -> dict[str, int]:
    """Parse a 'Group:n;Group:n' composition string."""
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


def load_hour_units() -> pd.DataFrame:
    """One row per dynamic social unit per hour, with centroid and support counts."""
    columns = [
        "animal_id", "window_start", "dynamic_social_unit", "is_observed",
        "longitude", "latitude", "temp_group_id", "temp_group_dynamic_counts",
    ]
    data = pd.read_parquet(MEMBERSHIP, columns=columns)
    data = data.dropna(subset=["dynamic_social_unit", "longitude", "latitude"]).copy()
    data["window_start"] = pd.to_datetime(data["window_start"])
    data["dynamic_social_unit"] = data["dynamic_social_unit"].astype(str)
    data["is_observed"] = data["is_observed"].fillna(False).astype(bool)

    hour_units = (
        data.groupby(["window_start", "dynamic_social_unit"], observed=True)
        .agg(
            n_animals=("animal_id", "nunique"),
            n_observed_animals=("is_observed", "sum"),
            centroid_latitude=("latitude", "mean"),
            centroid_longitude=("longitude", "mean"),
        )
        .reset_index()
    )
    hour_units["n_observed_animals"] = hour_units["n_observed_animals"].astype(int)
    return data, hour_units


def build_cluster_pairs(data: pd.DataFrame) -> pd.DataFrame:
    """Dyad-hours in which both units were members of the same spatial cluster."""
    clusters = (
        data.dropna(subset=["temp_group_id"])
        .drop_duplicates(["window_start", "temp_group_id"])
        [["window_start", "temp_group_id", "temp_group_dynamic_counts"]]
    )
    rows: list[dict[str, object]] = []
    for record in clusters.itertuples(index=False):
        counts = parse_counts(record.temp_group_dynamic_counts)
        if len(counts) < 2:
            continue
        for unit_a, unit_b in combinations(sorted(counts), 2):
            rows.append(
                {
                    "window_start": record.window_start,
                    "unit_a": unit_a,
                    "unit_b": unit_b,
                    "n_a_in_cluster": counts[unit_a],
                    "n_b_in_cluster": counts[unit_b],
                    "cluster_size": int(sum(counts.values())),
                }
            )
    if not rows:
        return pd.DataFrame(columns=["window_start", "unit_a", "unit_b", "n_a_in_cluster",
                                     "n_b_in_cluster", "cluster_size"])
    pairs = pd.DataFrame(rows)
    return (
        pairs.groupby(["window_start", "unit_a", "unit_b"], observed=True)
        .agg(
            n_a_in_cluster=("n_a_in_cluster", "max"),
            n_b_in_cluster=("n_b_in_cluster", "max"),
            cluster_size=("cluster_size", "max"),
        )
        .reset_index()
    )


def build_dyad_hours(hour_units: pd.DataFrame) -> pd.DataFrame:
    """All unordered unit pairs co-observed within the same hour."""
    left = hour_units.rename(columns={
        "dynamic_social_unit": "unit_a", "n_animals": "n_animals_a",
        "n_observed_animals": "n_observed_a", "centroid_latitude": "lat_a",
        "centroid_longitude": "lon_a",
    })
    right = hour_units.rename(columns={
        "dynamic_social_unit": "unit_b", "n_animals": "n_animals_b",
        "n_observed_animals": "n_observed_b", "centroid_latitude": "lat_b",
        "centroid_longitude": "lon_b",
    })
    merged = left.merge(right, on="window_start")
    merged = merged[merged["unit_a"] < merged["unit_b"]].copy()
    merged["centroid_distance_m"] = haversine_m(
        merged["lat_a"].to_numpy(), merged["lon_a"].to_numpy(),
        merged["lat_b"].to_numpy(), merged["lon_b"].to_numpy(),
    )
    merged["both_exact_support"] = (merged["n_observed_a"] > 0) & (merged["n_observed_b"] > 0)
    hour = merged["window_start"].dt.hour
    merged["is_daytime"] = (hour >= DAY_START_HOUR) & (hour < DAY_END_HOUR)
    return merged


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("loading membership ...")
    data, hour_units = load_hour_units()
    print(f"  hours={hour_units['window_start'].nunique():,}  unit-hours={len(hour_units):,}")

    print("building cluster co-membership ...")
    cluster_pairs = build_cluster_pairs(data)
    print(f"  cluster dyad-hours={len(cluster_pairs):,}")

    print("building dyad-hours ...")
    dyad_hours = build_dyad_hours(hour_units)
    print(f"  dyad-hours={len(dyad_hours):,}")

    dyad_hours = dyad_hours.merge(cluster_pairs, on=["window_start", "unit_a", "unit_b"], how="left")
    dyad_hours["is_cluster_hour"] = dyad_hours["n_a_in_cluster"].notna()
    dyad_hours["period_start"] = dyad_hours["window_start"].dt.floor("D")

    print("aggregating to dyad-days ...")
    grouped = dyad_hours.groupby(["period_start", "unit_a", "unit_b"], observed=True)
    days = grouped.agg(
        shared_hours=("window_start", "nunique"),
        shared_hours_exact=("both_exact_support", "sum"),
        shared_daytime_hours=("is_daytime", "sum"),
        encounter_hours=("is_cluster_hour", "sum"),
        min_centroid_distance_m=("centroid_distance_m", "min"),
        median_centroid_distance_m=("centroid_distance_m", "median"),
        max_animals_a=("n_animals_a", "max"),
        max_animals_b=("n_animals_b", "max"),
        max_observed_a=("n_observed_a", "max"),
        max_observed_b=("n_observed_b", "max"),
        max_n_a_in_cluster=("n_a_in_cluster", "max"),
        max_n_b_in_cluster=("n_b_in_cluster", "max"),
        max_cluster_size=("cluster_size", "max"),
    ).reset_index()

    for col in ["shared_hours_exact", "shared_daytime_hours", "encounter_hours"]:
        days[col] = days[col].astype(int)
    days["pair_key"] = [pair_key(a, b) for a, b in zip(days["unit_a"], days["unit_b"])]

    # -- fine-scale support ---------------------------------------------------
    fine = pd.read_csv(FINE_5M, usecols=["bin_2min", "pair_key"], parse_dates=["bin_2min"])
    fine["period_start"] = fine["bin_2min"].dt.floor("D")
    fine_days = (
        fine.groupby(["period_start", "pair_key"], observed=True)
        .agg(finescale_5m_bins=("bin_2min", "nunique"))
        .reset_index()
    )
    finescale_dyads = set(fine["pair_key"].unique())
    days = days.merge(fine_days, on=["period_start", "pair_key"], how="left")
    days["finescale_5m_bins"] = days["finescale_5m_bins"].fillna(0).astype(int)
    days["dyad_has_finescale_product"] = days["pair_key"].isin(finescale_dyads)

    bigmerge = pd.read_csv(BIGMERGE, usecols=["pair_key"])
    bigmerge_dyads = set(bigmerge["pair_key"].astype(str).unique())
    days["dyad_in_hurdle_input"] = days["pair_key"].isin(bigmerge_dyads)

    # -- state assignment -----------------------------------------------------
    state = pd.Series("observed_no_encounter", index=days.index, dtype=object)
    state[days["shared_hours_exact"].eq(0)] = "insufficient_support"
    state[days["encounter_hours"].gt(0)] = "detected_encounter"
    excluded = days["pair_key"].isin(UPSTREAM_EXCLUDED_DYADS)
    days["upstream_excluded_from_hurdle_input"] = excluded
    days["state"] = state
    days["state_as_audited_pipeline"] = np.where(excluded, "excluded", state)

    days["encounter_hour_fraction"] = np.where(
        days["shared_hours"].gt(0), days["encounter_hours"] / days["shared_hours"], np.nan
    )

    ordered = [
        "period_start", "pair_key", "unit_a", "unit_b", "state", "state_as_audited_pipeline",
        "shared_hours", "shared_hours_exact", "shared_daytime_hours", "encounter_hours",
        "encounter_hour_fraction", "min_centroid_distance_m", "median_centroid_distance_m",
        "max_animals_a", "max_animals_b", "max_observed_a", "max_observed_b",
        "max_n_a_in_cluster", "max_n_b_in_cluster", "max_cluster_size",
        "finescale_5m_bins", "dyad_has_finescale_product", "dyad_in_hurdle_input",
        "upstream_excluded_from_hurdle_input",
    ]
    days = days[ordered].sort_values(["pair_key", "period_start"]).reset_index(drop=True)
    days.to_csv(OUT_DIR / "opportunity_dyad_day.csv", index=False)

    # -- summaries ------------------------------------------------------------
    state_summary = (
        days.groupby("state", observed=True)
        .agg(
            dyad_days=("state", "size"),
            dyads=("pair_key", "nunique"),
            median_shared_hours=("shared_hours", "median"),
            median_min_distance_m=("min_centroid_distance_m", "median"),
            median_median_distance_m=("median_centroid_distance_m", "median"),
        )
        .reset_index()
        .sort_values("dyad_days", ascending=False)
    )
    state_summary["share_of_dyad_days"] = (
        state_summary["dyad_days"] / len(days)
    ).round(5)
    state_summary.to_csv(OUT_DIR / "opportunity_state_summary.csv", index=False)

    dyad_summary = (
        days.groupby("pair_key", observed=True)
        .agg(
            dyad_days=("state", "size"),
            detected_encounter_days=("state", lambda s: int((s == "detected_encounter").sum())),
            observed_no_encounter_days=("state", lambda s: int((s == "observed_no_encounter").sum())),
            insufficient_support_days=("state", lambda s: int((s == "insufficient_support").sum())),
            encounter_hours=("encounter_hours", "sum"),
            min_distance_m=("min_centroid_distance_m", "min"),
            median_distance_m=("median_centroid_distance_m", "median"),
            finescale_5m_bins=("finescale_5m_bins", "sum"),
            dyad_has_finescale_product=("dyad_has_finescale_product", "max"),
            dyad_in_hurdle_input=("dyad_in_hurdle_input", "max"),
        )
        .reset_index()
        .sort_values("detected_encounter_days", ascending=False)
    )
    dyad_summary.to_csv(OUT_DIR / "opportunity_dyad_summary.csv", index=False)

    # -- validation against the saved candidate table -------------------------
    saved = pd.read_csv(
        SAVED_CANDIDATES,
        usecols=["period_start", "pair_key", "any_interaction", "eligible_2min_bins",
                 "daily_centroid_distance_m"],
        parse_dates=["period_start"],
    )
    saved = saved.rename(columns={
        "any_interaction": "saved_any_interaction",
        "eligible_2min_bins": "saved_eligible_2min_bins",
        "daily_centroid_distance_m": "saved_daily_centroid_distance_m",
    })
    compare = days.merge(saved, on=["period_start", "pair_key"], how="outer", indicator=True)
    crosstab = (
        compare.assign(
            new_state=compare["state"].fillna("absent_from_phase1"),
            saved_call=np.where(
                compare["_merge"].eq("right_only"), "absent_from_phase1",
                np.where(compare["saved_any_interaction"].fillna(-1).eq(1), "saved_positive",
                         np.where(compare["saved_any_interaction"].fillna(-1).eq(0), "saved_zero",
                                  "absent_from_saved")),
            ),
        )
        .groupby(["new_state", "saved_call"], observed=True)
        .size()
        .rename("dyad_days")
        .reset_index()
        .sort_values("dyad_days", ascending=False)
    )
    crosstab.to_csv(OUT_DIR / "opportunity_vs_saved_candidates.csv", index=False)

    copper_lilac = compare[compare["pair_key"].eq("Copper - Lilac")]
    cl_states = copper_lilac["state"].value_counts(dropna=False).to_dict()
    cl_saved_positive = int(copper_lilac["saved_any_interaction"].fillna(0).sum())
    copper_lilac.sort_values("period_start").to_csv(OUT_DIR / "validation_copper_lilac_days.csv", index=False)

    metadata = {
        "phase": "1 - dyad-day opportunity table",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "membership_source": str(MEMBERSHIP),
        "row_unit": "unordered pair of dynamic social units on one UTC calendar day",
        "state_rules": {
            "detected_encounter": "both units in one spatial cluster in >= 1 shared hour",
            "observed_no_encounter": ">= 1 shared hour with exact-hour support for both units, never co-clustered",
            "insufficient_support": "no shared hour with exact-hour support for both units",
            "excluded": "dyad removed by an upstream rule of the audited pipeline",
        },
        "daytime_window_hours": [DAY_START_HOUR, DAY_END_HOUR],
        "upstream_excluded_dyads": sorted(UPSTREAM_EXCLUDED_DYADS),
        "dyad_days": int(len(days)),
        "dyads": int(days["pair_key"].nunique()),
        "social_units": int(pd.concat([days["unit_a"], days["unit_b"]]).nunique()),
        "days": int(days["period_start"].nunique()),
        "state_counts": days["state"].value_counts().to_dict(),
        "dyads_with_any_detected_encounter": int(
            days.loc[days["state"].eq("detected_encounter"), "pair_key"].nunique()
        ),
        "dyads_with_finescale_product": int(len(finescale_dyads)),
        "dyads_in_hurdle_input": int(len(bigmerge_dyads)),
        "copper_lilac_phase1_states": {str(k): int(v) for k, v in cl_states.items()},
        "copper_lilac_saved_positive_days": cl_saved_positive,
        "saved_candidate_rows": int(len(saved)),
        "outputs": {
            "dyad_day_table": str(OUT_DIR / "opportunity_dyad_day.csv"),
            "state_summary": str(OUT_DIR / "opportunity_state_summary.csv"),
            "dyad_summary": str(OUT_DIR / "opportunity_dyad_summary.csv"),
            "comparison_with_saved": str(OUT_DIR / "opportunity_vs_saved_candidates.csv"),
            "copper_lilac_validation": str(OUT_DIR / "validation_copper_lilac_days.csv"),
        },
    }
    (OUT_DIR / "phase1_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print()
    print(state_summary.to_string(index=False))
    print()
    print(crosstab.head(12).to_string(index=False))
    print()
    print(json.dumps({k: v for k, v in metadata.items() if k != "outputs"}, indent=2))


if __name__ == "__main__":
    main()
