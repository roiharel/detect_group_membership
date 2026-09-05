#!/usr/bin/env python
"""Create a log-log duration versus size plot for canonical social events."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_MEMBERSHIP = Path(
    r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_local_2h_support\canonical_hourly_membership.parquet"
)
DEFAULT_OUTPUT = Path(
    r"C:\Users\rharel\Documents\group_mebership\outputs\canonical_group_merge_scale_log_scatter"
)

EVENT_ORDER = [
    "single_animal_separation",
    "within_group_split",
    "small_subset_merge",
    "medium_partial_merge",
    "large_merge",
]
EVENT_COLORS = {
    "single_animal_separation": "#8d5a9e",
    "isolated": "#f2a7b5",
    "disperser": "#8d5a9e",
    "within_group_split": "#009e73",
    "small_subset_merge": "#9ecae1",
    "medium_partial_merge": "#4292c6",
    "large_merge": "#08519c",
}
EVENT_LABELS = {
    "single_animal_separation": "single-animal separation",
    "isolated": "isolated singleton",
    "disperser": "disperser / sustained outside-origin",
    "within_group_split": "within-group split",
    "small_subset_merge": "small subset (<50% on >=1 side)",
    "medium_partial_merge": "medium partial (50-80%)",
    "large_merge": "large/whole (>=80% both sides)",
}
PLOT_GROUPS = [
    ("single_animal_separation", ["single_animal_separation"], "single-animal separation", EVENT_COLORS["single_animal_separation"]),
    ("within_group_split", ["within_group_split"], "within-group split (2+ animals)", EVENT_COLORS["within_group_split"]),
    ("small_subset_merge", ["small_subset_merge"], "small subset merge", EVENT_COLORS["small_subset_merge"]),
    ("medium_partial_merge", ["medium_partial_merge"], "medium partial merge", EVENT_COLORS["medium_partial_merge"]),
    ("large_merge", ["large_merge"], "large/whole merge", EVENT_COLORS["large_merge"]),
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--membership-file", type=Path, default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--start-date", default="2025-01-01")
    parser.add_argument("--unit-col", default="dynamic_social_unit")
    parser.add_argument("--count-col", default="temp_group_dynamic_counts")
    parser.add_argument("--large-threshold", type=float, default=0.80)
    parser.add_argument("--small-threshold", type=float, default=0.50)
    parser.add_argument("--max-gap-hours", type=float, default=2.5)
    parser.add_argument(
        "--sample-per-type",
        type=int,
        default=0,
        help="Optional deterministic sample per plotted event type; 0 plots all events.",
    )
    return parser


def parse_counts(value: object) -> dict[str, int]:
    if not isinstance(value, str) or not value:
        return {}
    counts: dict[str, int] = {}
    for piece in value.split(";"):
        if ":" not in piece:
            continue
        group, raw_n = piece.rsplit(":", 1)
        try:
            counts[str(group)] = int(raw_n)
        except ValueError:
            continue
    return counts


def classify_scale(frac_a: float, frac_b: float, large_threshold: float, small_threshold: float) -> str:
    if frac_a >= large_threshold and frac_b >= large_threshold:
        return "large_merge"
    if frac_a < small_threshold or frac_b < small_threshold:
        return "small_subset_merge"
    return "medium_partial_merge"


def build_pair_rows(membership: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    unit_col = args.unit_col if args.unit_col in membership.columns else "origin_group"
    count_col = args.count_col if args.count_col in membership.columns else "temp_group_origin_counts"

    observed = (
        membership.groupby(["window_start", unit_col], observed=True)["animal_id"]
        .nunique()
        .rename("observed_n")
        .reset_index()
        .rename(columns={unit_col: "social_unit"})
    )
    observed_by_time = {
        timestamp: dict(zip(frame["social_unit"], frame["observed_n"]))
        for timestamp, frame in observed.groupby("window_start", observed=True)
    }

    clusters = (
        membership[membership[count_col].astype(str).str.contains(";", regex=False, na=False)]
        .drop_duplicates(["window_start", "temp_group_id"])
        .copy()
    )

    rows = []
    for record in clusters.itertuples(index=False):
        counts = parse_counts(getattr(record, count_col))
        if len(counts) < 2:
            continue
        observed_counts = observed_by_time.get(record.window_start, {})
        cluster_size = int(sum(counts.values()))
        for group_a, group_b in combinations(sorted(counts), 2):
            observed_a = observed_counts.get(group_a, np.nan)
            observed_b = observed_counts.get(group_b, np.nan)
            if not observed_a or not observed_b:
                continue
            frac_a = counts[group_a] / observed_a
            frac_b = counts[group_b] / observed_b
            rows.append(
                {
                    "window_start": record.window_start,
                    "temp_group_id": record.temp_group_id,
                    "group_a": group_a,
                    "group_b": group_b,
                    "pair": f"{group_a} - {group_b}",
                    "n_a_in_cluster": counts[group_a],
                    "n_b_in_cluster": counts[group_b],
                    "pair_event_size": counts[group_a] + counts[group_b],
                    "cluster_size": cluster_size,
                    "observed_a": observed_a,
                    "observed_b": observed_b,
                    "frac_a": frac_a,
                    "frac_b": frac_b,
                    "merge_scale": classify_scale(
                        frac_a,
                        frac_b,
                        args.large_threshold,
                        args.small_threshold,
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_group_merge_events(pair_rows: pd.DataFrame, max_gap_hours: float) -> pd.DataFrame:
    if pair_rows.empty:
        return pd.DataFrame()
    rows = pair_rows.sort_values(["pair", "merge_scale", "window_start"]).copy()
    gap_hours = (
        rows.groupby(["pair", "merge_scale"], observed=True)["window_start"]
        .diff()
        .dt.total_seconds()
        .div(3600)
    )
    starts_new = gap_hours.isna() | gap_hours.gt(max_gap_hours)
    rows["event_number"] = starts_new.groupby([rows["pair"], rows["merge_scale"]], observed=True).cumsum()
    rows["event_id"] = (
        rows["pair"].str.replace(" - ", "__", regex=False)
        + "__"
        + rows["merge_scale"]
        + "__"
        + rows["event_number"].astype(int).astype(str)
    )
    events = (
        rows.groupby("event_id", observed=True)
        .agg(
            event_type=("merge_scale", "first"),
            event_family=("merge_scale", lambda _: "group_merge"),
            label=("pair", "first"),
            pair=("pair", "first"),
            group_a=("group_a", "first"),
            group_b=("group_b", "first"),
            merge_scale=("merge_scale", "first"),
            start_time=("window_start", "min"),
            end_time=("window_start", "max"),
            observed_windows=("window_start", "nunique"),
            max_pair_event_size=("pair_event_size", "max"),
            mean_pair_event_size=("pair_event_size", "mean"),
            max_cluster_size=("cluster_size", "max"),
            max_group_a_n=("n_a_in_cluster", "max"),
            max_group_b_n=("n_b_in_cluster", "max"),
            max_group_a_fraction=("frac_a", "max"),
            max_group_b_fraction=("frac_b", "max"),
            median_group_a_fraction=("frac_a", "median"),
            median_group_b_fraction=("frac_b", "median"),
        )
        .reset_index()
    )
    events["duration_hours"] = (
        (events["end_time"] - events["start_time"]).dt.total_seconds().div(3600) + 1.0
    )
    events["duration_days"] = events["duration_hours"] / 24.0
    events["start_month"] = events["start_time"].dt.strftime("%Y-%m")
    return events


def add_contiguous_ids(
    frame: pd.DataFrame,
    key_cols: list[str],
    time_col: str,
    max_gap_hours: float,
    prefix: str,
) -> pd.Series:
    data = frame.sort_values(key_cols + [time_col]).copy()
    previous_time = data.groupby(key_cols, observed=True)[time_col].shift(1)
    gap_hours = (data[time_col] - previous_time).dt.total_seconds().div(3600)
    starts_new = previous_time.isna() | gap_hours.gt(max_gap_hours)
    event_number = starts_new.groupby([data[col] for col in key_cols], observed=True).cumsum().astype(int)
    key_text = data[key_cols].astype(str).agg("__".join, axis=1)
    return prefix + "__" + key_text + "__" + event_number.astype(str)


def summarize_individual_events(
    membership: pd.DataFrame,
    mask: pd.Series,
    event_type: str,
    max_gap_hours: float,
) -> pd.DataFrame:
    active = membership[mask.fillna(False)].copy()
    if active.empty:
        return pd.DataFrame()
    key_cols = ["animal_id"]
    if event_type == "disperser" and "dynamic_target_group" in active.columns:
        active["target_group_for_event"] = active["dynamic_target_group"].fillna(active["dynamic_social_unit"])
        key_cols = ["animal_id", "target_group_for_event"]
    active["event_id"] = add_contiguous_ids(active, key_cols, "window_start", max_gap_hours, event_type)
    label_col = "target_group_for_event" if "target_group_for_event" in active.columns else "origin_group"
    events = (
        active.groupby("event_id", observed=True)
        .agg(
            event_type=("animal_id", lambda _: event_type),
            event_family=("animal_id", lambda _: "individual"),
            label=(label_col, "first"),
            animal_id=("animal_id", "first"),
            origin_group=("origin_group", "first"),
            dynamic_social_unit=("dynamic_social_unit", "first"),
            start_time=("window_start", "min"),
            end_time=("window_start", "max"),
            observed_windows=("window_start", "nunique"),
            max_cluster_size=("temp_group_size", "max"),
            mean_cluster_size=("temp_group_size", "mean"),
            n_carried_night=("is_carried_night", "sum"),
            n_local_2h_supported=("is_local_2h_supported", "sum"),
        )
        .reset_index()
    )
    events["duration_hours"] = (
        (events["end_time"] - events["start_time"]).dt.total_seconds().div(3600) + 1.0
    )
    events["duration_days"] = events["duration_hours"] / 24.0
    events["start_month"] = events["start_time"].dt.strftime("%Y-%m")
    events["max_event_size"] = 1
    return events


def build_single_animal_separation_events(
    membership: pd.DataFrame,
    unit_col: str,
    max_gap_hours: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Union isolated, sustained disperser, and size-1 split rows without duplicates."""

    base = membership[
        membership["social_context"].eq("isolated")
        | membership["dynamic_assignment"].eq("sustained_non_origin_association")
    ].copy()
    if not base.empty:
        base["single_animal_reason"] = np.select(
            [
                base["dynamic_assignment"].eq("sustained_non_origin_association"),
                base["social_context"].eq("isolated"),
            ],
            ["sustained_non_origin", "isolated"],
            default="other_single_animal",
        )

    cluster_unit = (
        membership.groupby(["window_start", "temp_group_id", unit_col], observed=True)["animal_id"]
        .nunique()
        .rename("unit_cluster_n")
        .reset_index()
        .rename(columns={unit_col: "social_unit"})
    )
    group_hour = (
        cluster_unit.groupby(["window_start", "social_unit"], observed=True)
        .agg(
            n_temp_clusters=("temp_group_id", "nunique"),
            n_animals_observed=("unit_cluster_n", "sum"),
            largest_cluster_size=("unit_cluster_n", "max"),
        )
        .reset_index()
    )
    group_hour["split_subset_size"] = (
        group_hour["n_animals_observed"] - group_hour["largest_cluster_size"]
    ).clip(lower=0)
    single_split_hours = group_hour[
        group_hour["n_temp_clusters"].gt(1) & group_hour["split_subset_size"].eq(1)
    ][["window_start", "social_unit", "largest_cluster_size"]]
    split_components = cluster_unit.merge(single_split_hours, on=["window_start", "social_unit"], how="inner")
    split_components = split_components[
        split_components["largest_cluster_size"].eq(1)
        | split_components["unit_cluster_n"].lt(split_components["largest_cluster_size"])
    ].copy()
    split_rows = membership.merge(
        split_components[["window_start", "temp_group_id", "social_unit"]],
        left_on=["window_start", "temp_group_id", unit_col],
        right_on=["window_start", "temp_group_id", "social_unit"],
        how="inner",
    ).drop(columns=["social_unit"])
    if not split_rows.empty:
        split_rows["single_animal_reason"] = "size1_within_group_split"

    active = pd.concat([base, split_rows], ignore_index=True, sort=False)
    if active.empty:
        return pd.DataFrame(), active

    reason_rank = {
        "sustained_non_origin": 3,
        "isolated": 2,
        "size1_within_group_split": 1,
        "other_single_animal": 0,
    }
    active["reason_rank"] = active["single_animal_reason"].map(reason_rank).fillna(0)
    active = (
        active.sort_values(["animal_id", "window_start", "reason_rank"], ascending=[True, True, False])
        .drop_duplicates(["animal_id", "window_start"], keep="first")
        .copy()
    )
    active["event_id"] = add_contiguous_ids(
        active,
        ["animal_id"],
        "window_start",
        max_gap_hours,
        "single_animal_separation",
    )
    rows = (
        active.groupby("event_id", observed=True)
        .agg(
            event_type=("animal_id", lambda _: "single_animal_separation"),
            event_family=("animal_id", lambda _: "individual"),
            label=("origin_group", "first"),
            animal_id=("animal_id", "first"),
            origin_group=("origin_group", "first"),
            dynamic_social_unit=("dynamic_social_unit", "first"),
            start_time=("window_start", "min"),
            end_time=("window_start", "max"),
            observed_windows=("window_start", "nunique"),
            n_isolated_rows=("single_animal_reason", lambda x: int((x == "isolated").sum())),
            n_disperser_rows=("single_animal_reason", lambda x: int((x == "sustained_non_origin").sum())),
            n_size1_split_rows=("single_animal_reason", lambda x: int((x == "size1_within_group_split").sum())),
            max_cluster_size=("temp_group_size", "max"),
            mean_cluster_size=("temp_group_size", "mean"),
            n_carried_night=("is_carried_night", "sum"),
            n_local_2h_supported=("is_local_2h_supported", "sum"),
        )
        .reset_index()
    )
    rows["duration_hours"] = (
        (rows["end_time"] - rows["start_time"]).dt.total_seconds().div(3600) + 1.0
    )
    rows["duration_days"] = rows["duration_hours"] / 24.0
    rows["start_month"] = rows["start_time"].dt.strftime("%Y-%m")
    rows["max_event_size"] = 1
    return rows, active


def build_within_group_split_events(
    membership: pd.DataFrame,
    unit_col: str,
    max_gap_hours: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cluster_unit = (
        membership.groupby(["window_start", "temp_group_id", unit_col], observed=True)["animal_id"]
        .nunique()
        .rename("unit_cluster_n")
        .reset_index()
        .rename(columns={unit_col: "social_unit"})
    )
    group_hour = (
        cluster_unit.groupby(["window_start", "social_unit"], observed=True)
        .agg(
            n_temp_clusters=("temp_group_id", "nunique"),
            n_animals_observed=("unit_cluster_n", "sum"),
            largest_cluster_size=("unit_cluster_n", "max"),
        )
        .reset_index()
    )
    group_hour["split_subset_size"] = (
        group_hour["n_animals_observed"] - group_hour["largest_cluster_size"]
    ).clip(lower=0)
    group_hour["split_subset_fraction"] = (
        group_hour["split_subset_size"] / group_hour["n_animals_observed"].replace(0, np.nan)
    )
    active = group_hour[
        group_hour["n_temp_clusters"].gt(1) & group_hour["split_subset_size"].gt(0)
    ].copy()
    if active.empty:
        return pd.DataFrame(), group_hour
    active["event_id"] = add_contiguous_ids(
        active,
        ["social_unit"],
        "window_start",
        max_gap_hours,
        "within_group_split",
    )
    events = (
        active.groupby("event_id", observed=True)
        .agg(
            event_type=("social_unit", lambda _: "within_group_split"),
            event_family=("social_unit", lambda _: "within_group"),
            label=("social_unit", "first"),
            origin_group=("social_unit", "first"),
            start_time=("window_start", "min"),
            end_time=("window_start", "max"),
            observed_windows=("window_start", "nunique"),
            max_event_size=("split_subset_size", "max"),
            mean_event_size=("split_subset_size", "mean"),
            max_cluster_size=("n_animals_observed", "max"),
            max_n_temp_clusters=("n_temp_clusters", "max"),
            mean_n_temp_clusters=("n_temp_clusters", "mean"),
            max_split_fraction=("split_subset_fraction", "max"),
            mean_split_fraction=("split_subset_fraction", "mean"),
        )
        .reset_index()
    )
    events["duration_hours"] = (
        (events["end_time"] - events["start_time"]).dt.total_seconds().div(3600) + 1.0
    )
    events["duration_days"] = events["duration_hours"] / 24.0
    events["start_month"] = events["start_time"].dt.strftime("%Y-%m")
    return events, group_hour


def plot_log_scatter(events: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(12.5, 8.7))
    for _, event_types, _, color in PLOT_GROUPS:
        frame = events[events["event_type"].isin(event_types)]
        if event_types == ["within_group_split"]:
            frame = frame[frame["max_event_size"].ge(2)]
        if frame.empty:
            continue
        ax.scatter(
            frame["duration_hours"].clip(lower=1),
            frame["max_event_size"].clip(lower=1),
            s=24 if event_types == ["single_animal_separation"] else 32,
            alpha=0.22 if event_types == ["single_animal_separation"] else 0.34,
            edgecolor=color,
            facecolor=color,
            linewidth=0.6,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Duration (hours, log scale)")
    ax.set_ylabel("Max event size (animals, log scale)")
    ax.set_title("Canonical local-2h event size versus duration")
    ax.grid(True, which="major", color="#e2e2e2", linewidth=0.9)
    ax.grid(True, which="minor", color="#f0f0f0", linewidth=0.5)
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color=color,
            markerfacecolor=color,
            alpha=0.82,
            markersize=7,
            label=label,
        )
        for _, _, label, color in PLOT_GROUPS
    ]
    ax.legend(
        handles=handles,
        frameon=True,
        facecolor="white",
        edgecolor="#cccccc",
        loc="upper right",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(output_dir / "canonical_event_size_duration_log_scatter.png", dpi=230)
    plt.close(fig)


def build_plot_points(events: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for group_key, event_types, label, color in PLOT_GROUPS:
        frame = events[events["event_type"].isin(event_types)].copy()
        if event_types == ["within_group_split"]:
            frame = frame[frame["max_event_size"].ge(2)].copy()
        if frame.empty:
            continue
        frame["plot_group"] = group_key
        frame["plot_label"] = label
        frame["plot_color"] = color
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    plot_data = pd.concat(frames, ignore_index=True, sort=False)
    plot_data["duration_plot"] = plot_data["duration_hours"].clip(lower=1).round(6)
    plot_data["size_plot"] = plot_data["max_event_size"].clip(lower=1).round(6)
    points = (
        plot_data.groupby(["plot_group", "plot_label", "plot_color", "duration_plot", "size_plot"], observed=True)
        .agg(n_events=("event_id", "nunique"))
        .reset_index()
    )
    return points


def plot_log_scatter_weighted_jitter(events: pd.DataFrame, output_dir: Path) -> None:
    points = build_plot_points(events)
    if points.empty:
        return

    rng = np.random.default_rng(11)
    points["x_jitter"] = 10 ** (np.log10(points["duration_plot"]) + rng.normal(0, 0.018, len(points)))
    points["y_jitter"] = 10 ** (np.log10(points["size_plot"]) + rng.normal(0, 0.012, len(points)))
    points.loc[points["size_plot"].eq(1), "y_jitter"] = 1.0
    points["marker_size"] = 18 + 18 * np.sqrt(points["n_events"].clip(lower=1))
    points["marker_size"] = points["marker_size"].clip(upper=260)

    fig, ax = plt.subplots(figsize=(12.5, 8.7))
    for group_key, _, label, color in PLOT_GROUPS:
        frame = points[points["plot_group"].eq(group_key)]
        if frame.empty:
            continue
        ax.scatter(
            frame["x_jitter"],
            frame["y_jitter"],
            s=frame["marker_size"],
            alpha=0.28 if group_key == "single_animal_separation" else 0.42,
            edgecolor=color,
            facecolor=color,
            linewidth=0.7,
            label=label,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Duration (hours, log scale)")
    ax.set_ylabel("Max event size (animals, log scale)")
    ax.set_title("Canonical local-2h event size versus duration\njittered; circle area reflects repeated events at the same point")
    ax.grid(True, which="major", color="#e2e2e2", linewidth=0.9)
    ax.grid(True, which="minor", color="#f0f0f0", linewidth=0.5)

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color=color,
            markerfacecolor=color,
            alpha=0.82,
            markersize=7,
            label=label,
        )
        for _, _, label, color in PLOT_GROUPS
    ]
    legend1 = ax.legend(
        handles=handles,
        frameon=True,
        facecolor="white",
        edgecolor="#cccccc",
        loc="upper right",
        fontsize=9,
    )
    ax.add_artist(legend1)

    size_counts = np.array([1, 5, 25, 100])
    size_handles = [
        plt.scatter([], [], s=18 + 18 * np.sqrt(n), facecolor="#777777", edgecolor="#777777", alpha=0.35)
        for n in size_counts
    ]
    ax.legend(
        size_handles,
        [f"{n} event{'s' if n > 1 else ''}" for n in size_counts],
        title="Point weight",
        frameon=True,
        facecolor="white",
        edgecolor="#cccccc",
        loc="lower right",
        fontsize=8,
        title_fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(output_dir / "canonical_event_size_duration_log_scatter_weighted_jitter.png", dpi=230)
    plt.close(fig)
    points.to_csv(output_dir / "canonical_event_size_duration_weighted_plot_points.csv", index=False)


def plot_log_scatter_facets(events: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5), sharex=True, sharey=True)
    axes = axes.reshape(-1)
    for ax, event_type in zip(axes, EVENT_ORDER):
        frame = events[events["event_type"].eq(event_type)]
        title_suffix = ""
        if event_type == "within_group_split":
            frame = frame[frame["max_event_size"].ge(2)]
            title_suffix = "\nsubset size >=2"
        if not frame.empty:
            ax.scatter(
                frame["duration_hours"].clip(lower=1),
                frame["max_event_size"].clip(lower=1),
                s=22,
                alpha=0.34,
                color=EVENT_COLORS[event_type],
                edgecolor=EVENT_COLORS[event_type],
                linewidth=0.5,
            )
        ax.set_title(f"{EVENT_LABELS[event_type]}{title_suffix}\nn={len(frame):,}", fontsize=10)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="major", color="#e4e4e4")
        ax.grid(True, which="minor", color="#f1f1f1", linewidth=0.5)
    fig.supxlabel("Duration (hours, log scale)")
    fig.supylabel("Max event size (animals, log scale)")
    fig.suptitle("Canonical local-2h event size versus duration by event type")
    fig.tight_layout()
    fig.savefig(output_dir / "canonical_event_size_duration_log_scatter_facets.png", dpi=230)
    plt.close(fig)


def write_top_long_events(events: pd.DataFrame, output_dir: Path) -> None:
    top = events.sort_values("duration_hours", ascending=False).head(80)
    top.to_csv(output_dir / "canonical_event_size_duration_longest_events.csv", index=False)


def event_summary(events: pd.DataFrame) -> pd.DataFrame:
    if events.empty:
        return pd.DataFrame()
    rows = []
    for event_type, frame in events.groupby("event_type", observed=True):
        rows.append(
            {
                "event_type": event_type,
                "n_events": int(len(frame)),
                "duration_median_h": frame["duration_hours"].median(),
                "duration_p90_h": frame["duration_hours"].quantile(0.90),
                "duration_max_h": frame["duration_hours"].max(),
                "size_median": frame["max_event_size"].median(),
                "size_p90": frame["max_event_size"].quantile(0.90),
                "size_max": frame["max_event_size"].max(),
            }
        )
    return pd.DataFrame(rows).sort_values("event_type")


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    membership = pd.read_parquet(args.membership_file)
    membership["window_start"] = pd.to_datetime(membership["window_start"])
    if args.start_date:
        membership = membership[membership["window_start"].ge(pd.Timestamp(args.start_date))].copy()

    unit_col = args.unit_col if args.unit_col in membership.columns else "origin_group"
    pair_rows = build_pair_rows(membership, args)
    merge_events = build_group_merge_events(pair_rows, args.max_gap_hours)
    if not merge_events.empty:
        merge_events["max_event_size"] = merge_events["max_cluster_size"]

    isolated_mask = membership["social_context"].eq("isolated") & ~membership["dynamic_assignment"].eq(
        "sustained_non_origin_association"
    )
    disperser_mask = membership["dynamic_assignment"].eq("sustained_non_origin_association")
    isolated_events = summarize_individual_events(membership, isolated_mask, "isolated", args.max_gap_hours)
    disperser_events = summarize_individual_events(membership, disperser_mask, "disperser", args.max_gap_hours)
    single_events, single_rows = build_single_animal_separation_events(membership, unit_col, args.max_gap_hours)
    split_events, split_rows = build_within_group_split_events(membership, unit_col, args.max_gap_hours)

    events = pd.concat(
        [single_events, split_events, merge_events],
        ignore_index=True,
        sort=False,
    )
    if args.sample_per_type > 0 and not events.empty:
        events = (
            events.groupby("event_type", observed=True, group_keys=False)
            .apply(lambda x: x.sample(min(len(x), args.sample_per_type), random_state=7))
            .reset_index(drop=True)
        )

    pair_rows.to_csv(args.output_dir / "canonical_group_merge_scale_pair_rows.csv", index=False)
    merge_events.to_csv(args.output_dir / "canonical_group_merge_scale_events.csv", index=False)
    isolated_events.to_csv(args.output_dir / "canonical_isolated_events.csv", index=False)
    disperser_events.to_csv(args.output_dir / "canonical_disperser_events.csv", index=False)
    single_events.to_csv(args.output_dir / "canonical_single_animal_separation_events.csv", index=False)
    single_rows.to_csv(args.output_dir / "canonical_single_animal_separation_rows.csv", index=False)
    split_events.to_csv(args.output_dir / "canonical_within_group_split_events.csv", index=False)
    split_rows.to_csv(args.output_dir / "canonical_within_group_split_rows.csv", index=False)
    events.to_csv(args.output_dir / "canonical_event_size_duration_all_events.csv", index=False)
    summary_table = event_summary(events)
    summary_table.to_csv(args.output_dir / "canonical_event_size_duration_summary.csv", index=False)
    write_top_long_events(events, args.output_dir)
    plot_log_scatter(events, args.output_dir)
    plot_log_scatter_weighted_jitter(events, args.output_dir)
    plot_log_scatter_facets(events, args.output_dir)

    summary = {
        "membership_file": str(args.membership_file),
        "rows_used": int(len(membership)),
        "merge_pair_rows": int(len(pair_rows)),
        "events": int(len(events)),
        "start_date": args.start_date,
        "max_gap_hours_for_event_stitching": args.max_gap_hours,
        "event_counts_by_type": events["event_type"].value_counts().to_dict(),
        "event_definitions": {
            "single_animal_separation": (
                "deduplicated union of isolated rows, sustained non-origin rows, "
                "and size-1 within-group split components"
            ),
            "isolated": "social_context == isolated, excluding rows inside sustained non-origin dynamic segments",
            "disperser": "dynamic_assignment == sustained_non_origin_association; event size is one focal animal",
            "within_group_split": "a dynamic social unit appears in more than one temp_group_id at the same timestamp",
            "group_merge_scales": "mixed temp-group pairs split into small/medium/large by observed representation fractions",
        },
        "output_dir": str(args.output_dir),
    }
    with (args.output_dir / "canonical_group_merge_scale_log_scatter_metadata.json").open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
