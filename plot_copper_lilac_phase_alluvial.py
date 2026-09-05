#!/usr/bin/env python
"""Copper-Lilac phase plot with collared-individual alluvial summaries."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle


DEFAULT_MEMBERSHIP = Path(
    r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_local_2h_support\canonical_hourly_membership.parquet"
)
DEFAULT_PAIR_SUMMARY = Path(
    r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_local_2h_support\effect_plots\canonical_group_pair_merge_probabilities_monthly.csv"
)
DEFAULT_OUTPUT = Path(r"C:\Users\rharel\Documents\group_mebership\outputs\copper_lilac_phase_alluvial")

GROUPS = ("Copper", "Lilac")
CATEGORY_ORDER = ["Copper", "Copper-Lilac mixed", "Lilac", "single-animal separation", "other"]
COLORS = {
    "Copper": "#b87333",
    "Lilac": "#8e63a9",
    "Copper-Lilac mixed": "#67b7b0",
    "single-animal separation": "#b279a2",
    "other": "#bdbdbd",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--membership-file", type=Path, default=DEFAULT_MEMBERSHIP)
    parser.add_argument("--pair-summary", type=Path, default=DEFAULT_PAIR_SUMMARY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--start-date", default="2025-01-01")
    parser.add_argument("--min-phase-months", type=int, default=2)
    parser.add_argument("--n-phases", type=int, default=4)
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
            counts[group] = int(raw_n)
        except ValueError:
            continue
    return counts


def classify_row(row: pd.Series) -> str:
    counts = parse_counts(row.get("temp_group_dynamic_counts", ""))
    has_copper = counts.get("Copper", 0) > 0
    has_lilac = counts.get("Lilac", 0) > 0
    if bool(row.get("is_single_animal_separation", False)):
        return "single-animal separation"
    if has_copper and has_lilac:
        return "Copper-Lilac mixed"
    unit = str(row.get("dynamic_social_unit", row.get("origin_group", "")))
    origin = str(row.get("origin_group", ""))
    if unit in GROUPS:
        return unit
    if origin in GROUPS:
        return origin
    return "other"


def phase_breaks_dp(values: np.ndarray, n_phases: int, min_len: int) -> list[tuple[int, int]]:
    n = len(values)
    n_phases = min(n_phases, max(1, n // min_len))
    sse = np.full((n, n), np.inf)
    for i in range(n):
        for j in range(i + min_len - 1, n):
            segment = values[i : j + 1]
            sse[i, j] = float(((segment - segment.mean()) ** 2).sum())
    dp = np.full((n_phases + 1, n), np.inf)
    prev = np.full((n_phases + 1, n), -1, dtype=int)
    for j in range(min_len - 1, n):
        dp[1, j] = sse[0, j]
    for k in range(2, n_phases + 1):
        for j in range(k * min_len - 1, n):
            for cut in range((k - 1) * min_len - 1, j - min_len + 1):
                cost = dp[k - 1, cut] + sse[cut + 1, j]
                if cost < dp[k, j]:
                    dp[k, j] = cost
                    prev[k, j] = cut
    segments: list[tuple[int, int]] = []
    k = n_phases
    end = n - 1
    while k > 0:
        cut = prev[k, end]
        start = 0 if k == 1 else cut + 1
        segments.append((start, end))
        end = cut
        k -= 1
    return list(reversed(segments))


def compute_copper_lilac_monthly_merge(membership: pd.DataFrame) -> pd.DataFrame:
    anchor = membership[membership["window_start"].dt.hour.isin([11, 16])].copy()
    observed = (
        anchor[anchor["dynamic_social_unit"].isin(GROUPS)]
        .groupby(["window_start", "dynamic_social_unit"], observed=True)["animal_id"]
        .nunique()
        .rename("observed_n")
        .reset_index()
    )
    observed_by_time = {
        timestamp: dict(zip(frame["dynamic_social_unit"], frame["observed_n"]))
        for timestamp, frame in observed.groupby("window_start", observed=True)
    }
    opportunities = []
    large_windows = set()
    mixed_clusters = (
        anchor[anchor["temp_group_dynamic_counts"].astype(str).str.contains(";", regex=False, na=False)]
        .drop_duplicates(["window_start", "temp_group_id"])
        .copy()
    )
    for timestamp, observed_counts in observed_by_time.items():
        if all(observed_counts.get(group, 0) > 0 for group in GROUPS):
            opportunities.append(timestamp)
    for record in mixed_clusters.itertuples(index=False):
        counts = parse_counts(record.temp_group_dynamic_counts)
        if not all(counts.get(group, 0) > 0 for group in GROUPS):
            continue
        observed_counts = observed_by_time.get(record.window_start, {})
        if not all(observed_counts.get(group, 0) > 0 for group in GROUPS):
            continue
        fractions = [counts[group] / observed_counts[group] for group in GROUPS]
        if min(fractions) >= 0.80:
            large_windows.add(record.window_start)
    opportunity = pd.DataFrame({"window_start": sorted(set(opportunities))})
    if opportunity.empty:
        return pd.DataFrame(
            columns=["month_dt", "large_merge_hours", "opportunity_hours", "large_merge_probability"]
        )
    opportunity["month_dt"] = opportunity["window_start"].dt.to_period("M").dt.to_timestamp()
    monthly = (
        opportunity.groupby("month_dt", observed=True)["window_start"]
        .nunique()
        .rename("opportunity_hours")
        .reset_index()
    )
    large = pd.DataFrame({"window_start": sorted(large_windows)})
    if large.empty:
        monthly["large_merge_hours"] = 0
    else:
        large["month_dt"] = large["window_start"].dt.to_period("M").dt.to_timestamp()
        large_month = (
            large.groupby("month_dt", observed=True)["window_start"]
            .nunique()
            .rename("large_merge_hours")
            .reset_index()
        )
        monthly = monthly.merge(large_month, on="month_dt", how="left")
        monthly["large_merge_hours"] = monthly["large_merge_hours"].fillna(0)
    monthly["large_merge_probability"] = (
        monthly["large_merge_hours"] / monthly["opportunity_hours"].replace(0, np.nan)
    )
    return monthly.sort_values("month_dt").reset_index(drop=True)


def phase_label(mean_probability: float, phase_i: int) -> str:
    if mean_probability < 0.08:
        return "low contact"
    if mean_probability < 0.35:
        return "early mixing"
    if mean_probability < 0.75:
        return "transition"
    return "high merge"


def load_phase_table(args: argparse.Namespace, membership: pd.DataFrame) -> pd.DataFrame:
    monthly = compute_copper_lilac_monthly_merge(membership)
    monthly = (
        monthly[monthly["month_dt"].ge(pd.Timestamp(args.start_date))]
        .sort_values("month_dt")
        .reset_index(drop=True)
    )
    if monthly.empty:
        raise ValueError("No Copper-Lilac opportunities were found for the requested start date.")
    segments = phase_breaks_dp(
        monthly["large_merge_probability"].fillna(0).to_numpy(),
        args.n_phases,
        args.min_phase_months,
    )
    phase_rows = []
    for phase_i, (start_i, end_i) in enumerate(segments):
        phase_months = monthly.iloc[start_i : end_i + 1].copy()
        label = phase_label(float(phase_months["large_merge_probability"].mean()), phase_i)
        if label in [row["phase"] for row in phase_rows]:
            label = f"{label} {phase_i + 1}"
        phase_rows.append(
            {
                "phase": label,
                "phase_index": phase_i,
                "start_month": phase_months["month_dt"].min(),
                "end_month": phase_months["month_dt"].max(),
                "mean_large_merge_probability": phase_months["large_merge_probability"].mean(),
                "n_months": len(phase_months),
            }
        )
        monthly.loc[phase_months.index, "phase"] = label
        monthly.loc[phase_months.index, "phase_index"] = phase_i
    phases = pd.DataFrame(phase_rows)
    return monthly.merge(phases[["phase", "phase_index"]], on=["phase", "phase_index"], how="left"), phases


def load_membership(args: argparse.Namespace) -> pd.DataFrame:
    cols = [
        "window_start",
        "animal_id",
        "sex",
        "age",
        "origin_group",
        "dynamic_social_unit",
        "social_context",
        "dynamic_assignment",
        "temp_group_id",
        "temp_group_dynamic_counts",
    ]
    data = pd.read_parquet(args.membership_file, columns=cols)
    data["window_start"] = pd.to_datetime(data["window_start"])
    data = data[data["window_start"].ge(pd.Timestamp(args.start_date))].copy()
    data = data[data["origin_group"].isin(GROUPS) | data["dynamic_social_unit"].isin(GROUPS)].copy()
    data["month_dt"] = data["window_start"].dt.to_period("M").dt.to_timestamp()
    data["is_single_animal_separation"] = (
        data["social_context"].eq("isolated")
        | data["dynamic_assignment"].eq("sustained_non_origin_association")
    )
    return data


def add_single_split_flags(data: pd.DataFrame) -> pd.DataFrame:
    cluster_unit = (
        data.groupby(["window_start", "temp_group_id", "dynamic_social_unit"], observed=True)["animal_id"]
        .nunique()
        .rename("unit_cluster_n")
        .reset_index()
        .rename(columns={"dynamic_social_unit": "social_unit"})
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
    single_hours = group_hour[
        group_hour["n_temp_clusters"].gt(1) & group_hour["split_subset_size"].eq(1)
    ][["window_start", "social_unit", "largest_cluster_size"]]
    components = cluster_unit.merge(single_hours, on=["window_start", "social_unit"], how="inner")
    components = components[
        components["largest_cluster_size"].eq(1) | components["unit_cluster_n"].lt(components["largest_cluster_size"])
    ]
    keys = components[["window_start", "temp_group_id", "social_unit"]].drop_duplicates()
    flagged = data.merge(
        keys,
        left_on=["window_start", "temp_group_id", "dynamic_social_unit"],
        right_on=["window_start", "temp_group_id", "social_unit"],
        how="left",
    )
    flagged["is_single_animal_separation"] = flagged["is_single_animal_separation"] | flagged["social_unit"].notna()
    return flagged.drop(columns=["social_unit"])


def summarize_monthly_categories(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    data["category"] = data.apply(classify_row, axis=1)
    animal_month = (
        data.groupby(["month_dt", "animal_id", "sex", "age", "origin_group", "category"], observed=True)
        .size()
        .rename("n_rows")
        .reset_index()
    )
    idx = animal_month.groupby(["month_dt", "animal_id"], observed=True)["n_rows"].idxmax()
    dominant = animal_month.loc[idx].copy()
    counts = (
        dominant.groupby(["month_dt", "category"], observed=True)["animal_id"]
        .nunique()
        .rename("n_animals")
        .reset_index()
    )
    return dominant, counts


def summarize_phase_categories(dominant_month: pd.DataFrame, phases_by_month: pd.DataFrame) -> pd.DataFrame:
    phase_lookup = phases_by_month[["month_dt", "phase", "phase_index"]].drop_duplicates()
    frame = dominant_month.merge(phase_lookup, on="month_dt", how="inner")
    phase_cat = (
        frame.groupby(["phase", "phase_index", "animal_id", "sex", "age", "origin_group", "category"], observed=True)
        .size()
        .rename("n_months")
        .reset_index()
    )
    idx = phase_cat.groupby(["phase", "animal_id"], observed=True)["n_months"].idxmax()
    return phase_cat.loc[idx].copy()


def draw_flow(ax, x0: float, x1: float, y0a: float, y0b: float, y1a: float, y1b: float, color: str, alpha: float) -> None:
    dx = x1 - x0
    verts = [
        (x0, y0a),
        (x0 + dx * 0.45, y0a),
        (x1 - dx * 0.45, y1a),
        (x1, y1a),
        (x1, y1b),
        (x1 - dx * 0.45, y1b),
        (x0 + dx * 0.45, y0b),
        (x0, y0b),
        (x0, y0a),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1))


def plot_alluvial_panel(ax, phase_assignments: pd.DataFrame, phases: pd.DataFrame, sex: str) -> None:
    frame = phase_assignments[phase_assignments["sex"].astype(str).str.lower().str.startswith(sex)].copy()
    if frame.empty:
        ax.set_title(f"{sex} - no data")
        return
    phases = phases.sort_values("phase_index").reset_index(drop=True)
    bar_w = 0.22
    x_positions = dict(zip(phases["phase"], np.arange(len(phases))))
    segment_ranges: dict[tuple[str, str], tuple[float, float]] = {}

    for phase in phases["phase"]:
        counts = frame[frame["phase"].eq(phase)]["category"].value_counts().to_dict()
        y = 0.0
        for category in CATEGORY_ORDER:
            n = int(counts.get(category, 0))
            if n <= 0:
                continue
            segment_ranges[(phase, category)] = (y, y + n)
            x = x_positions[phase]
            ax.add_patch(
                Rectangle(
                    (x - bar_w / 2, y),
                    bar_w,
                    n,
                    facecolor=COLORS[category],
                    edgecolor="black",
                    linewidth=0.8,
                    zorder=3,
                )
            )
            if n >= 2:
                ax.text(x, y + n / 2, category.replace("Copper-Lilac mixed", "mixed"), ha="center", va="center", fontsize=7)
            y += n

    for phase_a, phase_b in zip(phases["phase"].iloc[:-1], phases["phase"].iloc[1:]):
        left = frame[frame["phase"].eq(phase_a)][["animal_id", "category"]].rename(columns={"category": "from_category"})
        right = frame[frame["phase"].eq(phase_b)][["animal_id", "category"]].rename(columns={"category": "to_category"})
        transitions = left.merge(right, on="animal_id", how="inner")
        if transitions.empty:
            continue
        trans_counts = (
            transitions.groupby(["from_category", "to_category"], observed=True)["animal_id"]
            .nunique()
            .rename("n")
            .reset_index()
        )
        source_offsets = {cat: segment_ranges.get((phase_a, cat), (0, 0))[0] for cat in CATEGORY_ORDER}
        target_offsets = {cat: segment_ranges.get((phase_b, cat), (0, 0))[0] for cat in CATEGORY_ORDER}
        for row in trans_counts.itertuples(index=False):
            n = int(row.n)
            if n <= 0:
                continue
            y0a = source_offsets[row.from_category]
            y0b = y0a + n
            source_offsets[row.from_category] = y0b
            y1a = target_offsets[row.to_category]
            y1b = y1a + n
            target_offsets[row.to_category] = y1b
            draw_flow(
                ax,
                x_positions[phase_a] + bar_w / 2,
                x_positions[phase_b] - bar_w / 2,
                y0a,
                y0b,
                y1a,
                y1b,
                COLORS[row.from_category],
                0.24,
            )

    ax.set_xticks(np.arange(len(phases)))
    ax.set_xticklabels([str(p).replace(" ", "\n") for p in phases["phase"]], rotation=0, ha="center")
    ax.set_ylabel("# collared adults/subadults")
    ax.set_title("females" if sex == "f" else "males", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-0.45, len(phases) - 0.55)
    ymax = max((rng[1] for rng in segment_ranges.values()), default=1)
    ax.set_ylim(0, ymax * 1.16)


def plot_figure(monthly_counts: pd.DataFrame, phases: pd.DataFrame, phase_assignments: pd.DataFrame, output_dir: Path) -> None:
    fig = plt.figure(figsize=(15, 9.8))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.35], hspace=0.62, wspace=0.30)
    ax_top = fig.add_subplot(grid[0, :])
    ax_f = fig.add_subplot(grid[1, 0])
    ax_m = fig.add_subplot(grid[1, 1])

    pivot = (
        monthly_counts.pivot_table(index="month_dt", columns="category", values="n_animals", fill_value=0)
        .reindex(columns=CATEGORY_ORDER, fill_value=0)
        .sort_index()
    )
    bottom = np.zeros(len(pivot))
    width_days = 23
    for category in CATEGORY_ORDER:
        values = pivot[category].to_numpy()
        ax_top.bar(
            pivot.index,
            values,
            bottom=bottom,
            width=width_days,
            color=COLORS[category],
            edgecolor=COLORS[category],
            alpha=0.72,
            label=category,
        )
        bottom += values

    for row in phases.itertuples(index=False):
        start = pd.Timestamp(row.start_month)
        end = pd.Timestamp(row.end_month) + pd.offsets.MonthEnd(1)
        ax_top.axvspan(start, end, color="#f2f2f2" if row.phase_index % 2 == 0 else "#e6e6e6", alpha=0.55, zorder=0)
        ax_top.text(
            start + (end - start) / 2,
            0.96,
            row.phase,
            ha="center",
            va="top",
            fontsize=10,
            transform=ax_top.get_xaxis_transform(),
        )

    ax_top.set_ylabel("# collared adults/subadults")
    ax_top.set_title("Copper-Lilac collared social states through data-derived merge phases", pad=18)
    n_months = len(pivot.index)
    ax_top.xaxis.set_major_locator(mdates.MonthLocator(interval=3 if n_months > 18 else 2))
    ax_top.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax_top.tick_params(axis="x", rotation=35, pad=4)
    if len(pivot.index):
        ax_top.set_xlim(pivot.index.min() - pd.Timedelta(days=25), pivot.index.max() + pd.Timedelta(days=45))
    ax_top.spines[["top", "right"]].set_visible(False)
    ax_top.legend(frameon=False, ncol=5, loc="upper left", bbox_to_anchor=(0, 1.28), fontsize=8)

    plot_alluvial_panel(ax_f, phase_assignments, phases, "f")
    plot_alluvial_panel(ax_m, phase_assignments, phases, "m")
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.10, top=0.90)
    fig.savefig(output_dir / "copper_lilac_phase_alluvial_collared.png", dpi=230)
    plt.close(fig)


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    membership = add_single_split_flags(load_membership(args))
    phase_months, phases = load_phase_table(args, membership)
    dominant_month, monthly_counts = summarize_monthly_categories(membership)
    phase_assignments = summarize_phase_categories(dominant_month, phase_months)

    phase_months.to_csv(args.output_dir / "copper_lilac_data_driven_phase_months.csv", index=False)
    phases.to_csv(args.output_dir / "copper_lilac_data_driven_phases.csv", index=False)
    dominant_month.to_csv(args.output_dir / "copper_lilac_animal_month_dominant_category.csv", index=False)
    phase_assignments.to_csv(args.output_dir / "copper_lilac_animal_phase_dominant_category.csv", index=False)
    monthly_counts.to_csv(args.output_dir / "copper_lilac_monthly_category_counts.csv", index=False)
    plot_figure(monthly_counts, phases, phase_assignments, args.output_dir)

    metadata = {
        "membership_file": str(args.membership_file),
        "pair_summary": str(args.pair_summary),
        "output_dir": str(args.output_dir),
        "groups": list(GROUPS),
        "category_order": CATEGORY_ORDER,
        "phase_definition": "contiguous dynamic-programming segments minimizing within-phase variance of monthly Copper-Lilac large/whole merge probability",
        "n_animals": int(phase_assignments["animal_id"].nunique()),
    }
    with (args.output_dir / "copper_lilac_phase_alluvial_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
