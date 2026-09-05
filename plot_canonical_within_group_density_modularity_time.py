from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_canonical_within_group_density_modularity_ridges import (
    ANCHOR_HOURS,
    MIN_WEEKLY_ANIMALS,
    MIN_WEEKLY_TIMESTAMPS,
    SCOPE_DEFS,
    SIZE_CLASS_COLORS,
    add_group_sizes,
    add_size_adjusted_metrics,
    compute_weekly_metrics,
    supplied_group_size_table,
)


BASE = Path(r"C:\Users\rharel\Documents")
CANONICAL = (
    BASE
    / "New project"
    / "outputs"
    / "canonical_robust_hourly_membership_local_2h_support_conservative_dynamic"
    / "canonical_hourly_membership.parquet"
)
OUT_DIR = (
    BASE
    / "group_mebership"
    / "outputs"
    / "canonical_within_group_density_modularity_time_conservative_14d"
)

START_DATE = "2024-03-01"
MIN_GROUP_WEEKS_FOR_TIME = 8
PLOT_GROUP_COLS = 4


def load_weekly_metrics() -> pd.DataFrame:
    columns = [
        "animal_id",
        "window_start",
        "dynamic_social_unit",
        "temp_group_id",
        "position_support_confidence",
        "membership_confidence",
        "is_carried_night",
        "is_local_2h_supported",
    ]
    df = pd.read_parquet(CANONICAL, columns=columns)
    df["window_start"] = pd.to_datetime(df["window_start"])
    df = df[df["window_start"].ge(pd.Timestamp(START_DATE))]
    df = df[df["window_start"].dt.hour.isin(ANCHOR_HOURS)]

    metrics = pd.concat(
        [
            compute_weekly_metrics(df, scope_name=scope_name, hour=hour)
            for scope_name, hour in SCOPE_DEFS.items()
        ],
        ignore_index=True,
    )
    metrics = add_group_sizes(metrics)
    metrics["valid_for_time"] = (
        metrics["n_animals_observed"].ge(MIN_WEEKLY_ANIMALS)
        & metrics["n_timestamps"].ge(MIN_WEEKLY_TIMESTAMPS)
        & metrics["association_density"].notna()
        & metrics["modularity"].notna()
        & metrics["group_size_class"].ne("unknown")
    )
    metrics = add_size_adjusted_metrics(metrics)
    return metrics


def valid_group_order(metrics: pd.DataFrame, scope: str) -> list[str]:
    use = metrics[metrics["scope"].eq(scope) & metrics["valid_for_time"]].copy()
    summary = (
        use.groupby("dynamic_social_unit", dropna=False)
        .agg(
            n_weeks=("period_start", "nunique"),
            nominal_size=("origin_group_total_size", "max"),
            median_tracked=("n_animals_observed", "median"),
        )
        .query("n_weeks >= @MIN_GROUP_WEEKS_FOR_TIME")
        .reset_index()
        .sort_values(["nominal_size", "dynamic_social_unit"], na_position="last")
    )
    return summary["dynamic_social_unit"].tolist()


def style_time_axis(ax: plt.Axes) -> None:
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(True, color="#e8e8e8", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_scope_facets(metrics: pd.DataFrame, scope: str, adjusted: bool, out_path: Path) -> None:
    groups = valid_group_order(metrics, scope)
    if not groups:
        return

    density_col = "association_density_size_adjusted" if adjusted else "association_density"
    modularity_col = "modularity_size_adjusted" if adjusted else "modularity"
    subtitle = "size adjusted to common tracked-animal count" if adjusted else "raw weekly values"

    ncols = PLOT_GROUP_COLS
    nrows = int(np.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.1 * ncols, 2.55 * nrows),
        sharex=True,
        sharey=False,
        squeeze=False,
    )
    for ax in axes.flat:
        ax.set_visible(False)

    for ax, group in zip(axes.flat, groups):
        ax.set_visible(True)
        gd = (
            metrics[
                metrics["scope"].eq(scope)
                & metrics["dynamic_social_unit"].eq(group)
                & metrics["valid_for_time"]
            ]
            .sort_values("period_start")
            .copy()
        )
        if gd.empty:
            continue
        size = gd["origin_group_total_size"].dropna()
        size_text = f"{int(size.iloc[0])}" if len(size) else "?"
        size_class = gd["group_size_class"].dropna().iloc[0] if gd["group_size_class"].notna().any() else "unknown"
        title_color = SIZE_CLASS_COLORS.get(size_class, "#7F7F7F")

        ax.plot(
            gd["period_start"],
            gd[density_col],
            color="#1b9e77",
            linewidth=1.6,
            marker="o",
            markersize=2.4,
            alpha=0.9,
            label="density",
        )
        ax.plot(
            gd["period_start"],
            gd[modularity_col],
            color="#6a3d9a",
            linewidth=1.6,
            marker="o",
            markersize=2.4,
            alpha=0.9,
            label="modularity",
        )
        ax.scatter(
            gd["period_start"],
            gd[density_col],
            s=np.clip(gd["n_animals_observed"].astype(float) * 4, 12, 90),
            color="#1b9e77",
            alpha=0.16,
            linewidth=0,
        )
        ax.scatter(
            gd["period_start"],
            gd[modularity_col],
            s=np.clip(gd["n_animals_observed"].astype(float) * 4, 12, 90),
            color="#6a3d9a",
            alpha=0.16,
            linewidth=0,
        )
        if len(gd) >= 4:
            smooth = gd.set_index("period_start")[[density_col, modularity_col]].rolling(
                4, min_periods=2
            ).mean()
            ax.plot(smooth.index, smooth[density_col], color="#0f6b51", linewidth=2.4, alpha=0.7)
            ax.plot(smooth.index, smooth[modularity_col], color="#3f1f63", linewidth=2.4, alpha=0.7)

        ax.set_title(f"{group} ({size_text})", fontsize=10, color=title_color, pad=5)
        style_time_axis(ax)
        if not adjusted:
            ax.set_ylim(-0.03, 1.03)

    handles = [
        plt.Line2D([0], [0], color="#1b9e77", lw=2.2, label="density"),
        plt.Line2D([0], [0], color="#6a3d9a", lw=2.2, label="modularity"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=10)
    fig.suptitle(
        f"Within-group density and modularity through time: {scope.replace('_', ' ')}\n"
        f"{subtitle}; point size scales with tracked animals",
        fontsize=15,
        y=0.995,
    )
    fig.supxlabel("Week", fontsize=11, y=0.035)
    fig.supylabel("Metric value", fontsize=11, x=0.006)
    fig.subplots_adjust(left=0.045, right=0.995, top=0.925, bottom=0.085, hspace=0.45, wspace=0.2)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_day_night_delta(metrics: pd.DataFrame, out_path: Path) -> None:
    use = metrics[metrics["scope"].isin(["day_1100", "night_1600"]) & metrics["valid_for_time"]].copy()
    day = use[use["scope"].eq("day_1100")][
        ["period_start", "dynamic_social_unit", "association_density", "modularity"]
    ].rename(
        columns={
            "association_density": "day_density",
            "modularity": "day_modularity",
        }
    )
    night = use[use["scope"].eq("night_1600")][
        ["period_start", "dynamic_social_unit", "association_density", "modularity"]
    ].rename(
        columns={
            "association_density": "night_density",
            "modularity": "night_modularity",
        }
    )
    delta = day.merge(night, on=["period_start", "dynamic_social_unit"], how="inner")
    delta = add_group_sizes(delta)
    delta["density_delta_night_minus_day"] = delta["night_density"] - delta["day_density"]
    delta["modularity_delta_night_minus_day"] = delta["night_modularity"] - delta["day_modularity"]
    delta.to_csv(OUT_DIR / "canonical_within_group_weekly_day_night_deltas.csv", index=False)

    groups = (
        delta.groupby("dynamic_social_unit")
        .agg(
            n_weeks=("period_start", "nunique"),
            nominal_size=("origin_group_total_size", "max"),
        )
        .query("n_weeks >= @MIN_GROUP_WEEKS_FOR_TIME")
        .reset_index()
        .sort_values(["nominal_size", "dynamic_social_unit"], na_position="last")[
            "dynamic_social_unit"
        ]
        .tolist()
    )
    if not groups:
        return

    ncols = PLOT_GROUP_COLS
    nrows = int(np.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.1 * ncols, 2.35 * nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for ax in axes.flat:
        ax.set_visible(False)

    for ax, group in zip(axes.flat, groups):
        ax.set_visible(True)
        gd = delta[delta["dynamic_social_unit"].eq(group)].sort_values("period_start")
        ax.axhline(0, color="#303030", linewidth=0.9)
        ax.plot(
            gd["period_start"],
            gd["density_delta_night_minus_day"],
            color="#1b9e77",
            marker="o",
            markersize=2.4,
            linewidth=1.5,
            label="density",
        )
        ax.plot(
            gd["period_start"],
            gd["modularity_delta_night_minus_day"],
            color="#6a3d9a",
            marker="o",
            markersize=2.4,
            linewidth=1.5,
            label="modularity",
        )
        size = gd["origin_group_total_size"].dropna()
        size_text = f"{int(size.iloc[0])}" if len(size) else "?"
        ax.set_title(f"{group} ({size_text})", fontsize=10, pad=5)
        style_time_axis(ax)
    fig.legend(
        [
            plt.Line2D([0], [0], color="#1b9e77", lw=2.2),
            plt.Line2D([0], [0], color="#6a3d9a", lw=2.2),
        ],
        ["density", "modularity"],
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=10,
    )
    fig.suptitle(
        "Night minus day within-group structure through time\npositive values mean higher at 16:00 than 11:00",
        fontsize=15,
        y=0.995,
    )
    fig.supxlabel("Week", fontsize=11, y=0.035)
    fig.supylabel("Night - day", fontsize=11, x=0.006)
    fig.subplots_adjust(left=0.045, right=0.995, top=0.92, bottom=0.09, hspace=0.45, wspace=0.2)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def write_time_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    valid = metrics[metrics["valid_for_time"]].copy()
    summary = (
        valid.groupby(["scope", "dynamic_social_unit", "group_size_class"], dropna=False)
        .agg(
            nominal_group_size=("origin_group_total_size", "max"),
            n_weeks=("period_start", "nunique"),
            first_week=("period_start", "min"),
            last_week=("period_start", "max"),
            median_tracked_animals=("n_animals_observed", "median"),
            median_density=("association_density", "median"),
            density_slope_per_year=(
                "association_density",
                lambda s: simple_slope_per_year(valid.loc[s.index, "period_start"], s),
            ),
            median_modularity=("modularity", "median"),
            modularity_slope_per_year=(
                "modularity",
                lambda s: simple_slope_per_year(valid.loc[s.index, "period_start"], s),
            ),
        )
        .reset_index()
        .sort_values(["scope", "nominal_group_size", "dynamic_social_unit"], na_position="last")
    )
    summary.to_csv(OUT_DIR / "canonical_within_group_density_modularity_time_summary.csv", index=False)
    return summary


def simple_slope_per_year(x_dates: pd.Series, y_values: pd.Series) -> float:
    y = y_values.astype(float).to_numpy()
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan
    x = pd.to_datetime(x_dates).map(pd.Timestamp.toordinal).astype(float).to_numpy()
    x = (x - np.nanmin(x)) / 365.25
    slope, _ = np.polyfit(x[ok], y[ok], deg=1)
    return float(slope)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    supplied_group_size_table().to_csv(
        OUT_DIR / "supplied_nominal_total_group_sizes_used.csv", index=False
    )

    metrics = load_weekly_metrics()
    metrics.to_csv(OUT_DIR / "canonical_within_group_weekly_network_metrics_time.csv", index=False)
    summary = write_time_summary(metrics)

    for scope in SCOPE_DEFS:
        plot_scope_facets(
            metrics,
            scope=scope,
            adjusted=False,
            out_path=OUT_DIR / f"{scope}_density_modularity_time_facets_raw.png",
        )
        plot_scope_facets(
            metrics,
            scope=scope,
            adjusted=True,
            out_path=OUT_DIR / f"{scope}_density_modularity_time_facets_size_adjusted.png",
        )

    plot_day_night_delta(
        metrics,
        OUT_DIR / "day_night_density_modularity_delta_time_facets.png",
    )

    print(f"Wrote weekly metrics: {OUT_DIR / 'canonical_within_group_weekly_network_metrics_time.csv'}")
    print(f"Wrote time summary: {OUT_DIR / 'canonical_within_group_density_modularity_time_summary.csv'}")
    print(summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
