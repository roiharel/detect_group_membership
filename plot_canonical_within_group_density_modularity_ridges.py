from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


BASE = Path(r"C:\Users\rharel\Documents")
CANONICAL = (
    BASE
    / "New project"
    / "outputs"
    / "canonical_robust_hourly_membership_local_2h_support"
    / "canonical_hourly_membership.parquet"
)
OUT_DIR = (
    BASE
    / "group_mebership"
    / "outputs"
    / "canonical_within_group_density_modularity_ridges"
)

START_DATE = "2024-03-01"
ANCHOR_HOURS = (11, 16)
MIN_WEEKLY_ANIMALS = 4
MIN_WEEKLY_TIMESTAMPS = 4
MIN_GROUP_WEEKS_FOR_RIDGE = 6

# User-supplied group sizes, lower/higher. The higher value is used as the
# nominal total group size for low/medium/high classes.
SUPPLIED_GROUP_SIZES = [
    ("Chartreuse", 24, 74),
    ("Emerald", 13, 31),
    ("Copper", 13, 50),
    ("Lilac", 13, 60),
    ("Magenta", 14, 56),
    ("Bronze", 11, 48),
    ("Lapis Splinter", 4, 23),
    ("Periwinkle", 8, 48),
    ("Lapis", 29, 74),
    ("Purple", 15, 62),
    ("Jade", 12, 108),
    ("Maroon", 9, 42),
    ("Sapphire", 5, 42),
    ("Tricky Teal", 1, 18),
    ("Sneaky Silver", 7, 31),
    ("Phantom West", 21, 78),
    ("Green", 1, 33),
    ("Ruby Runners", 14, 64),
    ("Fence Jumpers", 9, 61),
    ("Obsidian", 13, 75),
    ("Fire Opal", 13, 62),
    ("Turquoise", 13, 60),
    ("Pearl", 8, 44),
    ("Red", 3, 50),
    ("Gold", 9, 50),
    ("Mint", 9, 80),
    ("Beverly Blue", 9, 60),
]

SIZE_CLASS_COLORS = {
    "low (<40)": "#C7A76C",
    "medium (41-59)": "#59A89C",
    "high (>60)": "#8E63A4",
    "boundary (=40 or 60)": "#B8B8B8",
    "unknown": "#7F7F7F",
}

SCOPE_DEFS = {
    "combined_1100_1600": None,
    "day_1100": 11,
    "night_1600": 16,
}


def week_start(series: pd.Series) -> pd.Series:
    return series.dt.to_period("W-SUN").dt.start_time


def pair_key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def safe_greedy_modularity(graph: nx.Graph) -> tuple[float, int, float]:
    if graph.number_of_nodes() < 2 or graph.number_of_edges() < 1:
        return 0.0, graph.number_of_nodes(), 1.0 if graph.number_of_nodes() else np.nan
    try:
        communities = list(nx.algorithms.community.greedy_modularity_communities(graph, weight="weight"))
        modularity = nx.algorithms.community.modularity(graph, communities, weight="weight")
    except Exception:
        return np.nan, np.nan, np.nan
    largest = max((len(c) for c in communities), default=0)
    largest_fraction = largest / graph.number_of_nodes() if graph.number_of_nodes() else np.nan
    return float(modularity), len(communities), float(largest_fraction)


def network_metrics_one_group(g: pd.DataFrame) -> dict[str, float | int]:
    animals = sorted(g["animal_id"].dropna().unique())
    n_animals = len(animals)
    n_possible_dyads = comb(n_animals, 2) if n_animals >= 2 else 0

    co_observed: Counter[tuple[str, str]] = Counter()
    co_clustered: Counter[tuple[str, str]] = Counter()
    split_timestamps = 0
    multi_animal_split_timestamps = 0
    largest_cluster_fractions: list[float] = []

    for _, gt in g.groupby("window_start", sort=False):
        observed = sorted(gt["animal_id"].dropna().unique())
        if len(observed) >= 2:
            for a, b in combinations(observed, 2):
                co_observed[pair_key(a, b)] += 1

        cluster_sizes = (
            gt.groupby("temp_group_id")["animal_id"]
            .nunique()
            .sort_values(ascending=False)
            .astype(int)
            .tolist()
        )
        if cluster_sizes:
            largest_cluster_fractions.append(cluster_sizes[0] / max(len(observed), 1))
        if len(cluster_sizes) > 1:
            split_timestamps += 1
            if sum(size >= 2 for size in cluster_sizes) >= 2:
                multi_animal_split_timestamps += 1

        for _, gc in gt.groupby("temp_group_id", sort=False):
            members = sorted(gc["animal_id"].dropna().unique())
            if len(members) < 2:
                continue
            for a, b in combinations(members, 2):
                co_clustered[pair_key(a, b)] += 1

    associations = {
        pair: co_clustered.get(pair, 0) / observed_count
        for pair, observed_count in co_observed.items()
        if observed_count > 0
    }
    positive_associations = {pair: weight for pair, weight in associations.items() if weight > 0}

    graph = nx.Graph()
    graph.add_nodes_from(animals)
    graph.add_weighted_edges_from((a, b, w) for (a, b), w in positive_associations.items())
    modularity, n_communities, largest_community_fraction = safe_greedy_modularity(graph)

    n_timestamps = g["window_start"].nunique()
    n_observed_dyads = len(associations)
    association_density = (
        float(np.mean(list(associations.values()))) if associations else np.nan
    )
    positive_edge_density = (
        len(positive_associations) / n_observed_dyads if n_observed_dyads else np.nan
    )

    return {
        "n_animals_observed": n_animals,
        "n_timestamps": int(n_timestamps),
        "possible_dyads": int(n_possible_dyads),
        "observed_dyads": int(n_observed_dyads),
        "positive_edges": int(len(positive_associations)),
        "positive_edge_density": positive_edge_density,
        "association_density": association_density,
        "mean_positive_association": (
            float(np.mean(list(positive_associations.values())))
            if positive_associations
            else np.nan
        ),
        "total_co_observed_pair_timestamps": int(sum(co_observed.values())),
        "total_co_clustered_pair_timestamps": int(sum(co_clustered.values())),
        "modularity": modularity,
        "n_communities": n_communities,
        "largest_community_fraction": largest_community_fraction,
        "split_timestamp_fraction": split_timestamps / n_timestamps if n_timestamps else np.nan,
        "multi_animal_split_timestamp_fraction": (
            multi_animal_split_timestamps / n_timestamps if n_timestamps else np.nan
        ),
        "mean_largest_cluster_fraction": (
            float(np.mean(largest_cluster_fractions)) if largest_cluster_fractions else np.nan
        ),
    }


def compute_weekly_metrics(df: pd.DataFrame, scope_name: str, hour: int | None) -> pd.DataFrame:
    if hour is None:
        use = df[df["window_start"].dt.hour.isin(ANCHOR_HOURS)].copy()
    else:
        use = df[df["window_start"].dt.hour.eq(hour)].copy()

    use = use.dropna(subset=["dynamic_social_unit", "temp_group_id", "animal_id"])
    use["week_start"] = week_start(use["window_start"])

    records: list[dict[str, object]] = []
    grouped = use.groupby(["week_start", "dynamic_social_unit"], sort=True)
    for (period_start, group), g in grouped:
        metrics = network_metrics_one_group(g)
        records.append(
            {
                "scope": scope_name,
                "period_start": period_start,
                "dynamic_social_unit": group,
                **metrics,
            }
        )
    return pd.DataFrame.from_records(records)


def normalize_group_name(name: str) -> str:
    return "".join(str(name).split()).lower()


def nominal_size_class(total_size: int | float) -> str:
    if pd.isna(total_size):
        return "unknown"
    if total_size < 40:
        return "low (<40)"
    if total_size > 60:
        return "high (>60)"
    if 41 <= total_size <= 59:
        return "medium (41-59)"
    return "boundary (=40 or 60)"


def supplied_group_size_table() -> pd.DataFrame:
    sizes = pd.DataFrame(
        SUPPLIED_GROUP_SIZES,
        columns=["group_name_supplied", "supplied_lower_count", "origin_group_total_size"],
    )
    sizes["group_name_key"] = sizes["group_name_supplied"].map(normalize_group_name)
    sizes["group_size_class"] = sizes["origin_group_total_size"].map(nominal_size_class)
    return sizes


def add_group_sizes(metrics: pd.DataFrame) -> pd.DataFrame:
    sizes = supplied_group_size_table()
    out = metrics.copy()
    out["group_name_key"] = out["dynamic_social_unit"].map(normalize_group_name)
    out = out.merge(sizes, on="group_name_key", how="left")
    out["group_size_class"] = out["group_size_class"].fillna("unknown")
    out["origin_group_total_size"] = out["origin_group_total_size"].astype("Float64")
    out["supplied_lower_count"] = out["supplied_lower_count"].astype("Float64")
    return out


def add_size_adjusted_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    out = metrics.copy()
    for scope, scope_df in out.groupby("scope"):
        idx = scope_df.index
        for metric in ("association_density", "modularity"):
            valid = scope_df[metric].notna() & scope_df["n_animals_observed"].notna()
            if valid.sum() < 5:
                out.loc[idx, f"{metric}_size_adjusted"] = np.nan
                out.loc[idx, f"{metric}_size_slope"] = np.nan
                continue
            x = scope_df.loc[valid, "n_animals_observed"].astype(float).to_numpy()
            y = scope_df.loc[valid, metric].astype(float).to_numpy()
            slope, intercept = np.polyfit(x, y, deg=1)
            ref_n = float(np.nanmedian(x))
            adjusted = scope_df[metric].astype(float) - slope * (
                scope_df["n_animals_observed"].astype(float) - ref_n
            )
            out.loc[idx, f"{metric}_size_adjusted"] = adjusted
            out.loc[idx, f"{metric}_size_slope"] = slope
            out.loc[idx, f"{metric}_reference_n"] = ref_n
    return out


def kde_curve(values: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.zeros_like(x_grid)
    if values.size == 1:
        bw = max((x_grid.max() - x_grid.min()) / 35, 0.025)
    else:
        std = np.std(values, ddof=1)
        bw = 1.06 * std * values.size ** (-1 / 5)
        bw = max(bw, (x_grid.max() - x_grid.min()) / 60, 0.025)
    z = (x_grid[:, None] - values[None, :]) / bw
    dens = np.exp(-0.5 * z * z).sum(axis=1) / (values.size * bw * np.sqrt(2 * np.pi))
    return dens


def group_order_for_plot(metrics: pd.DataFrame, value_col: str) -> list[str]:
    valid = metrics[metrics["valid_for_ridge"]].copy()
    summary = (
        valid.groupby("dynamic_social_unit")
        .agg(
            median_metric=(value_col, "median"),
            n_weeks=(value_col, "count"),
            nominal_size=("origin_group_total_size", "max"),
            group_size_class=("group_size_class", "first"),
        )
        .reset_index()
    )
    summary = summary.sort_values(["nominal_size", "dynamic_social_unit"], na_position="first")
    return summary["dynamic_social_unit"].tolist()


def group_axis_labels(metrics: pd.DataFrame, groups: list[str]) -> list[str]:
    labels = []
    for group in groups:
        rows = metrics[metrics["dynamic_social_unit"].eq(group)]
        if rows.empty or rows["origin_group_total_size"].isna().all():
            labels.append(group)
            continue
        size = int(rows["origin_group_total_size"].dropna().iloc[0])
        labels.append(f"{group} ({size})")
    return labels


def draw_ridge_panel(
    ax: plt.Axes,
    metrics: pd.DataFrame,
    groups: list[str],
    value_col: str,
    title: str,
    x_label: str,
    xlim: tuple[float, float],
) -> None:
    x_grid = np.linspace(xlim[0], xlim[1], 400)
    y_positions = np.arange(len(groups))

    for y, group in zip(y_positions, groups):
        row = metrics[(metrics["dynamic_social_unit"] == group) & metrics["valid_for_ridge"]]
        vals = row[value_col].dropna().astype(float).to_numpy()
        size_class = row["group_size_class"].dropna().iloc[0] if len(row) else "unknown"
        color = SIZE_CLASS_COLORS.get(size_class, SIZE_CLASS_COLORS["unknown"])
        if vals.size == 0:
            continue

        dens = kde_curve(vals, x_grid)
        if dens.max() > 0:
            dens = dens / dens.max() * 0.72
        ax.fill_between(x_grid, y, y + dens, color=color, alpha=0.55, linewidth=0)
        ax.plot(x_grid, y + dens, color=color, linewidth=1.4, alpha=0.95)

        if vals.size > 80:
            rng = np.random.default_rng(42)
            rug_vals = rng.choice(vals, size=80, replace=False)
        else:
            rug_vals = vals
        ax.vlines(rug_vals, y - 0.13, y - 0.03, color="#303030", alpha=0.35, linewidth=0.8)
        ax.scatter(
            [np.nanmedian(vals)],
            [y + 0.10],
            s=22,
            color="#111111",
            edgecolor="white",
            linewidth=0.4,
            zorder=5,
        )

    ax.set_title(title, fontsize=13, pad=10)
    ax.set_xlabel(x_label, fontsize=11)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.45, len(groups) - 0.1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(group_axis_labels(metrics, groups), fontsize=9)
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_ridges(
    metrics: pd.DataFrame,
    scope: str,
    density_col: str,
    modularity_col: str,
    out_path: Path,
    adjusted: bool,
) -> None:
    use = metrics[metrics["scope"].eq(scope)].copy()
    order = group_order_for_plot(use, modularity_col)
    if not order:
        raise ValueError(f"No valid groups to plot for {scope}")

    if adjusted:
        xmins = [
            use.loc[use["valid_for_ridge"], density_col].min(),
            use.loc[use["valid_for_ridge"], modularity_col].min(),
        ]
        xmaxs = [
            use.loc[use["valid_for_ridge"], density_col].max(),
            use.loc[use["valid_for_ridge"], modularity_col].max(),
        ]
        density_xlim = (max(-0.15, float(xmins[0]) - 0.05), min(1.15, float(xmaxs[0]) + 0.05))
        modularity_xlim = (max(-0.15, float(xmins[1]) - 0.05), min(1.15, float(xmaxs[1]) + 0.05))
        subtitle = "size adjusted to median tracked-animal count"
    else:
        density_xlim = (0.0, 1.0)
        modularity_xlim = (0.0, 1.0)
        subtitle = "raw weekly values"

    height = max(7, 0.36 * len(order) + 2.2)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, height), sharey=True)
    draw_ridge_panel(
        axes[0],
        use,
        order,
        density_col,
        "Density",
        "Mean dyadic association",
        density_xlim,
    )
    draw_ridge_panel(
        axes[1],
        use,
        order,
        modularity_col,
        "Modularity",
        "Weighted modularity",
        modularity_xlim,
    )
    axes[1].tick_params(axis="y", labelleft=False)

    handles = [
        plt.Line2D(
            [0],
            [0],
            color=color,
            lw=8,
            alpha=0.55,
            label=label.replace(" (", "\n("),
        )
        for label, color in SIZE_CLASS_COLORS.items()
        if label != "unknown"
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.01),
        fontsize=9,
        title="Nominal total group-size class",
        title_fontsize=9,
    )
    fig.suptitle(
        f"Within-group social structure by group: {scope.replace('_', ' ')}\n"
        f"{subtitle}; weeks with >= {MIN_WEEKLY_ANIMALS} tracked animals and >= {MIN_WEEKLY_TIMESTAMPS} timestamps",
        fontsize=15,
        y=0.985,
    )
    fig.subplots_adjust(left=0.17, right=0.98, top=0.89, bottom=0.13, wspace=0.09)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_dendrogram_ridges(metrics: pd.DataFrame, scope: str, out_path: Path) -> None:
    use = metrics[(metrics["scope"].eq(scope)) & metrics["valid_for_ridge"]].copy()
    summary = (
        use.groupby("dynamic_social_unit")
        .agg(
            density=("association_density_size_adjusted", "median"),
            modularity=("modularity_size_adjusted", "median"),
            split_fraction=("split_timestamp_fraction", "median"),
            largest_cluster_fraction=("mean_largest_cluster_fraction", "median"),
            group_size_class=("group_size_class", "first"),
            nominal_size=("origin_group_total_size", "max"),
            n_weeks=("period_start", "nunique"),
        )
        .dropna(subset=["density", "modularity", "split_fraction", "largest_cluster_fraction"])
    )
    if len(summary) < 3:
        return

    features = summary[
        ["density", "modularity", "split_fraction", "largest_cluster_fraction"]
    ].astype(float)
    standardized = (features - features.mean()) / features.std(ddof=0).replace(0, 1)
    distances = pdist(standardized.to_numpy(), metric="euclidean")
    z = linkage(distances, method="average")

    fig = plt.figure(figsize=(15.5, max(7, 0.38 * len(summary) + 2.3)))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 3.5, 3.2], wspace=0.22)
    ax_tree = fig.add_subplot(gs[0, 0])
    dendro = dendrogram(
        z,
        labels=summary.index.tolist(),
        orientation="left",
        ax=ax_tree,
        no_labels=True,
        leaf_font_size=9,
        color_threshold=0,
        above_threshold_color="#3A3A3A",
    )
    ordered_groups = summary.index[dendro["leaves"]].tolist()
    ax_tree.set_xlabel("Distance", fontsize=10)
    ax_tree.tick_params(axis="x", labelsize=8)
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)
    ax_tree.spines["left"].set_visible(False)

    ax_density = fig.add_subplot(gs[0, 1])
    ax_mod = fig.add_subplot(gs[0, 2], sharey=ax_density)

    draw_ridge_panel(
        ax_density,
        use,
        ordered_groups,
        "association_density_size_adjusted",
        "Density",
        "Mean dyadic association\n(size adjusted)",
        (
            max(-0.15, float(use["association_density_size_adjusted"].min()) - 0.05),
            min(1.15, float(use["association_density_size_adjusted"].max()) + 0.05),
        ),
    )
    draw_ridge_panel(
        ax_mod,
        use,
        ordered_groups,
        "modularity_size_adjusted",
        "Modularity",
        "Weighted modularity\n(size adjusted)",
        (
            max(-0.15, float(use["modularity_size_adjusted"].min()) - 0.05),
            min(1.15, float(use["modularity_size_adjusted"].max()) + 0.05),
        ),
    )
    ax_mod.tick_params(axis="y", labelleft=False)

    handles = [
        plt.Line2D([0], [0], color=color, lw=8, alpha=0.55, label=label)
        for label, color in SIZE_CLASS_COLORS.items()
        if label != "unknown"
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.58, 0.01),
        fontsize=9,
        title="Nominal total group-size class",
        title_fontsize=9,
    )
    fig.suptitle(
        f"Within-group density and modularity across groups: {scope.replace('_', ' ')}\n"
        "dendrogram clusters groups by median size-adjusted density/modularity and split structure",
        fontsize=15,
        y=0.985,
    )
    fig.subplots_adjust(left=0.04, right=0.98, top=0.88, bottom=0.14)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def write_summary(metrics: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    valid = metrics[metrics["valid_for_ridge"]].copy()
    summary = (
        valid.groupby(["scope", "dynamic_social_unit", "group_size_class"], dropna=False)
        .agg(
            nominal_group_size=("origin_group_total_size", "max"),
            supplied_lower_count=("supplied_lower_count", "max"),
            n_weeks=("period_start", "nunique"),
            median_tracked_animals=("n_animals_observed", "median"),
            median_density=("association_density", "median"),
            iqr_density=("association_density", lambda x: x.quantile(0.75) - x.quantile(0.25)),
            median_modularity=("modularity", "median"),
            iqr_modularity=("modularity", lambda x: x.quantile(0.75) - x.quantile(0.25)),
            median_density_size_adjusted=("association_density_size_adjusted", "median"),
            median_modularity_size_adjusted=("modularity_size_adjusted", "median"),
            median_split_fraction=("split_timestamp_fraction", "median"),
            median_largest_cluster_fraction=("mean_largest_cluster_fraction", "median"),
        )
        .reset_index()
        .sort_values(["scope", "nominal_group_size", "dynamic_social_unit"], na_position="first")
    )
    summary.to_csv(out_dir / "canonical_within_group_density_modularity_summary_by_group.csv", index=False)
    return summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    supplied_group_size_table().to_csv(
        OUT_DIR / "supplied_nominal_total_group_sizes_used.csv", index=False
    )

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
    metrics["valid_for_ridge"] = (
        metrics["n_animals_observed"].ge(MIN_WEEKLY_ANIMALS)
        & metrics["n_timestamps"].ge(MIN_WEEKLY_TIMESTAMPS)
        & metrics["association_density"].notna()
        & metrics["modularity"].notna()
        & metrics["group_size_class"].ne("unknown")
    )
    group_week_counts = (
        metrics[metrics["valid_for_ridge"]]
        .groupby(["scope", "dynamic_social_unit"])["period_start"]
        .transform("nunique")
    )
    metrics.loc[metrics["valid_for_ridge"], "valid_group_weeks"] = group_week_counts
    metrics["valid_for_ridge"] = metrics["valid_for_ridge"] & metrics["valid_group_weeks"].fillna(0).ge(
        MIN_GROUP_WEEKS_FOR_RIDGE
    )
    metrics = add_size_adjusted_metrics(metrics)

    metrics.to_csv(OUT_DIR / "canonical_within_group_weekly_network_metrics.csv", index=False)
    summary = write_summary(metrics, OUT_DIR)

    for scope in SCOPE_DEFS:
        plot_ridges(
            metrics,
            scope=scope,
            density_col="association_density",
            modularity_col="modularity",
            out_path=OUT_DIR / f"{scope}_density_modularity_ridges_raw.png",
            adjusted=False,
        )
        plot_ridges(
            metrics,
            scope=scope,
            density_col="association_density_size_adjusted",
            modularity_col="modularity_size_adjusted",
            out_path=OUT_DIR / f"{scope}_density_modularity_ridges_size_adjusted.png",
            adjusted=True,
        )

    print(f"Wrote weekly metrics: {OUT_DIR / 'canonical_within_group_weekly_network_metrics.csv'}")
    print(f"Wrote group summary: {OUT_DIR / 'canonical_within_group_density_modularity_summary_by_group.csv'}")
    print(f"Wrote supplied group sizes: {OUT_DIR / 'supplied_nominal_total_group_sizes_used.csv'}")
    print(summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
