from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_METRICS = Path(
    r"C:\Users\rharel\Documents\New project\outputs"
    r"\canonical_robust_hourly_membership_shared_full_20260722"
    r"\canonical_2m_shared_history_shuffle_expectation"
    r"\canonical_5m_shuffle_expectation_2min_rows.csv"
)
DEFAULT_PANEL = Path(
    "outputs/chartreuse_social_reorganization/chartreuse_hourly_identity_panel.parquet"
)
DEFAULT_OUTPUT = Path("outputs/chartreuse_2m_integration_context")
CONTEXT_ORDER = ["splinter_present", "disperser_present"]
COLORS = {"splinter_present": "#e69f00", "disperser_present": "#7b3294"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Relate Chartreuse 2 m integration to hourly social context.")
    parser.add_argument("--metrics-file", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--panel-file", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--pair-key", default="Chartreuse - Purple")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-replicates", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260808)
    return parser.parse_args()


def summarize_context(rows: pd.DataFrame) -> pd.DataFrame:
    summary = (
        rows.groupby("social_context", observed=True)
        .agg(
            hours=("timestamp", "nunique"),
            days=("date", "nunique"),
            bins_2min=("bin_2min", "size"),
            positive_bins=("cross_edges", lambda x: int((x > 0).sum())),
            cross_edges=("cross_edges", "sum"),
            total_edges=("total_edges", "sum"),
            mean_shuffle_expected=("shuffle_expected_cross_edge_fraction", "mean"),
            mean_obs_minus_shuffle=("cross_edge_fraction_obs_minus_shuffle", "mean"),
            mean_modularity=("edge_modularity_q", "mean"),
            mean_composition_entropy=("composition_entropy_norm", "mean"),
            mean_pair_balance=("pair_balance", "mean"),
            mean_chartreuse_splinter_n=("n_splinter", "mean"),
            mean_chartreuse_isolated_n=("n_isolated", "mean"),
            mean_chartreuse_dispersed_n=("n_dispersed", "mean"),
        )
        .reset_index()
    )
    summary["weighted_observed_cross_edge_fraction"] = (
        summary["cross_edges"] / summary["total_edges"].replace(0, np.nan)
    )
    summary["positive_bin_fraction"] = summary["positive_bins"] / summary["bins_2min"]
    return summary


def cluster_bootstrap_difference(
    rows: pd.DataFrame, replicates: int, seed: int
) -> dict[str, float | int | str]:
    contexts = [c for c in CONTEXT_ORDER if c in set(rows["social_context"])]
    if len(contexts) != 2:
        return {"estimand": "unavailable", "reason": "fewer than two social contexts"}
    first, second = contexts
    by_day = (
        rows.groupby(["date", "social_context"], observed=True)
        .agg(cross_edges=("cross_edges", "sum"), total_edges=("total_edges", "sum"))
        .reset_index()
    )
    days = by_day["date"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed)

    def estimate(frame: pd.DataFrame) -> float:
        totals = frame.groupby("social_context", observed=True)[["cross_edges", "total_edges"]].sum()
        if first not in totals.index or second not in totals.index:
            return np.nan
        a = totals.loc[first, "cross_edges"] / totals.loc[first, "total_edges"]
        b = totals.loc[second, "cross_edges"] / totals.loc[second, "total_edges"]
        return float(a - b)

    observed = estimate(by_day)
    draws = []
    for _ in range(replicates):
        sampled = rng.choice(days, size=len(days), replace=True)
        pieces = [by_day[by_day["date"].eq(day)] for day in sampled]
        draws.append(estimate(pd.concat(pieces, ignore_index=True)))
    draws = np.asarray(draws, dtype=float)
    draws = draws[np.isfinite(draws)]
    return {
        "estimand": f"weighted 2m cross-edge fraction: {first} minus {second}",
        "estimate": observed,
        "cluster_unit": "calendar day",
        "days": int(len(days)),
        "valid_bootstrap_replicates": int(len(draws)),
        "ci_95_low": float(np.quantile(draws, 0.025)),
        "ci_95_high": float(np.quantile(draws, 0.975)),
    }


def make_plot(hourly: pd.DataFrame, summary: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
    use = summary.set_index("social_context").reindex(CONTEXT_ORDER).dropna(how="all")
    labels = [x.replace("_", " ") for x in use.index]
    colors = [COLORS[x] for x in use.index]
    axes[0].bar(labels, use["weighted_observed_cross_edge_fraction"], color=colors)
    for i, row in enumerate(use.itertuples()):
        axes[0].text(
            i,
            row.weighted_observed_cross_edge_fraction + 0.00015,
            f"{int(row.cross_edges)}/{int(row.total_edges)} edges",
            ha="center",
            fontsize=9,
        )
    axes[0].set_ylabel("Observed cross-group edge fraction")
    axes[0].set_title("Observed 2 m Chartreuse–Purple links")
    axes[0].tick_params(axis="x", rotation=15)

    rng = np.random.default_rng(117)
    for i, context in enumerate(CONTEXT_ORDER):
        values = hourly.loc[
            hourly["social_context"].eq(context), "mean_obs_minus_shuffle"
        ].dropna()
        if values.empty:
            continue
        x = i + rng.uniform(-0.14, 0.14, len(values))
        axes[1].scatter(x, values, s=19, alpha=0.55, color=COLORS[context])
        axes[1].plot([i - 0.23, i + 0.23], [values.median()] * 2, color="black", lw=2.2)
    axes[1].axhline(0, color="#555555", linestyle="--", lw=1)
    axes[1].set_xticks(range(len(CONTEXT_ORDER)))
    axes[1].set_xticklabels([x.replace("_", " ") for x in CONTEXT_ORDER], rotation=15)
    axes[1].set_ylabel("Observed minus shuffled cross-edge fraction")
    axes[1].set_title("Hourly integration relative to random mixing")
    fig.suptitle(
        "Chartreuse–Purple integration at 2 m by Chartreuse social context\n"
        "Only hours with a saved canonical merge and usable 2-minute fixes are represented",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(args.metrics_file, parse_dates=["timestamp", "bin_2min"])
    metrics = metrics[metrics["pair_key"].eq(args.pair_key)].copy()
    if metrics.empty:
        raise ValueError(f"No 2 m rows found for pair_key={args.pair_key!r}")

    panel = pd.read_parquet(args.panel_file)
    panel["window_start"] = pd.to_datetime(panel["window_start"])
    context = (
        panel.groupby("window_start", observed=True)
        .agg(
            n_chartreuse_animals=("animal_id", "nunique"),
            n_core=("social_state", lambda x: int((x == "with_chartreuse_core").sum())),
            n_splinter=("social_state", lambda x: int((x == "within_group_splinter").sum())),
            n_isolated=("social_state", lambda x: int((x == "isolated").sum())),
            n_dispersed=("social_state", lambda x: int((x == "dispersed_to_other_group").sum())),
        )
        .reset_index()
    )
    rows = metrics.merge(context, left_on="timestamp", right_on="window_start", how="left", validate="many_to_one")
    if rows["n_chartreuse_animals"].isna().any():
        raise ValueError("Some 2 m rows could not be matched to the Chartreuse hourly panel")
    rows["social_context"] = np.select(
        [rows["n_dispersed"].gt(0), rows["n_splinter"].gt(0), rows["n_isolated"].gt(0)],
        ["disperser_present", "splinter_present", "isolation_present"],
        default="core_only",
    )
    rows["date"] = rows["bin_2min"].dt.floor("D")
    hourly = (
        rows.groupby(["timestamp", "date", "social_context"], observed=True)
        .agg(
            bins_2min=("bin_2min", "size"),
            cross_edges=("cross_edges", "sum"),
            total_edges=("total_edges", "sum"),
            positive_bins=("cross_edges", lambda x: int((x > 0).sum())),
            mean_observed_cross_edge_fraction=("observed_cross_edge_fraction", "mean"),
            mean_shuffle_expected=("shuffle_expected_cross_edge_fraction", "mean"),
            mean_obs_minus_shuffle=("cross_edge_fraction_obs_minus_shuffle", "mean"),
            mean_modularity=("edge_modularity_q", "mean"),
            n_splinter=("n_splinter", "first"),
            n_isolated=("n_isolated", "first"),
            n_dispersed=("n_dispersed", "first"),
        )
        .reset_index()
    )
    hourly["weighted_observed_cross_edge_fraction"] = (
        hourly["cross_edges"] / hourly["total_edges"].replace(0, np.nan)
    )

    summary = summarize_context(rows)
    bootstrap = cluster_bootstrap_difference(rows, args.bootstrap_replicates, args.seed)
    rows.to_csv(args.output_dir / "chartreuse_purple_2m_rows_with_social_context.csv", index=False)
    hourly.to_csv(args.output_dir / "chartreuse_purple_2m_hourly_context.csv", index=False)
    summary.to_csv(args.output_dir / "chartreuse_purple_2m_context_summary.csv", index=False)
    make_plot(hourly, summary, args.output_dir / "chartreuse_purple_2m_integration_by_context.png")

    metadata = {
        "metrics_file": str(args.metrics_file.resolve()),
        "panel_file": str(args.panel_file.resolve()),
        "pair_key": args.pair_key,
        "edge_radius_m": 2,
        "two_minute_rows": len(rows),
        "hours": int(rows["timestamp"].nunique()),
        "days": int(rows["date"].nunique()),
        "positive_cross_edge_rows": int(rows["cross_edges"].gt(0).sum()),
        "contexts_observed": rows["social_context"].value_counts().to_dict(),
        "bootstrap_contrast": bootstrap,
        "limitations": [
            "Saved 2 m support for Chartreuse is available only for the Chartreuse-Purple dyad.",
            "The social-context categories describe Chartreuse at that hour, not the identity of each 2 m edge.",
            "All supported hours contained a splinter or a disperser; there is no cohesive-core control context.",
            "Only nine two-minute rows contain any observed cross-group 2 m edge, so results are descriptive.",
        ],
    }
    (args.output_dir / "chartreuse_2m_integration_context_metadata.json").write_text(
        json.dumps(metadata, indent=2, default=str), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, default=str))


if __name__ == "__main__":
    main()
