from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUT = Path(
    r"C:\Users\rharel\Documents\New project\outputs"
    r"\canonical_robust_hourly_membership_shared_full_20260722"
    r"\canonical_2m_shared_history_shuffle_expectation"
    r"\canonical_5m_shuffle_expectation_2min_rows.csv"
)
OUT = Path("outputs/all_supported_2m_group_merges")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize all group dyads with saved 2 m merge support.")
    p.add_argument("--input-file", type=Path, default=INPUT)
    p.add_argument("--output-dir", type=Path, default=OUT)
    p.add_argument("--bootstrap-replicates", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20260808)
    return p.parse_args()


def bootstrap_events(data: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    event = data.groupby(["pair_key", "event_number_for_pair"], observed=True).agg(
        cross_edges=("cross_edges", "sum"), total_edges=("total_edges", "sum"),
        mean_obs_minus_shuffle=("cross_edge_fraction_obs_minus_shuffle", "mean"),
        mean_modularity=("edge_modularity_q", "mean"), bins=("bin_2min", "size"),
    ).reset_index()
    rng = np.random.default_rng(seed)
    draws = []
    for pair, frame in event.groupby("pair_key", observed=True):
        records = frame.to_dict("records")
        for replicate in range(replicates):
            sample = [records[i] for i in rng.integers(0, len(records), len(records))]
            cross = sum(x["cross_edges"] for x in sample)
            total = sum(x["total_edges"] for x in sample)
            weights = np.array([x["bins"] for x in sample], float)
            draws.append({
                "pair_key": pair, "replicate": replicate,
                "observed_cross_edge_fraction": cross / total if total else np.nan,
                "mean_obs_minus_shuffle": np.average([x["mean_obs_minus_shuffle"] for x in sample], weights=weights),
                "mean_modularity": np.average([x["mean_modularity"] for x in sample], weights=weights),
            })
    return pd.DataFrame(draws)


def plot(summary: pd.DataFrame, draws: pd.DataFrame, output: Path) -> None:
    order = summary.sort_values("observed_cross_edge_fraction")["pair_key"].tolist()
    y = np.arange(len(order))
    fig, axes = plt.subplots(1, 2, figsize=(15, max(7, len(order) * 0.55)), sharey=True)
    specs = [
        ("observed_cross_edge_fraction", "Observed cross-group edge fraction", 0),
        ("mean_obs_minus_shuffle", "Observed minus shuffled expectation", 1),
    ]
    for metric, xlabel, axis_i in specs:
        ax = axes[axis_i]
        point = summary.set_index("pair_key").loc[order, metric]
        grouped = draws.groupby("pair_key")[metric]
        low = grouped.quantile(0.025).reindex(order)
        high = grouped.quantile(0.975).reindex(order)
        ax.errorbar(point, y, xerr=np.vstack([point - low, high - point]), fmt="o", color="#3b6fb6", ecolor="#777777", capsize=3)
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", color="#eeeeee")
        if metric == "mean_obs_minus_shuffle":
            ax.axvline(0, color="#444444", linestyle="--", lw=1)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([
        f"{pair}\n{int(summary.set_index('pair_key').loc[pair, 'events'])} events; "
        f"{int(summary.set_index('pair_key').loc[pair, 'bins_2min'])} bins"
        for pair in order
    ], fontsize=9)
    fig.suptitle("Integration at 2 m across all supported group-merge dyads\nEvent-bootstrap 95% confidence intervals")
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input_file, parse_dates=["timestamp", "bin_2min"])
    summary = data.groupby("pair_key", observed=True).agg(
        events=("event_number_for_pair", "nunique"), hours=("timestamp", "nunique"), days=("bin_2min", lambda x: x.dt.floor("D").nunique()),
        bins_2min=("bin_2min", "size"), positive_bins=("cross_edges", lambda x: int((x > 0).sum())),
        cross_edges=("cross_edges", "sum"), total_edges=("total_edges", "sum"),
        mean_obs_minus_shuffle=("cross_edge_fraction_obs_minus_shuffle", "mean"),
        mean_modularity=("edge_modularity_q", "mean"), mean_composition_entropy=("composition_entropy_norm", "mean"),
    ).reset_index()
    summary["observed_cross_edge_fraction"] = summary["cross_edges"] / summary["total_edges"]
    summary["positive_bin_fraction"] = summary["positive_bins"] / summary["bins_2min"]
    draws = bootstrap_events(data, args.bootstrap_replicates, args.seed)
    intervals = draws.groupby("pair_key").quantile([0.025, 0.5, 0.975]).reset_index()
    summary.to_csv(args.output_dir / "all_supported_2m_merge_dyad_summary.csv", index=False)
    intervals.to_csv(args.output_dir / "all_supported_2m_merge_event_bootstrap_intervals.csv", index=False)
    plot(summary, draws, args.output_dir / "all_supported_2m_group_merge_integration.png")
    metadata = {
        "input_file": str(args.input_file.resolve()), "edge_radius_m": 2,
        "supported_dyads": int(summary["pair_key"].nunique()), "events": int(data.groupby(["pair_key", "event_number_for_pair"]).ngroups),
        "hours": int(data["timestamp"].nunique()), "two_minute_rows": len(data),
        "bootstrap_cluster_unit": "event_number_for_pair within dyad",
        "interpretation": "negative observed-minus-shuffle means fewer 2 m cross-group links than expected from group composition",
    }
    (args.output_dir / "all_supported_2m_group_merges_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
