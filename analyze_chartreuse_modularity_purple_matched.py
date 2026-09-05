from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path(
    r"C:\Users\rharel\Documents\New project\outputs"
    r"\canonical_robust_hourly_membership_shared_full_20260722"
    r"\canonical_hourly_membership_with_association_events.parquet"
)
DEFAULT_OUTPUT = Path("outputs/chartreuse_modularity_purple_matched")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare internal Chartreuse modularity with Purple present to matched reference hours.")
    parser.add_argument("--input-file", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--min-collars", type=int, default=8)
    parser.add_argument("--max-date-gap-days", type=int, default=30)
    parser.add_argument("--max-distance-m", type=float, default=2000)
    parser.add_argument("--max-collar-difference", type=int, default=2)
    parser.add_argument("--bootstrap-replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260808)
    return parser.parse_args()


def haversine_m(lon1: float, lat1: float, lon2: pd.Series, lat2: pd.Series) -> np.ndarray:
    radius = 6_371_000.0
    lon1r, lat1r = np.radians(lon1), np.radians(lat1)
    lon2r, lat2r = np.radians(lon2.to_numpy()), np.radians(lat2.to_numpy())
    a = np.sin((lat2r - lat1r) / 2) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin((lon2r - lon1r) / 2) ** 2
    return radius * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def make_hour_table(rows: pd.DataFrame, min_collars: int) -> pd.DataFrame:
    rows = rows.copy()
    rows["purple_event_id"] = rows["association_event_id"].where(
        rows["association_id"].eq("MERGE_Chartreuse_Purple")
    )
    cluster_sizes = (
        rows.groupby(["window_start", "observed_cluster_id"], observed=True)["animal_id"]
        .nunique().rename("cluster_size").reset_index()
    )
    split = cluster_sizes.groupby("window_start", observed=True).agg(
        n_clusters=("observed_cluster_id", "nunique"),
        n_multi_animal_clusters=("cluster_size", lambda x: int((x >= 2).sum())),
        largest_cluster=("cluster_size", "max"),
    )
    hour = rows.groupby("window_start", observed=True).agg(
        n_collars=("animal_id", "nunique"),
        centroid_longitude=("longitude", "median"),
        centroid_latitude=("latitude", "median"),
        purple_present=("association_id", lambda x: bool((x == "MERGE_Chartreuse_Purple").any())),
        purple_event_id=("purple_event_id", lambda x: x.dropna().mode().iloc[0] if not x.dropna().empty else None),
    ).join(split).reset_index()
    hour["eligible"] = hour["n_collars"].ge(min_collars)
    hour["within_split"] = hour["eligible"] & hour["n_multi_animal_clusters"].ge(2)
    hour["largest_fraction"] = hour["largest_cluster"] / hour["n_collars"]
    return hour


def greedy_match(encounter: pd.DataFrame, reference: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    options: dict[pd.Timestamp, pd.DataFrame] = {}
    for rec in encounter.itertuples(index=False):
        candidate = reference[
            reference["window_start"].dt.hour.eq(rec.window_start.hour)
            & reference["window_start"].sub(rec.window_start).abs().le(pd.Timedelta(days=args.max_date_gap_days))
            & reference["n_collars"].sub(rec.n_collars).abs().le(args.max_collar_difference)
        ].copy()
        candidate["spatial_distance_m"] = haversine_m(
            rec.centroid_longitude, rec.centroid_latitude,
            candidate["centroid_longitude"], candidate["centroid_latitude"],
        )
        candidate = candidate[candidate["spatial_distance_m"].le(args.max_distance_m)].copy()
        candidate["date_gap_days"] = candidate["window_start"].sub(rec.window_start).abs().dt.total_seconds() / 86400
        candidate["collar_difference"] = candidate["n_collars"].sub(rec.n_collars).abs()
        candidate["match_score"] = (
            candidate["spatial_distance_m"] / 500
            + candidate["date_gap_days"] / 14
            + candidate["collar_difference"] / 2
        )
        options[rec.window_start] = candidate.sort_values("match_score")

    used: set[pd.Timestamp] = set()
    matched = []
    order = sorted(options, key=lambda t: len(options[t]))
    encounter_lookup = encounter.set_index("window_start")
    for encounter_time in order:
        available = options[encounter_time][~options[encounter_time]["window_start"].isin(used)]
        if available.empty:
            continue
        chosen = available.iloc[0]
        used.add(chosen["window_start"])
        event = encounter_lookup.loc[encounter_time]
        matched.append({
            "encounter_time": encounter_time,
            "encounter_event_id": event["purple_event_id"],
            "reference_time": chosen["window_start"],
            "hour_of_day": encounter_time.hour,
            "spatial_distance_m": chosen["spatial_distance_m"],
            "date_gap_days": chosen["date_gap_days"],
            "encounter_n_collars": event["n_collars"],
            "reference_n_collars": chosen["n_collars"],
            "encounter_within_split": event["within_split"],
            "reference_within_split": chosen["within_split"],
            "encounter_largest_fraction": event["largest_fraction"],
            "reference_largest_fraction": chosen["largest_fraction"],
            "match_score": chosen["match_score"],
        })
    return pd.DataFrame(matched).sort_values("encounter_time")


def network_metrics(rows: pd.DataFrame, hours: list[pd.Timestamp]) -> tuple[dict[str, float | int], pd.DataFrame]:
    work = rows[rows["window_start"].isin(set(hours))]
    hour_lookup = {time: frame for time, frame in work.groupby("window_start", observed=True, sort=False)}
    observed: dict[tuple[str, str], int] = {}
    together: dict[tuple[str, str], int] = {}
    for hour in hours:
        frame = hour_lookup.get(hour)
        if frame is None:
            continue
        animals = sorted(frame["animal_id"].astype(str).unique())
        for pair in combinations(animals, 2):
            observed[pair] = observed.get(pair, 0) + 1
        for _, cluster in frame.groupby("observed_cluster_id", observed=True):
            members = sorted(cluster["animal_id"].astype(str).unique())
            for pair in combinations(members, 2):
                together[pair] = together.get(pair, 0) + 1
    graph = nx.Graph()
    graph.add_nodes_from(sorted({animal for pair in observed for animal in pair}))
    edge_rows = []
    for pair, n_observed in sorted(observed.items()):
        n_together = together.get(pair, 0)
        weight = n_together / n_observed
        edge_rows.append({"animal_a": pair[0], "animal_b": pair[1], "observed_hours": n_observed, "same_cluster_hours": n_together, "coassociation_weight": weight})
        if n_together:
            graph.add_edge(*pair, weight=weight)
    if graph.number_of_edges() == 0 or graph.number_of_nodes() < 3:
        modularity, communities = np.nan, []
    else:
        communities = list(nx.algorithms.community.louvain_communities(graph, weight="weight", seed=42))
        modularity = float(nx.algorithms.community.modularity(graph, communities, weight="weight"))
    metrics = {
        "hours": len(hours),
        "animals": graph.number_of_nodes(),
        "observed_pairs": len(observed),
        "network_edges": graph.number_of_edges(),
        "edge_density": float(nx.density(graph)) if graph.number_of_nodes() > 1 else np.nan,
        "mean_coassociation_weight": float(np.mean([x["coassociation_weight"] for x in edge_rows])) if edge_rows else np.nan,
        "modularity": modularity,
        "n_communities": len(communities),
    }
    return metrics, pd.DataFrame(edge_rows)


def paired_bootstrap(rows: pd.DataFrame, matches: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    event_ids = matches["encounter_event_id"].drop_duplicates().to_numpy()
    draws = []
    for i in range(replicates):
        sampled_events = rng.choice(event_ids, size=len(event_ids), replace=True)
        sample = pd.concat(
            [matches[matches["encounter_event_id"].eq(event_id)] for event_id in sampled_events],
            ignore_index=True,
        )
        encounter_metrics, _ = network_metrics(rows, sample["encounter_time"].tolist())
        reference_metrics, _ = network_metrics(rows, sample["reference_time"].tolist())
        encounter_split = sample["encounter_within_split"].mean()
        reference_split = sample["reference_within_split"].mean()
        encounter_largest = sample["encounter_largest_fraction"].mean()
        reference_largest = sample["reference_largest_fraction"].mean()
        draws.append({
            "replicate": i,
            "purple_present_modularity": encounter_metrics["modularity"],
            "purple_absent_modularity": reference_metrics["modularity"],
            "modularity_difference": encounter_metrics["modularity"] - reference_metrics["modularity"],
            "purple_present_split_probability": encounter_split,
            "purple_absent_split_probability": reference_split,
            "split_probability_difference": encounter_split - reference_split,
            "purple_present_largest_fraction": encounter_largest,
            "purple_absent_largest_fraction": reference_largest,
            "largest_fraction_difference": encounter_largest - reference_largest,
        })
    return pd.DataFrame(draws)


def plot_results(summary: pd.DataFrame, bootstrap: pd.DataFrame, output: Path, n_events: int) -> None:
    contexts = ["Purple present", "Matched Purple absent"]
    colors = ["#7b2cbf", "#777777"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    metrics = [
        ("modularity", "Internal co-association modularity", "purple_present_modularity", "purple_absent_modularity"),
        ("split_probability", "Probability of a multi-animal split", "purple_present_split_probability", "purple_absent_split_probability"),
        ("mean_largest_fraction", "Mean fraction in largest component", "purple_present_largest_fraction", "purple_absent_largest_fraction"),
    ]
    for ax, (metric, title, present_col, absent_col) in zip(axes, metrics):
        values = summary.set_index("context").loc[contexts, metric]
        lows = np.array([bootstrap[present_col].quantile(0.025), bootstrap[absent_col].quantile(0.025)])
        highs = np.array([bootstrap[present_col].quantile(0.975), bootstrap[absent_col].quantile(0.975)])
        errors = np.vstack([values.to_numpy() - lows, highs - values.to_numpy()])
        ax.bar(contexts, values, color=colors, yerr=errors, capsize=5, error_kw={"lw": 1.3})
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=18)
        for i, value in enumerate(values):
            ax.text(i, value, f" {value:.3f}", ha="center", va="bottom", fontsize=9)
    fig.suptitle(
        f"Internal Chartreuse structure when Purple is present versus matched reference hours\n"
        f"{n_events} encounter events; error bars are event-clustered bootstrap 95% CIs"
    )
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    columns = ["animal_id", "window_start", "dynamic_social_unit", "longitude", "latitude", "observed_cluster_id", "association_id", "association_event_id"]
    data = pd.read_parquet(args.input_file, columns=columns)
    data["window_start"] = pd.to_datetime(data["window_start"])
    chartreuse = data[data["dynamic_social_unit"].eq("Chartreuse")].dropna(subset=["observed_cluster_id"]).copy()
    hour = make_hour_table(chartreuse, args.min_collars)
    encounter = hour[hour["eligible"] & hour["purple_present"]].copy()
    reference = hour[hour["eligible"] & ~hour["purple_present"]].copy()
    matches = greedy_match(encounter, reference, args)
    if len(matches) != len(encounter):
        raise ValueError(f"Only {len(matches)} of {len(encounter)} Purple-present hours could be matched")

    encounter_metrics, encounter_edges = network_metrics(chartreuse, matches["encounter_time"].tolist())
    reference_metrics, reference_edges = network_metrics(chartreuse, matches["reference_time"].tolist())
    summary = pd.DataFrame([
        {"context": "Purple present", **encounter_metrics,
         "split_probability": matches["encounter_within_split"].mean(),
         "mean_largest_fraction": matches["encounter_largest_fraction"].mean()},
        {"context": "Matched Purple absent", **reference_metrics,
         "split_probability": matches["reference_within_split"].mean(),
         "mean_largest_fraction": matches["reference_largest_fraction"].mean()},
    ])
    bootstrap = paired_bootstrap(chartreuse, matches, args.bootstrap_replicates, args.seed)
    intervals = bootstrap.drop(columns="replicate").quantile([0.025, 0.5, 0.975]).reset_index().rename(columns={"index": "quantile"})

    matches.to_csv(args.output_dir / "chartreuse_purple_matched_hours.csv", index=False)
    summary.to_csv(args.output_dir / "chartreuse_purple_matched_modularity_summary.csv", index=False)
    intervals.to_csv(args.output_dir / "chartreuse_purple_matched_bootstrap_intervals.csv", index=False)
    encounter_edges.assign(context="Purple present").to_csv(args.output_dir / "chartreuse_purple_present_network_edges.csv", index=False)
    reference_edges.assign(context="Matched Purple absent").to_csv(args.output_dir / "chartreuse_purple_reference_network_edges.csv", index=False)
    n_events = int(matches["encounter_event_id"].nunique())
    plot_results(summary, bootstrap, args.output_dir / "chartreuse_modularity_purple_matched.png", n_events)

    metadata = {
        "input_file": str(args.input_file.resolve()),
        "encounter_definition": "association_id == MERGE_Chartreuse_Purple for any Chartreuse animal in the hour",
        "reference_definition": "Purple absent; exact hour of day; within 30 days and 2 km; collar count within 2; greedy without replacement",
        "eligible_definition": f"at least {args.min_collars} Chartreuse collars",
        "matched_pairs": len(matches),
        "encounter_events": n_events,
        "median_spatial_distance_m": float(matches["spatial_distance_m"].median()),
        "median_date_gap_days": float(matches["date_gap_days"].median()),
        "exact_collar_count_fraction": float(matches["encounter_n_collars"].eq(matches["reference_n_collars"]).mean()),
        "modularity_definition": "Louvain modularity of weighted Chartreuse pair co-association network; weight is same-cluster hours / co-observed hours",
        "observed_difference_purple_minus_reference": {
            "modularity": encounter_metrics["modularity"] - reference_metrics["modularity"],
            "split_probability": float(matches["encounter_within_split"].mean() - matches["reference_within_split"].mean()),
            "largest_component_fraction": float(matches["encounter_largest_fraction"].mean() - matches["reference_largest_fraction"].mean()),
        },
        "bootstrap_95_intervals": intervals.to_dict(orient="records"),
        "bootstrap_cluster_unit": "Chartreuse-Purple association_event_id",
        "caution": "This is a matched observational comparison, not a causal estimate of Purple presence.",
    }
    (args.output_dir / "chartreuse_modularity_purple_matched_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
