from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
PANEL = Path("outputs/chartreuse_social_reorganization/chartreuse_hourly_identity_panel.parquet")
MATCHES = Path("outputs/chartreuse_modularity_purple_matched/chartreuse_purple_matched_hours.csv")
OUT = Path("outputs/chartreuse_finescale_clustering")
RADII = (2.0, 5.0, 10.0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fine-scale within-Chartreuse proximity clustering.")
    p.add_argument("--gps-file", type=Path, default=GPS)
    p.add_argument("--panel-file", type=Path, default=PANEL)
    p.add_argument("--matches-file", type=Path, default=MATCHES)
    p.add_argument("--output-dir", type=Path, default=OUT)
    p.add_argument("--bootstrap-replicates", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20260808)
    return p.parse_args()


def pairwise_m(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    radius = 6_371_000.0
    lonr, latr = np.radians(lon), np.radians(lat)
    dlon = lonr[:, None] - lonr[None, :]
    dlat = latr[:, None] - latr[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(latr[:, None]) * np.cos(latr[None, :]) * np.sin(dlon / 2) ** 2
    return radius * 2 * np.arctan2(np.sqrt(a), np.sqrt(np.maximum(0, 1 - a)))


def component_sizes(adjacency: np.ndarray) -> list[int]:
    n = len(adjacency)
    unseen = set(range(n))
    sizes = []
    while unseen:
        start = unseen.pop()
        stack = [start]
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            neighbors = set(np.flatnonzero(adjacency[node])) & unseen
            unseen -= neighbors
            stack.extend(neighbors)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def bin_metrics(frame: pd.DataFrame, radius_m: float) -> dict[str, float | int]:
    n = len(frame)
    distances = pairwise_m(frame["longitude"].to_numpy(), frame["latitude"].to_numpy())
    upper = np.triu_indices(n, 1)
    edges = int((distances[upper] <= radius_m).sum())
    candidates = n * (n - 1) // 2
    adjacency = (distances <= radius_m) & ~np.eye(n, dtype=bool)
    sizes = component_sizes(adjacency)
    multi = [size for size in sizes if size >= 2]
    return {
        "n_animals": n,
        "candidate_pairs": candidates,
        "proximity_edges": edges,
        "edge_density": edges / candidates if candidates else np.nan,
        "n_components": len(sizes),
        "n_multi_animal_components": len(multi),
        "largest_component_size": sizes[0],
        "largest_component_fraction": sizes[0] / n,
        "isolated_fraction": sum(size == 1 for size in sizes) / n,
        "fine_split": len(multi) >= 2,
    }


def bootstrap_events(rows: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    event = (
        rows.groupby(["encounter_event_id", "context", "radius_m"], observed=True)
        .agg(
            edge_density=("edge_density", "mean"),
            largest_component_fraction=("largest_component_fraction", "mean"),
            isolated_fraction=("isolated_fraction", "mean"),
            fine_split_probability=("fine_split", "mean"),
        )
        .reset_index()
    )
    ids = event["encounter_event_id"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed)
    draws = []
    for replicate in range(replicates):
        chosen = rng.choice(ids, len(ids), replace=True)
        sample = pd.concat([event[event["encounter_event_id"].eq(x)] for x in chosen], ignore_index=True)
        means = sample.groupby(["context", "radius_m"], observed=True)[
            ["edge_density", "largest_component_fraction", "isolated_fraction", "fine_split_probability"]
        ].mean()
        for radius in RADII:
            if not all((context, radius) in means.index for context in ["Purple present", "Matched Purple absent"]):
                continue
            present = means.loc[("Purple present", radius)]
            absent = means.loc[("Matched Purple absent", radius)]
            draws.append({
                "replicate": replicate,
                "radius_m": radius,
                **{f"present_{k}": v for k, v in present.items()},
                **{f"absent_{k}": v for k, v in absent.items()},
                **{f"difference_{k}": present[k] - absent[k] for k in present.index},
            })
    return pd.DataFrame(draws)


def plot(summary: pd.DataFrame, draws: pd.DataFrame, n_events: int, output: Path) -> None:
    metrics = [
        ("edge_density", "Pairwise proximity-edge density"),
        ("largest_component_fraction", "Fraction in largest proximity component"),
        ("fine_split_probability", "Fine-scale split probability"),
    ]
    contexts = ["Purple present", "Matched Purple absent"]
    colors = {"Purple present": "#7b2cbf", "Matched Purple absent": "#777777"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    x = np.arange(len(RADII))
    width = 0.36
    for ax, (metric, title) in zip(axes, metrics):
        for offset, context in zip([-width / 2, width / 2], contexts):
            values = summary[summary["context"].eq(context)].set_index("radius_m").reindex(RADII)[metric]
            prefix = "present" if context == "Purple present" else "absent"
            lows = draws.groupby("radius_m")[f"{prefix}_{metric}"].quantile(0.025).reindex(RADII)
            highs = draws.groupby("radius_m")[f"{prefix}_{metric}"].quantile(0.975).reindex(RADII)
            errors = np.vstack([values.to_numpy() - lows.to_numpy(), highs.to_numpy() - values.to_numpy()])
            ax.bar(x + offset, values, width, color=colors[context], label=context, yerr=errors, capsize=4)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{int(r)} m" for r in RADII])
        ax.set_title(title)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"Fine-scale within-Chartreuse proximity clustering\n{n_events} encounter events; event-bootstrap 95% CIs")
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matches = pd.read_csv(args.matches_file, parse_dates=["encounter_time", "reference_time"])
    panel = pd.read_parquet(args.panel_file, columns=["animal_id", "origin_group"])
    animals = sorted(panel.loc[panel["origin_group"].eq("Chartreuse"), "animal_id"].astype(str).unique())
    requested_hours = set(matches["encounter_time"]) | set(matches["reference_time"])
    start, end = min(requested_hours), max(requested_hours) + pd.Timedelta(hours=1)
    start_utc = start.tz_localize("UTC")
    end_utc = end.tz_localize("UTC")
    gps = pd.read_parquet(
        args.gps_file,
        columns=["animal_id", "timestamp", "location.long", "location.lat"],
        filters=[("animal_id", "in", animals), ("timestamp", ">=", start_utc), ("timestamp", "<", end_utc)],
    )
    gps["timestamp"] = pd.to_datetime(gps["timestamp"], utc=True).dt.tz_localize(None)
    gps["hour"] = gps["timestamp"].dt.floor("h")
    gps = gps[gps["hour"].isin(requested_hours)].rename(columns={"location.long": "longitude", "location.lat": "latitude"})
    gps["bin_2min"] = gps["timestamp"].dt.floor("2min")
    fixes = gps.groupby(["hour", "bin_2min", "animal_id"], observed=True, as_index=False).agg(
        longitude=("longitude", "median"), latitude=("latitude", "median")
    )
    hour_context = pd.concat([
        matches[["encounter_time", "encounter_event_id"]].rename(columns={"encounter_time": "hour"}).assign(context="Purple present"),
        matches[["reference_time", "encounter_event_id"]].rename(columns={"reference_time": "hour"}).assign(context="Matched Purple absent"),
    ], ignore_index=True)
    fixes = fixes.merge(hour_context, on="hour", how="inner", validate="many_to_one")

    rows = []
    for (hour, bin_time, context, event_id), frame in fixes.groupby(
        ["hour", "bin_2min", "context", "encounter_event_id"], observed=True
    ):
        if len(frame) < 2:
            continue
        for radius in RADII:
            rows.append({"hour": hour, "bin_2min": bin_time, "context": context, "encounter_event_id": event_id, "radius_m": radius, **bin_metrics(frame, radius)})
    rows = pd.DataFrame(rows)
    if rows.empty:
        raise ValueError("No matched fine-scale bins with at least two Chartreuse animals")
    summary = rows.groupby(["context", "radius_m"], observed=True).agg(
        encounter_events=("encounter_event_id", "nunique"), hours=("hour", "nunique"), bins_2min=("bin_2min", "size"),
        mean_animals=("n_animals", "mean"), total_edges=("proximity_edges", "sum"), total_candidate_pairs=("candidate_pairs", "sum"),
        edge_density=("edge_density", "mean"), largest_component_fraction=("largest_component_fraction", "mean"),
        isolated_fraction=("isolated_fraction", "mean"), fine_split_probability=("fine_split", "mean"),
    ).reset_index()
    draws = bootstrap_events(rows, args.bootstrap_replicates, args.seed)
    intervals = draws.groupby("radius_m").quantile([0.025, 0.5, 0.975]).reset_index()
    n_events = int(rows["encounter_event_id"].nunique())
    rows.to_parquet(args.output_dir / "chartreuse_finescale_bin_metrics.parquet", index=False)
    summary.to_csv(args.output_dir / "chartreuse_finescale_context_summary.csv", index=False)
    intervals.to_csv(args.output_dir / "chartreuse_finescale_bootstrap_intervals.csv", index=False)
    plot(summary, draws, n_events, args.output_dir / "chartreuse_finescale_clustering.png")
    metadata = {
        "gps_file": str(args.gps_file.resolve()), "matched_events_with_gps": n_events,
        "matched_hours_with_gps": {k: int(v) for k, v in rows.groupby("context")["hour"].nunique().items()},
        "radii_m": RADII, "fine_split_definition": "at least two connected components each containing at least two Chartreuse animals in a 2-minute bin",
        "edge_definition": "great-circle distance at or below the stated radius using one median location per animal per 2-minute bin",
        "bootstrap_cluster_unit": "Chartreuse-Purple encounter event",
        "caution": "GPS availability varies by bin; absence of an edge is evaluated only among co-observed animals.",
    }
    (args.output_dir / "chartreuse_finescale_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
