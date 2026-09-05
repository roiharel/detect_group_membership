from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BINS = Path("outputs/disperser_finescale_integration/disperser_2min_contact_rows.parquet")
EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/canonical_disperser_events.csv")
OUT = Path("outputs/disperser_integration_over_time")
RADII = (2.0, 5.0, 10.0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Describe recipient integration through dispersal-event time.")
    p.add_argument("--bins-file", type=Path, default=BINS)
    p.add_argument("--events-file", type=Path, default=EVENTS)
    p.add_argument("--output-dir", type=Path, default=OUT)
    p.add_argument("--progress-bins", type=int, default=10)
    p.add_argument("--bootstrap-replicates", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20260808)
    return p.parse_args()


def event_equal_summary(cells: pd.DataFrame) -> pd.DataFrame:
    return cells.groupby(["radius_m", "progress_bin"], observed=True).agg(
        events=("event_id", "nunique"),
        mean_recipient_contact_probability=("recipient_contact_probability", "mean"),
        mean_origin_contact_probability=("origin_contact_probability", "mean"),
        mean_nearest_recipient_m=("median_nearest_recipient_m", "mean"),
    ).reset_index()


def bootstrap(cells: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    ids = cells["event_id"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed)
    draws = []
    for replicate in range(replicates):
        chosen = rng.choice(ids, len(ids), replace=True)
        sample = pd.concat([cells[cells["event_id"].eq(x)] for x in chosen], ignore_index=True)
        means = sample.groupby(["radius_m", "progress_bin"], observed=True)[
            ["recipient_contact_probability", "origin_contact_probability"]
        ].mean().reset_index()
        means["replicate"] = replicate
        draws.append(means)
    return pd.concat(draws, ignore_index=True)


def slope_table(cells: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    ids = cells["event_id"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed + 19)

    def slopes(frame: pd.DataFrame) -> dict[float, float]:
        out = {}
        for radius, group in frame.groupby("radius_m", observed=True):
            x = group["progress_midpoint"].to_numpy(float)
            y = group["recipient_contact_probability"].to_numpy(float)
            out[float(radius)] = float(np.polyfit(x, y, 1)[0]) if len(np.unique(x)) > 1 else np.nan
        return out

    observed = slopes(cells)
    draw_rows = []
    for replicate in range(replicates):
        chosen = rng.choice(ids, len(ids), replace=True)
        sample = pd.concat([cells[cells["event_id"].eq(x)] for x in chosen], ignore_index=True)
        for radius, value in slopes(sample).items():
            draw_rows.append({"replicate": replicate, "radius_m": radius, "slope": value})
    draws = pd.DataFrame(draw_rows)
    rows = []
    for radius in RADII:
        values = draws.loc[draws["radius_m"].eq(radius), "slope"].dropna()
        rows.append({
            "radius_m": radius,
            "observed_change_start_to_end": observed.get(radius, np.nan),
            "ci_95_low": values.quantile(0.025),
            "ci_95_high": values.quantile(0.975),
        })
    return pd.DataFrame(rows)


def plot(summary: pd.DataFrame, intervals: pd.DataFrame, individual: pd.DataFrame, output: Path) -> None:
    colors = {2.0: "#4c78a8", 5.0: "#7b3294", 10.0: "#e45756"}
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for radius in RADII:
        frame = summary[summary["radius_m"].eq(radius)].sort_values("progress_bin")
        ci = intervals[intervals["radius_m"].eq(radius)].set_index("progress_bin").reindex(frame["progress_bin"])
        x = frame["progress_midpoint"] * 100
        y = frame["mean_recipient_contact_probability"]
        axes[0].plot(x, y, marker="o", color=colors[radius], label=f"{int(radius)} m")
        axes[0].fill_between(x, ci["ci_95_low"], ci["ci_95_high"], color=colors[radius], alpha=0.15)
    axes[0].set_xlabel("Progress through dispersal event (%)")
    axes[0].set_ylabel("Recipient contact probability\n(equal event weighting)")
    axes[0].set_title("Pooled recipient integration")
    axes[0].legend(frameon=False)
    five = individual[individual["radius_m"].eq(5)].copy()
    for (animal, target), frame in five.groupby(["animal_id", "recipient_group"], observed=True):
        frame = frame.sort_values("progress_bin")
        axes[1].plot(frame["progress_midpoint"] * 100, frame["recipient_contact_probability"], marker="o", alpha=0.8,
                     label=f"{animal} → {target}")
    axes[1].set_xlabel("Progress through dispersal event (%)")
    axes[1].set_ylabel("Recipient contact probability at 5 m")
    axes[1].set_title("Individual disperser trajectories")
    axes[1].legend(fontsize=7, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.suptitle("How fine-scale integration changes through dispersal events\n95% CIs resample complete events")
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    bins = pd.read_parquet(args.bins_file)
    events = pd.read_csv(args.events_file, parse_dates=["start_time", "end_time"])[["event_id", "start_time", "end_time"]]
    bins["bin_2min"] = pd.to_datetime(bins["bin_2min"])
    bins = bins.merge(events, on="event_id", how="left", validate="many_to_one")
    duration = (bins["end_time"] - bins["start_time"]).dt.total_seconds().clip(lower=3600)
    bins["event_progress"] = ((bins["bin_2min"] - bins["start_time"]).dt.total_seconds() / duration).clip(0, 0.999999)
    bins["progress_bin"] = np.floor(bins["event_progress"] * args.progress_bins).astype(int)
    bins["progress_midpoint"] = (bins["progress_bin"] + 0.5) / args.progress_bins
    cells = bins.groupby(
        ["event_id", "animal_id", "origin_group", "recipient_group", "radius_m", "progress_bin", "progress_midpoint"], observed=True
    ).agg(
        bins=("bin_2min", "size"), recipient_opportunities=("recipient_opportunity", "sum"),
        origin_opportunities=("origin_opportunity", "sum"), recipient_contacts=("recipient_contact", "sum"),
        origin_contacts=("origin_contact", "sum"), median_nearest_recipient_m=("nearest_recipient_m", "median"),
    ).reset_index()
    cells["recipient_contact_probability"] = cells["recipient_contacts"] / cells["recipient_opportunities"].replace(0, np.nan)
    cells["origin_contact_probability"] = cells["origin_contacts"] / cells["origin_opportunities"].replace(0, np.nan)
    valid = cells.dropna(subset=["recipient_contact_probability"])
    summary = event_equal_summary(valid)
    summary["progress_midpoint"] = (summary["progress_bin"] + 0.5) / args.progress_bins
    draws = bootstrap(valid, args.bootstrap_replicates, args.seed)
    intervals = draws.groupby(["radius_m", "progress_bin"])["recipient_contact_probability"].quantile([0.025, 0.5, 0.975]).unstack().reset_index().rename(columns={0.025: "ci_95_low", 0.5: "bootstrap_median", 0.975: "ci_95_high"})
    slopes = slope_table(valid, args.bootstrap_replicates, args.seed)
    individual = valid.groupby(["animal_id", "recipient_group", "radius_m", "progress_bin", "progress_midpoint"], observed=True).agg(
        events=("event_id", "nunique"), recipient_contacts=("recipient_contacts", "sum"), recipient_opportunities=("recipient_opportunities", "sum")
    ).reset_index()
    individual["recipient_contact_probability"] = individual["recipient_contacts"] / individual["recipient_opportunities"]
    cells.to_csv(args.output_dir / "disperser_event_progress_cells.csv", index=False)
    summary.to_csv(args.output_dir / "disperser_integration_progress_summary.csv", index=False)
    intervals.to_csv(args.output_dir / "disperser_integration_progress_intervals.csv", index=False)
    slopes.to_csv(args.output_dir / "disperser_integration_progress_slopes.csv", index=False)
    individual.to_csv(args.output_dir / "disperser_individual_progress_summary.csv", index=False)
    plot(summary, intervals, individual, args.output_dir / "disperser_integration_over_time.png")
    metadata = {"events": int(valid.event_id.nunique()), "animals": int(valid.animal_id.nunique()), "progress_bins": args.progress_bins,
                "weighting": "each event-progress cell receives equal weight in pooled curves", "bootstrap_cluster_unit": "dispersal event",
                "slope_scale": "change in contact probability from 0 to 100 percent event progress under a linear summary",
                "caution": "normalized event time compares trajectories of different duration; it is not calendar time"}
    (args.output_dir / "disperser_integration_over_time_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2)); print(slopes.to_string(index=False))


if __name__ == "__main__": main()
