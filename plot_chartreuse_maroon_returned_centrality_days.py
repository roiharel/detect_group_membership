from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_disperser_network_centrality import centralities


POSITIONS = Path("outputs/established_dispersal_centrality/event_member_positions_established_segments.parquet")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUTDIR = Path("outputs/established_dispersal_radius_samples")
MAX_GAP_MIN = 3
MIN_BOUT_BINS = 15
MAX_DAY = 41
METRICS = ["degree", "strength", "eigenvector", "harmonic_closeness", "betweenness"]


def main():
    events = pd.read_csv(EVENTS, parse_dates=["start_time"])
    route = events[(events.origin_group.eq("Chartreuse")) & events.recipient_group.eq("Maroon") &
                   events.established_outcome.eq("Returned to origin")]
    positions = pd.read_parquet(POSITIONS)
    positions = positions[positions.event_id.isin(route.event_id)].merge(
        route[["event_id", "start_time"]], on="event_id", how="left", suffixes=("", "_event"))
    focal = positions[positions.animal_id.astype(str).eq(positions.animal_id_focal.astype(str))]
    focal = focal[["event_id", "bin_2min"]].drop_duplicates().sort_values(["event_id", "bin_2min"])
    focal["gap_min"] = focal.groupby("event_id").bin_2min.diff().dt.total_seconds().div(60)
    focal["bout"] = focal.groupby("event_id").gap_min.transform(lambda x: x.gt(MAX_GAP_MIN).cumsum())
    sizes = focal.groupby(["event_id", "bout"], observed=True).size().rename("bout_bins").reset_index()
    valid = focal.merge(sizes, on=["event_id", "bout"])
    valid = valid[valid.bout_bins.ge(MIN_BOUT_BINS)][["event_id", "bin_2min"]]
    positions = positions.merge(valid, on=["event_id", "bin_2min"], how="inner")
    positions["day"] = np.floor((positions.bin_2min - positions.start_time_event).dt.total_seconds() / 86400).astype(int)
    positions = positions[positions.day.between(0, MAX_DAY)].copy()
    positions["week"] = positions.day
    lookup = route[["event_id", "animal_id", "recipient_group"]]
    cells = centralities(positions, lookup).rename(columns={"week": "day"})
    cells["segment"] = cells.event_id.str.extract(r"__(\d+)$").astype(int)
    cells.to_csv(OUTDIR / "chartreuse_to_maroon_returned_daily_centrality_cells.csv", index=False)
    mean = (cells.groupby(["metric", "day"], observed=True)
            .agg(mean_percentile=("focal_percentile", "mean"), segments=("event_id", "nunique")).reset_index())

    labels = {"degree": "Degree", "strength": "Weighted strength", "eigenvector": "Eigenvector",
              "harmonic_closeness": "Harmonic closeness", "betweenness": "Betweenness"}
    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(3, 2, figsize=(13, 13), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, metric in zip(axes, METRICS):
        z = cells[cells.metric.eq(metric)]
        for i, (_, event) in enumerate(z.groupby("event_id", observed=True)):
            event = event.sort_values("day")
            ax.plot(event.day, event.focal_percentile, color=colors[i], alpha=.35, linewidth=1.2,
                    label=f"Segment {int(event.segment.iloc[0])}")
        m = mean[mean.metric.eq(metric)].sort_values("day")
        ax.plot(m.day, m.mean_percentile, color="black", marker="o", markersize=3, linewidth=2.6,
                label="Event-equal mean")
        ax.axhline(.5, color="black", linestyle="--", linewidth=1)
        ax.set(title=labels[metric], xlabel="Days since Maroon-segment entry",
               ylabel="Percentile among Maroon members", ylim=(-.03, 1.03), xlim=(-1, MAX_DAY + 1))
        ax.grid(alpha=.2)
    handles, leglabels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, leglabels, frameon=False, ncol=2, fontsize=8)
    axes[-1].axis("off")
    fig.suptitle("24AC15_3F4G: Chartreuse -> Maroon daily 5 m centrality\n"
                 "Returned departures during verified 2-minute sampling bouts")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(OUTDIR / "chartreuse_to_maroon_returned_5m_centrality_by_day.png", dpi=220)
    plt.close(fig)
    print(cells.groupby("event_id").agg(days=("day", "nunique"), focal_bins=("focal_bins", "sum")).to_string())


if __name__ == "__main__":
    main()
