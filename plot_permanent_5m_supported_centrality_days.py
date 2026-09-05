from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_disperser_network_centrality import centralities


BINS = Path("outputs/established_dispersal_radius_samples/recipient_contact_2min_bins_1m_2m_5m.parquet")
POSITIONS = Path("outputs/established_dispersal_centrality/event_member_positions_established_segments.parquet")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUTDIR = Path("outputs/established_dispersal_radius_samples")
METRICS = ["degree", "strength", "eigenvector", "harmonic_closeness", "betweenness"]
MAX_DAY = 41


def main():
    bins = pd.read_parquet(BINS)
    perm = bins[bins.established_outcome.eq("Permanent dispersal")]
    support = perm.groupby("animal_id").contact_5m.sum().rename("positive_5m_bins")
    keep = support[support > 100].index.tolist()

    events = pd.read_csv(EVENTS, parse_dates=["start_time"])
    selected = events[(events.established_outcome.eq("Permanent dispersal")) & events.animal_id.isin(keep)]
    lookup = selected[["event_id", "animal_id", "recipient_group"]]
    positions = pd.read_parquet(POSITIONS)
    positions = positions[positions.event_id.isin(selected.event_id)].merge(
        selected[["event_id", "start_time"]], on="event_id", how="left", suffixes=("", "_event"))
    positions["day"] = np.floor((pd.to_datetime(positions.bin_2min) - positions.start_time_event).dt.total_seconds() / 86400).astype(int)
    positions = positions[positions.day.between(0, MAX_DAY)].copy()
    # Reuse the validated network calculation with day as its grouping index.
    positions["week"] = positions.day
    cells = centralities(positions, lookup).rename(columns={"week": "day"})
    cells = cells.merge(selected[["event_id", "animal_id"]], on="event_id", how="left")
    cells.to_csv(OUTDIR / "permanent_gt100_5m_daily_centrality_cells.csv", index=False)
    summary = (cells.groupby(["animal_id", "metric", "day"], observed=True)
               .agg(mean_percentile=("focal_percentile", "mean"), segments=("event_id", "nunique"))
               .reset_index().merge(support.reset_index(), on="animal_id", how="left"))
    summary.to_csv(OUTDIR / "permanent_gt100_5m_daily_centrality_summary.csv", index=False)

    labels = {"degree": "Degree", "strength": "Weighted strength", "eigenvector": "Eigenvector",
              "harmonic_closeness": "Harmonic closeness", "betweenness": "Betweenness"}
    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(3, 2, figsize=(13, 13), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, metric in zip(axes, METRICS):
        z = summary[summary.metric.eq(metric)]
        for i, animal in enumerate(keep):
            q = z[z.animal_id.eq(animal)].sort_values("day")
            if q.empty:
                continue
            ax.plot(q.day, q.mean_percentile, marker="o", markersize=3, linewidth=1.8,
                    color=colors[i], label=f"{animal} (5 m n={int(support.loc[animal])})")
        ax.axhline(.5, color="black", linestyle="--", linewidth=1)
        ax.set(title=labels[metric], xlabel="Days since recipient-segment entry",
               ylabel="Percentile among recipient members", ylim=(-.03, 1.03), xlim=(-1, MAX_DAY + 1))
        ax.grid(alpha=.2)
    axes[0].legend(frameon=False, fontsize=8)
    axes[-1].axis("off")
    fig.suptitle("Daily 5 m network centrality: permanent dispersers with >100 positive 5 m bins\n"
                 "Daily networks require at least 4 nodes and 30 focal 2-minute bins")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(OUTDIR / "permanent_gt100_5m_centrality_by_day.png", dpi=220)
    plt.close(fig)
    print(support[support > 100].sort_values(ascending=False).to_string())
    print(cells.groupby("animal_id").agg(segments=("event_id", "nunique"), days=("day", "nunique")).to_string())


if __name__ == "__main__":
    main()
