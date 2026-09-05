from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BINS = Path("outputs/established_dispersal_radius_samples/recipient_contact_2min_bins_1m_2m_5m.parquet")
CELLS = Path("outputs/established_dispersal_centrality/event_week_centrality_cells.csv")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUTDIR = Path("outputs/established_dispersal_radius_samples")
METRICS = ["degree", "strength", "eigenvector", "harmonic_closeness", "betweenness"]


def main():
    bins = pd.read_parquet(BINS)
    perm = bins[bins.established_outcome.eq("Permanent dispersal")]
    support = perm.groupby("animal_id").contact_5m.sum().rename("positive_5m_bins")
    keep = support[support > 100].index.tolist()

    cells = pd.read_csv(CELLS)
    events = pd.read_csv(EVENTS, usecols=["event_id", "animal_id"])
    cells = cells.merge(events, on="event_id", how="left")
    cells = cells[(cells.established_outcome.eq("Permanent dispersal")) & cells.animal_id.isin(keep)].copy()
    summary = (cells.groupby(["animal_id", "metric", "week"], observed=True)
               .agg(mean_percentile=("focal_percentile", "mean"), segments=("event_id", "nunique"))
               .reset_index())
    summary = summary.merge(support.reset_index(), on="animal_id", how="left")
    summary.to_csv(OUTDIR / "permanent_gt100_5m_centrality_weekly.csv", index=False)

    labels = {"degree": "Degree", "strength": "Weighted strength", "eigenvector": "Eigenvector",
              "harmonic_closeness": "Harmonic closeness", "betweenness": "Betweenness"}
    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(3, 2, figsize=(13, 13), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, metric in zip(axes, METRICS):
        z = summary[summary.metric.eq(metric)]
        for i, animal in enumerate(keep):
            q = z[z.animal_id.eq(animal)].sort_values("week")
            if q.empty:
                continue
            n5 = int(support.loc[animal])
            ax.plot(q.week, q.mean_percentile, marker="o", linewidth=2.2, color=colors[i],
                    label=f"{animal} (5 m n={n5})")
        ax.axhline(.5, color="black", linestyle="--", linewidth=1)
        ax.set(title=labels[metric], xlabel="Weeks since recipient-segment entry",
               ylabel="Percentile among recipient members", ylim=(-.03, 1.03))
        ax.set_xticks(range(6)); ax.grid(alpha=.2)
    axes[0].legend(frameon=False, fontsize=8)
    axes[-1].axis("off")
    fig.suptitle("5 m network centrality: permanent dispersers with >100 positive 5 m bins\n"
                 "Lines are event-equal weekly means; dashed line is recipient-group median")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(OUTDIR / "permanent_gt100_5m_centrality_over_time.png", dpi=220)
    plt.close(fig)
    print(support[support > 100].sort_values(ascending=False).to_string())
    print(cells.groupby("animal_id").event_id.nunique().rename("centrality_segments").to_string())


if __name__ == "__main__":
    main()
