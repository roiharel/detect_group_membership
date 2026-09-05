from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


CELLS = Path("outputs/established_dispersal_centrality/event_week_centrality_cells.csv")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUTDIR = Path("outputs/established_dispersal_radius_samples")
FOCAL = "24AC02_3C4D"
METRICS = ["degree", "strength", "eigenvector", "harmonic_closeness", "betweenness"]


def main():
    cells = pd.read_csv(CELLS)
    events = pd.read_csv(EVENTS, usecols=["event_id", "animal_id"])
    cells = cells.merge(events, on="event_id", how="left")
    cells = cells[(cells.animal_id.eq(FOCAL)) & cells.established_outcome.eq("Permanent dispersal")].copy()
    cells.to_csv(OUTDIR / "permanent_high_support_centrality_event_weeks.csv", index=False)

    fig, axes = plt.subplots(3, 2, figsize=(13, 13), sharex=True, sharey=True)
    axes = axes.ravel()
    labels = {"degree": "Degree", "strength": "Weighted strength", "eigenvector": "Eigenvector",
              "harmonic_closeness": "Harmonic closeness", "betweenness": "Betweenness"}
    for ax, metric in zip(axes, METRICS):
        z = cells[cells.metric.eq(metric)]
        for _, event in z.groupby("event_id", observed=True):
            ax.plot(event.week, event.focal_percentile, color="#7b3294", alpha=.20, linewidth=1)
        mean = z.groupby("week").agg(percentile=("focal_percentile", "mean"), segments=("event_id", "nunique"))
        ax.plot(mean.index, mean.percentile, color="#7b3294", marker="o", linewidth=2.7)
        for week, row in mean.iterrows():
            ax.annotate(f"n={int(row.segments)}", (week, row.percentile), xytext=(0, 8),
                        textcoords="offset points", ha="center", fontsize=8)
        ax.axhline(.5, color="black", linestyle="--", linewidth=1)
        ax.set(title=labels[metric], xlabel="Weeks since recipient-segment entry",
               ylabel="Percentile among recipient members", ylim=(-.03, 1.03))
        ax.set_xticks(range(6)); ax.grid(alpha=.2)
    axes[-1].axis("off")
    fig.suptitle(f"5 m network centrality over time: permanent disperser {FOCAL}\n"
                 "Thin lines: recipient segments; thick line: event-equal mean; dashed line: group median")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(OUTDIR / "permanent_high_support_5m_centrality_over_time.png", dpi=220)
    plt.close(fig)
    print(cells.groupby(["metric", "week"]).agg(segments=("event_id", "nunique"), mean=("focal_percentile", "mean")).to_string())


if __name__ == "__main__":
    main()
