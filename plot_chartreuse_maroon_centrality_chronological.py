from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd


CELLS = Path("outputs/established_dispersal_radius_samples/chartreuse_to_maroon_returned_daily_centrality_cells.csv")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUT = Path("outputs/established_dispersal_radius_samples/chartreuse_to_maroon_5m_centrality_chronological.png")
METRICS = ["degree", "strength", "eigenvector", "harmonic_closeness", "betweenness"]


def main():
    cells = pd.read_csv(CELLS)
    events = pd.read_csv(EVENTS, parse_dates=["start_time"])[["event_id", "start_time"]]
    cells = cells.merge(events, on="event_id", how="left")
    cells["date"] = cells.start_time.dt.normalize() + pd.to_timedelta(cells.day, unit="D")
    labels = {"degree": "Degree", "strength": "Weighted strength", "eigenvector": "Eigenvector",
              "harmonic_closeness": "Harmonic closeness", "betweenness": "Betweenness"}
    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(3, 2, figsize=(15, 13), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, metric in zip(axes, METRICS):
        z = cells[cells.metric.eq(metric)]
        for i, (_, event) in enumerate(z.groupby("event_id", observed=True)):
            event = event.sort_values("date")
            segment = int(event.segment.iloc[0])
            ax.plot(event.date, event.focal_percentile, color=colors[i], marker="o", markersize=2.5,
                    linewidth=1.6, label=f"Segment {segment}")
            ax.axvline(event.date.min(), color=colors[i], alpha=.18, linewidth=1)
        ax.axhline(.5, color="black", linestyle="--", linewidth=1)
        ax.set(title=labels[metric], xlabel="Calendar date", ylabel="Percentile among Maroon members",
               ylim=(-.03, 1.03))
        ax.grid(alpha=.2)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    axes[0].legend(frameon=False, ncol=3, fontsize=8)
    axes[-1].axis("off")
    fig.suptitle("24AC15_3F4G: Chartreuse -> Maroon 5 m centrality in chronological order\n"
                 "Returned-departure segments; verified 2-minute sampling bouts")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(OUT, dpi=220)
    plt.close(fig)
    print(events[events.event_id.isin(cells.event_id.unique())].sort_values("start_time").to_string(index=False))


if __name__ == "__main__":
    main()
