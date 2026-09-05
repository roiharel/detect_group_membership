from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BINS = Path("outputs/established_dispersal_radius_samples/recipient_contact_2min_bins_1m_2m_5m.parquet")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUTDIR = Path("outputs/established_dispersal_radius_samples")
MAX_DAY = 41


def main():
    bins = pd.read_parquet(BINS)
    bins = bins[bins.established_outcome.eq("Permanent dispersal")].copy()
    support = bins.groupby("animal_id").contact_5m.sum().rename("positive_5m_bins")
    keep = support[support > 100].index.tolist()
    bins = bins[bins.animal_id.isin(keep)]
    events = pd.read_csv(EVENTS, parse_dates=["start_time"])[["event_id", "start_time"]]
    bins = bins.merge(events, on="event_id", how="left")
    bins["day"] = np.floor((pd.to_datetime(bins.bin_2min) - bins.start_time).dt.total_seconds() / 86400).astype(int)
    bins = bins[bins.day.between(0, MAX_DAY)]

    event_day = (bins.groupby(["animal_id", "event_id", "day"], observed=True)
                 .agg(observed_bins=("bin_2min", "size"), contact_bins=("contact_5m", "sum"))
                 .reset_index())
    event_day["observed_minutes"] = 2 * event_day.observed_bins
    event_day["colocation_minutes"] = 2 * event_day.contact_bins
    event_day["colocation_proportion"] = event_day.contact_bins / event_day.observed_bins
    event_day.to_csv(OUTDIR / "permanent_gt100_5m_colocation_event_days.csv", index=False)
    daily = (event_day.groupby(["animal_id", "day"], observed=True)
             .agg(mean_minutes=("colocation_minutes", "mean"),
                  mean_proportion=("colocation_proportion", "mean"),
                  segments=("event_id", "nunique"), observed_minutes=("observed_minutes", "sum"))
             .reset_index())
    daily.to_csv(OUTDIR / "permanent_gt100_5m_colocation_daily_summary.csv", index=False)

    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
    for i, animal in enumerate(keep):
        z = daily[daily.animal_id.eq(animal)].sort_values("day")
        label = f"{animal} (5 m n={int(support.loc[animal])})"
        axes[0].plot(z.day, z.mean_minutes, marker="o", markersize=3, linewidth=1.8, color=colors[i], label=label)
        axes[1].plot(z.day, z.mean_proportion, marker="o", markersize=3, linewidth=1.8, color=colors[i])
    axes[0].set_ylabel("Mean sampled minutes within 5 m")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_ylabel("Proportion of observed 2-minute bins")
    axes[1].set_xlabel("Days since recipient-segment entry")
    for ax in axes:
        ax.set_xlim(-1, MAX_DAY + 1); ax.set_ylim(bottom=0); ax.grid(alpha=.2)
    fig.suptitle("Daily time within 5 m of recipient-group members\n"
                 "Permanent dispersers with >100 positive 5 m bins; event-equal daily means")
    fig.tight_layout(rect=[0, 0, 1, .94])
    fig.savefig(OUTDIR / "permanent_gt100_5m_colocation_time_by_day.png", dpi=220)
    plt.close(fig)
    print(event_day.groupby("animal_id").agg(segments=("event_id", "nunique"), days=("day", "nunique"),
          total_minutes=("colocation_minutes", "sum")).to_string())


if __name__ == "__main__":
    main()
