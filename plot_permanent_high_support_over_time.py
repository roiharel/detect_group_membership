from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BINS = Path("outputs/established_dispersal_radius_samples/recipient_contact_2min_bins_1m_2m_5m.parquet")
OUT = Path("outputs/established_dispersal_radius_samples/permanent_high_support_radius_over_time.png")
RADII = (1, 2, 5)


def main():
    bins = pd.read_parquet(BINS)
    bins = bins[bins.established_outcome.eq("Permanent dispersal")].copy()
    support = bins.groupby("animal_id").contact_1m.sum()
    keep = support[support > 50].index
    bins = bins[bins.animal_id.isin(keep)]
    weekly = (bins.groupby(["animal_id", "event_id", "week"], observed=True)
              .agg(bins=("bin_2min", "size"),
                   **{f"p{r}": (f"contact_{r}m", "mean") for r in RADII}).reset_index())

    fig, ax = plt.subplots(figsize=(9, 6))
    colors = {1: "#0072B2", 2: "#E69F00", 5: "#009E73"}
    for radius in RADII:
        for _, event in weekly.groupby("event_id", observed=True):
            ax.plot(event.week, event[f"p{radius}"], color=colors[radius], alpha=.18, linewidth=1)
        mean = weekly.groupby("week")[f"p{radius}"].mean()
        ax.plot(mean.index, mean.values, color=colors[radius], marker="o", linewidth=2.7, label=f"{radius} m")
    focal = keep[0]
    ax.set(title=f"Permanent dispersal: {focal} ({weekly.event_id.nunique()} recipient segments)",
           xlabel="Weeks since recipient-segment entry", ylabel="Contact proportion", ylim=(-.01, .40))
    ax.set_xticks(range(6)); ax.grid(alpha=.22); ax.legend(frameon=False, ncol=3)
    fig.suptitle("Only one permanent disperser has >50 positive 1 m bins\nThin lines: segments; thick lines: event-equal weekly means")
    fig.tight_layout(rect=[0, 0, 1, .91]); fig.savefig(OUT, dpi=220); plt.close(fig)
    print(support.sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
