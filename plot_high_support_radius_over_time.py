from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BINS = Path("outputs/established_dispersal_radius_samples/recipient_contact_2min_bins_1m_2m_5m.parquet")
OUTDIR = Path("outputs/established_dispersal_radius_samples")
RADII = (1, 2, 5)


def main():
    bins = pd.read_parquet(BINS)
    totals = bins.groupby("animal_id").contact_1m.sum()
    keep = totals[totals > 50].index
    bins = bins[bins.animal_id.isin(keep)].copy()

    event_week = (bins.groupby(["animal_id", "event_id", "established_outcome", "week"], observed=True)
                  .agg(bins=("bin_2min", "size"),
                       **{f"p_{r}m": (f"contact_{r}m", "mean") for r in RADII})
                  .reset_index())
    event_week.to_csv(OUTDIR / "high_support_event_week_integration_1m_2m_5m.csv", index=False)

    long = event_week.melt(
        id_vars=["animal_id", "event_id", "established_outcome", "week", "bins"],
        value_vars=[f"p_{r}m" for r in RADII], var_name="radius", value_name="integration")
    long["radius_m"] = long.radius.str.extract(r"(\d+)").astype(int)
    mean = (long.groupby(["animal_id", "week", "radius_m"], observed=True)
            .agg(mean_integration=("integration", "mean"), segments=("event_id", "nunique"), bins=("bins", "sum"))
            .reset_index())
    mean.to_csv(OUTDIR / "high_support_weekly_integration_1m_2m_5m.csv", index=False)

    animals = list(keep)
    fig, axes = plt.subplots(3, 2, figsize=(13, 13), sharex=True, sharey=True)
    axes = axes.ravel()
    colors = {1: "#0072B2", 2: "#E69F00", 5: "#009E73"}
    for ax, animal in zip(axes, animals):
        z = long[long.animal_id.eq(animal)]
        for radius in RADII:
            q = z[z.radius_m.eq(radius)]
            for _, event in q.groupby("event_id", observed=True):
                ax.plot(event.week, event.integration, color=colors[radius], alpha=.16, linewidth=1)
            m = mean[(mean.animal_id.eq(animal)) & mean.radius_m.eq(radius)].sort_values("week")
            ax.plot(m.week, m.mean_integration, color=colors[radius], marker="o", linewidth=2.5,
                    label=f"{radius} m")
        n_events = z.event_id.nunique()
        outcomes = ", ".join(sorted(z.established_outcome.str.replace(" to origin", "").unique()))
        ax.set_title(f"{animal} — {n_events} segments ({outcomes})")
        ax.set_ylim(-.01, .40)
        ax.set_xticks(range(6))
        ax.grid(alpha=.2)
        ax.set_xlabel("Weeks since recipient-segment entry")
        ax.set_ylabel("Contact proportion")
    axes[0].legend(frameon=False, ncol=3)
    for ax in axes[len(animals):]:
        ax.axis("off")
    fig.suptitle("Integration over time for individuals with >50 positive 1 m bins\n"
                 "Thin lines: recipient-association segments; thick lines: event-equal weekly means")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(OUTDIR / "high_support_radius_integration_over_time.png", dpi=220)
    plt.close(fig)
    print(mean.to_string(index=False))


if __name__ == "__main__":
    main()
