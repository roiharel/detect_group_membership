from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CONTACTS = Path("outputs/established_dispersal_radius_samples/recipient_contact_2min_bins_1m_2m_5m.parquet")
POSITIONS = Path("outputs/established_dispersal_centrality/event_member_positions_established_segments.parquet")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUTDIR = Path("outputs/established_dispersal_radius_samples")
MAX_GAP_MIN = 3
MIN_BOUT_BINS = 15
MAX_DAY = 41


def main():
    events = pd.read_csv(EVENTS, parse_dates=["start_time"])
    route = events[(events.origin_group.eq("Chartreuse")) & events.recipient_group.eq("Maroon") &
                   events.established_outcome.eq("Returned to origin")]
    positions = pd.read_parquet(POSITIONS, columns=["event_id", "bin_2min", "animal_id", "animal_id_focal"])
    focal = positions[(positions.event_id.isin(route.event_id)) &
                      positions.animal_id.astype(str).eq(positions.animal_id_focal.astype(str))]
    focal = focal[["event_id", "bin_2min"]].drop_duplicates().sort_values(["event_id", "bin_2min"])
    focal["gap_min"] = focal.groupby("event_id").bin_2min.diff().dt.total_seconds().div(60)
    focal["bout"] = focal.groupby("event_id").gap_min.transform(lambda x: x.gt(MAX_GAP_MIN).cumsum())
    sizes = focal.groupby(["event_id", "bout"], observed=True).size().rename("bout_bins").reset_index()
    valid = focal.merge(sizes, on=["event_id", "bout"])
    valid = valid[valid.bout_bins.ge(MIN_BOUT_BINS)][["event_id", "bin_2min"]]
    contacts = pd.read_parquet(CONTACTS)
    contacts = contacts[contacts.event_id.isin(route.event_id)].merge(valid, on=["event_id", "bin_2min"])
    contacts = contacts.merge(route[["event_id", "start_time"]], on="event_id", how="left")
    contacts["day"] = np.floor((contacts.bin_2min - contacts.start_time).dt.total_seconds() / 86400).astype(int)
    contacts = contacts[contacts.day.between(0, MAX_DAY)]
    daily = (contacts.groupby(["event_id", "day"], observed=True)
             .agg(eligible_bins=("bin_2min", "size"), contact_bins=("contact_5m", "sum")).reset_index())
    daily["minutes_5m"] = 2 * daily.contact_bins
    daily["proportion_5m"] = daily.contact_bins / daily.eligible_bins
    daily["segment"] = daily.event_id.str.extract(r"__(\d+)$").astype(int)
    daily.to_csv(OUTDIR / "chartreuse_to_maroon_returned_event_days_cadence_filtered.csv", index=False)
    mean = daily.groupby("day").agg(minutes_5m=("minutes_5m", "mean"), proportion_5m=("proportion_5m", "mean"),
                                     segments=("event_id", "nunique")).reset_index()

    fig, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
    colors = plt.get_cmap("tab10").colors
    for i, (event_id, z) in enumerate(daily.groupby("event_id", observed=True)):
        z = z.sort_values("day"); label = f"Segment {int(z.segment.iloc[0])}"
        axes[0].plot(z.day, z.minutes_5m, color=colors[i], alpha=.35, linewidth=1.2, label=label)
        axes[1].plot(z.day, z.proportion_5m, color=colors[i], alpha=.35, linewidth=1.2)
    axes[0].plot(mean.day, mean.minutes_5m, color="black", marker="o", markersize=3, linewidth=2.6, label="Event-equal mean")
    axes[1].plot(mean.day, mean.proportion_5m, color="black", marker="o", markersize=3, linewidth=2.6)
    axes[0].set_ylabel("Sampled minutes within 5 m")
    axes[0].legend(frameon=False, ncol=3, fontsize=8)
    axes[1].set_ylabel("Proportion of eligible 2-minute bins")
    axes[1].set_xlabel("Days since Maroon-segment entry")
    for ax in axes:
        ax.set_xlim(-1, MAX_DAY + 1); ax.set_ylim(bottom=0); ax.grid(alpha=.2)
    fig.suptitle("24AC15_3F4G: Chartreuse -> Maroon returned departures\n"
                 "Five segments during verified 2-minute sampling bouts")
    fig.tight_layout(rect=[0, 0, 1, .94])
    fig.savefig(OUTDIR / "chartreuse_to_maroon_returned_5m_by_day.png", dpi=220)
    plt.close(fig)
    print(daily.groupby("event_id").agg(days=("day", "nunique"), bins=("eligible_bins", "sum"),
          contacts=("contact_bins", "sum"), minutes=("minutes_5m", "sum")).to_string())


if __name__ == "__main__":
    main()
