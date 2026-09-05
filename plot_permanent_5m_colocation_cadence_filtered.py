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
    perm_events = events[events.established_outcome.eq("Permanent dispersal")]
    positions = pd.read_parquet(POSITIONS, columns=["event_id", "bin_2min", "animal_id", "animal_id_focal"])
    focal = positions[(positions.event_id.isin(perm_events.event_id)) &
                      (positions.animal_id.astype(str).eq(positions.animal_id_focal.astype(str)))]
    focal = focal[["event_id", "bin_2min"]].drop_duplicates().sort_values(["event_id", "bin_2min"])
    focal["gap_min"] = focal.groupby("event_id").bin_2min.diff().dt.total_seconds().div(60)
    focal["bout"] = focal.groupby("event_id").gap_min.transform(lambda x: x.gt(MAX_GAP_MIN).cumsum())
    sizes = focal.groupby(["event_id", "bout"], observed=True).size().rename("bout_bins").reset_index()
    focal = focal.merge(sizes, on=["event_id", "bout"], how="left")
    valid = focal[focal.bout_bins.ge(MIN_BOUT_BINS)][["event_id", "bin_2min", "bout", "bout_bins"]]

    contacts = pd.read_parquet(CONTACTS)
    contacts = contacts[contacts.established_outcome.eq("Permanent dispersal")]
    before = len(contacts)
    contacts = contacts.merge(valid, on=["event_id", "bin_2min"], how="inner")
    contacts = contacts.merge(
        perm_events[["event_id", "start_time", "origin_group", "recipient_group"]],
        on="event_id", how="left")
    contacts["day"] = np.floor((contacts.bin_2min - contacts.start_time).dt.total_seconds() / 86400).astype(int)
    contacts = contacts[contacts.day.between(0, MAX_DAY)]
    cadence_n = len(contacts)
    support = contacts.groupby("animal_id").contact_5m.sum().rename("cadence_filtered_positive_5m_bins")
    keep = support[support > 100].index.tolist()
    contacts = contacts[contacts.animal_id.isin(keep)]

    event_day = (contacts.groupby(
        ["animal_id", "origin_group", "recipient_group", "event_id", "day"], observed=True)
                 .agg(eligible_bins=("bin_2min", "size"), contact_bins=("contact_5m", "sum"), bouts=("bout", "nunique"))
                 .reset_index())
    event_day["eligible_minutes"] = 2 * event_day.eligible_bins
    event_day["colocation_minutes"] = 2 * event_day.contact_bins
    event_day["colocation_proportion"] = event_day.contact_bins / event_day.eligible_bins
    event_day.to_csv(OUTDIR / "permanent_gt100_5m_colocation_event_days_cadence_filtered.csv", index=False)
    daily = (event_day.groupby(["animal_id", "origin_group", "recipient_group", "day"], observed=True)
             .agg(mean_minutes=("colocation_minutes", "mean"), mean_proportion=("colocation_proportion", "mean"),
                  segments=("event_id", "nunique"), eligible_minutes=("eligible_minutes", "sum"))
             .reset_index())
    daily.to_csv(OUTDIR / "permanent_gt100_5m_colocation_daily_cadence_filtered.csv", index=False)

    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
    routes = daily[["animal_id", "origin_group", "recipient_group"]].drop_duplicates().itertuples(index=False)
    for i, route in enumerate(routes):
        z = daily[(daily.animal_id.eq(route.animal_id)) &
                  (daily.origin_group.eq(route.origin_group)) &
                  (daily.recipient_group.eq(route.recipient_group))].sort_values("day")
        label = (f"{route.animal_id}: {route.origin_group} -> {route.recipient_group} "
                 f"(animal 5 m n={int(support.loc[route.animal_id])})")
        axes[0].plot(z.day, z.mean_minutes, marker="o", markersize=3, linewidth=1.8, color=colors[i], label=label)
        axes[1].plot(z.day, z.mean_proportion, marker="o", markersize=3, linewidth=1.8, color=colors[i])
    axes[0].set_ylabel("Mean sampled minutes within 5 m")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_ylabel("Proportion of eligible 2-minute bins")
    axes[1].set_xlabel("Days since recipient-segment entry")
    for ax in axes:
        ax.set_xlim(-1, MAX_DAY + 1); ax.set_ylim(bottom=0); ax.grid(alpha=.2)
    fig.suptitle("Daily time within 5 m during verified 2-minute sampling bouts\n"
                 "Bouts: gaps <=3 min and at least 15 focal fixes; selection reapplied after filtering")
    fig.tight_layout(rect=[0, 0, 1, .94])
    fig.savefig(OUTDIR / "permanent_gt100_5m_colocation_by_day_cadence_filtered.png", dpi=220)
    plt.close(fig)
    print(f"contact opportunities before cadence filter={before:,}; after cadence filter={cadence_n:,}; after support selection={len(contacts):,}")
    print(support.sort_values(ascending=False).to_string())
    print(event_day.groupby("animal_id").agg(segments=("event_id", "nunique"), days=("day", "nunique"),
          minutes=("colocation_minutes", "sum")).to_string())


if __name__ == "__main__":
    main()
