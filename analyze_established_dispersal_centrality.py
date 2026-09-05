from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_disperser_network_centrality import GPS, MEMBERSHIP, centralities


REVERSIBLE = Path(r"C:\Users\rharel\Documents\New project\outputs\reversible_departure_transitions\reversible_departure_events.csv")
PERMANENT_PERIODS = Path(r"C:\Users\rharel\Documents\New project\outputs\latest_v1_cleaned_full_resolution_binary_origin_exclusion\full_resolution_binary_origin_exclusion_event_periods.csv")
CURRENT_EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/canonical_disperser_events.csv")
OUT = Path("outputs/established_dispersal_centrality")
MAX_WEEK = 5
N_BOOT = 2000
SEED = 20260810
METRICS = ["degree", "strength", "eigenvector", "betweenness", "harmonic_closeness"]


def established_events() -> pd.DataFrame:
    r = pd.read_csv(REVERSIBLE, parse_dates=["away_start", "return_time"])
    current = pd.read_csv(CURRENT_EVENTS, parse_dates=["start_time", "end_time"]).rename(columns={"dynamic_social_unit": "recipient_group"})
    keep = []
    for e in current.itertuples(index=False):
        candidates = r[(r.animal_id.astype(str).eq(str(e.animal_id))) & r.origin_group.eq(e.origin_group)]
        overlap = candidates.apply(lambda z: max(0, (min(e.end_time, z.return_time) - max(e.start_time, z.away_start)).total_seconds()), axis=1)
        if len(overlap) and overlap.max() > 0:
            row = e._asdict(); row["source_event_id"] = candidates.loc[overlap.idxmax(), "reversible_event_id"]; keep.append(row)
    returned = pd.DataFrame(keep)[["event_id", "animal_id", "origin_group", "recipient_group", "start_time", "end_time"]]
    returned["established_outcome"] = "Returned to origin"

    p = pd.read_csv(PERMANENT_PERIODS, parse_dates=["event_start_time", "event_end_time"])
    p = p[p.event_type.eq("WITH_OTHER_GROUP") & p.destination_group.notna()].copy()
    p["segment"] = p.groupby("episode_id").cumcount() + 1
    p["event_id"] = p.episode_id.astype(str) + "__" + p.destination_group.astype(str) + "__" + p.segment.astype(str)
    permanent = p.rename(columns={"destination_group": "recipient_group", "event_start_time": "start_time", "event_end_time": "end_time"})
    permanent = permanent[["event_id", "animal_id", "origin_group", "recipient_group", "start_time", "end_time"]]
    permanent["established_outcome"] = "Permanent dispersal"
    events = pd.concat([returned, permanent], ignore_index=True)
    events["animal_id"] = events.animal_id.astype(str)
    events["start_time"] = pd.to_datetime(events.start_time); events["end_time"] = pd.to_datetime(events.end_time)
    # A group that serves as an animal's origin in any established departure
    # episode is never treated as a recipient/new group for that animal, even
    # if it appears as a visited group in a different episode.
    origin_sets = events.groupby("animal_id").origin_group.agg(lambda x: set(x.dropna().astype(str)))
    events["recipient_is_any_origin"] = [
        str(recipient) in origin_sets.get(animal, set())
        for animal, recipient in zip(events.animal_id, events.recipient_group)
    ]
    excluded = events[events.recipient_is_any_origin].copy()
    excluded.to_csv(OUT / "excluded_origin_group_recipient_events.csv", index=False)
    events = events[~events.recipient_is_any_origin].drop(columns="recipient_is_any_origin").copy()
    events["analysis_end"] = np.minimum(events.end_time.values, (events.start_time + pd.Timedelta(weeks=MAX_WEEK + 1)).values)
    return events.sort_values(["established_outcome", "animal_id", "start_time"]).reset_index(drop=True)


def load_positions(events: pd.DataFrame) -> pd.DataFrame:
    hour_rows = []
    for e in events.itertuples(index=False):
        for hour in pd.date_range(e.start_time.floor("h"), e.analysis_end.ceil("h"), freq="h"):
            hour_rows.append((e.event_id, e.animal_id, e.recipient_group, e.start_time, e.analysis_end, hour))
    eh = pd.DataFrame(hour_rows, columns=["event_id", "animal_id_focal", "recipient_group", "start_time", "analysis_end", "hour"])
    membership = pd.read_parquet(MEMBERSHIP, columns=["animal_id", "window_start", "dynamic_social_unit"])
    membership["window_start"] = pd.to_datetime(membership.window_start)
    members = membership.merge(eh, left_on="window_start", right_on="hour", how="inner")
    members = members[members.dynamic_social_unit.eq(members.recipient_group)][["event_id", "hour", "animal_id", "animal_id_focal", "start_time", "analysis_end"]]
    focal = eh.assign(animal_id=lambda x: x.animal_id_focal)[["event_id", "hour", "animal_id", "animal_id_focal", "start_time", "analysis_end"]]
    members = pd.concat([members, focal], ignore_index=True).drop_duplicates(["event_id", "hour", "animal_id"])
    animals = sorted(members.animal_id.astype(str).unique())
    lo, hi = events.start_time.min().tz_localize("UTC"), events.analysis_end.max().tz_localize("UTC")
    gps = pd.read_parquet(GPS, columns=["animal_id", "timestamp", "location.long", "location.lat"],
                          filters=[("animal_id", "in", animals), ("timestamp", ">=", lo), ("timestamp", "<=", hi)])
    gps["timestamp"] = pd.to_datetime(gps.timestamp, utc=True).dt.tz_localize(None)
    gps["hour"] = gps.timestamp.dt.floor("h"); gps["bin_2min"] = gps.timestamp.dt.floor("2min")
    gps = (gps.groupby(["hour", "bin_2min", "animal_id"], observed=True)
           .agg(lon=("location.long", "median"), lat=("location.lat", "median")).reset_index())
    pos = gps.merge(members, on=["hour", "animal_id"], how="inner")
    pos = pos[pos.bin_2min.between(pos.start_time, pos.analysis_end)].copy()
    pos["week"] = np.floor((pos.bin_2min - pos.start_time).dt.total_seconds() / (7 * 86400)).astype(int)
    return pos[pos.week.between(0, MAX_WEEK)]


def summarize(cells: pd.DataFrame) -> pd.DataFrame:
    obs = (cells.groupby(["metric", "week", "established_outcome"], observed=True)
           .agg(events=("event_id", "nunique"), mean_percentile=("focal_percentile", "mean")).reset_index())
    rng = np.random.default_rng(SEED); rows = []
    for outcome, z in cells.groupby("established_outcome", observed=True):
        ids = z.event_id.unique()
        for rep in range(N_BOOT):
            chosen = rng.choice(ids, len(ids), replace=True)
            s = pd.concat([z[z.event_id.eq(e)] for e in chosen], ignore_index=True)
            q = s.groupby(["metric", "week"], observed=True).focal_percentile.mean().rename("value").reset_index()
            q["rep"] = rep; q["established_outcome"] = outcome; rows.append(q)
    draws = pd.concat(rows, ignore_index=True)
    ci = (draws.groupby(["metric", "week", "established_outcome"], observed=True).value.quantile([.025, .975]).unstack()
          .reset_index().rename(columns={.025: "low", .975: "high"}))
    return obs.merge(ci, on=["metric", "week", "established_outcome"]).query("events >= 3")


def plot(summary: pd.DataFrame) -> None:
    labels = {"degree": "Degree", "strength": "Weighted strength", "eigenvector": "Eigenvector",
              "betweenness": "Betweenness", "harmonic_closeness": "Harmonic closeness"}
    colors = {"Permanent dispersal": "#7b3294", "Returned to origin": "#e36c09"}
    fig, axes = plt.subplots(3, 2, figsize=(13, 13)); axes = axes.ravel()
    for ax, metric in zip(axes, METRICS):
        for outcome in ["Permanent dispersal", "Returned to origin"]:
            z = summary[(summary.metric.eq(metric)) & summary.established_outcome.eq(outcome) & summary.events.ge(5)].sort_values("week")
            ax.plot(z.week, z.mean_percentile, marker="o", color=colors[outcome], label=outcome)
            ax.fill_between(z.week, z.low, z.high, color=colors[outcome], alpha=.13)
        ax.axhline(.5, color="black", ls="--", lw=1); ax.set(title=labels[metric], xlabel="Weeks since recipient-group entry", ylabel="Focal percentile", ylim=(-.03, 1.03)); ax.grid(alpha=.2)
    axes[0].legend(frameon=False); axes[-1].axis("off")
    fig.suptitle("Recipient-group network integration under established dispersal outcomes\n"
                 "Returned: explicit return_time; permanent: binary origin-exclusion event")
    fig.tight_layout(rect=[0, 0, 1, .96]); fig.savefig(OUT / "established_outcome_centrality.png", dpi=220); plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True); events = established_events(); events.to_csv(OUT / "established_events.csv", index=False)
    cache = OUT / "event_member_positions_established_segments.parquet"
    positions = pd.read_parquet(cache) if cache.exists() else load_positions(events)
    if not cache.exists(): positions.to_parquet(cache, index=False)
    lookup = events[["event_id", "animal_id", "recipient_group"]]
    cells = centralities(positions, lookup).merge(events[["event_id", "established_outcome"]], on="event_id", how="left")
    cells.to_csv(OUT / "event_week_centrality_cells.csv", index=False)
    summary = summarize(cells); summary.to_csv(OUT / "established_outcome_centrality_summary.csv", index=False); plot(summary)
    excluded = pd.read_csv(OUT / "excluded_origin_group_recipient_events.csv")
    print(f"excluded recipient-is-origin events={len(excluded)}")
    print(events.established_outcome.value_counts().to_string()); print(f"analyzed events={cells.event_id.nunique()} positions={len(positions):,}")
    print(cells.groupby("established_outcome").event_id.nunique().to_string()); print(summary.to_string(index=False))


if __name__ == "__main__": main()
