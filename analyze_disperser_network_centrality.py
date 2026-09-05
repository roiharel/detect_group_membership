from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


BINS = Path("outputs/disperser_finescale_integration/disperser_2min_contact_rows.parquet")
EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/canonical_disperser_events.csv")
MEMBERSHIP = Path(r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_local_2h_support\canonical_hourly_membership.parquet")
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
OUT = Path("outputs/disperser_network_centrality")
RADIUS_M = 5.0
MAX_WEEK = 7
MIN_NODES = 4
MIN_FOCAL_BINS = 30
N_BOOT = 1000
SEED = 20260810


METRICS = ["degree", "strength", "eigenvector", "betweenness", "harmonic_closeness"]


def distance_matrix(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lo, la = np.radians(lon), np.radians(lat)
    dlo, dla = lo[:, None] - lo[None, :], la[:, None] - la[None, :]
    a = np.sin(dla / 2) ** 2 + np.cos(la[:, None]) * np.cos(la[None, :]) * np.sin(dlo / 2) ** 2
    return 2 * 6_371_000 * np.arctan2(np.sqrt(a), np.sqrt(np.maximum(0, 1 - a)))


def load_positions() -> tuple[pd.DataFrame, pd.DataFrame]:
    focal = pd.read_parquet(BINS, columns=["event_id", "animal_id", "recipient_group", "hour", "bin_2min"]).drop_duplicates()
    focal["hour"] = pd.to_datetime(focal.hour); focal["bin_2min"] = pd.to_datetime(focal.bin_2min)
    events = pd.read_csv(EVENTS, parse_dates=["start_time"])[["event_id", "start_time"]]
    focal = focal.merge(events, on="event_id", how="left")
    focal["week"] = np.floor((focal.bin_2min - focal.start_time).dt.total_seconds() / (7 * 86400)).astype(int)
    focal = focal[focal.week.between(0, MAX_WEEK)]
    event_hours = focal[["event_id", "animal_id", "recipient_group", "hour", "start_time"]].drop_duplicates()
    membership = pd.read_parquet(MEMBERSHIP, columns=["animal_id", "window_start", "dynamic_social_unit"])
    membership["window_start"] = pd.to_datetime(membership.window_start)
    members = membership.merge(event_hours, left_on="window_start", right_on="hour", suffixes=("", "_focal"), how="inner")
    members = members[members.dynamic_social_unit.eq(members.recipient_group)][["event_id", "hour", "animal_id", "animal_id_focal"]].drop_duplicates()
    # Guarantee that the focal node is retained in every event-hour even if its canonical assignment is intermittent.
    focal_members = event_hours.rename(columns={"animal_id": "animal_id_focal"}).assign(animal_id=lambda x: x.animal_id_focal)[["event_id", "hour", "animal_id", "animal_id_focal"]]
    members = pd.concat([members, focal_members], ignore_index=True).drop_duplicates()
    animals = sorted(members.animal_id.astype(str).unique())
    start = focal.hour.min().tz_localize("UTC"); end = (focal.hour.max() + pd.Timedelta(hours=1)).tz_localize("UTC")
    gps = pd.read_parquet(GPS, columns=["animal_id", "timestamp", "location.long", "location.lat"],
                          filters=[("animal_id", "in", animals), ("timestamp", ">=", start), ("timestamp", "<", end)])
    gps["timestamp"] = pd.to_datetime(gps.timestamp, utc=True).dt.tz_localize(None)
    gps["hour"] = gps.timestamp.dt.floor("h"); gps["bin_2min"] = gps.timestamp.dt.floor("2min")
    gps = (gps.groupby(["hour", "bin_2min", "animal_id"], observed=True)
           .agg(lon=("location.long", "median"), lat=("location.lat", "median")).reset_index())
    positions = gps.merge(members, on=["hour", "animal_id"], how="inner")
    valid_bins = focal[["event_id", "bin_2min", "week"]].drop_duplicates()
    positions = positions.merge(valid_bins, on=["event_id", "bin_2min"], how="inner")
    return positions, event_hours[["event_id", "animal_id", "recipient_group"]].drop_duplicates()


def centralities(positions: pd.DataFrame, focal_lookup: pd.DataFrame) -> pd.DataFrame:
    rows = []
    focal_map = focal_lookup.set_index("event_id").animal_id.astype(str).to_dict()
    for (event_id, week), ew in positions.groupby(["event_id", "week"], observed=True):
        focal = focal_map[event_id]; coobs, contact, node_bins = defaultdict(int), defaultdict(int), defaultdict(int)
        for _, b in ew.groupby("bin_2min", observed=True, sort=False):
            b = b.drop_duplicates("animal_id"); ids = b.animal_id.astype(str).to_numpy()
            for a in ids: node_bins[a] += 1
            if len(ids) < 2: continue
            dist = distance_matrix(b.lon.to_numpy(float), b.lat.to_numpy(float))
            for i, j in combinations(range(len(ids)), 2):
                pair = tuple(sorted((ids[i], ids[j]))); coobs[pair] += 1
                if dist[i, j] <= RADIUS_M: contact[pair] += 1
        nodes = sorted(node_bins); n = len(nodes)
        if n < MIN_NODES or node_bins.get(focal, 0) < MIN_FOCAL_BINS: continue
        g = nx.Graph(); g.add_nodes_from(nodes)
        for pair, denom in coobs.items():
            weight = contact[pair] / denom
            if weight > 0: g.add_edge(*pair, weight=weight, distance=1 / weight)
        degree = nx.degree_centrality(g)
        strength = {a: sum(x["weight"] for _, _, x in g.edges(a, data=True)) / max(n - 1, 1) for a in nodes}
        if g.number_of_edges():
            try:
                eigen = nx.eigenvector_centrality_numpy(g, weight="weight")
            except Exception:
                try:
                    eigen = nx.eigenvector_centrality(g, weight="weight", max_iter=5000, tol=1e-8)
                except nx.PowerIterationFailedConvergence:
                    # Daily networks can be disconnected and numerically degenerate.
                    # Keep the metric defined without aborting the other centralities.
                    eigen = {a: 0.0 for a in nodes}
            between = nx.betweenness_centrality(g, weight="distance", normalized=True)
            harmonic = {a: v / max(n - 1, 1) for a, v in nx.harmonic_centrality(g, distance="distance").items()}
        else: eigen = between = harmonic = {a: 0.0 for a in nodes}
        values = {"degree": degree, "strength": strength, "eigenvector": eigen, "betweenness": between, "harmonic_closeness": harmonic}
        for metric, vals in values.items():
            residents = [vals[a] for a in nodes if a != focal]
            if focal not in vals or not residents: continue
            rank = (sum(v < vals[focal] for v in residents) + .5 * sum(v == vals[focal] for v in residents)) / len(residents)
            rows.append({"event_id": event_id, "week": int(week), "metric": metric, "focal": vals[focal],
                         "resident_median": float(np.median(residents)), "focal_minus_resident": vals[focal] - np.median(residents),
                         "focal_percentile": rank, "n_nodes": n, "focal_bins": node_bins[focal]})
    return pd.DataFrame(rows)


def summarize(cells: pd.DataFrame) -> pd.DataFrame:
    obs = (cells.groupby(["metric", "week"], observed=True)
           .agg(events=("event_id", "nunique"), mean_percentile=("focal_percentile", "mean"), mean_gap=("focal_minus_resident", "mean")).reset_index())
    ids = cells.event_id.unique(); rng = np.random.default_rng(SEED); out = []
    for rep in range(N_BOOT):
        chosen = rng.choice(ids, len(ids), replace=True)
        s = pd.concat([cells[cells.event_id.eq(e)] for e in chosen], ignore_index=True)
        q = s.groupby(["metric", "week"], observed=True).focal_percentile.mean().rename("value").reset_index(); q["rep"] = rep; out.append(q)
    draws = pd.concat(out, ignore_index=True)
    ci = draws.groupby(["metric", "week"], observed=True).value.quantile([.025, .975]).unstack().reset_index().rename(columns={.025: "low", .975: "high"})
    return obs.merge(ci, on=["metric", "week"]).query("events >= 5")


def plot(summary: pd.DataFrame) -> None:
    labels = {"degree": "Degree centrality", "strength": "Weighted strength", "eigenvector": "Eigenvector centrality",
              "betweenness": "Betweenness", "harmonic_closeness": "Harmonic closeness"}
    fig, axes = plt.subplots(3, 2, figsize=(13, 13)); axes = axes.ravel()
    for ax, metric in zip(axes, METRICS):
        z = summary[summary.metric.eq(metric)].sort_values("week")
        ax.plot(z.week, z.mean_percentile, marker="o", color="#7b3294", lw=2); ax.fill_between(z.week, z.low, z.high, color="#7b3294", alpha=.16)
        ax.axhline(.5, color="black", ls="--", lw=1); ax.set(title=labels[metric], xlabel="Weeks since entry", ylabel="Focal percentile among group members", ylim=(-.03, 1.03)); ax.grid(alpha=.2)
    axes[-1].axis("off"); fig.suptitle("Disperser position in the recipient group’s 5 m network\n0.5 = median recipient member; 95% CIs resample dispersal events")
    fig.tight_layout(rect=[0, 0, 1, .96]); fig.savefig(OUT / "disperser_network_centrality.png", dpi=220); plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "event_member_positions.parquet"
    if cache.exists(): positions = pd.read_parquet(cache); focal_lookup = pd.read_csv(OUT / "event_focal_lookup.csv")
    else:
        positions, focal_lookup = load_positions(); positions.to_parquet(cache, index=False); focal_lookup.to_csv(OUT / "event_focal_lookup.csv", index=False)
    cells = centralities(positions, focal_lookup); cells.to_csv(OUT / "event_week_centrality_cells.csv", index=False)
    summary = summarize(cells); summary.to_csv(OUT / "centrality_summary.csv", index=False); plot(summary)
    print(f"positions={len(positions):,} events={cells.event_id.nunique()} event-weeks={cells[['event_id','week']].drop_duplicates().shape[0]}")
    print(summary.to_string(index=False))


if __name__ == "__main__": main()
