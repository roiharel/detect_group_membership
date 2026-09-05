from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import adjusted_rand_score, roc_auc_score


CANONICAL = Path(r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_shared_full_20260722\canonical_hourly_membership_with_association_events.parquet")
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
OUT = Path("outputs/presplit_cositting_prediction")
MIN_COLLARS = 8
MAX_PRE_HOURS = 6
MIN_SPLIT_HOURS = 2
MIN_PAIR_BINS = 10
THRESHOLD_M = 5.0
SEED = 20260808


def haversine_matrix(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon = np.radians(lon)
    lat = np.radians(lat)
    dlon = lon[:, None] - lon[None, :]
    dlat = lat[:, None] - lat[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
    return 2 * 6371000.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def make_events(c: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    c = c.copy()
    c["window_start"] = pd.to_datetime(c["window_start"], utc=True)
    sizes = (
        c.groupby(["dynamic_social_unit", "window_start", "observed_cluster_id"], observed=True)
        .animal_id.nunique().rename("cluster_n").reset_index()
    )
    hours = (
        c.groupby(["dynamic_social_unit", "window_start"], observed=True)
        .agg(n_animals=("animal_id", "nunique"), n_clusters=("observed_cluster_id", "nunique"))
        .reset_index()
    )
    multi = sizes[sizes.cluster_n.ge(2)].groupby(["dynamic_social_unit", "window_start"], observed=True).size().rename("n_multi").reset_index()
    hours = hours.merge(multi, how="left", on=["dynamic_social_unit", "window_start"])
    hours["n_multi"] = hours.n_multi.fillna(0).astype(int)
    hours["state"] = "other"
    hours.loc[(hours.n_animals.ge(MIN_COLLARS)) & hours.n_clusters.eq(1), "state"] = "cohesive"
    hours.loc[(hours.n_animals.ge(MIN_COLLARS)) & hours.n_multi.ge(2), "state"] = "split"
    hours = hours.sort_values(["dynamic_social_unit", "window_start"]).reset_index(drop=True)

    event_rows, pre_rows, outcome_rows = [], [], []
    for unit, h in hours.groupby("dynamic_social_unit", sort=True, observed=True):
        h = h.reset_index(drop=True)
        for i in range(1, len(h)):
            onset = h.iloc[i]
            prev = h.iloc[i - 1]
            if onset.state != "split" or prev.state != "cohesive" or onset.window_start - prev.window_start > pd.Timedelta("1.5h"):
                continue
            # Require a persistent split rather than a single-hour clustering flicker.
            persistent = 1
            j = i + 1
            while j < len(h) and h.iloc[j].state == "split" and h.iloc[j].window_start - h.iloc[j - 1].window_start <= pd.Timedelta("1.5h"):
                persistent += 1
                j += 1
            if persistent < MIN_SPLIT_HOURS:
                continue
            pre = []
            j = i - 1
            while j >= 0 and len(pre) < MAX_PRE_HOURS and h.iloc[j].state == "cohesive":
                if pre and pre[-1] - h.iloc[j].window_start > pd.Timedelta("1.5h"):
                    break
                pre.append(h.iloc[j].window_start)
                j -= 1
            if not pre:
                continue
            eid = f"{unit}_SPLIT_{onset.window_start.strftime('%Y%m%dT%H%M%S')}"
            event_rows.append({"event_id": eid, "dynamic_social_unit": unit, "split_onset": onset.window_start,
                               "n_pre_hours": len(pre), "persistent_split_hours": persistent,
                               "n_animals_onset": int(onset.n_animals), "n_clusters_onset": int(onset.n_clusters)})
            pre_rows.extend({"event_id": eid, "dynamic_social_unit": unit, "window_start": t} for t in pre)
            o = c[(c.dynamic_social_unit == unit) & (c.window_start == onset.window_start)].copy()
            osizes = o.groupby("observed_cluster_id").animal_id.transform("nunique")
            o = o[osizes.ge(2)]
            outcome_rows.extend({"event_id": eid, "animal_id": str(r.animal_id), "split_cluster": str(r.observed_cluster_id)} for r in o.itertuples())
    return pd.DataFrame(event_rows), pd.DataFrame(pre_rows), pd.DataFrame(outcome_rows)


def extract_positions(c: pd.DataFrame, pre_hours: pd.DataFrame) -> pd.DataFrame:
    keys = c.merge(pre_hours[["dynamic_social_unit", "window_start"]].drop_duplicates(), on=["dynamic_social_unit", "window_start"], how="inner")
    keys = keys[["animal_id", "dynamic_social_unit", "window_start"]].drop_duplicates()
    pieces = []
    pf = pq.ParquetFile(GPS)
    for batch_i, batch in enumerate(pf.iter_batches(columns=["animal_id", "timestamp", "location.long", "location.lat"], batch_size=750_000), 1):
        x = batch.to_pandas()
        x["animal_id"] = x.animal_id.astype(str)
        x["timestamp"] = pd.to_datetime(x.timestamp, utc=True)
        x["window_start"] = x.timestamp.dt.floor("h")
        x = x.merge(keys, on=["animal_id", "window_start"], how="inner")
        if len(x):
            pieces.append(x)
        if batch_i % 8 == 0:
            print(f"GPS batches {batch_i}/{pf.num_row_groups * 2}: retained {sum(map(len, pieces)):,} rows", flush=True)
    if not pieces:
        return pd.DataFrame()
    p = pd.concat(pieces, ignore_index=True)
    p["bin_2m"] = p.timestamp.dt.floor("2min")
    return (p.groupby(["dynamic_social_unit", "window_start", "bin_2m", "animal_id"], observed=True)
             .agg(lon=("location.long", "median"), lat=("location.lat", "median")).reset_index())


def hourly_dyads(pos: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (unit, hour, _), b in pos.groupby(["dynamic_social_unit", "window_start", "bin_2m"], sort=False, observed=True):
        if len(b) < 2:
            continue
        ids = b.animal_id.to_numpy(str)
        dist = haversine_matrix(b.lon.to_numpy(float), b.lat.to_numpy(float))
        for i, j in zip(*np.triu_indices(len(b), 1)):
            left, right = sorted((ids[i], ids[j]))
            rows.append((unit, hour, left, right, int(dist[i, j] <= THRESHOLD_M)))
    d = pd.DataFrame(rows, columns=["dynamic_social_unit", "window_start", "animal_a", "animal_b", "within5"])
    if d.empty:
        return d
    return (d.groupby(["dynamic_social_unit", "window_start", "animal_a", "animal_b"], observed=True)
             .agg(coobserved_bins=("within5", "size"), within5_bins=("within5", "sum")).reset_index())


def event_pairs(pre: pd.DataFrame, outcomes: pd.DataFrame, hourly: pd.DataFrame) -> pd.DataFrame:
    d = pre.merge(hourly, on=["dynamic_social_unit", "window_start"], how="inner")
    d = (d.groupby(["event_id", "dynamic_social_unit", "animal_a", "animal_b"], observed=True)
         .agg(coobserved_bins=("coobserved_bins", "sum"), within5_bins=("within5_bins", "sum")).reset_index())
    oa = outcomes.rename(columns={"animal_id": "animal_a", "split_cluster": "cluster_a"})
    ob = outcomes.rename(columns={"animal_id": "animal_b", "split_cluster": "cluster_b"})
    d = d.merge(oa, on=["event_id", "animal_a"]).merge(ob, on=["event_id", "animal_b"])
    d = d[d.coobserved_bins.ge(MIN_PAIR_BINS)].copy()
    d["cosit_rate_5m"] = d.within5_bins / d.coobserved_bins
    d["same_subgroup"] = (d.cluster_a == d.cluster_b).astype(int)
    return d


def event_metrics(pairs: pd.DataFrame, outcomes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (eid, unit), d in pairs.groupby(["event_id", "dynamic_social_unit"], observed=True):
        if d.same_subgroup.nunique() < 2:
            continue
        auc = roc_auc_score(d.same_subgroup, d.cosit_rate_5m)
        # Community recovery uses only animals represented in eligible dyads.
        g = nx.Graph()
        for r in d.itertuples():
            g.add_edge(r.animal_a, r.animal_b, weight=float(r.cosit_rate_5m) + 1e-9)
        comms = list(nx.algorithms.community.greedy_modularity_communities(g, weight="weight"))
        pred = {a: k for k, comm in enumerate(comms) for a in comm}
        truth = outcomes[outcomes.event_id.eq(eid)].set_index("animal_id").split_cluster.to_dict()
        common = sorted(set(pred) & set(truth))
        ari = adjusted_rand_score([truth[a] for a in common], [pred[a] for a in common]) if len(common) >= 4 else np.nan
        rows.append({"event_id": eid, "dynamic_social_unit": unit, "n_pairs": len(d),
                     "n_same_pairs": int(d.same_subgroup.sum()), "n_different_pairs": int((1-d.same_subgroup).sum()),
                     "auc": auc, "ari": ari,
                     "same_mean_5m": d.loc[d.same_subgroup.eq(1), "cosit_rate_5m"].mean(),
                     "different_mean_5m": d.loc[d.same_subgroup.eq(0), "cosit_rate_5m"].mean()})
    return pd.DataFrame(rows)


def bootstrap_and_permute(metrics: pd.DataFrame, pairs: pd.DataFrame, n_boot=2000, n_perm=500) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    aucs = metrics.auc.to_numpy()
    diffs = (metrics.same_mean_5m - metrics.different_mean_5m).to_numpy()
    b_auc = np.array([rng.choice(aucs, len(aucs), replace=True).mean() for _ in range(n_boot)])
    b_diff = np.array([rng.choice(diffs, len(diffs), replace=True).mean() for _ in range(n_boot)])
    null = []
    event_data = []
    for _, d in pairs[pairs.event_id.isin(metrics.event_id)].groupby("event_id"):
        animals = pd.Index(pd.unique(pd.concat([d.animal_a, d.animal_b], ignore_index=True)))
        cluster = pd.concat([
            d[["animal_a", "cluster_a"]].rename(columns={"animal_a": "animal", "cluster_a": "cluster"}),
            d[["animal_b", "cluster_b"]].rename(columns={"animal_b": "animal", "cluster_b": "cluster"}),
        ]).drop_duplicates("animal").set_index("animal").loc[animals, "cluster"].to_numpy()
        event_data.append((d.cosit_rate_5m.to_numpy(), animals.get_indexer(d.animal_a), animals.get_indexer(d.animal_b), cluster))
    for _ in range(n_perm):
        vals = []
        for scores, ia, ib, cluster in event_data:
            # Shuffle animals among the observed subgroups, preserving subgroup sizes
            # and the dependence among all dyads sharing an animal.
            shuffled = rng.permutation(cluster)
            y = (shuffled[ia] == shuffled[ib]).astype(int)
            if np.unique(y).size == 2:
                vals.append(roc_auc_score(y, scores))
        null.append(np.mean(vals))
    null = np.asarray(null)
    obs = aucs.mean()
    return pd.DataFrame([{"n_events": len(metrics), "n_groups": metrics.dynamic_social_unit.nunique(),
                          "mean_event_auc": obs, "auc_ci_low": np.quantile(b_auc, .025), "auc_ci_high": np.quantile(b_auc, .975),
                          "permutation_p_one_sided": (1 + np.sum(null >= obs)) / (1 + len(null)),
                          "mean_same_minus_different_5m": diffs.mean(), "difference_ci_low": np.quantile(b_diff, .025),
                          "difference_ci_high": np.quantile(b_diff, .975), "mean_event_ari": metrics.ari.mean()}])


def plots(metrics: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    ax[0].hist(metrics.auc, bins=np.linspace(0, 1, 21), color="#536dfe", alpha=.82)
    ax[0].axvline(.5, color="black", ls="--", lw=1)
    ax[0].axvline(metrics.auc.mean(), color="#d81b60", lw=2)
    ax[0].set(xlabel="Event-level AUC", ylabel="Split events", title="Does cohesive-period 5 m co-sitting predict split subgroup?")
    m = metrics[["same_mean_5m", "different_mean_5m"]].dropna()
    for r in m.itertuples(index=False):
        ax[1].plot([0, 1], [r.different_mean_5m, r.same_mean_5m], color="0.75", lw=.6, alpha=.5)
    ax[1].scatter(np.zeros(len(m)), m.different_mean_5m, s=12, color="#ef6c00", alpha=.65)
    ax[1].scatter(np.ones(len(m)), m.same_mean_5m, s=12, color="#00897b", alpha=.65)
    ax[1].plot([0, 1], [m.different_mean_5m.mean(), m.same_mean_5m.mean()], color="black", lw=3)
    ax[1].set(xticks=[0, 1], xticklabels=["Different subgroup", "Same subgroup"], ylabel="Pre-split proportion within 5 m", title="Event-level dyadic contrast")
    fig.tight_layout(); fig.savefig(OUT / "presplit_cositting_prediction.png", dpi=220); plt.close(fig)

    group = metrics.groupby("dynamic_social_unit").agg(n_events=("auc", "size"), mean_auc=("auc", "mean"), sd=("auc", "std")).reset_index()
    group = group[group.n_events.ge(3)].sort_values("mean_auc")
    group["se"] = group.sd / np.sqrt(group.n_events)
    fig, ax = plt.subplots(figsize=(8, max(4, .4 * len(group))))
    y = np.arange(len(group)); ax.errorbar(group.mean_auc, y, xerr=1.96*group.se, fmt="o", color="#3949ab", capsize=2)
    ax.axvline(.5, color="black", ls="--", lw=1); ax.set(yticks=y, yticklabels=[f"{u} (n={n})" for u,n in zip(group.dynamic_social_unit, group.n_events)], xlabel="Mean event AUC (approx. 95% CI)", title="Prediction by social unit")
    fig.tight_layout(); fig.savefig(OUT / "presplit_prediction_by_group.png", dpi=220); plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    c = pd.read_parquet(CANONICAL, columns=["window_start", "animal_id", "dynamic_social_unit", "observed_cluster_id"])
    c["animal_id"] = c.animal_id.astype(str)
    c["window_start"] = pd.to_datetime(c.window_start, utc=True)
    events, pre, outcomes = make_events(c)
    events.to_csv(OUT / "candidate_split_events.csv", index=False)
    print(f"Candidate persistent cohesive-to-split transitions: {len(events):,}", flush=True)
    pos_path = OUT / "cohesive_2m_positions.parquet"
    if pos_path.exists():
        pos = pd.read_parquet(pos_path)
    else:
        pos = extract_positions(c, pre)
        pos.to_parquet(pos_path, index=False)
    print(f"Retained 2-minute animal positions: {len(pos):,}", flush=True)
    hourly_path = OUT / "cohesive_hour_dyads_5m.parquet"
    if hourly_path.exists():
        hourly = pd.read_parquet(hourly_path)
    else:
        hourly = hourly_dyads(pos)
        hourly.to_parquet(hourly_path, index=False)
    pairs = event_pairs(pre, outcomes, hourly)
    pairs.to_parquet(OUT / "event_dyad_predictions.parquet", index=False)
    metrics = event_metrics(pairs, outcomes)
    metrics.to_csv(OUT / "event_prediction_metrics.csv", index=False)
    summary = bootstrap_and_permute(metrics, pairs)
    summary.to_csv(OUT / "prediction_summary.csv", index=False)
    plots(metrics)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
