"""Which clustering rule should define a social cluster? A comparison with criteria.

The canonical pipeline clusters ANIMALS BY THEIR MEDIAN POSITION IN AN HOUR, using an
adaptive kNN-scaled single linkage with edges clipped to 120-900 m. Three things about
that are worth testing rather than assuming:

  1. THE UPPER CLIP. 900 m is generous for a baboon group. Does 600 m change who counts
     as grouped, and in which direction?
  2. THE ALGORITHM. HDBSCAN disagrees sharply with every scale-fixing method (mean ARI
     0.42-0.45 in the earlier run). Is that HDBSCAN finding real substructure, or
     shattering groups?
  3. THE INPUT. Collars fix every ~2 minutes, up to 30 times an hour, and the fixes are
     NOT synchronised across animals. Collapsing that to one median position per hour
     throws away the whole within-hour record -- two animals can share a median position
     without ever having been close, and can be close all hour with medians 400 m apart.

There is no ground truth for "one group", so the methods are judged on four things that
do not require one:

  SEPARATION   silhouette of the hourly partition scored against the FINE-SCALE distance
               matrix (median dyadic distance over shared 2-minute bins). A partition
               that corresponds to a real gap in the fine-scale geometry scores high.
  STABILITY    adjusted Rand index between the partition from all fixes in the hour and
               the partition from a random half of them. A method reading noise is
               unstable.
  PERSISTENCE  adjusted Rand index between consecutive hours, on the animals present in
               both. Real groups persist; an artefact of the hour does not.
  COVERAGE BIAS  correlation between the number of clusters found and the number of
               animals collared that hour. A method whose cluster count tracks collar
               count is measuring the collars.

Methods compared, all on the same sampled hours:

  adaptive_900     median position, kNN single linkage k=2 factor 1.65, clip 120-900 m
                   -- the canonical pipeline
  adaptive_600     the same with the upper clip at 600 m
  adaptive_600_2m  the same rule, but on the FINE-SCALE distance matrix (median dyadic
                   distance across shared 2-minute bins) instead of median positions
  dbscan_500       DBSCAN eps 500 m, min_samples 1 -- the manuscript's stated setting
  dbscan_100       DBSCAN eps 100 m, min_samples 1
  hdbscan          HDBSCAN min_cluster_size 2, noise as singletons -- fixes no scale
  hdbscan_2m       HDBSCAN on the fine-scale distance matrix

Read-only. Outputs: outputs/general_structure_2026_09/phase5_clustering_options/
"""

from __future__ import annotations

import argparse
import json
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN, HDBSCAN
from sklearn.metrics import adjusted_rand_score, silhouette_score

# The canonical file is on the EAS share and is updated IN PLACE, so this points at a
# DATED local copy of it (2026-09-05: 30,143,804 fixes, 392 animals, to 2026-09-05).
# The previous local copy was undated and stale -- it carried 2,915 fixes for
# 25AA07_4S5T that the canonical cleaning has since retracted, which is exactly the kind
# of silent divergence a dated filename prevents.
GPS = Path("data/gps_v1_canonical_20260905.parquet")
OUT = Path("outputs/general_structure_2026_09/phase5_clustering_options")

WINDOW = ("2025-03-01", "2025-06-01")
N_HOURS = 240
MIN_ANIMALS = 12
BIN_MIN = 2                     # the fine-scale bin the pipeline already uses
ADAPT_K, ADAPT_FACTOR, ADAPT_MIN = 2, 1.65, 120.0
R_EARTH_M = 6371000.0
SEED = 20260905


# ------------------------------------------------------------------ geometry
def dist_matrix(lat, lon):
    """Pairwise metres on a local equirectangular projection."""
    lat0 = np.nanmean(lat)
    x = np.radians(lon) * R_EARTH_M * np.cos(np.radians(lat0))
    y = np.radians(lat) * R_EARTH_M
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    return np.sqrt(dx ** 2 + dy ** 2)


def components(adj):
    n = adj.shape[0]
    lab = np.zeros(n, dtype=int)
    cur = 0
    for s in range(n):
        if lab[s]:
            continue
        cur += 1
        stack = [s]
        lab[s] = cur
        while stack:
            v = stack.pop()
            for u in np.flatnonzero(adj[v]):
                if not lab[u]:
                    lab[u] = cur
                    stack.append(int(u))
    return lab


def adaptive_thresholds(D, max_edge, k=ADAPT_K, factor=ADAPT_FACTOR,
                        min_edge=ADAPT_MIN):
    """The per-dyad edge thresholds the adaptive rule would use."""
    n = D.shape[0]
    Dn = np.where(np.eye(n, dtype=bool), np.nan, D)
    srt = np.sort(Dn, axis=1)
    kth = min(max(k, 1), n - 1)
    scale = np.clip(srt[:, kth - 1], min_edge, max_edge)
    return np.minimum(max_edge,
                      np.maximum(min_edge,
                                 factor * np.minimum(scale[:, None], scale[None, :])))


def adaptive_labels(D, max_edge, k=ADAPT_K, factor=ADAPT_FACTOR, min_edge=ADAPT_MIN):
    """kNN-scaled single linkage, reimplemented from the canonical builder."""
    n = D.shape[0]
    if n < 2:
        return np.ones(n, dtype=int)
    Dn = np.where(np.eye(n, dtype=bool), np.nan, D)
    srt = np.sort(Dn, axis=1)
    kth = min(max(k, 1), n - 1)
    scale = np.clip(srt[:, kth - 1], min_edge, max_edge)
    thr = np.minimum(max_edge,
                     np.maximum(min_edge,
                                factor * np.minimum(scale[:, None], scale[None, :])))
    return components((D <= thr) & ~np.eye(n, dtype=bool))


def dbscan_labels(D, eps):
    return DBSCAN(eps=eps, min_samples=1,
                  metric="precomputed").fit_predict(np.array(D, copy=True)) + 1


def hdbscan_labels(D, allow_single=False):
    # allow_single_cluster matters a lot here: with it off, HDBSCAN can never return one
    # cluster, so an hour in which every animal really is together must be split. The
    # earlier agreement run used it ON, so both settings are compared.
    lab = HDBSCAN(min_cluster_size=2, allow_single_cluster=allow_single,
                  metric="precomputed").fit_predict(np.array(D, copy=True))
    out = lab.copy()
    nxt = (out.max() + 1) if out.max() >= 0 else 0
    for i in np.flatnonzero(out < 0):        # noise points become singletons
        out[i] = nxt
        nxt += 1
    return out + 1


METHODS = {
    "adaptive_900": lambda Dm, Df: adaptive_labels(Dm, 900.0),
    "adaptive_600": lambda Dm, Df: adaptive_labels(Dm, 600.0),
    "adaptive_600_2m": lambda Dm, Df: adaptive_labels(Df, 600.0),
    "dbscan_500": lambda Dm, Df: dbscan_labels(Dm, 500.0),
    # eps 300 m is what the parallel `New project` hourly-grouping pipeline actually
    # runs (build_hourly_grouping_table.py, --eps-meters 300, min_samples 1), so it is
    # the live alternative rather than a hypothetical one
    "dbscan_300": lambda Dm, Df: dbscan_labels(Dm, 300.0),
    "dbscan_300_2m": lambda Dm, Df: dbscan_labels(Df, 300.0),
    "dbscan_100": lambda Dm, Df: dbscan_labels(Dm, 100.0),
    "hdbscan": lambda Dm, Df: hdbscan_labels(Dm),
    "hdbscan_1ok": lambda Dm, Df: hdbscan_labels(Dm, allow_single=True),
    "hdbscan_2m": lambda Dm, Df: hdbscan_labels(Df),
    "hdbscan_2m_1ok": lambda Dm, Df: hdbscan_labels(Df, allow_single=True),
}
FINE_INPUT = {"adaptive_600_2m", "dbscan_300_2m", "hdbscan_2m", "hdbscan_2m_1ok"}


# ------------------------------------------------------------------ per hour
def hour_matrices(g, rng, half=False, floor=MIN_ANIMALS):
    """Median-position and fine-scale distance matrices for one hour.

    The fine-scale matrix is the MEDIAN dyadic distance over the 2-minute bins in which
    both animals have a fix -- the honest hourly summary of how far apart two animals
    actually were, as opposed to how far apart their median positions are. Dyads with no
    shared bin fall back to the median-position distance, which is all there is for them.
    """
    if half:
        keep = rng.random(len(g)) < 0.5
        g = g[keep]
        if g["animal_id"].nunique() < floor:
            return None
    animals = sorted(g["animal_id"].unique())
    ai = {a: i for i, a in enumerate(animals)}
    n = len(animals)

    med = g.groupby("animal_id")[["location.lat", "location.long"]].median()
    med = med.loc[animals]
    Dm = dist_matrix(med["location.lat"].to_numpy(), med["location.long"].to_numpy())

    bins = sorted(g["bin"].unique())
    lat = np.full((n, len(bins)), np.nan)
    lon = np.full((n, len(bins)), np.nan)
    bi = {b: j for j, b in enumerate(bins)}
    per = g.groupby(["animal_id", "bin"])[["location.lat", "location.long"]].mean()
    for (a, b), r in per.iterrows():
        lat[ai[a], bi[b]] = r["location.lat"]
        lon[ai[a], bi[b]] = r["location.long"]

    lat0 = np.nanmean(lat)
    x = np.radians(lon) * R_EARTH_M * np.cos(np.radians(lat0))
    y = np.radians(lat) * R_EARTH_M
    acc = np.full((n, n, len(bins)), np.nan)
    for j in range(len(bins)):
        dx = x[:, j][:, None] - x[:, j][None, :]
        dy = y[:, j][:, None] - y[:, j][None, :]
        acc[:, :, j] = np.sqrt(dx ** 2 + dy ** 2)
    shared = np.isfinite(acc).sum(axis=2)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Df = np.nanmedian(acc, axis=2)
    # dyads that never share a bin keep the median-position distance, which is all the
    # evidence there is for them
    Df = np.where(np.isfinite(Df), Df, Dm).astype(float)
    Df = (Df + Df.T) / 2.0                 # float noise can break exact symmetry
    np.fill_diagonal(Df, 0.0)
    np.fill_diagonal(Dm, 0.0)
    return animals, Dm, Df, shared


def score_hour(animals, Dm, Df, labels_by_method, groups):
    """Separation, cluster counts and group-pair calls for one hour."""
    out = {}
    for m, lab in labels_by_method.items():
        k = len(set(lab))
        sil = np.nan
        if 1 < k < len(lab):
            sil = float(silhouette_score(Df, lab, metric="precomputed"))
        calls = 0
        pairs = 0
        by_lab = {}
        for a, l in zip(animals, lab):
            by_lab.setdefault(l, set()).add(groups.get(a))
        present = sorted({groups.get(a) for a in animals if groups.get(a)})
        for u, v in combinations(present, 2):
            pairs += 1
            calls += int(any(u in s and v in s for s in by_lab.values()))
        out[m] = {"k": k, "silhouette": sil, "pair_calls": calls, "pairs": pairs}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    ap.add_argument("--n-hours", type=int, default=N_HOURS)
    # the window is a parameter because the conclusions were first measured on one
    # three-month stretch with 83 collared animals, which is not enough to commit to
    ap.add_argument("--start", default=WINDOW[0])
    ap.add_argument("--end", default=WINDOW[1])
    ap.add_argument("--tag", default="", help="suffix for the output filenames")
    # collar coverage rises ~120-fold across the record, and a kNN-scaled rule is exactly
    # the kind that can misbehave when an animal's 2nd-nearest neighbour is kilometres
    # away, so the sparse end has to be testable rather than filtered out
    ap.add_argument("--min-animals", type=int, default=MIN_ANIMALS)
    args = ap.parse_args()
    window = (args.start, args.end)
    tag = ("_" + args.tag) if args.tag else ""
    floor = args.min_animals
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    import pyarrow.parquet as pq
    cols = ["animal_id", "timestamp", "location.long", "location.lat", "group_id"]
    d = pq.read_table(GPS, columns=cols, filters=[
        ("timestamp", ">=", pd.Timestamp(window[0], tz="UTC")),
        ("timestamp", "<", pd.Timestamp(window[1], tz="UTC"))]).to_pandas()
    d["animal_id"] = d["animal_id"].astype(str)
    d["group_id"] = d["group_id"].astype(str)
    d = d[d["location.lat"].notna() & d["location.long"].notna()]
    d["hour"] = d["timestamp"].dt.floor("h")
    d["bin"] = d["timestamp"].dt.floor("%dmin" % BIN_MIN)
    groups = dict(d.groupby("animal_id")["group_id"].agg(
        lambda x: x.value_counts().index[0]))
    print("window %s to %s: %d fixes, %d animals"
          % (window[0], window[1], len(d), d["animal_id"].nunique()))

    n_per_hour = d.groupby("hour")["animal_id"].nunique()
    ok = sorted(n_per_hour[n_per_hour >= args.min_animals].index)
    pick = sorted(rng.choice(len(ok), size=min(args.n_hours, len(ok)),
                             replace=False).tolist())
    hours = [ok[i] for i in pick]
    hour_set = set(hours) | {h + pd.Timedelta(hours=1) for h in hours}
    d = d[d["hour"].isin(hour_set)]
    by_hour = {h: g for h, g in d.groupby("hour")}
    print("scoring %d hours (%d+ animals each)" % (len(hours), args.min_animals))

    rows, ari_rows, stab_rows, pers_rows, clip_rows = [], [], [], [], []
    for n_done, h in enumerate(hours, 1):
        got = hour_matrices(by_hour[h], rng, floor=floor)
        if got is None:
            continue
        animals, Dm, Df, shared = got
        labs = {m: f(Dm, Df) for m, f in METHODS.items()}
        sc = score_hour(animals, Dm, Df, labs, groups)
        thr = adaptive_thresholds(Dm, 900.0)
        off = ~np.eye(len(animals), dtype=bool)
        for m, s in sc.items():
            rows.append({"hour": h, "n_animals": len(animals), "method": m, **s,
                         "median_shared_bins": float(np.median(shared[shared > 0]))})
        clip_rows.append({
            "hour": h, "n_animals": len(animals),
            "thr_median": float(np.median(thr[off])),
            "thr_p95": float(np.percentile(thr[off], 95)),
            "thr_max": float(thr[off].max()),
            "frac_thr_over_600": float((thr[off] > 600.0).mean()),
            "frac_at_min_120": float((thr[off] <= 120.0 + 1e-9).mean()),
            "frac_at_max_900": float((thr[off] >= 900.0 - 1e-9).mean()),
            "dyads_600_vs_900": int(((thr[off] > 600.0)
                                     & (Dm[off] > 600.0)
                                     & (Dm[off] <= thr[off])).sum())})
        for a, b in combinations(METHODS, 2):
            ari_rows.append({"hour": h, "a": a, "b": b,
                             "ari": adjusted_rand_score(labs[a], labs[b])})
        # stability: the same hour from a random half of its fixes
        halfgot = hour_matrices(by_hour[h], rng, half=True, floor=floor)
        if halfgot is not None:
            ha, hDm, hDf, _ = halfgot
            keep = [i for i, a in enumerate(animals) if a in set(ha)]
            idx = {a: i for i, a in enumerate(ha)}
            for m, f in METHODS.items():
                l_half = f(hDm, hDf)
                l_full = [labs[m][i] for i in keep]
                l_sub = [l_half[idx[animals[i]]] for i in keep]
                if len(set(l_full)) > 1 or len(set(l_sub)) > 1:
                    stab_rows.append({"hour": h, "method": m,
                                      "ari": adjusted_rand_score(l_full, l_sub)})
        # persistence: the next hour, on the animals present in both
        nxt = by_hour.get(h + pd.Timedelta(hours=1))
        if nxt is not None and nxt["animal_id"].nunique() >= floor:
            got2 = hour_matrices(nxt, rng, floor=floor)
            if got2 is not None:
                a2, Dm2, Df2, _ = got2
                both = [a for a in animals if a in set(a2)]
                if len(both) >= floor:
                    i1 = {a: i for i, a in enumerate(animals)}
                    i2 = {a: i for i, a in enumerate(a2)}
                    for m, f in METHODS.items():
                        l2 = f(Dm2, Df2)
                        pers_rows.append({"hour": h, "method": m, "n_both": len(both),
                                          "ari": adjusted_rand_score(
                                              [labs[m][i1[a]] for a in both],
                                              [l2[i2[a]] for a in both])})
        if n_done % 40 == 0:
            print("  %d/%d hours" % (n_done, len(hours)))

    per_hour = pd.DataFrame(rows)
    ari = pd.DataFrame(ari_rows)
    stab = pd.DataFrame(stab_rows)
    pers = pd.DataFrame(pers_rows)
    per_hour.to_csv(args.output_dir / ("per_hour_by_method%s.csv" % tag), index=False)
    ari.to_csv(args.output_dir / ("pairwise_ari%s.csv" % tag), index=False)
    stab.to_csv(args.output_dir / ("stability_half_fixes%s.csv" % tag), index=False)
    pers.to_csv(args.output_dir / ("persistence_next_hour%s.csv" % tag), index=False)
    clip = pd.DataFrame(clip_rows)
    clip.to_csv(args.output_dir / ("adaptive_threshold_clip%s.csv" % tag), index=False)

    summary = {}
    for m in METHODS:
        p = per_hour[per_hour["method"].eq(m)]
        s = stab[stab["method"].eq(m)]["ari"]
        q = pers[pers["method"].eq(m)]["ari"]
        summary[m] = {
            "input": "fine-scale 2 min" if m in FINE_INPUT else "hourly median",
            "clusters_median": float(p["k"].median()),
            "clusters_mean": round(float(p["k"].mean()), 2),
            "silhouette_median": round(float(p["silhouette"].median()), 3),
            "stability_ari": round(float(s.median()), 3) if len(s) else None,
            "persistence_ari": round(float(q.median()), 3) if len(q) else None,
            "pair_call_rate": round(float(p["pair_calls"].sum() / p["pairs"].sum()), 4),
            "k_vs_n_animals_r": round(float(np.corrcoef(p["n_animals"], p["k"])[0, 1]), 3),
        }
    report = {
        "window": list(window), "gps_file": str(GPS),
        "hours_scored": int(per_hour["hour"].nunique()),
        "min_animals": args.min_animals, "fine_bin_minutes": BIN_MIN,
        "animals_per_hour_median": float(per_hour.groupby("hour")["n_animals"]
                                         .first().median()),
        "one_cluster_hours": {m: round(float(
            (per_hour[per_hour["method"].eq(m)]["k"] == 1).mean()), 3)
            for m in METHODS},
        "median_shared_bins_per_dyad": round(
            float(per_hour["median_shared_bins"].median()), 1),
        "criteria": {
            "silhouette": "hourly partition scored against the fine-scale (2 min) "
                          "distance matrix; higher means the partition matches a real "
                          "gap in the geometry",
            "stability_ari": "same hour, random half of its fixes; higher is better",
            "persistence_ari": "next hour, animals present in both; higher is better",
            "k_vs_n_animals_r": "cluster count against collar count; nearer zero is "
                                "better, a high value means the method measures collars",
        },
        "upper_clip": {
            "threshold_median_m": round(float(clip["thr_median"].median()), 1),
            "threshold_p95_m": round(float(clip["thr_p95"].median()), 1),
            "frac_dyads_over_600m": round(float(clip["frac_thr_over_600"].mean()), 4),
            "frac_dyads_pinned_at_120m": round(float(clip["frac_at_min_120"].mean()), 4),
            "frac_dyads_pinned_at_900m": round(float(clip["frac_at_max_900"].mean()), 4),
            "dyad_hours_where_600_would_cut": int(clip["dyads_600_vs_900"].sum()),
        },
        "methods": summary,
        "pairwise_ari_median": {("%s__%s" % (a, b)): round(float(
            ari[(ari["a"].eq(a)) & (ari["b"].eq(b))]["ari"].median()), 3)
            for a, b in combinations(METHODS, 2)},
    }
    with open(args.output_dir / ("clustering_options_report%s.json" % tag), "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("=" * 78)
    print("CLUSTERING OPTIONS")
    print("=" * 78)
    print("hours scored %d, median shared 2-min bins per dyad %.1f"
          % (report["hours_scored"], report["median_shared_bins_per_dyad"]))
    hdr = ("method", "input", "k", "silh", "stab", "pers", "callrate", "k~n")
    print("\n%-16s %-16s %5s %6s %6s %6s %9s %6s" % hdr)
    for m, s in summary.items():
        print("%-16s %-16s %5.1f %6s %6s %6s %9.4f %6s"
              % (m, s["input"], s["clusters_median"], s["silhouette_median"],
                 s["stability_ari"], s["persistence_ari"], s["pair_call_rate"],
                 s["k_vs_n_animals_r"]))
    print("\nmedian pairwise ARI:")
    for k, v in report["pairwise_ari_median"].items():
        print("    %-34s %s" % (k, v))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
