"""Does the choice of clustering algorithm change who counts as grouped?

The manuscript uses DBSCAN with eps = 500 m. The canonical pipeline uses an
adaptive kNN-scaled single-linkage rule (120-900 m). Neither is a hierarchical
method, and both fix a scale - the first explicitly, the second through its edge
bounds. HDBSCAN fixes no scale: it extracts clusters from a condensed hierarchy
at locally-varying density.

If the three agree, the choice of radius is a detail. If they disagree, then
"are these animals one group?" has no algorithm-free answer, and that is the
paper's thesis in methodological form.

The comparison runs at the level that matters for the paper's conclusions:
per hour, for each pair of origin groups present, does the method place them in
the same cluster? Agreement is then measured on that binary call, not only on
the raw partition.

Methods compared
  fixed_500   DBSCAN(eps=500 m, min_samples=1)      - the manuscript
  fixed_100   DBSCAN(eps=100 m, min_samples=1)      - a finer scale, for contrast
  adaptive    kNN-scaled single linkage, k=2, factor 1.65, clipped 120-900 m
              - reimplemented from build_canonical_robust_hourly_membership.py
  hdbscan     HDBSCAN(min_cluster_size=2), noise points treated as singletons

Read-only.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq
from sklearn.cluster import DBSCAN, HDBSCAN
from sklearn.metrics import adjusted_rand_score

PROJECT = Path(__file__).resolve().parent
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
R_EARTH_M = 6371000.0
# adaptive parameters, from the saved canonical run metadata
ADAPT_K, ADAPT_FACTOR, ADAPT_MIN, ADAPT_MAX = 2, 1.65, 120.0, 900.0


def haversine_matrix(lat, lon):
    la = np.radians(lat)[:, None]
    lo = np.radians(lon)[:, None]
    dla = la - la.T
    dlo = lo - lo.T
    a = np.sin(dla / 2) ** 2 + np.cos(la) * np.cos(la.T) * np.sin(dlo / 2) ** 2
    return 2 * R_EARTH_M * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def connected_components(adj: np.ndarray) -> np.ndarray:
    n = adj.shape[0]
    lab = np.full(n, -1)
    cur = 0
    for s in range(n):
        if lab[s] != -1:
            continue
        stack, lab[s] = [s], cur
        while stack:
            v = stack.pop()
            for w in np.flatnonzero(adj[v] & (lab == -1)):
                lab[w] = cur
                stack.append(w)
        cur += 1
    return lab


def adaptive_labels(lat, lon):
    n = len(lat)
    if n < 2:
        return np.zeros(n, dtype=int)
    d = haversine_matrix(lat, lon)
    srt = np.sort(np.where(np.eye(n, dtype=bool), np.nan, d), axis=1)
    kth = min(max(ADAPT_K, 1), n - 1)
    local = np.clip(srt[:, kth - 1], ADAPT_MIN, ADAPT_MAX)
    thr = np.minimum(ADAPT_MAX, np.maximum(ADAPT_MIN,
                                           ADAPT_FACTOR * np.minimum(local[:, None], local[None, :])))
    adj = (d <= thr) & ~np.eye(n, dtype=bool)
    return connected_components(adj)


def dbscan_labels(xy, eps):
    return DBSCAN(eps=eps, min_samples=1).fit_predict(xy)


def hdbscan_labels(xy):
    if len(xy) < 3:
        return np.arange(len(xy))
    lab = HDBSCAN(min_cluster_size=2, allow_single_cluster=True).fit_predict(xy)
    # noise (-1) becomes its own singleton, so partitions are comparable to
    # DBSCAN with min_samples=1, which never emits noise
    out, nxt = lab.copy(), lab.max() + 1
    for i in np.flatnonzero(lab == -1):
        out[i] = nxt
        nxt += 1
    return out


def load_hourly(path: Path, start: str, end: str) -> pd.DataFrame:
    f = pq.ParquetFile(path)
    cols = ["animal_id", "timestamp", "location.long", "location.lat", "group_id"]
    keep = []
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=cols).to_pandas()
        t["timestamp"] = pd.to_datetime(t["timestamp"], utc=True).dt.tz_localize(None)
        t = t[(t.timestamp >= start) & (t.timestamp < end)]
        if len(t):
            keep.append(t)
    d = pd.concat(keep, ignore_index=True).rename(
        columns={"location.long": "lon", "location.lat": "lat", "group_id": "group"})
    for c in ("animal_id", "group"):
        d[c] = d[c].astype(str)
    d = d.dropna(subset=["lat", "lon"])
    d["hour"] = d.timestamp.dt.floor("h")
    return (d.groupby(["hour", "animal_id", "group"], as_index=False)
             .agg(lat=("lat", "median"), lon=("lon", "median")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2025-01-01")
    ap.add_argument("--end", default="2025-07-01")
    ap.add_argument("--sample-hours", type=int, default=1200)
    ap.add_argument("--min-animals", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260904)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "clustering_method_agreement_2026-09-04")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    print(f"loading hourly positions {a.start} to {a.end} ...")
    h = load_hourly(GPS, a.start, a.end)
    per_hour = h.groupby("hour").animal_id.nunique()
    usable = per_hour[per_hour >= a.min_animals].index
    print(f"  {len(h):,} animal-hours | {h.animal_id.nunique()} animals | "
          f"{len(per_hour):,} hours, {len(usable):,} with >={a.min_animals} animals")
    hours = rng.choice(usable, size=min(a.sample_hours, len(usable)), replace=False)
    print(f"  sampling {len(hours):,} hours")

    METHODS = ["fixed_500", "fixed_100", "adaptive", "hdbscan"]
    part_rows, pair_rows = [], []
    for hr in hours:
        s = h[h.hour == hr]
        if len(s) < a.min_animals:
            continue
        lat, lon = s.lat.to_numpy(), s.lon.to_numpy()
        lat0 = lat.mean()
        xy = np.column_stack([lon * 111320.0 * np.cos(np.radians(lat0)), lat * 110540.0])
        labs = {
            "fixed_500": dbscan_labels(xy, 500.0),
            "fixed_100": dbscan_labels(xy, 100.0),
            "adaptive": adaptive_labels(lat, lon),
            "hdbscan": hdbscan_labels(xy),
        }
        rec = {"hour": hr, "n_animals": len(s), "n_groups": s.group.nunique()}
        for m in METHODS:
            rec[f"n_clusters_{m}"] = int(len(np.unique(labs[m])))
        for m1, m2 in combinations(METHODS, 2):
            rec[f"ari_{m1}__{m2}"] = float(adjusted_rand_score(labs[m1], labs[m2]))
        part_rows.append(rec)

        # per group-pair: does each method call them co-clustered this hour?
        groups = sorted(s.group.unique())
        gi = {g: np.flatnonzero((s.group == g).to_numpy()) for g in groups}
        for ga, gb in combinations(groups, 2):
            r = {"hour": hr, "group_a": ga, "group_b": gb,
                 "n_a": len(gi[ga]), "n_b": len(gi[gb])}
            for m in METHODS:
                la, lb = labs[m][gi[ga]], labs[m][gi[gb]]
                r[m] = bool(len(np.intersect1d(la, lb)) > 0)
            pair_rows.append(r)

    part = pd.DataFrame(part_rows)
    pair = pd.DataFrame(pair_rows)
    part.to_csv(a.output_dir / "partition_agreement_by_hour.csv", index=False)
    pair.to_csv(a.output_dir / "group_pair_calls_by_hour.csv", index=False)

    print(f"\n{'=' * 84}\nCLUSTER COUNTS per hour (median across {len(part):,} hours)\n{'=' * 84}")
    for m in METHODS:
        v = part[f"n_clusters_{m}"]
        print(f"  {m:<12} median {v.median():>5.1f}   mean {v.mean():>6.2f}   "
              f"range {v.min()}-{v.max()}")

    print(f"\n{'=' * 84}\nPARTITION AGREEMENT (adjusted Rand, median over hours)\n{'=' * 84}")
    for m1, m2 in combinations(METHODS, 2):
        v = part[f"ari_{m1}__{m2}"]
        print(f"  {m1:<12} vs {m2:<12} ARI median {v.median():>6.3f}   "
              f"IQR {v.quantile(.25):.3f}-{v.quantile(.75):.3f}")

    print(f"\n{'=' * 84}\nIS GROUP A WITH GROUP B? ({len(pair):,} group-pair-hours)\n{'=' * 84}")
    print(f"  {'method':<12}{'says together':>15}{'rate':>9}")
    for m in METHODS:
        print(f"  {m:<12}{int(pair[m].sum()):>15,}{pair[m].mean():>9.3%}")

    print(f"\n  pairwise agreement on that call:")
    agree = {}
    for m1, m2 in combinations(METHODS, 2):
        same = (pair[m1] == pair[m2]).mean()
        both = (pair[m1] & pair[m2]).sum()
        either = (pair[m1] | pair[m2]).sum()
        j = both / either if either else np.nan
        agree[f"{m1}__{m2}"] = {"agreement": float(same), "jaccard": float(j)}
        print(f"    {m1:<12} vs {m2:<12} agree {same:>7.2%}   "
              f"Jaccard on 'together' {j:>6.3f}")

    print(f"\n{'=' * 84}\nWHERE THEY DISAGREE\n{'=' * 84}")
    d = pair[pair.fixed_500 != pair.hdbscan]
    print(f"  fixed_500 vs hdbscan disagree on {len(d):,} of {len(pair):,} "
          f"({len(d)/len(pair):.1%}) group-pair-hours")
    print(f"    fixed_500 says together, hdbscan does not: "
          f"{int((pair.fixed_500 & ~pair.hdbscan).sum()):,}")
    print(f"    hdbscan says together, fixed_500 does not: "
          f"{int((~pair.fixed_500 & pair.hdbscan).sum()):,}")
    top = (d.groupby(["group_a", "group_b"]).size().sort_values(ascending=False).head(8))
    print("\n  dyads where the call differs most often:")
    for (ga, gb), n in top.items():
        tot = int(((pair.group_a == ga) & (pair.group_b == gb)).sum())
        print(f"    {ga} - {gb:<18} {n:>5} of {tot:>5} hours ({n/tot:.1%})")

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "window": [a.start, a.end], "hours_sampled": int(len(part)),
        "min_animals_per_hour": a.min_animals, "seed": a.seed,
        "methods": {
            "fixed_500": "DBSCAN eps=500 m, min_samples=1 - the manuscript's setting",
            "fixed_100": "DBSCAN eps=100 m, min_samples=1",
            "adaptive": f"kNN single linkage, k={ADAPT_K}, factor={ADAPT_FACTOR}, "
                        f"clipped {ADAPT_MIN}-{ADAPT_MAX} m - the canonical pipeline",
            "hdbscan": "HDBSCAN min_cluster_size=2, noise as singletons - no fixed scale",
        },
        "pair_call_rates": {m: float(pair[m].mean()) for m in METHODS},
        "pairwise_agreement": agree,
        "caveat": "Agreement is measured on hourly median positions from the raw GPS, not on "
                  "the pipeline's own temporal-support logic, so this isolates the clustering "
                  "step from everything downstream of it.",
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
