"""Where HDBSCAN and the scale-fixing methods disagree, drawn on the ground.

The agreement numbers say HDBSCAN partitions the same hour very differently (median ARI
0.42 against the adaptive rule). A number cannot say WHY, so this draws the worst
disagreements: the hours with the lowest ARI between HDBSCAN and the fine-scale adaptive
rule, with the same animals in the same positions under each method.

Reading the panels: a dot is an animal at its median position for that hour, with its
2-minute track behind it, and dot colour is the animal's ORIGIN GROUP -- fixed across all
panels, so the points do not move between columns. Only the outlines change: each
outline is one cluster as that method defines it. Disagreement is therefore visible as
different carving of identical points.

Output: docs/framing_2026-09-04/clustering_disagreement.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import adjusted_rand_score

from analyze_clustering_options import (GPS, MIN_ANIMALS, R_EARTH_M, SEED, WINDOW,
                                        adaptive_labels, dbscan_labels, hdbscan_labels,
                                        hour_matrices)
from build_level2_figure import CSS

OUT = Path("docs/framing_2026-09-04/clustering_disagreement.html")

N_SCAN = 150            # hours scanned for the worst disagreements
N_SHOW = 3              # examples drawn
CELL = 196.0
GAP = 22.0
SHOW = [("adaptive_600_2m", "adaptive, fine-scale"),
        ("dbscan_500", "DBSCAN 500 m"),
        ("hdbscan", "HDBSCAN")]
REF = "adaptive_600_2m"
PALETTE = ("#8a5a1f", "#1f6f8b", "#0f5f57", "#6b4a7a", "#b0448c", "#8a9a2b",
           "#b03a4a", "#2a5fa8", "#c46a3a", "#4a4a52")


def labels_for(name, Dm, Df):
    if name == "adaptive_600_2m":
        return adaptive_labels(Df, 600.0)
    if name == "adaptive_900":
        return adaptive_labels(Dm, 900.0)
    if name == "dbscan_500":
        return dbscan_labels(Dm, 500.0)
    return hdbscan_labels(Dm)


def cell(g, animals, labels, groups, gcol, x0, y0, title, sub):
    """One method's partition of one hour, on fixed points."""
    med = g.groupby("animal_id")[["location.lat", "location.long"]].median()
    lat0 = float(med["location.lat"].mean())
    kx = R_EARTH_M * math.cos(math.radians(lat0)) * math.pi / 180.0
    ky = R_EARTH_M * math.pi / 180.0
    tr = g[["animal_id", "bin", "location.lat", "location.long"]].copy()
    tr["mx"] = (tr["location.long"] - float(med["location.long"].mean())) * kx
    tr["my"] = (tr["location.lat"] - lat0) * ky
    pad = 60.0
    half = max(tr["mx"].abs().max(), tr["my"].abs().max()) + pad
    sc = (CELL - 26.0) / (2 * half)

    def sx(v):
        return x0 + CELL / 2 + v * sc

    def sy(v):
        return y0 + CELL / 2 - v * sc

    o = ['<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--paper)" '
         'stroke="var(--rule-soft)" stroke-width=".9"/>' % (x0, y0, CELL, CELL)]
    pos = {a: (float(med.loc[a, "location.long"] - med["location.long"].mean()) * kx,
               float(med.loc[a, "location.lat"] - lat0) * ky) for a in animals}
    lab = dict(zip(animals, labels))
    # one outline per cluster: a circle round its members, which is enough to see the
    # carving without implying a boundary the method never drew
    for cl in sorted(set(labels)):
        mem = [a for a in animals if lab[a] == cl]
        if len(mem) < 2:
            continue
        px = [sx(pos[a][0]) for a in mem]
        py = [sy(pos[a][1]) for a in mem]
        cx, cy = float(np.mean(px)), float(np.mean(py))
        rad = max(6.0, max(math.hypot(p - cx, q - cy) for p, q in zip(px, py)) + 5.0)
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="var(--ink-2)" '
                 'opacity=".07" stroke="var(--ink-2)" stroke-width=".9" '
                 'stroke-dasharray="3 2"/>' % (cx, cy, rad))
    for a, t in tr.sort_values("bin").groupby("animal_id"):
        if len(t) < 2:
            continue
        pts = " ".join("%.1f,%.1f" % (sx(v), sy(q))
                       for v, q in zip(t["mx"], t["my"]))
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width=".7" '
                 'opacity=".45"/>' % (pts, gcol[groups.get(a, "?")]))
    for a in animals:
        o.append('<circle cx="%.1f" cy="%.1f" r="2.2" fill="%s"/>'
                 % (sx(pos[a][0]), sy(pos[a][1]), gcol[groups.get(a, "?")]))
    o.append('<text class="ph" x="%.1f" y="%.1f">%s</text>' % (x0, y0 - 15, title))
    o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (x0, y0 - 5, sub))
    return chr(10).join(o)


def build():
    rng = np.random.default_rng(SEED)
    cols = ["animal_id", "timestamp", "location.long", "location.lat", "group_id"]
    d = pq.read_table(GPS, columns=cols, filters=[
        ("timestamp", ">=", pd.Timestamp(WINDOW[0], tz="UTC")),
        ("timestamp", "<", pd.Timestamp(WINDOW[1], tz="UTC"))]).to_pandas()
    d["animal_id"] = d["animal_id"].astype(str)
    d["group_id"] = d["group_id"].astype(str)
    d = d[d["location.lat"].notna() & d["location.long"].notna()]
    d["hour"] = d["timestamp"].dt.floor("h")
    d["bin"] = d["timestamp"].dt.floor("2min")
    groups = dict(d.groupby("animal_id")["group_id"].agg(
        lambda x: x.value_counts().index[0]))

    n_per = d.groupby("hour")["animal_id"].nunique()
    ok = sorted(n_per[n_per >= MIN_ANIMALS].index)
    pick = sorted(rng.choice(len(ok), size=min(N_SCAN, len(ok)),
                             replace=False).tolist())
    hours = [ok[i] for i in pick]
    by = {h: g for h, g in d[d["hour"].isin(set(hours))].groupby("hour")}

    scored = []
    for h in hours:
        got = hour_matrices(by[h], rng)
        if got is None:
            continue
        animals, Dm, Df, _ = got
        labs = {name: labels_for(name, Dm, Df) for name, _ in SHOW}
        scored.append((adjusted_rand_score(labs[REF], labs["hdbscan"]), h,
                       animals, labs))
    scored.sort(key=lambda t: t[0])
    picks = scored[:N_SHOW]
    print("worst disagreements (ARI %s vs hdbscan):" % REF)
    for ari, h, animals, labs in picks:
        print("  %s  ARI %.3f  n=%d  k: %s" % (h, ari, len(animals),
                                               {k: len(set(v)) for k, v in labs.items()}))

    units = sorted({groups[a] for _, _, animals, _ in picks for a in animals})
    gcol = {u: PALETTE[i % len(PALETTE)] for i, u in enumerate(units)}

    x0, y0 = 34.0, 74.0
    o, rows = [], []
    for r, (ari, h, animals, labs) in enumerate(picks):
        yy = y0 + r * (CELL + 52)
        o.append('<text class="pl" x="0" y="%.1f">%s</text>'
                 % (yy + CELL / 2, "abc"[r]))
        o.append('<text class="ts" x="0" y="%.1f">%s</text>'
                 % (yy + CELL / 2 + 12, h.strftime("%d %b")))
        for c, (name, nice) in enumerate(SHOW):
            xx = x0 + c * (CELL + GAP)
            k = len(set(labs[name]))
            a2 = adjusted_rand_score(labs[REF], labs[name])
            sub = "%d clusters%s" % (k, "" if name == REF else ", ARI %.2f vs left" % a2)
            o.append(cell(by[h], animals, labs[name], groups, gcol, xx, yy, nice, sub))
        rows.append({"hour": str(h), "ari": round(float(ari), 3),
                     "n_animals": len(animals),
                     "k": {kk: len(set(vv)) for kk, vv in labs.items()}})
    H = y0 + N_SHOW * (CELL + 52) + 44
    lx = x0
    for u in units:
        o.append('<circle cx="%.1f" cy="%.1f" r="2.6" fill="%s"/>' % (lx, H - 26, gcol[u]))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (lx + 7, H - 23, u))
        lx += 22 + len(u) * 4.4
    head = ('<svg viewBox="0 0 700 %.0f" role="img" aria-label="Three hours in which '
            'HDBSCAN partitions the animals very differently from the adaptive and '
            'DBSCAN rules, with the same animals in the same positions under each '
            'method and only the cluster outlines changing.">' % H)
    return head + chr(10) + chr(10).join(o) + chr(10) + "</svg>", rows


def page():
    fig, rows = build()
    worst = rows[0]
    return f"""<title>Clustering Disagreement</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Clustering methods &middot; diagnostic</div>
  <h1>Where HDBSCAN disagrees, and what it is doing</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      The three hours with the lowest adjusted Rand index between HDBSCAN and the
      fine-scale adaptive rule, out of {N_SCAN} hours sampled from {WINDOW[0]} to
      {WINDOW[1]}. <b>The points are identical across each row</b> &mdash; a dot is an
      animal at its median position for the hour, with its 2-minute track behind it, and
      colour is the animal's origin group. Only the dashed outlines change: each is one
      cluster as that column's method defines it, and singletons are left unringed. In
      the worst hour ({worst["hour"][:10]}, {worst["n_animals"]} animals, ARI
      {worst["ari"]:.2f}) the adaptive rule finds {worst["k"]["adaptive_600_2m"]} clusters
      and HDBSCAN finds {worst["k"]["hdbscan"]}. <b>The pattern is consistent: HDBSCAN
      splits single tight aggregations into several pieces</b>, because it reads the
      density gradient inside a group as a cluster boundary rather than as the ordinary
      internal structure of a foraging group. That is the right behaviour for finding
      substructure and the wrong behaviour for deciding which animals are together.
    </p>
  </div>

  <p class="note">
    Distances are metres on a local equirectangular projection about each hour's centroid.
    The adaptive column uses the fine-scale distance matrix &mdash; median dyadic distance
    over the 2-minute bins in which both animals have a fix &mdash; and the other two use
    hourly median positions, which is how each performs best. HDBSCAN is
    <code>min_cluster_size=2</code> with noise points treated as singletons;
    <code>allow_single_cluster</code> makes no difference to its partitions here (ARI
    1.00 between the two settings). Generated by
    <code>build_clustering_disagreement_figure.py</code> from
    <code>analyze_clustering_options.py</code>; seed {SEED}.
  </p>
</div>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(page(), encoding="utf-8")
    print("wrote %s (%d bytes)" % (args.output, args.output.stat().st_size))


if __name__ == "__main__":
    main()
