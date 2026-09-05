"""Level-2 figure: the structures the level-1 events build up to.

Level 1 is the event -- a split, a merge, an animal alone. Level 2 is what those events
aggregate into: modularity inside a unit, an encounter network between units, and a
transfer network of animals moving between them. Extent and variation only, no causes
or consequences.

Four panels on one grid, every plotting box the same size:

  (a) ONE GROUP, WEEK BY WEEK. Lilac, the unit that carries most of the modularity in
      this population -- because a per-unit average hides the thing worth seeing, which
      is that the group separates and re-fuses repeatedly.
  (b) BOTH BETWEEN-UNIT NETWORKS ON ONE LAYOUT. Group encounters and animal transfers
      share their nodes, so drawing them apart forces a comparison by eye. The layout is
      an ellipse rather than a circle so the panel fills the same box as the others.
  (c) CONCENTRATION, ALL ON ONE AXIS. Four quantities, one diagonal, each Gini with a
      cluster-bootstrap interval printed in the empty half of the plot.
  (d) DURATION, ALL ON ONE AXIS, IN HOURS. Encounters last hours to a day, modular
      phases weeks, excursions a night to permanent -- with bootstrap bands.

Nothing in the figure carries a sentence of prose: numbers live in the caption, the
panels carry axes, keys and marks only.

A note on what counts as a between-group encounter. A detected encounter fires when two
units share an association cluster, and 22.3% of those dyad-days have a SINGLE animal on
one side -- that is one animal visiting, which is the individual axis, not the group one.
Panel (b) therefore requires two or more animals on both sides, which is what leaves the
blue encounter chords and the plum transfer chords measuring different things.

Durations in (d) are ELAPSED hours in the state, on one convention for all three: 24 h
per encounter-day, 168 h per well-covered modular week, 24 h per away-night. Nothing is
therefore shorter than a day, which is a statement about the resolution of a nightly
state assignment rather than about the animals.

Every interval is a cluster bootstrap on the natural unit of replication -- dyad for
encounters, unit for modularity, animal for excursions -- never on the row.

Output: docs/framing_2026-09-04/level2_figure.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("outputs/general_structure_2026_09")
FROZEN = BASE / "phase4d_axis_b_frozen/weekly_network_metrics_frozen.csv"
OPP = BASE / "phase1_opportunity/opportunity_dyad_day.csv"
EDGES = BASE / "phase4g_excursion_destinations/transfer_edges_dominant.csv"
DESTS = BASE / "phase4g_excursion_destinations/excursion_destinations_dominant.csv"
EXC = BASE / "phase4b_individual_axis/excursions_dominant_gap0.csv"
INPUTS = BASE / "phase4h_level2_inputs"
BOUTS_H = INPUTS / "encounter_bouts_hourly.csv"
SPLITS_H = INPUTS / "split_bouts_hourly.csv"
EXC_H = INPUTS / "excursion_bouts_hourly.csv"
OUT = Path("docs/framing_2026-09-04/level2_figure.html")

WELL_COVERED = 12
MIN_WEEKS = 20
MOD_EPS = 0.001
EXAMPLE_UNIT = "Lilac"
N_BOOT = 1200
DUR_MAX_H = 13000.0
DAY_H = 24.0
WEEK_H = 7 * DAY_H
MONTH_H = 30 * DAY_H
YEAR_H = 365 * DAY_H
SEED = 20260904
FLICKER_H = 2.0
COMM_COLOUR = ("var(--c4)", "var(--c1)", "var(--c6)", "var(--c5)")

BOX_W, BOX_H = 262.0, 194.0
COL_X = (58.0, 388.0)
ROW_Y = (78.0, 380.0)


# ------------------------------------------------------------------ statistics
def lorenz(values):
    v = np.sort(np.asarray(values, dtype=float))[::-1]
    v = v[v > 0]
    if len(v) == 0:
        return np.array([0.0]), np.array([0.0])
    x = np.arange(1, len(v) + 1) / len(v)
    y = np.cumsum(v) / v.sum()
    return np.concatenate([[0.0], x]), np.concatenate([[0.0], y])


def gini(values):
    v = np.sort(np.asarray(values, dtype=float))
    v = v[v > 0]
    n = len(v)
    if n < 2:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(v) / (n * v.sum()))


def by_cluster(values, clusters):
    """Group values by their cluster label, as a list of arrays."""
    v = np.asarray(values, dtype=float)
    c = np.asarray(clusters)
    return [v[c == k] for k in pd.unique(c)]


def boot_gini(values, clusters, rng, n_boot=N_BOOT):
    """Cluster bootstrap CI for a Gini: resample clusters, not rows."""
    blocks = by_cluster(values, clusters)
    k = len(blocks)
    if k < 3:
        return None
    out = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, k, k)
        out[i] = gini(np.concatenate([blocks[j] for j in pick]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def boot_cdf(values, clusters, grid, rng, n_boot=N_BOOT):
    """Cluster bootstrap band for an empirical CDF evaluated on `grid`."""
    blocks = by_cluster(values, clusters)
    k = len(blocks)
    if k < 3:
        return None
    reps = np.empty((n_boot, len(grid)))
    for i in range(n_boot):
        pick = rng.integers(0, k, k)
        s = np.sort(np.concatenate([blocks[j] for j in pick]))
        reps[i] = np.searchsorted(s, grid, side="right") / len(s)
    return np.percentile(reps, 2.5, axis=0), np.percentile(reps, 97.5, axis=0)


def boot_median(values, clusters, rng, n_boot=N_BOOT):
    """Median CI. Not drawn any more, but still reported in the caption."""
    blocks = by_cluster(values, clusters)
    k = len(blocks)
    if k < 3:
        return None
    out = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, k, k)
        out[i] = np.median(np.concatenate([blocks[j] for j in pick]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def runs_of(flags, weekno):
    """(first, last) week numbers of each consecutive-week run where flags is true."""
    out, cur = [], None
    for f, w in zip(flags, weekno):
        if f and cur is not None and w == cur[1] + 1:
            cur = (cur[0], w)
        elif f:
            if cur:
                out.append(cur)
            cur = (w, w)
        else:
            if cur:
                out.append(cur)
            cur = None
    if cur:
        out.append(cur)
    return out


# ------------------------------------------------------------------ drawing
def axes(x0, y0, w, h, xticks, yticks, xlab, ylab=None, xticks2=None):
    """Shared axis furniture, so every panel is read the same way.

    `xticks2` is a second, quieter row of tick labels -- used on the duration axis to
    name the hour / day / week / month / year landmarks a bare log scale hides.
    """
    o = []
    for v, lab in yticks:
        yy = y0 + h - v * h
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0, yy, x0 + w, yy))
        if lab:
            o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                     % (x0 - 5, yy + 3, lab))
    for v, lab in xticks:
        xx = x0 + v * w
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (xx, y0, xx, y0 + h))
        if lab:
            o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                     % (xx, y0 + h + 11, lab))
    lab_y = y0 + h + 24
    if xticks2:
        for v, lab in xticks2:
            o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" '
                     'fill="var(--muted)">%s</text>' % (x0 + v * w, y0 + h + 21, lab))
        lab_y = y0 + h + 34
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + h, x0 + w, y0 + h))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0, x0, y0 + h))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
             % (x0 + w / 2, lab_y, xlab))
    if ylab:
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" '
                 'transform="rotate(-90 %.1f %.1f)">%s</text>'
                 % (x0 - 28, y0 + h / 2, x0 - 28, y0 + h / 2, ylab))
    return "\n".join(o)


def swatch(kind, x, y, colour, dash=None):
    """One key mark. Line, bar or dot, matching how the series is actually drawn."""
    if kind == "bar":
        return ('<rect x="%.1f" y="%.1f" width="3" height="9" fill="%s"/>'
                % (x + 5, y - 8, colour))
    if kind == "dot":
        return ('<circle cx="%.1f" cy="%.1f" r="3.4" fill="%s"/>'
                % (x + 6.5, y - 3.5, colour))
    if kind == "band":
        return ('<rect x="%.1f" y="%.1f" width="13" height="8" fill="%s" '
                'opacity=".16"/>' % (x, y - 7.5, colour))
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
            'stroke-width="1.9"%s/>'
            % (x, y - 3.5, x + 13, y - 3.5, colour,
               "" if not dash else ' stroke-dasharray="%s"' % dash))


def key_row(items, x, y, colw, ncol=2, step=13.0):
    """A compact key: mark plus name, nothing else. The numbers are in the caption."""
    o = []
    for i, it in enumerate(items):
        label, kind, colour = it[0], it[1], it[2]
        dash = it[3] if len(it) > 3 else None
        cx, cy = x + (i % ncol) * colw, y + (i // ncol) * step
        o.append(swatch(kind, cx, cy, colour, dash))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>'
                 % (cx + 18, cy, label))
    return "\n".join(o)


def ellipse_network(nodes, edges, x0, y0, w, h, cov, node_col, edge_col,
                    label_all=True, nodes_off=False):
    """Units on an ellipse, edge width by weight, node area by collar coverage.

    An ellipse rather than a circle so this panel fills the same box as the others.
    `nodes_off` draws chords only, for overlaying a second edge type on a layout whose
    nodes are already on the canvas; those chords bow harder so the two stay separable.
    """
    o = []
    cx, cy = x0 + w / 2, y0 + h / 2
    rx, ry = w / 2 - 52, h / 2 - 34
    pos = {}
    for i, u in enumerate(nodes):
        a = -math.pi / 2 + 2 * math.pi * i / len(nodes)
        pos[u] = (cx + rx * math.cos(a), cy + ry * math.sin(a), a)
    wmax = max((wt for _, _, wt in edges), default=1.0)
    for u, v, wt in sorted(edges, key=lambda e: e[2]):
        if u not in pos or v not in pos:
            continue
        x1, y1, _ = pos[u]
        x2, y2, _ = pos[v]
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        pull = 0.62 if nodes_off else 0.32
        qx, qy = cx + (mx - cx) * pull, cy + (my - cy) * pull
        sw = 0.5 + 2.9 * (wt / wmax) ** 0.55
        o.append('<path d="M %.1f %.1f Q %.1f %.1f %.1f %.1f" fill="none" stroke="%s" '
                 'stroke-width="%.2f" opacity=".55"/>'
                 % (x1, y1, qx, qy, x2, y2, edge_col, sw))
    if nodes_off:
        return "\n".join(o)
    cmax = max(cov.values()) if cov else 1
    for u in nodes:
        x, y, a = pos[u]
        r = 2.3 + 4.0 * math.sqrt(max(cov.get(u, 0), 0) / cmax)
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s"/>'
                 % (x, y, r, node_col))
        if label_all:
            left = math.cos(a) < 0
            d = r + 4.0
            lx, ly = x + d * math.cos(a), y + d * math.sin(a) + 2.2
            rot = math.degrees(a) + (180 if left else 0)
            o.append('<text class="nl" x="%.1f" y="%.1f" text-anchor="%s" '
                     'transform="rotate(%.1f %.1f %.1f)">%s</text>'
                     % (lx, ly, "end" if left else "start", rot, lx, ly, u))
    return "\n".join(o)


# ------------------------------------------------------------------ data
def load_modularity():
    wk = pd.read_csv(FROZEN, parse_dates=["period_start"])
    wc = wk[wk["n_animals_observed"].ge(WELL_COVERED)].copy()
    counts = wc.groupby("dynamic_social_unit").size()
    units = sorted(counts[counts >= MIN_WEEKS].index)
    sub = wc[wc["dynamic_social_unit"].isin(units)].copy()
    origin = sub["period_start"].min()
    sub["weekno"] = ((sub["period_start"] - origin).dt.days // 7).astype(int)
    return sub, units


def load_networks():
    wk = pd.read_csv(FROZEN)
    cov = wk.groupby("dynamic_social_unit")["n_animals_observed"].median().to_dict()
    opp = pd.read_csv(OPP, parse_dates=["period_start"])
    enc = opp[opp["state"].eq("detected_encounter")].copy()
    enc["min_side"] = enc[["max_n_a_in_cluster", "max_n_b_in_cluster"]].min(axis=1)
    grp = enc[enc["min_side"].ge(2)].copy()
    solo = enc[enc["min_side"].eq(1)].copy()
    # keep the table's own pair_key rather than rebuilding one, so the solo and group
    # dyad sets can actually be differenced
    dy = (grp.groupby(["unit_a", "unit_b", "pair_key"]).size()
          .rename("days").reset_index())
    ed = pd.read_csv(EDGES)
    deg = {}
    for r in dy.itertuples():
        deg[r.unit_a] = deg.get(r.unit_a, 0) + 1
        deg[r.unit_b] = deg.get(r.unit_b, 0) + 1
    nodes = sorted(set(deg) | set(ed["origin_group"]) | set(ed["destination"]))
    nodes = sorted(nodes, key=lambda u: (-deg.get(u, 0), u))
    return nodes, cov, grp, solo, dy, ed, deg


def encounter_bouts():
    """Group-encounter bouts at hourly resolution, from derive_level2_inputs.py.

    A dyad-day table floors every encounter at 24 h, which is not a fact about baboons
    but about the table: two units sharing a cluster is an hourly observation, and 95%
    of these bouts turn out to be shorter than a day. Elapsed hours, gaps of up to two
    hours bridged, both units contributing two or more animals throughout.
    """
    b = pd.read_csv(BOUTS_H, parse_dates=["start", "end"])
    b["pair_key"] = ["|".join(sorted((x, y)))
                     for x, y in zip(b["unit_a"], b["unit_b"])]
    return b


def modular_phases(sub, units):
    """Every run of consecutive modular weeks, with the unit that owns it."""
    rows = []
    for u in units:
        g = sub[sub["dynamic_social_unit"].eq(u)].sort_values("period_start")
        flags = (g["modularity"] > MOD_EPS).to_numpy()
        for a, b in runs_of(flags, g["weekno"].to_numpy()):
            rows.append({"unit": u, "weeks": b - a + 1,
                         "hours": WEEK_H * (b - a + 1)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ panels
def panel_lilac(sub, x0, y0, w, h):
    """One group's modularity through time -- the dynamic a per-unit number hides."""
    g = sub[sub["dynamic_social_unit"].eq(EXAMPLE_UNIT)].sort_values("period_start")
    lo, hi = g["period_start"].min(), g["period_start"].max()
    span = max(1.0, (hi - lo).days)
    top = max(0.05, float(g["modularity"].max()))
    o = [axes(x0, y0, w, h, [(0.0, ""), (0.5, ""), (1.0, "")],
              [(0.0, "0"), (0.25, ""), (0.5, "%.2f" % (top / 2)), (0.75, ""),
               (1.0, "%.2f" % top)],
              "week", "modularity")]
    for tv, tl, ta in ((0.0, str(lo.date()), "start"), (1.0, str(hi.date()), "end")):
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="%s">%s</text>'
                 % (x0 + tv * w, y0 + h + 11, ta, tl))
    mod = (g["modularity"] > MOD_EPS).to_numpy()
    wkno = g["weekno"].to_numpy()
    wk2date = dict(zip(wkno, g["period_start"]))
    phases = [(a, b) for a, b in runs_of(mod, wkno) if b - a + 1 >= 2]
    for a, b in phases:
        xa = x0 + ((wk2date[a] - lo).days / span) * w
        xb = x0 + ((wk2date[b] - lo).days / span) * w + 2.6
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--c4)" '
                 'opacity=".11"/>' % (xa - 1.3, y0, max(2.8, xb - xa), h))
    pts = []
    for r in g.itertuples():
        xx = x0 + ((r.period_start - lo).days / span) * w
        pts.append("%.1f,%.1f" % (xx, y0 + h - float(r.largest_community_fraction) * h))
        hgt = (r.modularity / top) * h
        if hgt <= 0.5:
            o.append('<circle cx="%.1f" cy="%.1f" r="1" fill="var(--muted)"/>'
                     % (xx, y0 + h))
        else:
            o.append('<rect x="%.1f" y="%.1f" width="2.6" height="%.1f" '
                     'fill="var(--c4)"/>' % (xx - 1.3, y0 + h - hgt, hgt))
    o.append('<polyline points="%s" fill="none" stroke="var(--c6)" stroke-width="1.4" '
             'opacity=".9"/>' % " ".join(pts))
    # the teal line runs on its own 0-to-1 scale, so it gets its own axis rather than
    # borrowing the modularity one and quietly meaning something else
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
             'stroke="var(--c6)" opacity=".55"/>' % (x0 + w, y0, x0 + w, y0 + h))
    for tv, tl in ((0.0, "0"), (0.5, "0.5"), (1.0, "1")):
        o.append('<text class="ts" x="%.1f" y="%.1f" fill="var(--c6)">%s</text>'
                 % (x0 + w + 4, y0 + h - tv * h + 3, tl))
    o.append(key_row([("weekly modularity", "bar", "var(--c4)"),
                      ("modular phase, 2+ weeks", "band", "var(--c4)"),
                      ("largest sub-community", "line", "var(--c6)")],
                     x0 - 28, y0 + h + 46, colw=w / 2 + 22, ncol=2))
    longest = max((b - a + 1 for a, b in phases), default=0)
    nc = g["n_communities"]
    return "\n".join(o), {
        "unit": EXAMPLE_UNIT, "weeks": int(len(g)), "modular_weeks": int(mod.sum()),
        "max_mod": round(top, 3), "phases": len(phases),
        "longest_phase": int(longest), "max_communities": int(nc.max()),
        "min_largest": round(float(g["largest_community_fraction"].min()), 2),
        "window": [str(lo.date()), str(hi.date())]}


def panel_network(nodes, cov, dy, ed, deg, solo, grp, with_d, x0, y0, w, h):
    """Both between-unit networks on one layout: encounters and individual transfers.

    They share nodes, so two circles force a comparison the reader has to make by eye.
    Overlaid -- and with the single-animal dyad-days taken out of the encounter side, so
    the two are not partly the same data -- a dyad with a blue chord and no plum one is a
    pair that meets without exchanging animals, and the reverse is a pair that exchanges
    animals without their groups ever being recorded together.
    """
    enc_edges = [(r.unit_a, r.unit_b, float(r.hours)) for r in dy.itertuples()]
    tr_edges = [(r.origin_group, r.destination, float(r.excursions))
                for r in ed.itertuples()]
    net = ellipse_network(nodes, enc_edges, x0, y0, w, h, cov,
                          "var(--ink-2)", "var(--c1)")
    net += "\n" + ellipse_network(nodes, tr_edges, x0, y0, w, h, cov,
                                  "var(--ink-2)", "var(--c5)", nodes_off=True)
    net += "\n" + key_row([("group encounter", "line", "var(--c1)"),
                           ("animal transfer", "line", "var(--c5)"),
                           ("node area = collars", "dot", "var(--ink-2)")],
                          x0 - 28, y0 + h + 46, colw=w / 2 + 22, ncol=2)
    enc_pairs = {frozenset((r.unit_a, r.unit_b)) for r in dy.itertuples()}
    tr_pairs = {frozenset((r.origin_group, r.destination)) for r in ed.itertuples()}
    recip = sum(1 for r in ed.itertuples()
                if ((ed["origin_group"] == r.destination)
                    & (ed["destination"] == r.origin_group)).any())
    hrs = dy["hours"].to_numpy(dtype=float)
    return net, {
        "units": len(nodes), "dyads": int(len(dy)),
        "degree_median": float(np.median(list(deg.values()))),
        "degree_max": int(max(deg.values())),
        "hours_median": float(np.median(hrs)), "hours_max": int(hrs.max()),
        "bouts": int(dy["bouts"].sum()),
        "top5_share": round(float(np.sort(hrs)[::-1][:5].sum() / hrs.sum()), 3),
        "solo_dyad_days": int(len(solo)),
        "solo_share": round(float(len(solo) / (len(solo) + len(grp))), 3),
        "solo_dyads": int(solo["pair_key"].nunique()),
        "solo_only_dyads": int(len(set(solo["pair_key"]) - set(dy["pair_key"]))),
        "solo_distance": int(round(float(solo["min_centroid_distance_m"].median()))),
        "grp_distance": int(round(float(grp["min_centroid_distance_m"].median()))),
        "grp_dyad_days": int(len(grp)),
        "transfer_edges": int(len(ed)), "reciprocal_edges": int(recip),
        "meet_only": int(len(enc_pairs - tr_pairs)),
        "exchange_only": int(len(tr_pairs - enc_pairs)),
        "both": int(len(enc_pairs & tr_pairs)),
        "animals": int(with_d["animal_id"].nunique()),
        "excursions": int(len(with_d)),
        "max_edge": int(ed["excursions"].max()),
        "top_edge": "%s to %s" % (ed.loc[ed["excursions"].idxmax(), "origin_group"],
                                  ed.loc[ed["excursions"].idxmax(), "destination"])}


def panel_conc(series, x0, y0, w, h, rng):
    """Every concentration curve on one axis, so the Ginis are directly comparable.

    The Ginis and their intervals sit in the lower right of the plot, which no Lorenz
    curve can ever reach, rather than in a separate block of text.
    """
    o = [axes(x0, y0, w, h, [(0.0, "0"), (0.5, "0.5"), (1.0, "1")],
              [(0.0, "0"), (0.5, "0.5"), (1.0, "1")],
              "share of units / dyads / edges / animals, ranked",
              "share of the total")]
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--muted)" '
             'stroke-width=".9" stroke-dasharray="4 3"/>' % (x0, y0 + h, x0 + w, y0))
    stats, items, ty = {}, [], y0 + h - 44.0
    for label, key, vals, clus, colour, dash in series:
        v = np.asarray(vals, dtype=float)
        xs, ys = lorenz(v)
        pts = " ".join("%.1f,%.1f" % (x0 + a * w, y0 + h - b * h)
                       for a, b in zip(xs, ys))
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.9"%s/>'
                 % (pts, colour, "" if not dash else ' stroke-dasharray="%s"' % dash))
        ci = boot_gini(v, clus, rng)
        stats[key] = {"gini": round(gini(v), 2), "n": int((v > 0).sum()),
                      "ci": None if ci is None else [round(c, 2) for c in ci]}
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end" fill="%s">%s</text>'
                 % (x0 + w - 4, ty, colour,
                    "%.2f [%.2f, %.2f]" % (stats[key]["gini"], ci[0], ci[1]) if ci
                    else "%.2f" % stats[key]["gini"]))
        ty += 10.5
        items.append((label, "line", colour, dash))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end" '
             'fill="var(--muted)">Gini [95%% CI]</text>' % (x0 + w - 4, y0 + h - 55))
    o.append(key_row(items, x0 - 28, y0 + h + 46, colw=w / 2 + 22, ncol=2))
    return "\n".join(o), stats


def panel_duration(curves, x0, y0, w, h, rng):
    """One duration axis for all three structures, in hours. Log x: the tails matter."""
    lmax = math.log10(DUR_MAX_H)

    def sx(v):
        return x0 + (math.log10(max(1.0, float(v))) / lmax) * w

    decades = [(10 ** k, format(10 ** k, ",")) for k in range(5)]
    marks = [(DAY_H, "day"), (WEEK_H, "week"), (MONTH_H, "month"), (YEAR_H, "year")]
    o = [axes(x0, y0, w, h,
              [(math.log10(t) / lmax, lab) for t, lab in decades],
              [(0.0, "0"), (0.25, ""), (0.5, "0.5"), (0.75, ""), (1.0, "1")],
              "hours in the state (log scale)", "share at or below",
              xticks2=[(math.log10(t) / lmax, lab) for t, lab in marks])]
    # minor ticks make a log axis legible as a log axis; the calendar landmarks get a
    # ruled line each so the word underneath attaches to a position, not to a gap
    for k in range(5):
        for m in range(2, 10):
            v = m * 10 ** k
            if v <= DUR_MAX_H:
                o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                         'opacity=".55"/>' % (sx(v), y0 + h, sx(v), y0 + h + 3))
    for t, _ in marks:
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--muted)" '
                 'stroke-width=".8" stroke-dasharray="2 3" opacity=".8"/>'
                 % (sx(t), y0, sx(t), y0 + h + 14))
    grid = np.logspace(0, lmax, 100)
    # a one-hour state is not a one-hour event: half the split and excursion bouts sit
    # at exactly 1 h, which is the hourly assignment flickering, not an animal leaving
    # and returning. Shade the floor rather than filter it away silently.
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--muted)" '
             'opacity=".07"/>' % (x0, y0, sx(FLICKER_H) - x0, h))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" '
             'fill="var(--muted)" transform="rotate(-90 %.1f %.1f)">flicker</text>'
             % ((x0 + sx(FLICKER_H)) / 2, y0 + h * 0.62,
                (x0 + sx(FLICKER_H)) / 2, y0 + h * 0.62))
    stats, items, bands, lines = {}, [], [], []
    for label, key, noun, vals, clus, colour in curves:
        v = np.asarray(vals, dtype=float)
        v = np.sort(v[v >= 1])
        n = len(v)
        band = boot_cdf(v, clus, grid, rng)
        if band is not None:
            up = ["%.1f,%.1f" % (sx(g), y0 + h - b * h) for g, b in zip(grid, band[1])]
            dn = ["%.1f,%.1f" % (sx(g), y0 + h - b * h)
                  for g, b in zip(grid[::-1], band[0][::-1])]
            bands.append('<polygon points="%s" fill="%s" opacity=".15" stroke="none"/>'
                         % (" ".join(up + dn), colour))
        pts = ["%.1f,%.1f" % (sx(v[0]), y0 + h)]
        for i, x in enumerate(v):
            pts.append("%.1f,%.1f" % (sx(x), y0 + h - (i / n) * h))
            pts.append("%.1f,%.1f" % (sx(x), y0 + h - ((i + 1) / n) * h))
        lines.append('<polyline points="%s" fill="none" stroke="%s" '
                     'stroke-width="1.7"/>' % (" ".join(pts), colour))
        mci = boot_median(v, clus, rng)
        stats[key] = {"n": n, "clusters": int(len(pd.unique(np.asarray(clus)))),
                      "noun": noun, "median": float(np.median(v)),
                      "median_ci": None if mci is None else [round(c, 1) for c in mci],
                      "p90": float(np.percentile(v, 90)), "max": float(v.max()),
                      "under_day": round(float((v < DAY_H).mean()), 3),
                      "under_6h": round(float((v < 6).mean()), 3),
                      "at_1h": round(float((v <= 1).mean()), 3),
                      "over_day": round(float((v >= 24).mean()), 3),
                      "over_week": round(float((v >= WEEK_H).mean()), 3),
                      "over_month": round(float((v >= 720).mean()), 3)}
        items.append((label, "line", colour))
    o.extend(bands)
    o.extend(lines)
    o.append(key_row(items, x0 - 28, y0 + h + 46, colw=w / 2 + 22, ncol=2))
    return "\n".join(o), stats


def build():
    rng = np.random.default_rng(SEED)
    sub, units = load_modularity()
    nodes, cov, grp, solo, dy, ed, deg = load_networks()
    dests = pd.read_csv(DESTS)
    with_d = dests[dests["destination"].notna()]
    phases = modular_phases(sub, units)
    bouts = encounter_bouts()
    dyh = (bouts.groupby(["unit_a", "unit_b", "pair_key"])
           .agg(hours=("hours", "sum"), bouts=("hours", "size")).reset_index())
    exc = pd.read_csv(EXC)
    per_unit = sub.groupby("dynamic_social_unit")["modularity"].apply(
        lambda x: float((x > MOD_EPS).sum()))
    per_animal = with_d.groupby("animal_id").size()

    pa, sa = panel_lilac(sub, COL_X[0], ROW_Y[0], BOX_W, BOX_H)
    pb, sb = panel_network(nodes, cov, dyh, ed, deg, solo, grp, with_d,
                           COL_X[1], ROW_Y[0], BOX_W, BOX_H)
    pc, gs = panel_conc(
        [("modular weeks per unit", "mod", per_unit.to_numpy(),
          per_unit.index.to_numpy(), "var(--c4)", None),
         ("encounter-hours per dyad", "enc", dyh["hours"].to_numpy(dtype=float),
          dyh["pair_key"].to_numpy(), "var(--c1)", None),
         ("excursions per transfer edge", "edge", ed["excursions"].to_numpy(dtype=float),
          (ed["origin_group"] + ">" + ed["destination"]).to_numpy(), "var(--c5)", None),
         ("excursions per animal", "animal", per_animal.to_numpy(dtype=float),
          per_animal.index.to_numpy(), "var(--c5)", "4 2.5")],
        COL_X[0], ROW_Y[1], BOX_W, BOX_H, rng)
    # the excursion curve is held to the audited axis-C cohort, so panel (d) counts the
    # same animals the individual axis does everywhere else in the paper
    exch = pd.read_csv(EXC_H)
    exch = exch[exch["animal_id"].isin(set(exc["animal_id"]))]
    splits = pd.read_csv(SPLITS_H)
    pdur, ds = panel_duration(
        [("group encounters", "enc", "dyads", bouts["hours"].to_numpy(dtype=float),
          bouts["pair_key"].to_numpy(), "var(--c1)"),
         ("within-group splits", "mod", "units", splits["hours"].to_numpy(dtype=float),
          splits["unit"].to_numpy(), "var(--c4)"),
         ("individual excursions", "exc", "animals",
          exch["hours"].to_numpy(dtype=float),
          exch["animal_id"].to_numpy(), "var(--c5)")],
        COL_X[1], ROW_Y[1], BOX_W, BOX_H, rng)
    H = ROW_Y[1] + BOX_H + 86

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Four panels on one grid: '
         'weekly modularity for one well-covered group, the between-unit encounter and '
         'transfer networks overlaid on one node layout, four concentration curves on '
         'shared axes with bootstrap intervals, and the duration distributions of all '
         'three structures on one logarithmic hours axis with bootstrap bands.">' % H]
    for letter, ci, ri in (("a", 0, 0), ("b", 1, 0), ("c", 0, 1), ("d", 1, 1)):
        o.append('<text class="pl" x="%.0f" y="%.1f">%s</text>'
                 % (COL_X[ci] - 44, ROW_Y[ri] - 16, letter))
    o.extend([pa, pb, pc, pdur])
    o.append("</svg>")
    return "\n".join(o), sa, sb, gs, ds


CSS = """
:root{--ground:#f4f5f3;--paper:#fff;--ink:#15191a;--ink-2:#3a4244;--muted:#6b7476;
--rule:#d5dad7;--rule-soft:#e6eae7;
--c1:#1f6f8b;--c4:#8a5a1f;--c5:#6b4a7a;--c6:#0f5f57;
--sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;
--serif:"Newsreader",Georgia,serif;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0f1213;
--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;--muted:#8f9998;--rule:#2c3435;
--rule-soft:#242c2d;--c1:#74b6d0;--c4:#d3a061;--c5:#b494c4;--c6:#5fc0b0;}}
:root[data-theme="dark"]{--ground:#0f1213;--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;
--muted:#8f9998;--rule:#2c3435;--rule-soft:#242c2d;--c1:#74b6d0;--c4:#d3a061;
--c5:#b494c4;--c6:#5fc0b0;}
*{box-sizing:border-box}
body{background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:15px;
line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:800px;margin:0 auto;padding:48px 26px 80px}
.kicker{font-family:var(--mono);font-size:10px;letter-spacing:.16em;
text-transform:uppercase;color:var(--muted);margin-bottom:10px}
h1{font-family:var(--serif);font-size:30px;font-weight:500;line-height:1.15;
letter-spacing:-.01em;margin:0 0 26px;text-wrap:balance}
.plate{background:var(--paper);border:1px solid var(--rule-soft);padding:22px 20px 8px}
.plate svg{display:block;width:100%;height:auto;overflow:visible}
.scroll{overflow-x:auto}.scroll>svg{min-width:660px}
.cap{font-size:12.5px;line-height:1.62;color:var(--ink-2);margin:16px 0 0;
padding-top:14px;border-top:1px solid var(--rule)}
.cap b{font-weight:600;color:var(--ink)}.cap .fignum{font-weight:600;color:var(--ink)}
.cap code{font-family:var(--mono);font-size:.9em}
.note{font-size:11.5px;line-height:1.6;color:var(--muted);margin-top:22px;
padding-top:14px;border-top:1px solid var(--rule-soft)}
.note code{font-family:var(--mono);font-size:.92em}
text{font-family:var(--sans)}
.nl{font-size:7.6px;fill:var(--ink-2)}
.ts{font-size:8.4px;fill:var(--muted)}
.ph{font-size:12px;fill:var(--ink);font-weight:600}
.pl{font-size:13px;font-weight:700;fill:var(--ink)}
.g{stroke:var(--rule-soft);stroke-width:.8}
.ax{stroke:var(--rule);stroke-width:1}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""



def page():
    fig, sa, sb, gs, ds = build()
    return f"""<title>Level Two Structures</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 2 &middot; draft</div>
  <h1>What the events add up to</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 2.</span> The three structures the level-one events
      aggregate into &mdash; modularity within a unit, encounters between units, animals
      moving between them &mdash; as extent (top row) and as variation (bottom row).
      Colour is constant throughout: <b>ochre</b> within-group modularity, <b>blue</b>
      between-group encounters, <b>plum</b> individual transfer, <b>teal</b> the
      largest-sub-community trace in (a). All intervals are 95% cluster bootstraps
      ({N_BOOT} replicates) on the natural unit of replication &mdash; dyad, unit or
      animal &mdash; never on the row.
      <b>(a)</b> {sa["unit"]}, the unit carrying most of the modularity in this
      population, across its {sa["weeks"]} weeks at {WELL_COVERED} or more collars
      ({sa["window"][0]} to {sa["window"][1]}). Bars are weekly modularity, to a maximum
      of {sa["max_mod"]:.2f}; {sa["modular_weeks"]} weeks are modular, falling in
      {sa["phases"]} shaded phases of two weeks or more, the longest
      {sa["longest_phase"]} consecutive weeks. Up to {sa["max_communities"]} communities
      appear, and the teal trace &mdash; the largest of them as a share of the group, on
      the right-hand axis &mdash; falls to {sa["min_largest"]:.0%}. The group does not
      drift apart and stay apart: it separates and re-fuses, on a scale of weeks.
      <b>(b)</b> Both between-unit networks on one node layout: {sb["units"]} units,
      node area by median collar count so a five-collar node is not read as a
      five-animal group. Blue chords are group encounters ({sb["dyads"]} dyads,
      {sb["bouts"]:,} bouts), plum chords animals moving
      ({sb["transfer_edges"]} directed edges, {sb["excursions"]} excursions by
      {sb["animals"]} animals). The encounter side requires <b>two or more animals on
      both sides</b>: {sb["solo_dyad_days"]} dyad-days ({sb["solo_share"]:.0%} of all
      detected encounters, across {sb["solo_dyads"]} dyads, {sb["solo_only_dyads"]} of
      which qualify no other way) had a single animal on one side, and those are one
      animal visiting, not two groups meeting &mdash; their group centroids sit a median
      {sb["solo_distance"]:,} m apart against {sb["grp_distance"]} m for a group
      encounter. Separating them is what makes the overlay informative:
      {sb["meet_only"]} dyads meet without ever exchanging an animal,
      {sb["exchange_only"]} exchange without a recorded group encounter, {sb["both"]} do
      both. {sb["reciprocal_edges"]} transfer edges are reciprocal; the heaviest,
      {sb["top_edge"]}, carries {sb["max_edge"]} excursions.
      <b>(c)</b> Concentration, all four quantities against one diagonal, with Gini and
      its interval printed in the panel. Modular weeks per unit is the most concentrated
      (Gini {gs["mod"]["gini"]:.2f}) and encounter-hours per dyad next
      ({gs["enc"]["gini"]:.2f}, across {gs["enc"]["n"]} dyads); excursions spread more
      evenly, over {gs["edge"]["n"]} transfer edges ({gs["edge"]["gini"]:.2f}) and
      {gs["animal"]["n"]} animals ({gs["animal"]["gini"]:.2f}, dashed). The interval on
      modularity is wide because only {gs["mod"]["n"]} units are above zero, which is
      itself the honest reading: seven units cannot pin down a Gini.
      <b>(d)</b> Duration, all three states at hourly resolution on one logarithmic
      axis with bootstrap bands; decade ticks carry the numbers and the dashed rules mark
      day, week, month and year. This is where a nightly or weekly table misleads:
      <b>{ds["enc"]["under_day"]:.0%} of group encounters end within a day</b> and
      {ds["enc"]["under_6h"]:.0%} within six hours, median {ds["enc"]["median"]:.0f} h
      [{ds["enc"]["median_ci"][0]:.0f}, {ds["enc"]["median_ci"][1]:.0f}]; the longest,
      {ds["enc"]["max"]:,.0f} h, is one pair of units that stopped separating.
      Within-group splits are shorter still ({ds["mod"]["n"]:,} bouts in
      {ds["mod"]["clusters"]} units, median {ds["mod"]["median"]:.0f} h,
      {ds["mod"]["under_day"]:.0%} inside a day) &mdash; the hourly face of the weekly
      modularity in (a), and the reason a group can be modular for seven consecutive
      weeks without ever being apart for a whole day. Excursions still carry the longest
      tail (median {ds["exc"]["median"]:.0f} h [{ds["exc"]["median_ci"][0]:.0f},
      {ds["exc"]["median_ci"][1]:.0f}], {ds["exc"]["over_week"]:.0%} past a week, out to
      {ds["exc"]["max"]:,.0f} h), but read this curve for what it is: an uninterrupted
      away-state run, not a whole dispersal. A months-long absence is cut here every time
      the animal&rsquo;s hourly context returns to its origin, so the audited nightly view
      &mdash; 338 excursions, 22% past a week, out to 500 nights &mdash; remains the
      measure of how long leaving lasts. What the hourly view adds is the bottom of the
      distribution, which a nightly rule cannot see at all.
      The shaded strip below {FLICKER_H:.0f} h is the resolution floor:
      {ds["mod"]["at_1h"]:.0%} of splits and {ds["exc"]["at_1h"]:.0%} of excursions are
      single hours, against {ds["enc"]["at_1h"]:.0%} of encounters, because one
      animal&rsquo;s cluster membership flickers where a whole group&rsquo;s
      co-occurrence does not.
      <b>The three differ not only in how often they happen but in what kind of thing
      they are:</b> a split is a within-day rearrangement, an encounter an event of
      hours, and an excursion either an event or a permanent change of group.
    </p>
  </div>

  <p class="note">
    Panel (a) is restricted to well-covered unit-weeks, because below {WELL_COVERED}
    collars a modularity estimator returns a single community regardless of the animals.
    All three curves in (d) are derived at hourly resolution by
    <code>derive_level2_inputs.py</code> from the same export, gaps of up to two hours
    bridged because 02:00 is a known coverage hole: an encounter bout is contiguous hours
    in which both units had two or more animals in one association cluster
    ({ds["enc"]["n"]:,} bouts, {ds["enc"]["clusters"]} dyads); a split is contiguous
    hours in which a unit with {WELL_COVERED} or more animals present occupied two or
    more clusters of two or more animals each ({ds["mod"]["n"]:,} bouts,
    {ds["mod"]["clusters"]} units); an excursion is contiguous hours of a non-origin
    social context ({ds["exc"]["n"]:,} bouts), held to the {ds["exc"]["clusters"]}
    animals of the audited axis-C cohort so that the panel counts the same animals the
    individual axis does elsewhere &mdash; unrestricted, the hourly excursion state fires
    for 319 animals, which is the sparse-collar population the located-reference-cluster
    rule exists to exclude. Panels (b) and (c) are not
    coverage-restricted, which is why node area carries coverage instead. Transfer
    destinations are resolved per away-night as the modal other unit sharing the focal
    animal&rsquo;s association clusters; 86.3% of joined away-nights resolve and the rest
    are left unassigned rather than guessed, so the plum network is a lower bound. All
    panels come from the frozen 1,924,104-row export. Generated by
    <code>build_level2_figure.py</code>; seed {SEED}.
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
