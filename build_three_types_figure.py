"""Manuscript figure: three event types, each with one continuous depth.

REPLACES the six-panel taxonomy. Two reasons the six-way scheme should go.

  1. The merge-scale subdivision (large / medium partial / small subset) rests on a
     classifier that keys on the MAXIMUM participation a unit ever reaches, not its
     typical state, so events labelled small-subset can have maximum fractions near 1.0.
  2. Figure 4 showed that subdivision does not distinguish groups: ICC 0.05-0.13 with
     permutation p = 0.09-0.26, and well-sampled groups traverse most of its range
     within their own record. It is a gradient, measured on a faulty rule, being
     reported as three categories.

So: three types, defined by WHO crosses, and one continuous depth per type, all on the
same 0-1 scale.

  BETWEEN-GROUP   two units share a cluster    depth = mutual participation,
                                               min(median frac A, median frac B)
  WITHIN-GROUP    one unit holds >1 cluster     depth = fraction of the unit outside
                                               its own largest cluster
  INDIVIDUAL      one animal leaves its unit    depth = fraction of away-nights spent
                                               with another unit

Each panel pairs a real worked example from the record with the depth distribution over
every event of that type, so the type reads as a category and the depth as a continuum.

Output: docs/framing_2026-09-04/three_types_figure.html
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
STAGE1 = Path("outputs/general_structure_2026_09/phase2_two_stage_events/"
              "stage1_events_with_stage2_mixing.csv")
EXC = Path("outputs/general_structure_2026_09/phase4b_individual_axis/"
           "excursions_dominant_gap0.csv")
OUT = Path("docs/framing_2026-09-04/three_types_figure.html")

MAX_ROWS = 20
MAX_COLS = 26
MIN_OBS_FOR_SPLIT = 8
WITH_ORIGIN = {"with_origin", "mixed_with_origin_present"}
AWAY_WITH_OTHERS = {"other", "mixed_without_origin_unit"}

FILL = {"own": "var(--n2)", "mixed": "var(--c1)", "alone": "var(--c6)",
        "main": "var(--n2)", "secondary": "var(--c4)",
        "origin": "var(--n2)", "other": "var(--c5)"}


# ------------------------------------------------------------------ depths
def between_depth(s1: pd.DataFrame) -> np.ndarray:
    v = s1[["median_frac_a", "median_frac_b"]].min(axis=1).dropna()
    return v[(v > 0) & (v <= 1)].to_numpy()


def within_depth(d: pd.DataFrame) -> np.ndarray:
    """Per unit-hour with more than one cluster: share of the unit outside the largest."""
    obs = d[d["is_observed"]]
    sizes = (obs.groupby(["window_start", "dynamic_social_unit",
                          "association_event_id"]).size().rename("k").reset_index())
    tot = sizes.groupby(["window_start", "dynamic_social_unit"]).agg(
        total=("k", "sum"), biggest=("k", "max"), n_clusters=("k", "size")).reset_index()
    tot = tot[(tot["n_clusters"] > 1) & (tot["total"] >= MIN_OBS_FOR_SPLIT)]
    return (1.0 - tot["biggest"] / tot["total"]).to_numpy()


def individual_depth(exc: pd.DataFrame) -> np.ndarray:
    f = exc[exc["away_nights"] > 0]
    return (f["other_nights"] / f["away_nights"]).to_numpy()


# ------------------------------------------------------------------ examples
def ex_between(d, s1):
    s = s1.copy()
    s["med_min"] = s[["median_frac_a", "median_frac_b"]].min(axis=1)
    s["obs_min"] = s[["mean_observed_a", "mean_observed_b"]].min(axis=1)
    c = s[s["structural_span_hours"].between(8, MAX_COLS)
          & s["observed_hour_fraction"].ge(0.8) & s["obs_min"].ge(6)
          & s["med_min"].between(0.35, 0.85)]
    if c.empty:
        c = s[s["obs_min"].ge(6)]
    r = c.sort_values("obs_min", ascending=False).iloc[0]
    ua, ub = r["unit_a"], r["unit_b"]
    t0 = pd.Timestamp(r["start_hour"])
    t1 = t0 + pd.Timedelta(hours=MAX_COLS - 1)
    w = d[d["window_start"].between(t0, t1)
          & d["dynamic_social_unit"].isin([ua, ub])]
    cols = sorted(w["window_start"].unique())[:MAX_COLS]
    w = w[w["window_start"].isin(cols)]
    mixed = set()
    for ts, g in w[w["is_observed"]].groupby("window_start"):
        for ev, gg in g.groupby("association_event_id"):
            if gg["dynamic_social_unit"].nunique() > 1:
                mixed.add((ts, ev))
    rows = []
    for u in (ua, ub):
        ids = (w[w["dynamic_social_unit"].eq(u) & w["is_observed"]]
               .groupby("animal_id").size().sort_values(ascending=False))
        rows += [(u, a) for a in list(ids.index)[:MAX_ROWS // 2]]
    cells = {}
    for t in w.itertuples():
        k = (t.dynamic_social_unit, t.animal_id)
        if k not in rows or not t.is_observed:
            continue
        cells[(k, t.window_start)] = ("mixed" if (t.window_start, t.association_event_id)
                                      in mixed else
                                      ("alone" if t.temp_group_size <= 1 else "own"))
    return {"cols": cols, "rows": rows, "cells": cells,
            "legend": [("own", "own unit only"), ("mixed", "both units"),
                       ("alone", "alone")],
            "meta": "%s, %s; depth %.2f" % (r["pair_key"], t0.date(), r["med_min"]),
            "xlab": "%d hours" % len(cols)}


def ex_within(d):
    obs = d[d["is_observed"] & d["window_start"].dt.year.ge(2025)]
    sizes = (obs.groupby(["window_start", "dynamic_social_unit",
                          "association_event_id"]).size().rename("k").reset_index())
    tot = sizes.groupby(["window_start", "dynamic_social_unit"]).agg(
        total=("k", "sum"), biggest=("k", "max"), n_clusters=("k", "size")).reset_index()
    tot["depth"] = 1.0 - tot["biggest"] / tot["total"]
    cand = tot[(tot["n_clusters"] > 1) & (tot["total"] >= 12)
               & tot["depth"].between(0.25, 0.6)]
    r = cand.sort_values("total", ascending=False).iloc[0]
    unit, centre = r["dynamic_social_unit"], pd.Timestamp(r["window_start"])
    t0 = centre - pd.Timedelta(hours=MAX_COLS // 2)
    t1 = t0 + pd.Timedelta(hours=MAX_COLS - 1)
    w = d[d["window_start"].between(t0, t1) & d["dynamic_social_unit"].eq(unit)]
    cols = sorted(w["window_start"].unique())[:MAX_COLS]
    w = w[w["window_start"].isin(cols)]
    biggest = {}
    for ts, g in w[w["is_observed"]].groupby("window_start"):
        biggest[ts] = g.groupby("association_event_id").size().idxmax()
    ids = w[w["is_observed"]].groupby("animal_id").size().sort_values(ascending=False)
    rows = [(unit, a) for a in list(ids.index)[:MAX_ROWS]]
    cells = {}
    for t in w.itertuples():
        k = (unit, t.animal_id)
        if k not in rows or not t.is_observed:
            continue
        cells[(k, t.window_start)] = ("alone" if t.temp_group_size <= 1 else
                                      ("main" if biggest.get(t.window_start)
                                       == t.association_event_id else "secondary"))
    return {"cols": cols, "rows": rows, "cells": cells,
            "legend": [("main", "largest cluster"), ("secondary", "secondary cluster"),
                       ("alone", "alone")],
            "meta": "%s, %s; depth %.2f" % (unit, centre.date(), r["depth"]),
            "xlab": "%d hours" % len(cols)}


def ex_individual(d, exc):
    c = exc[exc["away_nights"].between(5, 12)
            & (exc["other_nights"] / exc["away_nights"]).between(0.3, 0.9)]
    if c.empty:
        c = exc[exc["depth"].eq("alone_and_joined")]
    r = c.sort_values("away_nights").iloc[len(c) // 2]
    animal, unit = r["animal_id"], r["origin_group"]
    t0 = pd.Timestamp(r["start_night"]) - pd.Timedelta(days=3)
    t1 = pd.Timestamp(r["end_night"]) + pd.Timedelta(days=4)
    w = d[d["window_start"].between(t0, t1)
          & d["dynamic_social_unit"].eq(unit)].copy()
    w["night"] = (w["window_start"] - pd.Timedelta(hours=10)).dt.normalize()
    cols = sorted(w["night"].unique())[:MAX_COLS]
    w = w[w["night"].isin(cols)]
    ids = w[w["is_observed"]].groupby("animal_id").size().sort_values(ascending=False)
    others = [a for a in ids.index if a != animal][:MAX_ROWS - 1]
    rows = [(unit, animal)] + [(unit, a) for a in others]
    cells = {}
    for (a, n), g in w.groupby(["animal_id", "night"]):
        k = (unit, a)
        if k not in rows:
            continue
        g = g[g["is_observed"]]
        if g.empty:
            continue
        sc = g["social_context"]
        v = ("other" if sc.isin(AWAY_WITH_OTHERS).any() else
             ("alone" if sc.eq("isolated").any() else
              ("origin" if sc.isin(WITH_ORIGIN).any() else None)))
        if v:
            cells[(k, n)] = v
    dep = float(r["other_nights"]) / float(r["away_nights"])
    return {"cols": cols, "rows": rows, "cells": cells, "focal": (unit, animal),
            "legend": [("origin", "with origin group"), ("alone", "alone"),
                       ("other", "another unit")],
            "meta": "%s of %s, %s; depth %.2f"
                    % (animal.split("_")[0], unit,
                       pd.Timestamp(r["start_night"]).date(), dep),
            "xlab": "%d nights" % len(cols)}


# ------------------------------------------------------------------ drawing
def raster(g, x0, y0, w):
    cols, rows, cells = g["cols"], g["rows"], g["cells"]
    lab_w = 52.0
    cw = min(8.0, (w - lab_w) / max(1, len(cols)))
    ch = 6.6
    o, prev = [], None
    for i, key in enumerate(rows):
        unit, animal = key
        yy = y0 + i * ch
        if prev is not None and unit != prev:
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
                     'stroke-width=".8"/>' % (x0, yy - 1.3,
                                              x0 + lab_w + cw * len(cols), yy - 1.3))
        prev = unit
        cls = "rlf" if g.get("focal") == key else "rl"
        o.append('<text class="%s" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (cls, x0 + lab_w - 4, yy + ch - 1.5, animal.split("_")[0]))
        for j, c in enumerate(cols):
            v = cells.get((key, c))
            xx = x0 + lab_w + j * cw
            if v is None:
                o.append('<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" fill="none" '
                         'stroke="var(--rule-soft)" stroke-width=".4"/>'
                         % (xx, yy, cw - .6, ch - .9))
            else:
                o.append('<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" fill="%s"/>'
                         % (xx, yy, cw - .6, ch - .9, FILL[v]))
    ybot = y0 + len(rows) * ch
    seen = []
    for i, (u, _) in enumerate(rows):
        if u not in [s[0] for s in seen]:
            seen.append((u, i))
    for u, i in seen:
        n = sum(1 for uu, _ in rows if uu == u)
        o.append('<text class="ul" transform="translate(%.1f %.1f) rotate(-90)" '
                 'text-anchor="middle">%s</text>' % (x0 - 4, y0 + (i + n / 2) * ch, u))
    right = x0 + lab_w + cw * len(cols)
    lx, ly = x0 + lab_w, ybot + 7
    for v, lab in list(g["legend"]) + [(None, "not observed")]:
        iw = 9 + len(lab) * 4.3 + 12
        if lx + iw > right and lx > x0 + lab_w:
            lx, ly = x0 + lab_w, ly + 10
        if v is None:
            o.append('<rect x="%.1f" y="%.1f" width="6" height="6" fill="none" '
                     'stroke="var(--rule-soft)" stroke-width=".6"/>' % (lx, ly))
        else:
            o.append('<rect x="%.1f" y="%.1f" width="6" height="6" fill="%s"/>'
                     % (lx, ly, FILL[v]))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (lx + 9, ly + 5.5, lab))
        lx += iw
    o.append('<text class="ts" x="%.1f" y="%.1f">%s &middot; %s</text>'
             % (x0 + lab_w, ly + 17, g["xlab"], g["meta"]))
    return "\n".join(o), (ly + 22) - y0


def depth_hist(vals, x0, y0, w, h, colour, xlab, rng):
    """Histogram of the depth measure with a cluster-free bootstrap band on the median."""
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    nb = 20
    edges = np.linspace(0, 1, nb + 1)
    cnt, _ = np.histogram(v, bins=edges)
    frac = cnt / cnt.sum()
    o = []
    top = frac.max() * 1.12
    for t in (0, 0.5, 1.0):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0 + t * w, y0, x0 + t * w, y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (x0 + t * w, y0 + h + 11, "0" if t == 0 else
                    ("1" if t == 1 else "0.5")))
    bw = w / nb
    for i, f in enumerate(frac):
        if f <= 0:
            continue
        bh = f / top * h
        o.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s" '
                 'opacity=".8"/>' % (x0 + i * bw, y0 + h - bh, bw - .7, bh, colour))
    med = float(np.median(v))
    boot = [np.median(rng.choice(v, len(v), replace=True)) for _ in range(1000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="1.1" stroke-dasharray="3 2"/>'
             % (x0 + med * w, y0 - 4, x0 + med * w, y0 + h))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="2"/>' % (x0 + lo * w, y0 - 4, x0 + hi * w, y0 - 4))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" fill="var(--ink)">'
             'median %.2f</text>' % (x0 + med * w, y0 - 8, med))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + h, x0 + w, y0 + h))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
             % (x0 + w / 2, y0 + h + 23, xlab))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">n = %s</text>'
             % (x0 + w / 2, y0 + h + 34, format(len(v), ",")))
    return "\n".join(o)


TYPES = [
    ("Between-group", "two units share one cluster",
     "mutual participation: min(share A, share B)", "var(--c1)"),
    ("Within-group", "one unit holds more than one cluster",
     "share of the unit outside its largest cluster", "var(--c4)"),
    ("Individual", "one animal leaves its unit",
     "share of away-nights spent with another unit", "var(--c5)"),
]


def build():
    cols = ["window_start", "animal_id", "origin_group", "dynamic_social_unit",
            "social_context", "association_event_id", "temp_group_size", "is_observed"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
    s1 = pd.read_csv(STAGE1, parse_dates=["start_hour"])
    exc = pd.read_csv(EXC, parse_dates=["start_night", "end_night"])
    rng = np.random.default_rng(20260904)

    examples = [ex_between(d, s1), ex_within(d), ex_individual(d, exc)]
    depths = [between_depth(s1), within_depth(d), individual_depth(exc)]

    o, y = [], 68.0
    rw, hx, hw, hh = 396.0, 476.0, 196.0, 92.0
    stats = []
    for i, ((name, defn, dlab, col), g, v) in enumerate(zip(TYPES, examples, depths)):
        o.append('<text class="pl" x="0" y="%.1f">%s</text>' % (y - 20, "abc"[i]))
        o.append('<text class="ph" x="16" y="%.1f">%s</text>' % (y - 20, name))
        o.append('<text class="ts" x="16" y="%.1f">%s</text>' % (y - 10, defn))
        rs, rh = raster(g, 22.0, y, rw)
        o.append(rs)
        o.append(depth_hist(v, hx, y + 8, hw, hh, col, dlab, rng))
        stats.append({"name": name, "n": int(len(v)),
                      "median": float(np.median(v)),
                      "q1": float(np.quantile(v, .25)),
                      "q3": float(np.quantile(v, .75)),
                      "at_zero": float((v <= 0.001).mean()),
                      "at_one": float((v >= 0.999).mean()),
                      "below_half": float((v < 0.5).mean()),
                      "meta": g["meta"]})
        y += max(rh, hh + 46) + 42
    H = y - 20
    head = ('<svg viewBox="0 0 700 %.0f" role="img" aria-label="Three social event '
            'types, each with a worked example raster from the tracking record and the '
            'distribution of its continuous depth measure on a zero to one scale.">' % H)
    return head + "\n" + "\n".join(o) + "\n</svg>", stats


CSS = """
:root{--ground:#f4f5f3;--paper:#fff;--ink:#15191a;--ink-2:#3a4244;--muted:#6b7476;
--rule:#d5dad7;--rule-soft:#e6eae7;--n2:#c3cbc8;
--c1:#1f6f8b;--c4:#8a5a1f;--c5:#6b4a7a;--c6:#0f5f57;
--sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;
--serif:"Newsreader",Georgia,serif;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0f1213;
--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;--muted:#8f9998;--rule:#2c3435;
--rule-soft:#242c2d;--n2:#3f4a4b;--c1:#74b6d0;--c4:#d3a061;--c5:#b494c4;--c6:#5fc0b0;}}
:root[data-theme="dark"]{--ground:#0f1213;--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;
--muted:#8f9998;--rule:#2c3435;--rule-soft:#242c2d;--n2:#3f4a4b;--c1:#74b6d0;
--c4:#d3a061;--c5:#b494c4;--c6:#5fc0b0;}
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
.note{font-size:11.5px;line-height:1.6;color:var(--muted);margin-top:22px;
padding-top:14px;border-top:1px solid var(--rule-soft)}
.note code{font-family:var(--mono);font-size:.92em}
text{font-family:var(--sans)}
.rl{font-size:6.5px;fill:var(--muted);font-family:var(--mono)}
.rlf{font-size:6.5px;fill:var(--ink);font-weight:700;font-family:var(--mono)}
.ul{font-size:7.4px;fill:var(--ink-2);font-weight:600}
.ts{font-size:8.2px;fill:var(--muted)}
.ph{font-size:12px;fill:var(--ink);font-weight:600}
.pl{font-size:13px;font-weight:700;fill:var(--ink)}
.g{stroke:var(--rule-soft);stroke-width:.8}
.ax{stroke:var(--rule);stroke-width:1}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def page():
    fig, st = build()
    a, b, c = st
    return f"""<title>Three Event Types</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;700&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 1 &middot; draft</div>
  <h1>Three ways a group boundary opens, each with one continuous depth</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 1.</span> The three social event types, defined by who
      crosses. <b>Left of each panel:</b> a worked example from the tracking record &mdash; rows
      are collared animals, columns hours (nights in c), cell colour is the animal's social
      context, and unobserved animal-intervals are left <b>blank</b> so observation and behaviour
      stay separable. <b>Right:</b> the distribution of that type's depth measure over every
      event, all three on a common 0&ndash;1 scale, with the median and its bootstrap interval.
      The three depth distributions have three characteristically different shapes, and that is
      the point of the panel.
      <b>(a)</b> Between-group: two units share one cluster; depth is mutual participation.
      Median {a["median"]:.2f}, IQR {a["q1"]:.2f}&ndash;{a["q3"]:.2f}, n&nbsp;=&nbsp;{a["n"]:,}
      encounters &mdash; and <b>{100 * a["at_one"]:.0f}% sit at exactly 1.0</b> against
      {100 * a["below_half"]:.0f}% below half. Between-group opening is largely all-or-nothing:
      when two units meet, most of both units are usually in it.
      <b>(b)</b> Within-group: one unit holds more than one cluster; depth is the share of the
      unit outside its own largest cluster. Median {b["median"]:.2f}, IQR
      {b["q1"]:.2f}&ndash;{b["q3"]:.2f}, n&nbsp;=&nbsp;{b["n"]:,} unit-hours &mdash; low and
      continuous, so a unit that divides usually detaches a small minority rather than halving.
      <b>(c)</b> Individual: one animal leaves its unit; depth is the share of away-nights spent
      with another unit. Median {c["median"]:.2f}, n&nbsp;=&nbsp;{c["n"]:,} excursions, and
      <b>strongly bimodal</b> &mdash; {100 * c["at_zero"]:.0f}% at 0 and
      {100 * c["at_one"]:.0f}% at 1. An animal is either alone for the whole trip or with another
      unit for the whole trip; passing through both is rare.
    </p>
  </div>

  <p class="note">
    <b>This replaces a six-type taxonomy</b> that split merges into large, medium-partial and
    small-subset. Two reasons. The saved <code>merge_scale</code> label is assigned from the
    <em>maximum</em> participation a unit ever reaches rather than its typical state, so events
    labelled small-subset can have maximum fractions near 1.0. And the subdivision does not
    distinguish groups: between-group ICC 0.05&ndash;0.13 with permutation p = 0.09&ndash;0.26,
    and well-sampled groups traverse most of its range within their own record. It is a gradient
    measured on a faulty rule and reported as three categories, so it is better carried as the
    continuous depth in panel (a). Examples: {a["meta"]}; {b["meta"]}; {c["meta"]}. All panels
    and depth distributions come from the frozen 1,924,104-row export (2024-03-01 to
    2026-07-22); within-group depth is computed over unit-hours with at least
    {MIN_OBS_FOR_SPLIT} observed animals. Generated by
    <code>build_three_types_figure.py</code>.
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
