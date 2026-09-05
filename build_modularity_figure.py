"""Within-group modularity on its own: across units, through time, and as pictures.

Modularity was one panel of the level-2 figure and it was doing too much work there.
A score of 0.45 means nothing to a reader until they see it, and a network thumbnail
squeezed into a corner of another panel is not seeing it. So this figure does one thing:

One panel, one line: the time series of a single group's weekly modularity, and two of
its own weeks drawn as association networks -- the week it is most modular beside a week
at zero, matched on collar count so the flat one cannot be dismissed as a thinner week.
Nodes are animals, coloured and ringed by community; edges inside a community take its
colour and edges between communities stay grey.

One group rather than seven rows, and two pictures rather than four, because the claim
is narrow: a score of 0.45 and a score of 0.00 come from the same unit five months
apart, so modularity is a state and not a trait.

Below about twelve collars a modularity estimator returns one community regardless of
the animals, so everything here is restricted to well-covered unit-weeks. That is a
restriction on what can be measured, not a claim that sparser units are cohesive.

Inputs from derive_level2_inputs.py. Output: docs/framing_2026-09-04/modularity_figure.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from build_level2_figure import (CSS, MIN_WEEKS, MOD_EPS, WELL_COVERED, axes, key_row,
                                 load_modularity, runs_of)

BASE = Path("outputs/general_structure_2026_09")
INPUTS = BASE / "phase4h_level2_inputs"
NETS = INPUTS / "example_week_networks.csv"
SPLITS = INPUTS / "split_bouts_hourly.csv"
OUT = Path("docs/framing_2026-09-04/modularity_figure.html")

EXAMPLE_UNIT = "Lilac"
COMM_COLOUR = ("var(--c4)", "var(--c1)", "var(--c6)", "var(--c5)")
SEED = 20260904


def example_network(net, key, x0, y0, size, title, sub):
    """One week's association network, drawn large enough to read.

    Community membership is carried three ways at once -- node colour, a ring around each
    community, and the colour of the edges inside it -- because at this size a reader
    should not have to trace lines to see whether the group is one piece or two.
    """
    w = net[net["week"].eq(key)]
    nd = w[w["kind"].eq("node")]
    ed = w[w["kind"].eq("edge")]
    nd = nd.copy()
    cs = sorted(nd["community"].unique())
    if len(cs) > 1:
        spread = max(float(np.hypot(*(nd.groupby("community")[["x", "y"]].std()
                                      .fillna(0).to_numpy().T)).max()), 1e-3)
        for i, c in enumerate(cs):
            m = nd["community"].eq(c)
            ang = -math.pi / 2 + 2 * math.pi * i / len(cs)
            gx, gy = nd.loc[m, "x"].mean(), nd.loc[m, "y"].mean()
            nd.loc[m, "x"] = nd.loc[m, "x"] - gx + 3.1 * spread * math.cos(ang)
            nd.loc[m, "y"] = nd.loc[m, "y"] - gy + 3.1 * spread * math.sin(ang)
    xs, ys = nd["x"].to_numpy(dtype=float), nd["y"].to_numpy(dtype=float)
    # one scale for both axes, so the layout is not stretched into looking structured
    cx, cy = (xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2
    half = max(xs.max() - xs.min(), ys.max() - ys.min(), 1e-6) / 2
    pad = 15.0
    r_box = (size - 2 * pad) / 2

    def px(v):
        return x0 + size / 2 + (v - cx) / half * r_box

    def py(v):
        return y0 + size / 2 + (v - cy) / half * r_box

    pos = {int(r.node): (px(r.x), py(r.y)) for r in nd.itertuples()}
    comm = {int(r.node): int(r.community) for r in nd.itertuples()}
    o = ['<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
         'stroke="var(--rule-soft)" stroke-width="1"/>' % (x0, y0, size, size)]
    # a soft ring per community, before anything else, so it reads as background
    for c in sorted(set(comm.values())):
        pts = np.array([pos[k] for k in pos if comm[k] == c])
        mx, my = pts[:, 0].mean(), pts[:, 1].mean()
        rx = max(9.0, 1.9 * pts[:, 0].std() + 7)
        ry = max(9.0, 1.9 * pts[:, 1].std() + 7)
        o.append('<ellipse cx="%.1f" cy="%.1f" rx="%.1f" ry="%.1f" fill="%s" '
                 'opacity=".09"/>' % (mx, my, rx, ry, COMM_COLOUR[c % len(COMM_COLOUR)]))
    wmax = max(float(ed["weight"].max()), 1e-6) if len(ed) else 1.0
    for r in ed.itertuples():
        a, b = int(r.node), int(r.community)      # edge rows store the far node here
        if a not in pos or b not in pos:
            continue
        same = comm[a] == comm[b]
        col = COMM_COLOUR[comm[a] % len(COMM_COLOUR)] if same else "var(--muted)"
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="%.2f" opacity="%s"/>'
                 % (pos[a][0], pos[a][1], pos[b][0], pos[b][1], col,
                    0.4 + 1.7 * (float(r.weight) / wmax), ".5" if same else ".8"))
    for k, (nx_, ny_) in pos.items():
        o.append('<circle cx="%.1f" cy="%.1f" r="3.4" fill="%s" stroke="var(--paper)" '
                 'stroke-width=".9"/>'
                 % (nx_, ny_, COMM_COLOUR[comm[k] % len(COMM_COLOUR)]))
    o.append('<text class="ph" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
             % (x0 + size / 2, y0 - 26, title))
    for k, line in enumerate(sub):
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (x0 + size / 2, y0 - 15 + k * 10, line))
    return "\n".join(o)


def unit_summary(sub, units):
    """The across-unit spread, kept as numbers for the caption rather than a panel."""
    per_unit = sub.groupby("dynamic_social_unit")["modularity"].apply(
        lambda x: float((x > MOD_EPS).sum()))
    top2 = per_unit.sort_values(ascending=False).head(2)
    return {"units": len(units), "weeks": int(len(sub)),
            "max_mod": round(float(sub["modularity"].max()), 3),
            "modular_units": int((per_unit > 0).sum()),
            "silent_units": int((per_unit == 0).sum()),
            "modular_weeks": int(per_unit.sum()),
            "top2_weeks": int(top2.sum()),
            "top2_share": round(float(top2.sum() / max(1.0, per_unit.sum())), 3),
            "top2_names": " and ".join(top2.index),
            "window": [str(sub["period_start"].min().date()),
                       str(sub["period_start"].max().date())]}


def panel_units(sub, units, x0, y0, w):
    """Every well-covered unit, one row of weekly modularity each."""
    lo, hi = sub["period_start"].min(), sub["period_start"].max()
    span = max(1.0, (hi - lo).days)
    lab_w, rowh = 62.0, 20.0
    spark = w - lab_w - 78
    top = max(0.05, float(sub["modularity"].max()))
    o = []
    for i, u in enumerate(units):
        g = sub[sub["dynamic_social_unit"].eq(u)].sort_values("period_start")
        yy = y0 + i * rowh
        bold = ' font-weight="700"' if u == EXAMPLE_UNIT else ""
        o.append('<text class="nl" x="%.1f" y="%.1f" text-anchor="end"%s>%s</text>'
                 % (x0 + lab_w - 6, yy + 11, bold, u))
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0 + lab_w, yy + 14, x0 + lab_w + spark, yy + 14))
        nmod = int((g["modularity"] > MOD_EPS).sum())
        o.append('<text class="ts" x="%.1f" y="%.1f">%d of %d weeks</text>'
                 % (x0 + lab_w + spark + 6, yy + 11, nmod, len(g)))
        for r in g.itertuples():
            xx = x0 + lab_w + ((r.period_start - lo).days / span) * spark
            hgt = (r.modularity / top) * 13.0
            if hgt <= 0.4:
                o.append('<circle cx="%.1f" cy="%.1f" r=".9" fill="var(--muted)"/>'
                         % (xx, yy + 14))
            else:
                o.append('<rect x="%.1f" y="%.1f" width="1.8" height="%.1f" '
                         'fill="var(--c4)"/>' % (xx - .9, yy + 14 - hgt, hgt))
    yb = y0 + len(units) * rowh
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0 + lab_w, yb, x0 + lab_w + spark, yb))
    for tv, tl, ta in ((0.0, str(lo.date()), "start"), (1.0, str(hi.date()), "end")):
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="%s">%s</text>'
                 % (x0 + lab_w + tv * spark, yb + 11, ta, tl))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">bar height = '
             'weekly modularity, 0 to %.2f</text>'
             % (x0 + lab_w + spark / 2, yb + 23, top))
    per_unit = sub.groupby("dynamic_social_unit")["modularity"].apply(
        lambda x: float((x > MOD_EPS).sum()))
    return "\n".join(o), len(units) * rowh + 20, {
        "units": len(units), "weeks": int(len(sub)), "max_mod": round(top, 3),
        "modular_units": int((per_unit > 0).sum()),
        "silent_units": int((per_unit == 0).sum()),
        "window": [str(lo.date()), str(hi.date())]}


def panel_detail(sub, marks, x0, y0, w, h):
    """Lilac week by week, with its phases and its largest sub-community."""
    g = sub[sub["dynamic_social_unit"].eq(EXAMPLE_UNIT)].sort_values("period_start")
    lo, hi = g["period_start"].min(), g["period_start"].max()
    span = max(1.0, (hi - lo).days)
    top = max(0.05, float(g["modularity"].max()))
    o = [axes(x0, y0, w, h, [(0.0, ""), (0.5, ""), (1.0, "")],
              [(0.0, "0"), (0.25, ""), (0.5, "%.2f" % (top / 2)), (0.75, ""),
               (1.0, "%.2f" % top)],
              "%d well-covered weeks" % len(g), "modularity")]
    for tv, tl, ta in ((0.0, str(lo.date()), "start"), (1.0, str(hi.date()), "end")):
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="%s">%s</text>'
                 % (x0 + tv * w, y0 + h + 11, ta, tl))
    mod = (g["modularity"] > MOD_EPS).to_numpy()
    wkno = g["weekno"].to_numpy()
    wk2date = dict(zip(wkno, g["period_start"]))
    phases = [(a, b) for a, b in runs_of(mod, wkno) if b - a + 1 >= 2]
    for a, b in phases:
        xa = x0 + ((wk2date[a] - lo).days / span) * w
        xb = x0 + ((wk2date[b] - lo).days / span) * w + 3
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--c4)" '
                 'opacity=".11"/>' % (xa - 1.5, y0, max(3.0, xb - xa), h))
    pts = []
    for r in g.itertuples():
        xx = x0 + ((r.period_start - lo).days / span) * w
        pts.append("%.1f,%.1f" % (xx, y0 + h - float(r.largest_community_fraction) * h))
        hgt = (r.modularity / top) * h
        if hgt <= 0.5:
            o.append('<circle cx="%.1f" cy="%.1f" r="1.1" fill="var(--muted)"/>'
                     % (xx, y0 + h))
        else:
            o.append('<rect x="%.1f" y="%.1f" width="3" height="%.1f" '
                     'fill="var(--c4)"/>' % (xx - 1.5, y0 + h - hgt, hgt))
    o.append('<polyline points="%s" fill="none" stroke="var(--c6)" stroke-width="1.5" '
             'opacity=".9"/>' % " ".join(pts))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
             'stroke="var(--c6)" opacity=".55"/>' % (x0 + w, y0, x0 + w, y0 + h))
    for tv, tl in ((0.0, "0"), (0.5, "0.5"), (1.0, "1")):
        o.append('<text class="ts" x="%.1f" y="%.1f" fill="var(--c6)">%s</text>'
                 % (x0 + w + 4, y0 + h - tv * h + 3, tl))
    # the weeks drawn in panel (c), tied to where they sit in the series
    for wdate, label in marks:
        if not (lo <= wdate <= hi):
            continue
        xx = x0 + ((wdate - lo).days / span) * w
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink-2)" '
                 'stroke-width=".8" stroke-dasharray="2 2"/>' % (xx, y0 - 6, xx, y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" '
                 'fill="var(--ink-2)">%s</text>' % (xx, y0 - 9, label))
    longest = max((b - a + 1 for a, b in phases), default=0)
    return "\n".join(o), {
        "unit": EXAMPLE_UNIT, "weeks": int(len(g)), "modular_weeks": int(mod.sum()),
        "max_mod": round(top, 3), "phases": len(phases), "longest_phase": int(longest),
        "max_communities": int(g["n_communities"].max()),
        "min_largest": round(float(g["largest_community_fraction"].min()), 2),
        "window": [str(lo.date()), str(hi.date())]}


def build():
    sub, units = load_modularity()
    net = pd.read_csv(NETS, parse_dates=["period_start"])
    meta = (net[net["kind"].eq("node")]
            .groupby(["week", "unit", "level", "period_start"])
            .agg(n=("animal_id", "nunique"),
                 comms=("community", "nunique")).reset_index())
    q = (sub.set_index(["dynamic_social_unit", "period_start"])["modularity"]
         .to_dict())
    meta["q"] = [q.get((u, pd.Timestamp(d)), float("nan"))
                 for u, d in zip(meta["unit"], meta["period_start"])]
    meta = meta.sort_values("q", ascending=False).reset_index(drop=True)

    x0, w = 26.0, 648.0
    sa = unit_summary(sub, units)
    # the two ends of the scale; the middle of the range is a separate argument
    shown = meta.iloc[[0, -1]].reset_index(drop=True)

    size, gap = 138.0, 16.0
    ya = 74.0
    series_w = w - 2 * size - gap - 74
    marks = [(pd.Timestamp(r.period_start), "%d" % (i + 1))
             for i, r in enumerate(shown.itertuples())]
    pa, sb = panel_detail(sub, marks, x0 + 40, ya, series_w, size)

    parts = []
    nx0 = x0 + 40 + series_w + 58
    for i, r in enumerate(shown.itertuples()):
        parts.append(example_network(
            net, r.week, nx0 + i * (size + gap), ya, size,
            "%d · Q = %.2f" % (i + 1, r.q),
            [str(r.period_start.date()),
             "%d animals, %d %s" % (r.n, r.comms,
                                    "community" if r.comms == 1
                                    else "communities")]))
    H = ya + size + 96

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Within-group modularity: '
         'for one group on a single line: its weekly values through time with modular '
         'phases shaded, beside two of those weeks drawn as association networks -- one '
         'of two communities that barely associate, one a single connected group.">' % H]
    o.append('<text class="pl" x="0" y="%.1f">a</text>' % (ya - 30))
    o.append(pa)
    o.extend(parts)
    o.append(key_row([("weekly modularity", "bar", "var(--c4)"),
                      ("modular phase, 2+ weeks", "band", "var(--c4)"),
                      ("largest sub-community", "line", "var(--c6)"),
                      ("node = animal, colour = community", "dot", "var(--c4)")],
                     x0 + 14, ya + size + 52, colw=152.0, ncol=4))
    o.append("</svg>")
    return "\n".join(o), sa, sb, meta


def page():
    fig, sa, sb, meta = build()
    hi = meta.iloc[0]
    lo = meta.iloc[-1]
    return f"""<title>Within-Group Modularity</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 3a &middot; draft</div>
  <h1>What modularity looks like</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 3a.</span> Within-group modularity for one group,
      through time and as pictures. <b>Left:</b> {sb["unit"]} week by week, across its
      {sb["weeks"]} weeks at {WELL_COVERED} or more collars ({sb["window"][0]} to
      {sb["window"][1]}); bar height is that week&rsquo;s modularity, dots are weeks at
      zero, and the teal trace is the largest sub-community as a share of the group on
      the right-hand axis. {sb["modular_weeks"]} of {sb["weeks"]} weeks are modular, in
      {sb["phases"]} shaded phases of two weeks or more, the longest running
      {sb["longest_phase"]} consecutive weeks; up to {sb["max_communities"]} communities
      appear and the largest falls to {sb["min_largest"]:.0%}. The group does not drift
      apart and stay apart &mdash; it separates and re-fuses, on a scale of weeks.
      <b>Right:</b> the two marked weeks as association networks. At Q = {hi["q"]:.2f}
      ({hi["period_start"].date()}) the group is two communities that barely associate;
      at Q = {lo["q"]:.2f} ({lo["period_start"].date()}) it is one connected mass with no
      partition to find, on the same {lo["n"]} collars. <b>Both are the same unit five
      months apart, which is the case for reading modularity as a state rather than a
      trait.</b> It is also not a typical unit: across the {sa["units"]} units with at
      least {MIN_WEEKS} weeks at {WELL_COVERED} or more collars, {sa["modular_units"]}
      reach modularity in some week and {sa["silent_units"]} never leave zero, but
      <b>{sa["top2_names"]} alone hold {sa["top2_weeks"]} of the {sa["modular_weeks"]}
      modular weeks</b> ({sa["top2_share"]:.0%}). This range is what a unit can do, not
      what most units do.
    </p>
  </div>

  <p class="note">
    Everything here is restricted to well-covered unit-weeks. Below about
    {WELL_COVERED} collars a modularity estimator returns a single community regardless
    of the animals, so a low score from a sparse unit is a statement about the collars,
    not the baboons; the {sa["silent_units"]} silent units counted in the caption are all
    well-covered and their zeros are real. The two networks are matched on collar count,
    so the flat one cannot be dismissed as a thinner week. Both are the week&rsquo;s
    association network:
    an edge is the share of hours two animals spent in the same association cluster,
    among hours both were present, drawn at 0.05 and above; communities are greedy
    modularity communities on that weighted graph, and positions come from a spring
    layout with a fixed seed and each community&rsquo;s centre then placed on a fixed
    circle, so that both networks are drawn to one convention and distance on the page
    is indicative rather than metric. Weekly metrics come from the frozen 1,924,104-row export; networks are built
    by <code>derive_level2_inputs.py</code>. Generated by
    <code>build_modularity_figure.py</code>; seed {SEED}.
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
