"""Figure 3: the three routes across a group boundary, each shown and then compared.

Row by row, the figure moves from what each structure looks like to how they differ:

  (a) WITHIN-GROUP MODULARITY. One group's weekly modularity, with two of its own weeks
      drawn as association networks inset where the line leaves room, each ruled back to
      the point it illustrates.
  (b) MEETING AND DISPERSAL ON ONE LAYOUT. Blue chords are group meetings, width by
      the hours two units spent together with both sides contributing two or more
      animals; plum chords are single animals moving between the same units.
  (c) CONCENTRATION. All four quantities as Lorenz curves against one diagonal, each
      Gini with a cluster-bootstrap interval.
  (d) DURATION. All three states at hourly resolution on one logarithmic axis.

The two chord types in (b) share a layout on purpose -- the comparison is the point --
but they are mostly, not entirely, independent evidence: 16 dyads carry both, and on
Copper-Lilac and Maroon-Sapphire most of the encounter record falls inside a period when
one of their animals was living in the other unit. A doubled chord on those two dyads is
one finding, not two, and the caption says so.

Everything is built from the frozen export; hourly states come from
derive_level2_inputs.py. Output: docs/framing_2026-09-04/figure3.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from build_level2_figure import (CSS, COMM_COLOUR, DESTS, EXC, EXC_H, MOD_EPS, N_BOOT,
                                 SPLITS_H, WELL_COVERED, axes, ellipse_network,
                                 encounter_bouts, key_row, load_modularity,
                                 load_networks, panel_conc, panel_duration, runs_of)
from build_modularity_figure import NETS, unit_summary

OUT = Path("docs/framing_2026-09-04/figure3.html")

EXAMPLE_UNIT = "Lilac"
SEED = 20260904
X_END = "2026-01-31"    # panel (a) axis end; the tail beyond it is reported

BOX_W, BOX_H = 262.0, 194.0
COL_X = (58.0, 388.0)
ROW_Y = (78.0, 380.0)


NEUTRAL_COMM = ("var(--ink-2)", "var(--c6)", "var(--muted)", "var(--c1)")


def inset_network(net, key, x0, y0, size, label, palette=None):
    """A small association network, sized to sit inside another panel.

    No title and no node labels: at this size the inset carries one fact, which is
    whether the group is one piece or two, and the caption carries the rest.
    """
    pal = palette or COMM_COLOUR
    w = net[net["week"].eq(key)]
    nd = w[w["kind"].eq("node")].copy()
    ed = w[w["kind"].eq("edge")]
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
    cx, cy = (xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2
    half = max(xs.max() - xs.min(), ys.max() - ys.min(), 1e-6) / 2
    r_box = (size - 16.0) / 2
    pos = {int(r.node): (x0 + size / 2 + (r.x - cx) / half * r_box,
                         y0 + size / 2 + (r.y - cy) / half * r_box)
           for r in nd.itertuples()}
    comm = {int(r.node): int(r.community) for r in nd.itertuples()}
    o = ['<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--paper)" '
         'stroke="var(--rule)" stroke-width=".8"/>' % (x0, y0, size, size + 11)]
    wmax = max(float(ed["weight"].max()), 1e-6) if len(ed) else 1.0
    for r in ed.itertuples():
        a, b = int(r.node), int(r.community)      # edge rows store the far node here
        if a not in pos or b not in pos:
            continue
        same = comm[a] == comm[b]
        col = pal[comm[a] % len(pal)] if same else "var(--c4)"
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="%.2f" opacity="%s"/>'
                 % (pos[a][0], pos[a][1], pos[b][0], pos[b][1], col,
                    (0.3 if same else 0.45) + 1.1 * (float(r.weight) / wmax),
                    ".45" if same else ".8"))
    for k, (nx_, ny_) in pos.items():
        o.append('<circle cx="%.1f" cy="%.1f" r="2.2" fill="%s"/>'
                 % (nx_, ny_, pal[comm[k] % len(pal)]))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
             % (x0 + size / 2, y0 + size + 8, label))
    return "\n".join(o)


def free_slot(g, lo, span, x0, w, size, top, h, taken, near_x):
    """Where an inset of `size` fits over the series without covering the line.

    Slides a window across the plot, measures how high the data reaches inside it, and
    keeps the windows where the inset would clear the peak. Among those it takes the one
    nearest `near_x`, so each picture sits close to the week it came from rather than
    wherever happens to be empty -- subject to staying clear of an inset already placed,
    because two boxes that touch read as one. Falls back to the emptiest window if
    nothing clears.
    """
    xs = np.array([(r.period_start - lo).days / span for r in g.itertuples()])
    qs = np.array([r.modularity / top for r in g.itertuples()])
    need = (size + 13.0) / h          # fraction of the box height the inset occupies
    best, fallback = None, None
    for frac in np.linspace(0.0, 1.0, 61):
        ix = x0 + frac * (w - size)
        a, b = (ix - x0) / w, (ix - x0 + size) / w
        inside = qs[(xs >= a - 0.02) & (xs <= b + 0.02)]
        ceiling = float(inside.max()) if len(inside) else 0.0
        # keep real air between them, or two insets read as one block
        if any(abs(ix - t) < size + 44 for t in taken):
            continue
        cand = (abs(ix + size / 2 - near_x), ix)
        if ceiling <= 1.0 - need:
            if best is None or cand < best[0:2]:
                best = (cand[0], ix, ceiling)
        if fallback is None or ceiling < fallback[1]:
            fallback = (ix, ceiling)
    return best[1] if best else (fallback[0] if fallback else x0)


def panel_modularity(sub, net, shown, x0, y0, w, h, size=62.0):
    """Modularity through time, with two of the group's own weeks inset in the gaps.

    The line spends most of the panel near zero, so the pictures go where the data is
    not: each inset is placed in the emptiest window nearest its own week and ruled back
    to the point on the line it illustrates. One line only -- ochre is modularity here
    and nothing else, which is why the communities take neutral colours.
    """
    allg = sub[sub["dynamic_social_unit"].eq(EXAMPLE_UNIT)].sort_values("period_start")
    end = pd.Timestamp(X_END)
    # the axis stops at X_END; the weeks past it are flat and would squeeze the part of
    # the record that has anything in it into the left half of an already small panel.
    # They are not dropped from the argument -- the caption states what they contain.
    g = allg[allg["period_start"].le(end)]
    tail = allg[allg["period_start"].gt(end)]
    lo, hi = g["period_start"].min(), end
    span = max(1.0, (hi - lo).days)
    top = max(0.05, float(allg["modularity"].max()))
    months = pd.date_range(lo.normalize().replace(day=1) + pd.offsets.MonthBegin(1),
                           hi, freq="MS")
    o = [axes(x0, y0, w, h,
              [(float((m - lo).days) / span, "") for m in months],
              [(0.0, "0"), (0.25, ""), (0.5, "%.2f" % (top / 2)), (0.75, ""),
               (1.0, "%.2f" % top)], "date", "modularity")]
    for tv, tl, ta in ((0.0, str(lo.date()), "start"), (1.0, str(hi.date()), "end")):
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="%s">%s</text>'
                 % (x0 + tv * w, y0 + h + 11, ta, tl))
    mod = (g["modularity"] > MOD_EPS).to_numpy()
    wk2date = dict(zip(g["weekno"].to_numpy(), g["period_start"]))
    phases = [(a, b) for a, b in runs_of(mod, g["weekno"].to_numpy())
              if b - a + 1 >= 2]
    for a, b in phases:
        xa = x0 + ((wk2date[a] - lo).days / span) * w
        xb = x0 + ((wk2date[b] - lo).days / span) * w + 3
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--c4)" '
                 'opacity=".11"/>' % (xa - 1.5, y0, max(3.0, xb - xa), h))
    pts = ["%.1f,%.1f" % (x0 + ((r.period_start - lo).days / span) * w,
                          y0 + h - (r.modularity / top) * h) for r in g.itertuples()]
    o.append('<polyline points="%s" fill="none" stroke="var(--c4)" '
             'stroke-width="1.7"/>' % " ".join(pts))
    taken = []
    for r in shown.itertuples():
        wdate = pd.Timestamp(r.period_start)
        wx = x0 + max(0.0, min(1.0, (wdate - lo).days / span)) * w
        wy = y0 + h - (max(0.0, float(r.q)) / top) * h
        ix = free_slot(g, lo, span, x0, w, size, top, h, taken, wx)
        taken.append(ix)
        iy = y0 + 3
        o.append('<path d="M %.1f %.1f L %.1f %.1f" fill="none" stroke="var(--ink-2)" '
                 'stroke-width=".8" stroke-dasharray="2 2"/>'
                 % (ix + size / 2, iy + size + 11, wx, wy))
        o.append('<circle cx="%.1f" cy="%.1f" r="2" fill="var(--c4)"/>' % (wx, wy))
        o.append(inset_network(net, r.week, ix, iy, size, "Q = %.2f" % r.q,
                               palette=NEUTRAL_COMM))
    o.append(key_row([("weekly modularity", "line", "var(--c4)"),
                      ("modular phase, 2+ weeks", "band", "var(--c4)"),
                      ("insets: ochre edges bridge the two communities", "dot",
                       "var(--ink-2)")],
                     x0 - 28, y0 + h + 46, colw=w / 2 + 22, ncol=2))
    return "\n".join(o), {
        "unit": EXAMPLE_UNIT, "weeks": int(len(g)), "modular_weeks": int(mod.sum()),
        "all_weeks": int(len(allg)),
        "all_modular": int((allg["modularity"] > MOD_EPS).sum()),
        "tail_weeks": int(len(tail)),
        "tail_max": round(float(tail["modularity"].max()), 3) if len(tail) else 0.0,
        "tail_end": str(allg["period_start"].max().date()),
        "max_mod": round(top, 3), "phases": len(phases),
        "longest_phase": max((b - a + 1 for a, b in phases), default=0),
        "max_communities": int(g["n_communities"].max()),
        "min_largest": round(float(g["largest_community_fraction"].min()), 2),
        "window": [str(lo.date()), str(hi.date())]}


def build():
    rng = np.random.default_rng(SEED)
    sub, units = load_modularity()
    nodes, cov, grp, solo, dyh_raw, ed, deg = load_networks()
    bouts = encounter_bouts()
    dyh = (bouts.groupby(["unit_a", "unit_b", "pair_key"])
           .agg(hours=("hours", "sum"), bouts=("hours", "size")).reset_index())
    dests = pd.read_csv(DESTS)
    with_d = dests[dests["destination"].notna()]
    per_unit = sub.groupby("dynamic_social_unit")["modularity"].apply(
        lambda x: float((x > MOD_EPS).sum()))
    per_animal = with_d.groupby("animal_id").size()
    exc = pd.read_csv(EXC)
    exch = pd.read_csv(EXC_H)
    exch = exch[exch["animal_id"].isin(set(exc["animal_id"]))]
    splits = pd.read_csv(SPLITS_H)

    net = pd.read_csv(NETS, parse_dates=["period_start"])
    meta = (net[net["kind"].eq("node")]
            .groupby(["week", "unit", "level", "period_start"])
            .agg(n=("animal_id", "nunique"), comms=("community", "nunique"))
            .reset_index())
    q = sub.set_index(["dynamic_social_unit", "period_start"])["modularity"].to_dict()
    meta["q"] = [q.get((u, pd.Timestamp(d)), float("nan"))
                 for u, d in zip(meta["unit"], meta["period_start"])]
    meta = meta.sort_values("q", ascending=False).reset_index(drop=True)

    sa = unit_summary(sub, units)
    shown = meta.iloc[[0, 2]].reset_index(drop=True)
    pa, sb = panel_modularity(sub, net, shown, COL_X[0], ROW_Y[0], BOX_W, BOX_H)

    # one layout, both edge types: the two networks share their nodes, so drawing them
    # apart makes the reader do the comparison by eye across two ellipses
    enc_edges = [(r.unit_a, r.unit_b, float(r.hours)) for r in dyh.itertuples()]
    tr_edges = [(r.origin_group, r.destination, float(r.excursions))
                for r in ed.itertuples()]
    pb = ellipse_network(nodes, enc_edges, COL_X[1], ROW_Y[0], BOX_W, BOX_H, cov,
                         "var(--ink-2)", "var(--c1)")
    pb += "\n" + ellipse_network(nodes, tr_edges, COL_X[1], ROW_Y[0], BOX_W, BOX_H,
                                 cov, "var(--ink-2)", "var(--c5)", nodes_off=True)
    pb += "\n" + key_row([("group meeting, width = hours", "line", "var(--c1)"),
                          ("animal transfer, width = excursions", "line", "var(--c5)"),
                          ("node area = median collars", "dot", "var(--ink-2)")],
                         COL_X[1] - 28, ROW_Y[0] + BOX_H + 46,
                         colw=BOX_W / 2 + 22, ncol=2)

    yd = ROW_Y[1]
    pd_svg, gs = panel_conc(
        [("modular weeks per unit", "mod", per_unit.to_numpy(),
          per_unit.index.to_numpy(), "var(--c4)", None),
         ("encounter-hours per dyad", "enc", dyh["hours"].to_numpy(dtype=float),
          dyh["pair_key"].to_numpy(), "var(--c1)", None),
         ("excursions per transfer edge", "edge", ed["excursions"].to_numpy(dtype=float),
          (ed["origin_group"] + ">" + ed["destination"]).to_numpy(), "var(--c5)", None),
         ("excursions per animal", "animal", per_animal.to_numpy(dtype=float),
          per_animal.index.to_numpy(), "var(--c5)", "4 2.5")],
        COL_X[0], yd, BOX_W, BOX_H, rng)
    pe, ds = panel_duration(
        [("group meetings", "enc", "dyads", bouts["hours"].to_numpy(dtype=float),
          bouts["pair_key"].to_numpy(), "var(--c1)"),
         ("within-group splits", "mod", "units", splits["hours"].to_numpy(dtype=float),
          splits["unit"].to_numpy(), "var(--c4)"),
         ("individual excursions", "exc", "animals",
          exch["hours"].to_numpy(dtype=float), exch["animal_id"].to_numpy(),
          "var(--c5)")],
        COL_X[1], yd, BOX_W, BOX_H, rng)
    H = yd + BOX_H + 88

    enc_pairs = {frozenset((r.unit_a, r.unit_b)) for r in dyh.itertuples()}
    tr_pairs = {frozenset((r.origin_group, r.destination)) for r in ed.itertuples()}
    hrs = dyh["hours"].to_numpy(dtype=float)
    sc = {"units": len(nodes), "dyads": int(len(dyh)), "bouts": int(dyh["bouts"].sum()),
          "hours_max": int(hrs.max()),
          "top5_share": round(float(np.sort(hrs)[::-1][:5].sum() / hrs.sum()), 3),
          "degree_max": int(max(deg.values())),
          "solo_dyad_days": int(len(solo)),
          "solo_share": round(float(len(solo) / (len(solo) + len(grp))), 3),
          "solo_dyads": int(solo["pair_key"].nunique()),
          "solo_distance": int(round(float(solo["min_centroid_distance_m"].median()))),
          "grp_distance": int(round(float(grp["min_centroid_distance_m"].median()))),
          "edges": int(len(ed)), "animals": int(with_d["animal_id"].nunique()),
          "excursions": int(len(with_d)),
          "reciprocal": sum(1 for r in ed.itertuples()
                            if ((ed["origin_group"] == r.destination)
                                & (ed["destination"] == r.origin_group)).any()),
          "meet_only": int(len(enc_pairs - tr_pairs)),
          "exchange_only": int(len(tr_pairs - enc_pairs)),
          "both": int(len(enc_pairs & tr_pairs)),
          "max_edge": int(ed["excursions"].max()),
          "top_edge": "%s to %s" % (ed.loc[ed["excursions"].idxmax(), "origin_group"],
                                    ed.loc[ed["excursions"].idxmax(), "destination"])}

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Four panels: one '
         'group&#39;s weekly modularity through time with two of its weeks inset as '
         'association networks; the group-meeting and animal-transfer networks overlaid '
         'on one node layout; four concentration curves on shared axes; and the duration '
         'of all three states on one logarithmic hourly axis.">' % H]
    for letter, ci, ri in (("a", 0, 0), ("b", 1, 0), ("c", 0, 1), ("d", 1, 1)):
        o.append('<text class="pl" x="%.0f" y="%.1f">%s</text>'
                 % (COL_X[ci] - 44, ROW_Y[ri] - 16, letter))
    o.append(pa)
    o.append(pb)
    o.append(pd_svg)
    o.append(pe)
    o.append("</svg>")
    return "\n".join(o), sa, sb, sc, gs, ds, shown


def page():
    fig, sa, sb, sc, gs, ds, shown = build()
    hi, lo = shown.iloc[0], shown.iloc[1]
    return f"""<title>Three Routes Across a Boundary</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 3 &middot; draft</div>
  <h1>Three routes across a group boundary</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 3.</span> The three ways a group boundary becomes
      fluid, shown (a&ndash;b) and then compared (c&ndash;d). Colour is constant:
      <b>ochre</b> within-group modularity and splits, <b>blue</b> group meetings,
      <b>plum</b> individual transfer. Intervals are 95% cluster bootstraps ({N_BOOT}
      replicates) on the natural unit of replication &mdash; unit, dyad or animal &mdash;
      never on the row.
      <b>(a)</b> {sb["unit"]} week by week, {sb["window"][0]} to {sb["window"][1]}, over
      the weeks it carried {WELL_COVERED} or more collars. {sb["modular_weeks"]} of the
      {sb["weeks"]} weeks shown are modular, in {sb["phases"]} shaded phases of two
      weeks or more, the longest {sb["longest_phase"]} consecutive weeks, reaching
      {sb["max_mod"]:.2f} at most and dividing into as many as
      {sb["max_communities"]} communities, the largest holding as little as
      {sb["min_largest"]:.0%} of the group. The two insets are that group&rsquo;s own
      weeks, placed where the series leaves room and ruled back to the point each one
      illustrates: at
      Q = {hi["q"]:.2f} ({hi["n"]} animals) it is two communities that barely associate;
      at Q = {lo["q"]:.2f} the same two are still there but so heavily bridged that the
      group is nearly one piece &mdash; <b>which is what the lower half of this scale
      means</b>, and it is not the same thing as a group with no partition at all. The
      axis stops at {sb["window"][1]} for legibility; the record runs on to
      {sb["tail_end"]}, and in those further {sb["tail_weeks"]} weeks modularity never
      exceeds {sb["tail_max"]:.2f} ({sb["all_modular"]} modular weeks of
      {sb["all_weeks"]} across the whole record). Across the {sa["units"]} well-covered units,
      {sa["modular_units"]} reach modularity in some week and {sa["silent_units"]} never
      leave zero, but <b>{sa["top2_names"]} hold {sa["top2_weeks"]} of the
      {sa["modular_weeks"]} modular weeks</b> ({sa["top2_share"]:.0%}).
      <b>(b)</b> Both between-unit networks on one node layout, node area by median
      collar count. Blue chords are group meetings ({sc["dyads"]} dyads over
      {sc["units"]} units, {sc["bouts"]:,} bouts, width by hours together); plum chords
      are animals moving ({sc["edges"]} directed edges, {sc["excursions"]} excursions by
      {sc["animals"]} animals, {sc["reciprocal"]} of the edges reciprocal, the heaviest
      {sc["top_edge"]} carrying {sc["max_edge"]}). A meeting requires <b>two or more
      animals on both sides</b>: {sc["solo_dyad_days"]} dyad-days
      ({sc["solo_share"]:.0%} of all detected encounters, {sc["solo_dyads"]} dyads) had a
      single animal on one side and count as transfer, not meeting &mdash; their group
      centroids sit a median {sc["solo_distance"]:,} m apart against
      {sc["grp_distance"]} m for a real meeting. That separation is what makes the
      overlay readable: <b>{sc["meet_only"]} dyads meet without ever exchanging an
      animal, {sc["exchange_only"]} exchange without a recorded meeting, and
      {sc["both"]} do both</b> &mdash; meeting and exchanging are largely different
      relationships between the same units.
      <b>(c)</b> Concentration, four quantities against one diagonal with Gini and its
      interval in the panel. Encounter-hours per dyad is the most unequal
      ({gs["enc"]["gini"]:.2f} across {gs["enc"]["n"]} dyads), modular weeks per unit
      next ({gs["mod"]["gini"]:.2f}); excursions spread more evenly, over
      {gs["edge"]["n"]} transfer edges ({gs["edge"]["gini"]:.2f}) and
      {gs["animal"]["n"]} animals ({gs["animal"]["gini"]:.2f}, dashed).
      <b>(d)</b> Duration, all three at hourly resolution.
      {ds["enc"]["under_day"]:.0%} of meetings end within a day and
      {ds["enc"]["under_6h"]:.0%} within six hours (median {ds["enc"]["median"]:.0f} h);
      splits are shorter still ({ds["mod"]["median"]:.0f} h,
      {ds["mod"]["under_day"]:.0%} inside a day), which is how a group can be modular for
      seven consecutive weeks without being apart for a whole day. The excursion curve is
      an uninterrupted away-state run rather than a whole dispersal &mdash; a long
      absence is cut every time the hourly context returns to origin &mdash; so the
      audited nightly view (338 excursions, 22% past a week, out to 500 nights) remains
      the measure of how long leaving lasts. <b>The shaded strip below 2 h is a
      resolution floor</b>: {ds["mod"]["at_1h"]:.0%} of splits and
      {ds["exc"]["at_1h"]:.0%} of excursions are single hours against
      {ds["enc"]["at_1h"]:.0%} of meetings, because one animal's cluster membership
      flickers where a whole group's co-occurrence does not.
    </p>
  </div>

  <p class="note">
    Panel (a) is restricted to well-covered unit-weeks: below about {WELL_COVERED}
    collars a modularity estimator returns a single community regardless of the animals.
    Its two networks are matched on collar count, so the flat one cannot be dismissed as
    a thinner week; an edge is the share of hours two animals shared an association
    cluster among hours both were present, drawn at 0.05 and above, communities are
    greedy modularity communities, and positions come from a spring layout with each
    community&rsquo;s centre placed on a fixed circle, so distance on the page is
    indicative rather than metric. Panels (b) and (c) are not coverage-restricted, which
    is why node area carries coverage instead. The two chord types in (b) are mostly
    but not entirely independent evidence: 16 dyads carry both, and on
    Copper&ndash;Lilac (84% of encounter hours) and
    Maroon&ndash;Sapphire (91%) most of the meeting record falls inside a period when one
    of their animals was living in the other unit, so those two should not be cited as
    evidence that meeting and exchanging are separate channels. Transfer destinations are
    the modal other unit sharing the focal animal&rsquo;s clusters on that night; 86.3% of
    joined away-nights resolve and the rest are left unassigned, so (c) is a lower bound.
    All three curves in (e) come from <code>derive_level2_inputs.py</code> at hourly
    resolution with gaps of up to two hours bridged, the excursion curve held to the
    {ds["exc"]["clusters"]} animals of the audited axis-C cohort. Weekly metrics come from
    the frozen 1,924,104-row export. Generated by <code>build_figure3.py</code>;
    seed {SEED}.
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
