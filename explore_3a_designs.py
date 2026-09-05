"""Four designs for panel 3a, drawn at the size they would actually occupy.

Panel 3a has to do two jobs in a quarter of a figure: show one group's modularity through
time, and show what a modularity score looks like. The four options trade those jobs off
differently, so they are drawn here side by side at the real 262x194 box rather than
described, because the only question that matters is which one survives being small.

  1  LINE WITH CORNER INSETS (the current build). Series across the full box, two example
     networks dropped into the top corners over the data, numbered to their weeks.
     Cheapest in space; the insets sit on top of the series.
  2  SERIES ABOVE, NETWORKS BELOW. The same series, with the two networks under the axis
     directly beneath the week they come from. Nothing overlaps and the link between a
     week and its picture is positional rather than numbered. Costs ~70px of height.
  3  COMPOSITION RIBBON. Drops the modularity line as the primary mark and shows, per
     week, the share of the group in its largest community as a filled area, with a dot
     row for the number of communities. Modularity rides on top as a thin line. Answers
     "how split was it" directly rather than through a score; the two insets become
     illustrations rather than the only picture of structure.
  4  ILLUSTRATED AXIS. The four example weeks stacked as a column beside the modularity
     axis at their own Q values, so the y-axis itself is annotated: any point on the
     series can be read across to a picture. Uses all four examples instead of two.

Output: docs/framing_2026-09-04/figure3a_options.html
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from build_figure3 import inset_network
from build_level2_figure import CSS, MOD_EPS, axes, key_row, load_modularity, runs_of
from build_modularity_figure import NETS

OUT = Path("docs/framing_2026-09-04/figure3a_options.html")
EXAMPLE_UNIT = "Lilac"
BOX_W, BOX_H = 262.0, 194.0
SEED = 20260904


def series_data(sub):
    g = sub[sub["dynamic_social_unit"].eq(EXAMPLE_UNIT)].sort_values("period_start")
    lo, hi = g["period_start"].min(), g["period_start"].max()
    return g, lo, hi, max(1.0, (hi - lo).days)


def frame(g, lo, span, x0, y0, w, h, ylab="modularity", ytop=None):
    """The shared furniture: axes, date ends, and the modular phases shaded."""
    top = ytop or max(0.05, float(g["modularity"].max()))
    o = [axes(x0, y0, w, h, [(0.0, ""), (0.5, ""), (1.0, "")],
              [(0.0, "0"), (0.5, "%.2f" % (top / 2)), (1.0, "%.2f" % top)],
              "week", ylab)]
    for tv, tl, ta in ((0.0, str(lo.date()), "start"), (1.0, str(g["period_start"].max()
                                                                .date()), "end")):
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="%s">%s</text>'
                 % (x0 + tv * w, y0 + h + 11, ta, tl))
    mod = (g["modularity"] > MOD_EPS).to_numpy()
    wk2date = dict(zip(g["weekno"].to_numpy(), g["period_start"]))
    for a, b in [(a, b) for a, b in runs_of(mod, g["weekno"].to_numpy())
                 if b - a + 1 >= 2]:
        xa = x0 + ((wk2date[a] - lo).days / span) * w
        xb = x0 + ((wk2date[b] - lo).days / span) * w + 3
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--c4)" '
                 'opacity=".11"/>' % (xa - 1.5, y0, max(3.0, xb - xa), h))
    return o, top


def q_line(g, lo, span, x0, y0, w, h, top, colour="var(--c4)", width=1.7):
    pts = ["%.1f,%.1f" % (x0 + ((r.period_start - lo).days / span) * w,
                          y0 + h - (r.modularity / top) * h) for r in g.itertuples()]
    return ('<polyline points="%s" fill="none" stroke="%s" stroke-width="%.1f"/>'
            % (" ".join(pts), colour, width))


def largest_line(g, lo, span, x0, y0, w, h):
    pts = ["%.1f,%.1f" % (x0 + ((r.period_start - lo).days / span) * w,
                          y0 + h - float(r.largest_community_fraction) * h)
           for r in g.itertuples()]
    return ('<polyline points="%s" fill="none" stroke="var(--c6)" stroke-width="1.2" '
            'opacity=".75"/>' % " ".join(pts))


def right_axis(x0, y0, w, h):
    o = ['<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--c6)" '
         'opacity=".55"/>' % (x0 + w, y0, x0 + w, y0 + h)]
    for tv, tl in ((0.0, "0"), (0.5, "0.5"), (1.0, "1")):
        o.append('<text class="ts" x="%.1f" y="%.1f" fill="var(--c6)">%s</text>'
                 % (x0 + w + 4, y0 + h - tv * h + 3, tl))
    return o


# ------------------------------------------------------------------ the options
def option1(g, lo, span, net, shown, x0, y0):
    """Line with the two example networks inset in the top corners."""
    w, h = BOX_W, BOX_H
    o, top = frame(g, lo, span, x0, y0, w, h)
    o.append(largest_line(g, lo, span, x0, y0, w, h))
    o.append(q_line(g, lo, span, x0, y0, w, h, top))
    o.extend(right_axis(x0, y0, w, h))
    size = 62.0
    for i, r in enumerate(shown.itertuples()):
        wx = x0 + max(0.0, min(1.0, (pd.Timestamp(r.period_start) - lo).days / span)) * w
        ix = x0 + 3 if i == 0 else x0 + w - size - 3
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink-2)" '
                 'stroke-width=".8" stroke-dasharray="2 2"/>' % (wx, y0, wx, y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" '
                 'fill="var(--ink-2)">%d</text>' % (wx, y0 - 3, i + 1))
        o.append(inset_network(net, r.week, ix, y0 + 4, size,
                               "%d · Q = %.2f" % (i + 1, r.q)))
    return "\n".join(o), h + 34


def option2(g, lo, span, net, shown, x0, y0):
    """Series above, the two networks below the axis under their own weeks."""
    w, h = BOX_W, 118.0
    o, top = frame(g, lo, span, x0, y0, w, h)
    o.append(largest_line(g, lo, span, x0, y0, w, h))
    o.append(q_line(g, lo, span, x0, y0, w, h, top))
    o.extend(right_axis(x0, y0, w, h))
    size, ty = 62.0, y0 + h + 30
    for i, r in enumerate(shown.itertuples()):
        wx = x0 + max(0.0, min(1.0, (pd.Timestamp(r.period_start) - lo).days / span)) * w
        ix = min(max(x0, wx - size / 2), x0 + w - size)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink-2)" '
                 'stroke-width=".8" stroke-dasharray="2 2"/>' % (wx, y0, wx, y0 + h))
        o.append('<path d="M %.1f %.1f L %.1f %.1f L %.1f %.1f" fill="none" '
                 'stroke="var(--ink-2)" stroke-width=".8" stroke-dasharray="2 2"/>'
                 % (wx, y0 + h, wx, ty - 8, ix + size / 2, ty))
        o.append(inset_network(net, r.week, ix, ty, size, "Q = %.2f" % r.q))
    return "\n".join(o), (ty + size + 22) - y0


def option3(g, lo, span, net, shown, x0, y0):
    """Composition ribbon: the share of the group in its largest community."""
    w, h = BOX_W, BOX_H
    o = [axes(x0, y0, w, h, [(0.0, ""), (0.5, ""), (1.0, "")],
              [(0.0, "0"), (0.5, "0.5"), (1.0, "1")], "week",
              "share in largest community")]
    for tv, tl, ta in ((0.0, str(lo.date()), "start"),
                       (1.0, str(g["period_start"].max().date()), "end")):
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="%s">%s</text>'
                 % (x0 + tv * w, y0 + h + 11, ta, tl))
    xs = [x0 + ((r.period_start - lo).days / span) * w for r in g.itertuples()]
    ys = [y0 + h - float(r.largest_community_fraction) * h for r in g.itertuples()]
    # everything above the line is the rest of the group, so fill it: the filled band IS
    # the fraction of the group outside its main body, which is the quantity of interest
    poly = (["%.1f,%.1f" % (xs[0], y0)]
            + ["%.1f,%.1f" % (a, b) for a, b in zip(xs, ys)]
            + ["%.1f,%.1f" % (xs[-1], y0)])
    o.append('<polygon points="%s" fill="var(--c6)" opacity=".20" stroke="none"/>'
             % " ".join(poly))
    o.append('<polyline points="%s" fill="none" stroke="var(--c6)" '
             'stroke-width="1.5"/>' % " ".join("%.1f,%.1f" % (a, b)
                                               for a, b in zip(xs, ys)))
    top = max(0.05, float(g["modularity"].max()))
    o.append(q_line(g, lo, span, x0, y0, w, h, top, colour="var(--c4)", width=1.2))
    # a dot row for the number of communities, since the ribbon cannot show it
    for x, r in zip(xs, g.itertuples()):
        for k in range(int(r.n_communities)):
            o.append('<circle cx="%.1f" cy="%.1f" r="1.1" fill="var(--ink-2)" '
                     'opacity=".65"/>' % (x, y0 + h + 27 + k * 3.4))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end">communities</text>'
             % (x0 - 5, y0 + h + 31))
    size = 58.0
    for i, r in enumerate(shown.itertuples()):
        ix = x0 + 3 if i == 0 else x0 + w - size - 3
        o.append(inset_network(net, r.week, ix, y0 + h - size - 14, size,
                               "Q = %.2f" % r.q))
    return "\n".join(o), h + 52


def option4(g, lo, span, net, meta, x0, y0):
    """The four example weeks stacked beside the axis at their own Q values."""
    col, gap = 52.0, 4.0
    w, h = BOX_W - col - 30, BOX_H
    sx = x0 + col + 30
    o, top = frame(g, lo, span, sx, y0, w, h)
    o.append(largest_line(g, lo, span, sx, y0, w, h))
    o.append(q_line(g, lo, span, sx, y0, w, h, top))
    o.extend(right_axis(sx, y0, w, h))
    n = len(meta)
    for i, r in enumerate(meta.itertuples()):
        iy = y0 + i * (col + gap + 8)
        o.append(inset_network(net, r.week, x0, iy, col, "Q = %.2f" % r.q))
        # tie each picture to the height on the axis that it represents
        ty = y0 + h - (max(0.0, r.q) / top) * h
        o.append('<path d="M %.1f %.1f L %.1f %.1f" fill="none" stroke="var(--ink-2)" '
                 'stroke-width=".7" stroke-dasharray="2 2"/>'
                 % (x0 + col + 1, iy + col / 2, sx - 1, ty))
    return "\n".join(o), max(h + 34, n * (col + gap + 8) + 10)


def build():
    sub, _ = load_modularity()
    g, lo, hi, span = series_data(sub)
    net = pd.read_csv(NETS, parse_dates=["period_start"])
    meta = (net[net["kind"].eq("node")]
            .groupby(["week", "unit", "period_start"])
            .agg(n=("animal_id", "nunique"), comms=("community", "nunique"))
            .reset_index())
    q = sub.set_index(["dynamic_social_unit", "period_start"])["modularity"].to_dict()
    meta["q"] = [q.get((u, pd.Timestamp(d)), float("nan"))
                 for u, d in zip(meta["unit"], meta["period_start"])]
    meta = meta.sort_values("q", ascending=False).reset_index(drop=True)
    shown = meta.iloc[[0, -1]].reset_index(drop=True)

    titles = ["1 · line with corner insets",
              "2 · series above, networks below",
              "3 · composition ribbon",
              "4 · illustrated axis"]
    notes = ["current build; cheapest in space, insets sit on the data",
             "nothing overlaps; costs about 70px of height",
             "shows how split, not just how modular; ribbon = share outside the "
             "largest community",
             "all four examples; the y-axis itself is annotated"]
    xs = (58.0, 388.0)
    y1 = 76.0
    p1, h1 = option1(g, lo, span, net, shown, xs[0], y1)
    p2, h2 = option2(g, lo, span, net, shown, xs[1], y1)
    y2 = y1 + max(h1, h2) + 78
    p3, h3 = option3(g, lo, span, net, shown, xs[0], y2)
    p4, h4 = option4(g, lo, span, net, meta, xs[1], y2)
    H = y2 + max(h3, h4) + 54

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Four alternative designs '
         'for the modularity panel, each drawn at the size it would occupy in the '
         'figure.">' % H]
    for i, (t, note) in enumerate(zip(titles, notes)):
        x = xs[i % 2] - 44
        y = (y1 if i < 2 else y2) - 34
        o.append('<text class="ph" x="%.0f" y="%.1f">%s</text>' % (x, y, t))
        o.append('<text class="ts" x="%.0f" y="%.1f">%s</text>' % (x, y + 11, note))
    o.extend([p1, p2, p3, p4])
    o.append("</svg>")
    return "\n".join(o)


def page():
    return f"""<title>Panel 3a Design Options</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Panel 3a &middot; options</div>
  <h1>Four ways to draw one group's modularity</h1>

  <div class="plate">
    <div class="scroll">{build()}</div>
    <p class="cap">
      Each option is drawn at the {BOX_W:.0f}&times;{BOX_H:.0f} box it would occupy as a
      quarter of Figure 3, because the only question worth asking is which survives being
      small. <b>1</b> is the current build. <b>2</b> removes every overlap and makes the
      week-to-picture link positional instead of numbered, at the cost of roughly seventy
      pixels of height, which Figure 3 would have to find somewhere. <b>3</b> changes the
      question: the filled band is the share of the group sitting <i>outside</i> its
      largest community, so it answers &ldquo;how split was it&rdquo; directly rather
      than through a score, with modularity as a thin line and a dot row for the number
      of communities. <b>4</b> spends the horizontal space on all four example weeks
      stacked beside the axis at their own Q values, turning the y-axis into a legend a
      reader can read across to &mdash; the series gets narrower in exchange.
    </p>
  </div>

  <p class="note">
    Same data throughout: {EXAMPLE_UNIT}&rsquo;s weekly modularity over its well-covered
    weeks, and the example networks built by <code>derive_level2_inputs.py</code>. Ochre
    is modularity, teal the largest sub-community&rsquo;s share of the group, shading a
    modular phase of two weeks or more. Generated by <code>explore_3a_designs.py</code>;
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
