"""Two ways to make Figure 1's panels share a scale, drawn side by side.

Panel (a) and panel (b) currently disagree on two things at once, and only one of them
can be fixed without cost:

  A  ONE TIME SCALE. Panel (a)'s stage inset and panel (b)'s raster are both hourly
     strips, drawn at wildly different hours-per-pixel: the inset compresses fourteen
     days into about 309 px, the raster spreads two days across 648. Matching them means
     the inset shows ONE DAY, because 309 px at the raster's pixels-per-hour buys no
     more. The claim that the dynamic unit is stable for a fortnight then lives in the
     caption rather than in the picture.

  B  ONE PANEL WIDTH. The map occupies 292 px of a 648 px panel, so (a) reads as half a
     panel next to (b). Growing the map to 430 px gives the two panels equal weight and
     puts more space on the thing the panel is about -- at the cost of squeezing the
     legend and the inset into the remaining strip, and leaving the time scales as
     different as they are now.

Both are rendered at full size below with panel (b) underneath, because the mismatch is
only visible when the two are on the same page.

Output: docs/framing_2026-09-04/fig1_scaling_options.html
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import build_event_examples_figure as F

OUT = Path("docs/framing_2026-09-04/fig1_scaling_options.html")


def variant(d, got, x0, y, label, side, inset_days):
    """One layout: the map at `side` px, its inset spanning `inset_days`."""
    info = {}

    def draw_inset(ix, iy, iw):
        head = ('<text class="ts" x="%.1f" y="%.1f">one animal, %d day%s, hourly</text>'
                % (ix, iy - 6, inset_days, "" if inset_days == 1 else "s"))
        svg, h, i = F.stages_panel(d, ix, iy, iw, compact=True, days=inset_days)
        info.update(i)
        return head + chr(10) + svg, h

    o = ['<text class="pl" x="%.0f" y="%.1f">a</text>' % (x0 - 20, y - 12)]
    gsvg, gh, ginfo = F.gps_panel(d, x0, y, 648.0, inset=draw_inset, side=side)
    o.append(gsvg)
    y2 = y + gh + 40
    o.append('<text class="pl" x="%.0f" y="%.1f">b</text>' % (x0 - 20, y2 - 12))
    rs, rh, _ = F.raster(got, x0, y2, 648.0)
    o.append(rs)
    return chr(10).join(o), (y2 + rh) - y, ginfo, info


def build():
    cols = ["window_start", "animal_id", "origin_group", "dynamic_social_unit",
            "social_context", "association_event_id", "temp_group_size", "is_observed",
            "is_carried_night", "assigned_social_unit", "longitude", "latitude"]
    d = pd.read_csv(F.NARROW, usecols=cols, parse_dates=["window_start"])
    s1 = pd.read_csv(F.STAGE1, parse_dates=["start_hour"])
    idx = F.build_index(d)
    got = F.pick_all_patterns(d, s1, idx)
    if got is None:
        raise SystemExit("no gapless window found")

    x0 = 26.0
    ya = 96.0
    pa, ha, ga, ia = variant(d, got, x0, ya, "A", side=292.0, inset_days=1)
    yb = ya + ha + 96
    pb, hb, gb, ib = variant(d, got, x0, yb, "B", side=430.0, inset_days=F.STAGE_DAYS)
    H = yb + hb + 40

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Two scaling options for '
         'Figure 1, each drawn in full: one matching the hourly time scale between the '
         'panels, one matching the panel widths.">' % H]
    for title, note, y in (
            ("A · one time scale",
             "inset redrawn at the raster's pixels-per-hour, so it shows one day",
             ya - 42),
            ("B · one panel width",
             "map grown to 430 px; the inset keeps its fortnight and its own scale",
             yb - 42)):
        o.append('<text class="ph" x="%.0f" y="%.1f">%s</text>' % (x0 - 20, y, title))
        o.append('<text class="ts" x="%.0f" y="%.1f">%s</text>'
                 % (x0 - 20, y + 12, note))
    o.extend([pa, pb])
    o.append("</svg>")
    return chr(10).join(o), ga, ia, ib


def page():
    fig, ga, ia, ib = build()
    return f"""<title>Figure 1 Scaling Options</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>{F.CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 1 &middot; scaling options</div>
  <h1>Two ways to make the panels agree</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      Both options use the same frame &mdash; {ga["hour"]}, {ga["animals"]} animals in
      {ga["clusters"]} clusters, {ga["off_origin"]} away from its origin group &mdash; and
      the same raster underneath, so the only difference is the geometry.
      <b>A</b> matches the <b>time</b> scale: the stage inset is redrawn at the raster's
      pixels-per-hour, which buys {ia["days"]} day rather than {F.STAGE_DAYS}. The two
      hourly strips can then be compared directly, but the inset no longer shows the
      dynamic unit holding across a fortnight &mdash; that claim moves to the caption.
      <b>B</b> matches the <b>panel</b>: the map grows from 292 to 430 px so (a) and (b)
      carry equal visual weight and the spatial structure is easier to read, while the
      inset keeps its {ib["days"]} days at its own compressed scale and the legend gets a
      narrower column. The time-scale mismatch stays, and has to be stated.
      <b>The two cannot both be had:</b> 309 px of inset at the raster's resolution is a
      day, and a fortnight at the raster's resolution is 4,700 px.
    </p>
  </div>

  <p class="note">
    Rendered by <code>explore_fig1_scaling.py</code>, which imports the panel functions
    from <code>build_event_examples_figure.py</code> unchanged, so what is drawn here is
    what the figure would draw. Seed {F.SEED}.
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
