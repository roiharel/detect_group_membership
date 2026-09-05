"""Methods figure: how GPS positions become dynamic social units.

Three panels.

  (a) The chain, with the real cardinality at every stage. Positions cluster into an
      hourly spatial cluster; clusters stitch across consecutive hours into an
      association event; origin labels turn a cluster into a social context; contexts
      resolve into an immediate assigned unit; and a sustained-association rule
      collapses those into a dynamic unit.

  (b) What the sustained rule actually does, as spell lengths. An assigned-unit spell
      has a median of 7 hours; a dynamic-unit spell 4,323. The rule discards everything
      transient -- 11.5% of animal-hours sit in an assigned spell shorter than a day,
      against 0.0% for dynamic.

  (c) The same resolution for one animal over two weeks at hourly resolution: cluster
      size moves, social context moves, the assigned label flickers between its own unit
      and `mixed:` composites, and the dynamic unit does not move at all.

Output: docs/framing_2026-09-04/pipeline_figure.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
OUT = Path("docs/framing_2026-09-04/pipeline_figure.html")

FOCAL_ANIMAL = "24AA01_5O8B"
FOCAL_START = pd.Timestamp("2025-06-20")
FOCAL_DAYS = 14
N_BOOT = 600
SEED = 20260904

CONTEXT_COLOUR = {
    "with_origin": "var(--n2)",
    "mixed_with_origin_present": "var(--c1)",
    "other": "var(--c5)",
    "mixed_without_origin_unit": "var(--c5b)",
    "isolated": "var(--alone)",
    "mixed_tie_or_unclear": "var(--muted)",
}
CONTEXT_LABEL = {
    "with_origin": "with origin unit",
    "mixed_with_origin_present": "merged, origin unit present",
    "other": "with another unit",
    "mixed_without_origin_unit": "merged, origin unit absent",
    "isolated": "alone",
    "mixed_tie_or_unclear": "unclear",
}
UNIT_COLOUR = {
    "Bronze": "#a8703a", "Emerald": "#2e8b6b", "Chartreuse": "#8a9a2b",
    "Purple": "#7a5aa0", "Copper": "#b5702f", "Lilac": "#9b8ac4",
    "Jade": "#2f8f7a", "Lapis": "#2a5fa8", "Magenta": "#b0448c",
    "Maroon": "#9c3a4a", "Periwinkle": "#6d82c4", "Sapphire": "#2a6fb0",
    "Turquoise": "#2a9aa5", "RubyRunners": "#b03a4a", "FireOpal": "#c46a3a",
    "Obsidian": "#4a4a52", "Pearl": "#9aa0a8", "Red": "#b83a3a",
    "Green": "#3f8f45", "Gold": "#b8912a", "Mint": "#4fb08a",
    "SneakySilver": "#8f97a0", "PhantomWest": "#6a7480",
    "FenceJumpers": "#7a8a5a", "TrickyTeal": "#2a8a95",
    "LapisSplinter": "#4a7fc0",
}


# ------------------------------------------------------------------ data
def load():
    cols = ["window_start", "animal_id", "origin_group", "assigned_social_unit",
            "dynamic_social_unit", "social_context", "association_event_id",
            "temp_group_size", "is_observed", "is_carried_night",
            "is_local_2h_supported"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
# A row's state is KNOWN if it has a position from a real fix -- exactly at the
# hour, borrowed from a fix within 60 min (`local_2h`), or carried across the
# night. Omitting local_2h while accepting carried_night was indefensible and
# became visible only when coverage improved: at 17:00 the 2026-09-05 build
# supplies 106,147 local_2h rows where the frozen build had carried_night, so the
# old predicate silently discarded almost the whole hour. Interpolated rows are
# deliberately NOT known: their position is inferred, not observed.
    d["known"] = (d["is_observed"] | d["is_carried_night"]
                  | d["is_local_2h_supported"])
    return d


def stage_counts(d: pd.DataFrame) -> list[dict]:
    origins = set(d["origin_group"].unique())
    labels = d["assigned_social_unit"].dropna().unique()
    mixed = [v for v in labels if str(v).startswith("mixed:")]
    return [
        dict(stage="GPS fixes, hourly",
             rule="one row per animal per hour",
             count="%s animal-hours" % f"{len(d):,}",
             detail="%d animals x %s hours" % (d["animal_id"].nunique(),
                                               f"{d['window_start'].nunique():,}")),
        dict(stage="Spatial cluster",
             rule="adaptive clustering, 120-900 m edge range",
             count="median %d animals" % d["temp_group_size"].median(),
             detail="one cluster in one hour (temp_group_id); max %d"
                    % d["temp_group_size"].max()),
        dict(stage="Association event",
             rule="clusters stitched across consecutive hours",
             count="%s events" % f"{d['association_event_id'].nunique():,}",
             detail="the Stage-1 encounter identity"),
        dict(stage="Social context",
             rule="cluster composition read against origin labels",
             count="%d contexts" % d["social_context"].nunique(),
             detail="with origin / merged / another unit / alone / unclear"),
        dict(stage="Assigned unit",
             rule="immediate membership, recomputed every hour",
             count="%d labels" % len(labels),
             detail="%d units + ISOLATED + %d `mixed:` composites"
                    % (len(origins), len(mixed))),
        dict(stage="Dynamic unit",
             rule="reassign only on 7 days of sustained association",
             count="%d units" % d["dynamic_social_unit"].nunique(),
             detail="the membership used for analysis"),
    ]


def spells(d: pd.DataFrame, col: str) -> tuple[np.ndarray, np.ndarray]:
    """Length in hours of each run of a constant label, and its animal."""
    x = d.sort_values(["animal_id", "window_start"])
    change = (x[col] != x.groupby("animal_id")[col].shift()) | (
        x["animal_id"] != x["animal_id"].shift())
    x = x.assign(spell=change.cumsum())
    g = x.groupby("spell").agg(n=("window_start", "size"),
                               animal=("animal_id", "first"))
    return g["n"].to_numpy(dtype=float), g["animal"].to_numpy()


# ------------------------------------------------------------------ panel a
def panel_a(rows, y0):
    x0, w, rowh = 8.0, 660.0, 40.0
    o = []
    for i, r in enumerate(rows):
        y = y0 + i * rowh
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="30" fill="var(--paper)" '
                 'stroke="var(--rule)" stroke-width=".9"/>' % (x0, y, w))
        o.append('<rect x="%.1f" y="%.1f" width="3.2" height="30" fill="%s"/>'
                 % (x0, y, "var(--c1)" if i < 3 else "var(--c5)"))
        o.append('<text class="sn" x="%.1f" y="%.1f">%s</text>'
                 % (x0 + 12, y + 13, r["stage"]))
        o.append('<text class="ss" x="%.1f" y="%.1f">%s</text>'
                 % (x0 + 12, y + 24, r["rule"]))
        o.append('<text class="sc" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 + w - 10, y + 13, r["count"]))
        o.append('<text class="ss" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 + w - 10, y + 24, r["detail"]))
        if i < len(rows) - 1:
            xm = x0 + 26
            o.append('<path d="M %.1f %.1f L %.1f %.1f" stroke="var(--muted)" '
                     'stroke-width="1"/>' % (xm, y + 30, xm, y + rowh))
            o.append('<path d="M %.1f %.1f l -2.6 -3.4 l 5.2 0 z" fill="var(--muted)"/>'
                     % (xm, y + rowh + 0.6))
    return "\n".join(o), rowh * len(rows) + 4


# ------------------------------------------------------------------ panel b
def panel_b(d, y0, rng):
    x0, w, h = 62.0, 402.0, 128.0
    lo, hi = 1.0, 20000.0
    xs = np.geomspace(lo, hi, 80)

    def px(v):
        v = min(max(v, lo), hi)
        return x0 + (math.log10(v) - math.log10(lo)) / (
            math.log10(hi) - math.log10(lo)) * w

    def py(p):
        return y0 + h - p * h

    series = [("Assigned unit, hourly", "assigned_social_unit", "var(--c5)"),
              ("Dynamic unit, 7-day rule", "dynamic_social_unit", "var(--c1)")]
    o, stats = [], {}
    for t, lab in ((1, "1 h"), (24, "1 d"), (168, "1 wk"), (730, "1 mo"),
                   (8760, "1 yr")):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (px(t), y0, px(t), y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (px(t), y0 + h + 12, lab))
    for p in (0, 0.5, 1.0):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0, py(p), x0 + w, py(p)))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 - 6, py(p) + 3, "0" if p == 0 else ("1" if p == 1 else "0.5")))
    for label, col, colour in series:
        v, who = spells(d, col)
        keys, inv = np.unique(who, return_inverse=True)
        buckets = [v[inv == i] for i in range(len(keys))]
        point = np.array([(v > x).sum() / len(v) for x in xs])
        draws = np.empty((N_BOOT, len(xs)))
        for b in range(N_BOOT):
            pick = rng.integers(0, len(keys), len(keys))
            sam = np.concatenate([buckets[i] for i in pick])
            sam.sort()
            draws[b] = 1.0 - np.searchsorted(sam, xs, side="right") / len(sam)
        blo = np.percentile(draws, 2.5, axis=0)
        bhi = np.percentile(draws, 97.5, axis=0)
        fwd = " ".join("%.1f,%.1f" % (px(x), py(p)) for x, p in zip(xs, bhi))
        revp = " ".join("%.1f,%.1f" % (px(x), py(p))
                        for x, p in zip(xs[::-1], blo[::-1]))
        o.append('<polygon points="%s %s" fill="%s" opacity=".16"/>'
                 % (fwd, revp, colour))
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.9"/>'
                 % (" ".join("%.1f,%.1f" % (px(x), py(p))
                             for x, p in zip(xs, point)), colour))
        med = float(np.median(v))
        stats[label] = {"spells": int(len(v)), "median": med,
                        "p90": float(np.quantile(v, .9)),
                        "under_day": float((v < 24).mean())}
        o.append('<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>'
                 % (px(med), py(0.5), colour))
    ly = y0 + 12
    for label, col, colour in series:
        st = stats[label]
        o.append('<rect x="%.1f" y="%.1f" width="8" height="8" fill="%s"/>'
                 % (x0 + w + 14, ly - 7, colour))
        o.append('<text class="sn" x="%.1f" y="%.1f">%s</text>'
                 % (x0 + w + 26, ly, label))
        o.append('<text class="ss" x="%.1f" y="%.1f">%s spells, median %s h</text>'
                 % (x0 + w + 26, ly + 11, f"{st['spells']:,}",
                    f"{st['median']:,.0f}"))
        ly += 30
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + h, x0 + w, y0 + h))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0, x0, y0 + h))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">spell length, '
             'log scale</text>' % (x0 + w / 2, y0 + h + 28))
    o.append('<text class="ts" transform="translate(%.1f %.1f) rotate(-90)" '
             'text-anchor="middle">P(spell &gt; x)</text>' % (x0 - 30, y0 + h / 2))
    return "\n".join(o), h + 36, stats


# ------------------------------------------------------------------ panel c
def panel_c(d, y0):
    t0 = FOCAL_START
    t1 = t0 + pd.Timedelta(days=FOCAL_DAYS)
    w = d[d["animal_id"].eq(FOCAL_ANIMAL) & d["window_start"].between(t0, t1)
          & d["known"]].sort_values("window_start")
    if w.empty:
        return "", 0.0, {}
    hours = list(w["window_start"])
    x0, width = 96.0, 566.0
    cw = width / len(hours)
    unit = w["dynamic_social_unit"].iloc[0]
    origin = w["origin_group"].iloc[0]

    o = []
    y = y0
    # strip 1: cluster size as a thin area
    sh = 30.0
    sizes = w["temp_group_size"].to_numpy(dtype=float)
    top = max(1.0, sizes.max())
    pts = " ".join("%.2f,%.2f" % (x0 + i * cw, y + sh - (v / top) * sh)
                   for i, v in enumerate(sizes))
    o.append('<polygon points="%.2f,%.2f %s %.2f,%.2f" fill="var(--c1)" opacity=".3"/>'
             % (x0, y + sh, pts, x0 + width, y + sh))
    o.append('<polyline points="%s" fill="none" stroke="var(--c1)" stroke-width="1"/>'
             % pts)
    o.append('<text class="sn" x="%.1f" y="%.1f" text-anchor="end">cluster size</text>'
             % (x0 - 8, y + 12))
    o.append('<text class="ss" x="%.1f" y="%.1f" text-anchor="end">1 &ndash; %d</text>'
             % (x0 - 8, y + 23, int(top)))
    y += sh + 10

    def strip(values, label, colour_of, note):
        nonlocal y
        h = 15.0
        run_start, run_val = 0, values[0]
        for i in range(1, len(values) + 1):
            if i == len(values) or values[i] != run_val:
                o.append('<rect x="%.2f" y="%.1f" width="%.2f" height="%.1f" '
                         'fill="%s"/>' % (x0 + run_start * cw, y,
                                          max(0.4, (i - run_start) * cw), h,
                                          colour_of(run_val)))
                if i < len(values):
                    run_start, run_val = i, values[i]
        o.append('<text class="sn" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 - 8, y + 7, label))
        o.append('<text class="ss" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 - 8, y + 17, note))
        y += h + 12

    ctx = list(w["social_context"])
    strip(ctx, "social context",
          lambda v: CONTEXT_COLOUR.get(v, "var(--muted)"),
          "%d distinct" % len(set(ctx)))
    asg = list(w["assigned_social_unit"])
    n_mixed = sum(1 for v in asg if str(v).startswith("mixed:"))

    def asg_colour(v):
        v = str(v)
        if v.startswith("mixed:"):
            return "var(--c1)"
        if v == "ISOLATED":
            return "var(--alone)"
        return UNIT_COLOUR.get(v, "var(--n2)")

    strip(asg, "assigned unit", asg_colour,
          "%d labels, %d h mixed" % (len(set(asg)), n_mixed))
    dyn = list(w["dynamic_social_unit"])
    strip(dyn, "dynamic unit",
          lambda v: UNIT_COLOUR.get(str(v), "var(--c1)"),
          "%d label" % len(set(dyn)))

    for k in range(0, FOCAL_DAYS + 1, 2):
        xx = x0 + (k / FOCAL_DAYS) * width
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (xx, y0, xx, y - 12))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (xx, y + 2, (t0 + pd.Timedelta(days=k)).strftime("%d %b")))
    y += 10
    used = []
    for v in ctx:
        if v not in used:
            used.append(v)
    lx, ly = x0, y + 8
    for v in used:
        lab = CONTEXT_LABEL.get(v, v)
        iw = 9 + len(lab) * 4.3 + 12
        if lx + iw > x0 + width:
            lx, ly = x0, ly + 11
        o.append('<rect x="%.1f" y="%.1f" width="6" height="6" fill="%s"/>'
                 % (lx, ly, CONTEXT_COLOUR.get(v, "var(--muted)")))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (lx + 9, ly + 5.5, lab))
        lx += iw
    return "\n".join(o), (ly + 16) - y0, {
        "animal": FOCAL_ANIMAL.split("_")[0], "unit": unit, "origin": origin,
        "hours": len(hours), "contexts": len(set(ctx)),
        "assigned_labels": len(set(asg)), "mixed_hours": n_mixed,
        "dynamic_labels": len(set(dyn)),
    }


def build():
    d = load()
    rng = np.random.default_rng(SEED)
    rows = stage_counts(d)

    ya = 62.0
    pa, ha = panel_a(rows, ya)
    yb = ya + ha + 56
    pb, hb, sb = panel_b(d, yb, rng)
    yc = yb + hb + 62
    pc, hc, sc = panel_c(d, yc)
    H = yc + hc + 22

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Methods figure: the chain '
         'from hourly GPS fixes through spatial clusters, association events and social '
         'contexts to assigned and dynamic social units; the spell-length distributions '
         'of the assigned and dynamic labels; and the same resolution for one animal '
         'over two weeks.">' % H]
    o.append('<text class="pl" x="0" y="%.1f">a</text>' % (ya - 18))
    o.append('<text class="ph" x="16" y="%.1f">From fixes to units</text>' % (ya - 18))
    o.append(pa)
    o.append('<text class="pl" x="0" y="%.1f">b</text>' % (yb - 18))
    o.append('<text class="ph" x="16" y="%.1f">What the 7-day rule discards</text>'
             % (yb - 18))
    o.append(pb)
    o.append('<text class="pl" x="0" y="%.1f">c</text>' % (yc - 18))
    o.append('<text class="ph" x="16" y="%.1f">One animal, two weeks, every stage</text>'
             % (yc - 18))
    o.append(pc)
    o.append("</svg>")
    return "\n".join(o), rows, sb, sc


CSS = """
:root{--ground:#f4f5f3;--paper:#fff;--ink:#15191a;--ink-2:#3a4244;--muted:#6b7476;
--rule:#d5dad7;--rule-soft:#e6eae7;--n2:#c3cbc8;--alone:#4a5356;
--c1:#1f6f8b;--c4:#8a5a1f;--c5:#6b4a7a;--c5b:#a58cba;
--sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;
--serif:"Newsreader",Georgia,serif;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0f1213;
--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;--muted:#8f9998;--rule:#2c3435;
--rule-soft:#242c2d;--n2:#3f4a4b;--alone:#a7b0b2;--c1:#74b6d0;--c4:#d3a061;
--c5:#b494c4;--c5b:#7d6291;}}
:root[data-theme="dark"]{--ground:#0f1213;--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;
--muted:#8f9998;--rule:#2c3435;--rule-soft:#242c2d;--n2:#3f4a4b;--alone:#a7b0b2;
--c1:#74b6d0;--c4:#d3a061;--c5:#b494c4;--c5b:#7d6291;}
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
.sn{font-size:9.6px;fill:var(--ink);font-weight:600}
.ss{font-size:8.4px;fill:var(--muted)}
.sc{font-size:9.6px;fill:var(--ink);font-weight:600;font-family:var(--mono)}
.ts{font-size:8.4px;fill:var(--muted)}
.ph{font-size:12px;fill:var(--ink);font-weight:600}
.pl{font-size:13px;font-weight:700;fill:var(--ink)}
.g{stroke:var(--rule-soft);stroke-width:.8}
.ax{stroke:var(--rule);stroke-width:1}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def page():
    fig, rows, sb, sc = build()
    a = sb["Assigned unit, hourly"]
    dy = sb["Dynamic unit, 7-day rule"]
    ratio = dy["median"] / a["median"]
    return f"""<title>From Fixes To Social Units</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Methods figure &middot; draft</div>
  <h1>How GPS fixes become dynamic social units</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Methods figure.</span> The clustering chain and what each step
      keeps.
      <b>(a)</b> Positions are clustered within each hour by adaptive clustering over a
      120&ndash;900 m edge range, giving one spatial cluster per hour
      (<code>temp_group_id</code>, median {rows[1]["count"].split()[1]} animals). Consecutive
      hours of the same cluster stitch into an association event
      (<code>association_event_id</code>, {rows[2]["count"]}), which is the Stage-1 encounter
      identity. Reading a cluster's composition against the animals' origin labels gives a
      social context; contexts resolve into an <em>assigned</em> unit recomputed every hour;
      and a sustained-association rule collapses those into the <em>dynamic</em> unit used for
      analysis. The assigned label carries {rows[4]["detail"]} &mdash; a merge shows up as a
      <code>mixed:</code> composite rather than as one of the units.
      <b>(b)</b> Spell lengths, log axis, with 95% cluster-bootstrap bands resampling animals.
      An assigned-unit spell has a median of <b>{a["median"]:,.0f} h</b> across
      {a["spells"]:,} spells; a dynamic-unit spell <b>{dy["median"]:,.0f} h</b> across
      {dy["spells"]:,} &mdash; a <b>{ratio:,.0f}-fold</b> difference.
      {100 * a["under_day"]:.1f}% of animal-hours sit in an assigned spell shorter than a day,
      against {100 * dy["under_day"]:.1f}% for dynamic. The rule is not a smoothing choice at
      the margin; it discards essentially all of the hourly structure.
      <b>(c)</b> The same resolution for one animal &mdash; {sc["animal"]} of {sc["origin"]}
      &mdash; over {FOCAL_DAYS} days at hourly resolution. Cluster size moves between 1 and
      {rows[1]["detail"].split()[-1]}; social context takes {sc["contexts"]} distinct values;
      the assigned label takes {sc["assigned_labels"]}, spending {sc["mixed_hours"]} h in
      <code>mixed:</code> composites; and the dynamic unit takes
      {sc["dynamic_labels"]} value for the whole window.
    </p>
  </div>

  <p class="note">
    Panel (b) counts a spell as a run of constant label for one animal, so a spell ends when the
    label changes or the animal's record does; spells are right-censored at the ends of the
    record, which affects the dynamic curve far more than the assigned one and makes the
    {ratio:,.0f}-fold ratio a lower bound. Panel (c) uses hours whose state is fixed or carried
    overnight. Choosing between these identifiers is not cosmetic: <code>temp_group_id</code>
    cannot persist across hours, <code>association_cluster_id</code> is the join key for
    fine-scale 2-minute rows, and <code>association_event_id</code> is the only one that counts
    events &mdash; the full mapping is in
    <code>phase0_cohort_ledger/event_id_crosswalk.csv</code>. Frozen export, 1,924,104
    animal-hours, 2024-03-01 to 2026-07-22. Generated by
    <code>build_pipeline_figure.py</code>; seed {SEED}.
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
