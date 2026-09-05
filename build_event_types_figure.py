"""Manuscript figure: the event taxonomy, each type with a real worked example.

Every panel is an animal x hour raster built from the frozen hourly export: rows are
collared animals, columns are hours, and each cell is coloured by that animal's social
context in that hour. Unobserved animal-hours are left blank rather than filled, so the
reader can see the observation as well as the behaviour.

Three encodings, matching the three axes:

  BETWEEN-GROUP   with own unit only / in a cluster containing both units / alone
  WITHIN-GROUP    in the unit's largest cluster / in a secondary cluster / alone
  INDIVIDUAL      with origin group / alone / with another unit

Examples are chosen on the MEDIAN participation fraction over the event, not the max.
The saved `merge_scale` label is assigned from the maximum a unit ever reaches, so an
event that begins as a subset merge and grows into a large one carries the large label
while being mostly a subset -- picking on the median gives examples that actually
illustrate the category they are labelled with.

Output: docs/framing_2026-09-04/event_types_figure.html
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
OUT = Path("docs/framing_2026-09-04/event_types_figure.html")

MAX_ROWS = 22
MAX_HOURS = 30
WITH_ORIGIN = {"with_origin", "mixed_with_origin_present"}
AWAY_WITH_OTHERS = {"other", "mixed_without_origin_unit"}


# ------------------------------------------------------------------ picking
def pick_merges(s1: pd.DataFrame) -> list[dict]:
    s1 = s1.copy()
    s1["med_min"] = s1[["median_frac_a", "median_frac_b"]].min(axis=1)
    s1["med_max"] = s1[["median_frac_a", "median_frac_b"]].max(axis=1)
    s1["obs_min"] = s1[["mean_observed_a", "mean_observed_b"]].min(axis=1)
    base = s1[s1["structural_span_hours"].between(6, MAX_HOURS)
              & s1["observed_hour_fraction"].ge(0.75)
              & s1["obs_min"].ge(5)]
    out = []
    # a large merge: most of BOTH units, typically
    c = base[base["med_min"].ge(0.7)].sort_values("obs_min", ascending=False)
    if len(c):
        out.append(dict(row=c.iloc[0], kind="between", label="Large merge",
                        desc="most of both units in one cluster"))
    # a partial merge: most of one unit, part of the other
    c = base[base["med_max"].ge(0.6) & base["med_min"].between(0.2, 0.55)]
    c = c.sort_values("obs_min", ascending=False)
    if len(c):
        out.append(dict(row=c.iloc[0], kind="between", label="Partial merge",
                        desc="most of one unit, part of the other"))
    # a small subset merge: a few animals from each side
    c = base[base["med_max"].le(0.45)].sort_values("obs_min", ascending=False)
    if len(c):
        out.append(dict(row=c.iloc[0], kind="between", label="Small subset merge",
                        desc="a few animals from each side"))
    return out


def pick_split(d: pd.DataFrame) -> dict | None:
    """A run of hours where one well-covered unit holds two clusters of >= 2 animals."""
    obs = d[d["is_observed"]]
    best = None
    for (ts, unit), g in obs.groupby(["window_start", "dynamic_social_unit"], sort=False):
        if len(g) < 10:
            continue
        sizes = g.groupby("association_event_id").size().sort_values(ascending=False)
        if len(sizes) >= 2 and sizes.iloc[1] >= 3:
            score = (sizes.iloc[1], len(g))
            if best is None or score > best[0]:
                best = (score, ts, unit)
    if best is None:
        return None
    _, ts, unit = best
    return dict(kind="within", label="Within-group split",
                desc="one unit in two or more clusters",
                unit=unit, centre=ts)


def pick_individual(exc: pd.DataFrame) -> list[dict]:
    out = []
    a = exc[exc["depth"].eq("alone_only") & exc["away_nights"].between(3, 8)]
    if len(a):
        r = a.sort_values("away_nights").iloc[len(a) // 2]
        out.append(dict(kind="individual", label="Individual separation",
                        desc="one animal alone, then back",
                        animal=r["animal_id"], unit=r["origin_group"],
                        start=pd.Timestamp(r["start_night"]),
                        end=pd.Timestamp(r["end_night"])))
    b = exc[exc["depth"].eq("joined_only") & exc["away_nights"].between(4, 10)]
    if len(b):
        r = b.sort_values("away_nights").iloc[len(b) // 2]
        out.append(dict(kind="individual", label="Excursion to another unit",
                        desc="one animal leaves and joins elsewhere",
                        animal=r["animal_id"], unit=r["origin_group"],
                        start=pd.Timestamp(r["start_night"]),
                        end=pd.Timestamp(r["end_night"])))
    return out


# ------------------------------------------------------------------ grids
def between_grid(d: pd.DataFrame, unit_a: str, unit_b: str,
                 t0: pd.Timestamp, t1: pd.Timestamp) -> dict:
    w = d[d["window_start"].between(t0, t1)
          & d["dynamic_social_unit"].isin([unit_a, unit_b])].copy()
    hours = sorted(w["window_start"].unique())[:MAX_HOURS]
    w = w[w["window_start"].isin(hours)]
    # which clusters hold animals from both units, per hour
    mixed = set()
    for ts, g in w[w["is_observed"]].groupby("window_start"):
        for ev, gg in g.groupby("association_event_id"):
            if gg["dynamic_social_unit"].nunique() > 1:
                mixed.add((ts, ev))
    cells, rows = {}, []
    for unit in (unit_a, unit_b):
        ids = (w[w["dynamic_social_unit"].eq(unit) & w["is_observed"]]
               .groupby("animal_id").size().sort_values(ascending=False))
        for a in list(ids.index)[:MAX_ROWS // 2]:
            rows.append((unit, a))
    for r in w.itertuples():
        key = (r.dynamic_social_unit, r.animal_id)
        if key not in rows or not r.is_observed:
            continue
        if (r.window_start, r.association_event_id) in mixed:
            v = "mixed"
        elif r.temp_group_size <= 1:
            v = "alone"
        else:
            v = "own"
        cells[(key, r.window_start)] = v
    return {"hours": hours, "rows": rows, "cells": cells,
            "legend": [("own", "with own unit"), ("mixed", "cluster holds both units"),
                       ("alone", "alone")]}


def within_grid(d: pd.DataFrame, unit: str, centre: pd.Timestamp) -> dict:
    t0, t1 = centre - pd.Timedelta(hours=12), centre + pd.Timedelta(hours=17)
    w = d[d["window_start"].between(t0, t1) & d["dynamic_social_unit"].eq(unit)].copy()
    hours = sorted(w["window_start"].unique())[:MAX_HOURS]
    w = w[w["window_start"].isin(hours)]
    biggest = {}
    for ts, g in w[w["is_observed"]].groupby("window_start"):
        sizes = g.groupby("association_event_id").size()
        biggest[ts] = sizes.idxmax()
    ids = (w[w["is_observed"]].groupby("animal_id").size()
           .sort_values(ascending=False))
    rows = [(unit, a) for a in list(ids.index)[:MAX_ROWS]]
    cells = {}
    for r in w.itertuples():
        key = (unit, r.animal_id)
        if key not in rows or not r.is_observed:
            continue
        if r.temp_group_size <= 1:
            v = "alone"
        elif biggest.get(r.window_start) == r.association_event_id:
            v = "main"
        else:
            v = "secondary"
        cells[(key, r.window_start)] = v
    return {"hours": hours, "rows": rows, "cells": cells,
            "legend": [("main", "largest cluster"), ("secondary", "secondary cluster"),
                       ("alone", "alone")]}


def individual_grid(d: pd.DataFrame, animal: str, unit: str,
                    start: pd.Timestamp, end: pd.Timestamp) -> dict:
    t0 = start - pd.Timedelta(days=2)
    t1 = end + pd.Timedelta(days=3)
    w = d[d["window_start"].between(t0, t1)
          & d["dynamic_social_unit"].eq(unit)].copy()
    # one column per night, to keep the panel readable over days
    w["night"] = (w["window_start"] - pd.Timedelta(hours=10)).dt.normalize()
    nights = sorted(w["night"].unique())[:MAX_HOURS]
    w = w[w["night"].isin(nights)]

    def state(g):
        g = g[g["is_observed"]]
        if g.empty:
            return None
        sc = g["social_context"]
        if (sc.isin(AWAY_WITH_OTHERS)).any():
            return "other"
        if (sc.eq("isolated")).any():
            return "alone"
        if (sc.isin(WITH_ORIGIN)).any():
            return "origin"
        return None

    ids = w[w["is_observed"]].groupby("animal_id").size().sort_values(ascending=False)
    others = [a for a in ids.index if a != animal][:MAX_ROWS - 1]
    rows = [(unit, animal)] + [(unit, a) for a in others]
    cells = {}
    for (a, n), g in w.groupby(["animal_id", "night"]):
        key = (unit, a)
        if key not in rows:
            continue
        v = state(g)
        if v:
            cells[(key, n)] = v
    return {"hours": nights, "rows": rows, "cells": cells, "focal": (unit, animal),
            "legend": [("origin", "with origin group"), ("alone", "alone"),
                       ("other", "with another unit")]}


# ------------------------------------------------------------------ drawing
FILL = {
    "own": "var(--n2)", "mixed": "var(--c1)", "alone": "var(--c6)",
    "main": "var(--n2)", "secondary": "var(--c4)",
    "origin": "var(--n2)", "other": "var(--c5)",
}


def panel(g: dict, spec: dict, x0: float, y0: float, w: float) -> tuple[str, float]:
    hours, rows, cells = g["hours"], g["rows"], g["cells"]
    if not hours or not rows:
        return "", 0.0
    lab_w = 58.0
    cw = min(7.0, (w - lab_w) / max(1, len(hours)))
    ch = 6.4
    o = []
    o.append('<text class="ph" x="%.1f" y="%.1f">%s</text>' % (x0, y0 - 15, spec["label"]))
    o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (x0, y0 - 5, spec["desc"]))
    prev_unit = None
    for i, key in enumerate(rows):
        unit, animal = key
        yy = y0 + i * ch
        if prev_unit is not None and unit != prev_unit:
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                     'stroke="var(--ink)" stroke-width=".8"/>'
                     % (x0, yy - 1.2, x0 + lab_w + cw * len(hours), yy - 1.2))
        prev_unit = unit
        focal = g.get("focal") == key
        o.append('<text class="%s" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % ("rlf" if focal else "rl", x0 + lab_w - 4, yy + ch - 1.4,
                    animal.split("_")[0]))
        for j, h in enumerate(hours):
            v = cells.get((key, h))
            xx = x0 + lab_w + j * cw
            if v is None:
                o.append('<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" '
                         'fill="none" stroke="var(--rule-soft)" stroke-width=".4"/>'
                         % (xx, yy, cw - .6, ch - .9))
            else:
                o.append('<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" '
                         'fill="%s"/>' % (xx, yy, cw - .6, ch - .9, FILL[v]))
    ybot = y0 + len(rows) * ch
    # unit labels down the left
    seen = []
    for i, (unit, _) in enumerate(rows):
        if unit not in [s[0] for s in seen]:
            seen.append((unit, i))
    for unit, i in seen:
        n = sum(1 for u, _ in rows if u == unit)
        o.append('<text class="ul" transform="translate(%.1f %.1f) rotate(-90)" '
                 'text-anchor="middle">%s</text>'
                 % (x0 - 5, y0 + (i + n / 2) * ch, unit))
    # legend, wrapping at the panel edge so nothing spills into the next column
    items = list(g["legend"]) + [(None, "not observed")]
    right = x0 + lab_w + cw * len(hours)
    lx, ly = x0 + lab_w, ybot + 6
    lines = 1
    for v, lab in items:
        item_w = 10 + len(lab) * 4.5 + 14
        if lx + item_w > right and lx > x0 + lab_w:
            lx, ly = x0 + lab_w, ly + 11
            lines += 1
        if v is None:
            o.append('<rect x="%.1f" y="%.1f" width="7" height="7" fill="none" '
                     'stroke="var(--rule-soft)" stroke-width=".6"/>' % (lx, ly))
        else:
            o.append('<rect x="%.1f" y="%.1f" width="7" height="7" fill="%s"/>'
                     % (lx, ly, FILL[v]))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>'
                 % (lx + 10, ly + 6, lab))
        lx += item_w
    unit_lab = "nights" if spec["kind"] == "individual" else "hours"
    o.append('<text class="ts" x="%.1f" y="%.1f">%d %s</text>'
             % (x0 + lab_w, ly + 18, len(hours), unit_lab))
    return "\n".join(o), (ly + 24) - y0


def build():
    cols = ["window_start", "animal_id", "origin_group", "dynamic_social_unit",
            "social_context", "association_event_id", "temp_group_size", "is_observed"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
    s1 = pd.read_csv(STAGE1, parse_dates=["start_hour", "end_hour"])
    exc = pd.read_csv(EXC, parse_dates=["start_night", "end_night"])

    specs, grids = [], []
    for m in pick_merges(s1):
        r = m["row"]
        g = between_grid(d, r["unit_a"], r["unit_b"],
                         pd.Timestamp(r["start_hour"]),
                         pd.Timestamp(r["start_hour"]) + pd.Timedelta(hours=MAX_HOURS))
        m["meta"] = ("%s, %s; median participation %.2f / %.2f"
                     % (r["pair_key"], pd.Timestamp(r["start_hour"]).date(),
                        r["median_frac_a"], r["median_frac_b"]))
        specs.append(m)
        grids.append(g)

    sp = pick_split(d[d["window_start"].dt.year.ge(2025)])
    if sp:
        g = within_grid(d, sp["unit"], sp["centre"])
        sp["meta"] = "%s, %s" % (sp["unit"], pd.Timestamp(sp["centre"]).date())
        specs.append(sp)
        grids.append(g)

    for iv in pick_individual(exc):
        g = individual_grid(d, iv["animal"], iv["unit"], iv["start"], iv["end"])
        iv["meta"] = ("%s of %s, %s to %s"
                      % (iv["animal"].split("_")[0], iv["unit"],
                         iv["start"].date(), iv["end"].date()))
        specs.append(iv)
        grids.append(g)

    # two columns
    colw, gap = 336.0, 28.0
    y = [66.0, 66.0]
    out = []
    for i, (spec, g) in enumerate(zip(specs, grids)):
        c = i % 2
        x0 = 34.0 + c * (colw + gap)
        svg, h = panel(g, spec, x0, y[c], colw)
        out.append('<text class="pl" x="%.1f" y="%.1f">%s</text>'
                   % (x0 - 26, y[c] - 15, "abcdefgh"[i]))
        out.append(svg)
        out.append('<text class="ts" x="%.1f" y="%.1f">%s</text>'
                   % (x0, y[c] + h - 2, spec["meta"]))
        y[c] += h + 40
    H = max(y) + 8
    head = ('<svg viewBox="0 0 734 %.0f" role="img" aria-label="Six worked examples of '
            'social event types, each an animal-by-hour raster from the tracking record: '
            'large merge, partial merge, small subset merge, within-group split, '
            'individual separation and an excursion to another unit.">' % H)
    return head + "\n" + "\n".join(out) + "\n</svg>", specs


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
.wrap{max-width:830px;margin:0 auto;padding:48px 24px 80px}
.kicker{font-family:var(--mono);font-size:10px;letter-spacing:.16em;
text-transform:uppercase;color:var(--muted);margin-bottom:10px}
h1{font-family:var(--serif);font-size:30px;font-weight:500;line-height:1.15;
letter-spacing:-.01em;margin:0 0 26px;text-wrap:balance}
.plate{background:var(--paper);border:1px solid var(--rule-soft);padding:20px 18px 8px}
.plate svg{display:block;width:100%;height:auto;overflow:visible}
.scroll{overflow-x:auto}.scroll>svg{min-width:700px}
.cap{font-size:12.5px;line-height:1.62;color:var(--ink-2);margin:16px 0 0;
padding-top:14px;border-top:1px solid var(--rule)}
.cap b{font-weight:600;color:var(--ink)}.cap .fignum{font-weight:600;color:var(--ink)}
.note{font-size:11.5px;line-height:1.6;color:var(--muted);margin-top:22px;
padding-top:14px;border-top:1px solid var(--rule-soft)}
.note code{font-family:var(--mono);font-size:.92em}
text{font-family:var(--sans)}
.rl{font-size:6.4px;fill:var(--muted);font-family:var(--mono)}
.rlf{font-size:6.4px;fill:var(--ink);font-weight:700;font-family:var(--mono)}
.ul{font-size:7.4px;fill:var(--ink-2);font-weight:600}
.ts{font-size:8.2px;fill:var(--muted)}
.ph{font-size:11px;fill:var(--ink);font-weight:600}
.pl{font-size:13px;font-weight:700;fill:var(--ink)}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def page():
    fig, specs = build()
    labels = ", ".join("(%s) %s" % ("abcdefgh"[i], s["label"].lower())
                       for i, s in enumerate(specs))
    return f"""<title>Event Types With Examples</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;700&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 1 &middot; draft</div>
  <h1>The six social events, as they appear in the record</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 1.</span> Worked examples of each social event type, drawn
      directly from the tracking record. In every panel rows are collared animals, columns are
      hours (nights in the two individual panels), and cell colour is the animal's social context
      in that interval; <b>unobserved animal-intervals are left blank</b>, so observation and
      behaviour are separable by eye. Panels: {labels}.
      <b>(a&ndash;c)</b> Between-group events, coloured by whether the animal's cluster contains
      only its own unit or both units. The horizontal rule separates the two units. The three
      panels differ only in how much of each unit participates &mdash; the same underlying
      structure at three depths.
      <b>(d)</b> A within-group split: one unit occupying two or more clusters at once, coloured
      by whether the animal sits in the unit's largest cluster or a secondary one.
      <b>(e&ndash;f)</b> Individual events at nightly resolution, with the focal animal in bold
      at the top and its unit-mates beneath for reference: one animal alone and then back, and
      one animal leaving to join another unit.
    </p>
  </div>

  <p class="note">
    Examples are selected on the <b>median</b> participation fraction over the event, not the
    maximum. The saved <code>merge_scale</code> label is assigned from the highest fraction a
    unit ever reaches, so an event that begins as a subset merge and grows carries the large
    label while being mostly a subset &mdash; several events labelled small-subset have maximum
    participation fractions near 1.0. Selecting on the median gives examples that illustrate
    the category they are named for, and the discrepancy is worth resolving in the event
    classifier. All panels are from the frozen 1,924,104-row export (2024-03-01 to 2026-07-22).
    Generated by <code>build_event_types_figure.py</code>.
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
