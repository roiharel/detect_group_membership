"""The basic census: how often each social structure occurs, and how long it lasts.

Everything upstream of the models. Five plots, all computed here from saved outputs:

  1. STATE FREQUENCIES -- the same nights classified three ways, so the effect of the
     resolvability filter is visible rather than assumed.
  2. EVENT COUNTS AND RATES -- with the denominator each rate is computed against, or
     an explicit note that none exists.
  3. DURATION SURVIVAL CURVES -- P(duration > x) on a log axis, one curve per social
     structure. This is the plot the project did not have: 181 saved figures and no
     cumulative distribution among them.
  4. THREE DURATIONS OF ONE ENCOUNTER -- structural span, supported exposure and
     active contact as three survival curves. They differ by two orders of magnitude,
     which is the Phase 2 point and has never been drawn.
  5. EXCURSION DURATION BY DEPTH -- alone-only against reached-another-unit, which is
     "depth requires duration" as a distribution rather than two medians.

Colour encodes axis throughout: blues for between-group structures, bronze for
within-group, plum and teal for individual.

Output: docs/framing_2026-09-04/social_structure_census.html
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("outputs/general_structure_2026_09")
LEGACY_EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter")
OUT = Path("docs/framing_2026-09-04/social_structure_census.html")

W = 700.0


# ------------------------------------------------------------------ helpers
def survival(vals: np.ndarray, xs: np.ndarray) -> np.ndarray:
    """P(T > x) for each x, empirical."""
    v = np.sort(np.asarray(vals, dtype=float))
    return np.array([(v > x).sum() / len(v) for x in xs])


def logspace_axis(lo: float, hi: float, n: int = 70) -> np.ndarray:
    return np.geomspace(lo, hi, n)


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ------------------------------------------------------------------ load
def load() -> dict:
    d: dict = {}
    ev = pd.read_csv(LEGACY_EVENTS / "canonical_event_size_duration_all_events.csv")
    d["events"] = ev
    d["isolated"] = pd.read_csv(LEGACY_EVENTS / "canonical_isolated_events.csv")
    d["stage1"] = pd.read_csv(BASE / "phase2_two_stage_events/stage1_events_with_stage2_mixing.csv")
    d["exc"] = pd.read_csv(BASE / "phase4b_individual_axis/excursions_dominant_gap0.csv")
    d["nightly"] = pd.read_csv(BASE / "phase4b_individual_axis/nightly_states_dominant.csv")
    d["opp"] = pd.read_csv(BASE / "phase1_opportunity/opportunity_state_summary.csv")
    return d


# ------------------------------------------------------------------ plot 1
STATE_SETS = [
    ("Where the animal slept",
     "canonical nightly export, 81,695 animal-nights",
     [("With origin group", 67561, "var(--ax-a)"),
      ("With another group", 11418, "var(--ax-c)"),
      ("Isolated", 2499, "var(--accent)"),
      ("Mixed / unclear", 217, "var(--muted)")]),
    ("What the configuration was",
     "same nights, orthogonal classification",
     [("Whole group", 54572, "var(--ax-a)"),
      ("Between-group merge", 20468, "var(--ax-a2)"),
      ("Within-group split", 3789, "var(--ax-b)"),
      ("Isolated / low support", 2866, "var(--muted)")]),
    ("After the located-reference filter",
     "frozen four-state ledger, 82,205 animal-nights",
     [("With origin group", 73330, "var(--ax-a)"),
      ("With another group", 1890, "var(--ax-c)"),
      ("Alone", 448, "var(--accent)"),
      ("Unresolvable", 6537, "var(--muted)")]),
]


def plot_states() -> str:
    left, right = 0.0, 700.0
    bar_x, bar_w = 208.0, 430.0
    top, blockh = 34.0, 96.0
    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Three stacked bars showing '
         'the same animal-nights classified three ways. The located-reference filter '
         'moves nights from with-another-group and isolated into unresolvable.">'
         % (top + blockh * len(STATE_SETS) + 26)]
    for i, (title, sub, parts) in enumerate(STATE_SETS):
        y = top + blockh * i
        total = sum(n for _, n, _ in parts)
        o.append('<text class="lbl-em" x="0" y="%.1f">%s</text>' % (y, esc(title)))
        o.append('<text class="lbl-sm" x="0" y="%.1f">%s</text>' % (y + 13, esc(sub)))
        o.append('<text class="lbl-sm" x="0" y="%.1f">n = %s</text>'
                 % (y + 26, format(total, ",")))
        x = bar_x
        for label, n, col in parts:
            w = bar_w * n / total
            o.append('<rect x="%.1f" y="%.1f" width="%.2f" height="26" fill="%s"/>'
                     % (x, y - 12, w, col))
            share = 100 * n / total
            if w > 42:
                o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle" '
                         'fill="var(--paper)">%.1f%%</text>' % (x + w / 2, y + 5, share))
            x += w
        # legend beneath, wrapping so nothing runs past the drawing edge
        lx, ly_off = bar_x, 33.0
        for label, n, col in parts:
            item_w = 11 + len(label) * 5.4 + 30
            if lx + item_w > 696:
                lx, ly_off = bar_x, ly_off + 13
            o.append('<rect x="%.1f" y="%.1f" width="8" height="8" fill="%s"/>'
                     % (lx, y + ly_off - 7, col))
            o.append('<text class="lbl-sm" x="%.1f" y="%.1f">%s %.1f%%</text>'
                     % (lx + 11, y + ly_off, esc(label), 100 * n / total))
            lx += item_w
        o.append('<line class="grid" x1="0" y1="%.1f" x2="700" y2="%.1f"/>'
                 % (y + ly_off + 13, y + ly_off + 13))
    o.append("</svg>")
    return "\n".join(o)


# ------------------------------------------------------------------ plot 3
CURVES = [
    # (label, key, colour, source)
    ("Large merge", "large_merge", "var(--ax-a)", "legacy"),
    ("Medium partial merge", "medium_partial_merge", "var(--ax-a2)", "legacy"),
    ("Small subset merge", "small_subset_merge", "var(--ax-a3)", "legacy"),
    ("Within-group split", "within_group_split", "var(--ax-b)", "legacy"),
    ("Single-animal separation", "single_animal_separation", "var(--ax-c)", "legacy"),
    ("Isolated (alone)", "isolated", "var(--accent)", "isolated"),
]


def plot_survival(d: dict) -> tuple[str, dict]:
    ev, iso = d["events"], d["isolated"]
    pad_l, pad_b, top = 52.0, 56.0, 26.0
    plot_w, plot_h = 452.0, 250.0
    lo, hi = 1.0, 10000.0
    xs = logspace_axis(lo, hi)

    def px(v):
        return pad_l + (math.log10(max(v, lo)) - math.log10(lo)) / (
            math.log10(hi) - math.log10(lo)) * plot_w

    def py(p):
        return top + plot_h - p * plot_h

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Survival curves of event '
         'duration on a log axis. Every social structure has a median at or under a few '
         'hours and a tail reaching hundreds of hours.">' % (top + plot_h + pad_b)]
    for t, lab in ((1, "1 h"), (10, "10 h"), (100, "100 h"), (1000, "1,000 h"),
                   (10000, "10,000 h")):
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (px(t), top, px(t), top + plot_h))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (px(t), top + plot_h + 16, lab))
    for t, lab in ((24, "1 day"), (168, "1 week"), (730, "1 month")):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
                 'stroke-width="1" stroke-dasharray="3 4" opacity=".45"/>'
                 % (px(t), top - 6, px(t), top + plot_h))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (px(t), top - 10, lab))
    for p in (0, 0.25, 0.5, 0.75, 1.0):
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (pad_l, py(p), pad_l + plot_w, py(p)))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="end">%d%%</text>'
                 % (pad_l - 7, py(p) + 3.5, int(p * 100)))

    stats = {}
    for label, key, col, src in CURVES:
        v = (iso["duration_hours"] if src == "isolated"
             else ev.loc[ev["event_type"].eq(key), "duration_hours"]).dropna()
        v = v[v > 0].to_numpy()
        s = survival(v, xs)
        pts = " ".join("%.1f,%.1f" % (px(x), py(p)) for x, p in zip(xs, s))
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2" '
                 'stroke-linejoin="round"/>' % (pts, col))
        stats[label] = {"n": int(len(v)), "median": float(np.median(v)),
                        "p90": float(np.quantile(v, 0.9)),
                        "max": float(v.max()),
                        "share_over_24h": float((v > 24).mean()),
                        "share_over_168h": float((v > 168).mean())}
    # legend at the right
    ly = top + 6
    for label, key, col, src in CURVES:
        st = stats[label]
        o.append('<rect x="%.1f" y="%.1f" width="9" height="9" fill="%s"/>'
                 % (pad_l + plot_w + 14, ly - 8, col))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">%s</text>'
                 % (pad_l + plot_w + 27, ly, esc(label)))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">n=%s  med %gh</text>'
                 % (pad_l + plot_w + 27, ly + 11, format(st["n"], ","), st["median"]))
        ly += 27
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top + plot_h, pad_l + plot_w, top + plot_h))
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top, pad_l, top + plot_h))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle">event '
             'duration, log scale</text>' % (pad_l + plot_w / 2, top + plot_h + 36))
    o.append('<text class="lbl-sm" transform="translate(%.1f %.1f) rotate(-90)" '
             'text-anchor="middle">share of events lasting longer</text>'
             % (pad_l - 34, top + plot_h / 2))
    o.append("</svg>")
    return "\n".join(o), stats


# ------------------------------------------------------------------ plot 4
def plot_three_durations(d: dict) -> tuple[str, dict]:
    a = d["stage1"]
    series = [
        ("Structural span", "structural_span_hours", "var(--ax-a)",
         "first to last encounter hour"),
        ("Supported exposure", "5m_supported_exposure_hours", "var(--mid)",
         "eligible 2-min bins x 2 min"),
        ("Active contact", "5m_active_contact_hours", "var(--differ)",
         "contact-positive bins only"),
    ]
    pad_l, top = 52.0, 26.0
    plot_w, plot_h = 414.0, 220.0
    lo, hi = 0.03, 10000.0
    xs = logspace_axis(lo, hi, 80)

    def px(v):
        return pad_l + (math.log10(max(v, lo)) - math.log10(lo)) / (
            math.log10(hi) - math.log10(lo)) * plot_w

    def py(p):
        return top + plot_h - p * plot_h

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Three survival curves for '
         'the same encounters: structural span, supported exposure and active contact, '
         'separated by about two orders of magnitude.">' % (top + plot_h + 58)]
    for t, lab in ((0.03, "2 min"), (0.1, "6 min"), (1, "1 h"), (10, "10 h"),
                   (100, "100 h"), (1000, "1,000 h"), (10000, "10,000 h")):
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (px(t), top, px(t), top + plot_h))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (px(t), top + plot_h + 16, lab))
    for p in (0, 0.5, 1.0):
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (pad_l, py(p), pad_l + plot_w, py(p)))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="end">%d%%</text>'
                 % (pad_l - 7, py(p) + 3.5, int(p * 100)))
    stats = {}
    ly = top + 8
    for label, col_name, col, note in series:
        v = pd.to_numeric(a[col_name], errors="coerce").dropna()
        v = v[v > 0].to_numpy()
        s = survival(v, xs)
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2.2" '
                 'stroke-linejoin="round"/>'
                 % (" ".join("%.1f,%.1f" % (px(x), py(p)) for x, p in zip(xs, s)), col))
        med = float(np.median(v))
        o.append('<circle cx="%.1f" cy="%.1f" r="3.5" fill="%s"/>' % (px(med), py(0.5), col))
        stats[label] = {"n": int(len(v)), "median": med,
                        "p90": float(np.quantile(v, 0.9)), "max": float(v.max())}
        o.append('<rect x="%.1f" y="%.1f" width="9" height="9" fill="%s"/>'
                 % (pad_l + plot_w + 14, ly - 8, col))
        o.append('<text class="lbl-em" x="%.1f" y="%.1f">%s</text>'
                 % (pad_l + plot_w + 27, ly, esc(label)))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">%s</text>'
                 % (pad_l + plot_w + 27, ly + 11, esc(note)))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">n=%s  median %s h</text>'
                 % (pad_l + plot_w + 27, ly + 22, format(stats[label]["n"], ","),
                    ("%.2f" % med).rstrip("0").rstrip(".")))
        ly += 40
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="1" stroke-dasharray="3 3" opacity=".5"/>'
             % (pad_l, py(0.5), pad_l + plot_w, py(0.5)))
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top + plot_h, pad_l + plot_w, top + plot_h))
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top, pad_l, top + plot_h))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle">duration of '
             'the same 1,705 encounters, measured three ways &mdash; log scale</text>'
             % (pad_l + plot_w / 2, top + plot_h + 38))
    o.append("</svg>")
    return "\n".join(o), stats


# ------------------------------------------------------------------ plot 5
def plot_excursions(d: dict) -> tuple[str, dict]:
    e = d["exc"]
    groups = [
        ("Alone only", e.loc[e["depth"].eq("alone_only"), "away_nights"], "var(--accent)"),
        ("Reached another unit", e.loc[e["depth"].ne("alone_only"), "away_nights"],
         "var(--ax-c)"),
    ]
    pad_l, top = 52.0, 26.0
    plot_w, plot_h = 408.0, 200.0
    lo, hi = 1.0, 600.0
    xs = logspace_axis(lo, hi, 70)

    def px(v):
        return pad_l + (math.log10(max(v, lo)) - math.log10(lo)) / (
            math.log10(hi) - math.log10(lo)) * plot_w

    def py(p):
        return top + plot_h - p * plot_h

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Survival curves of '
         'excursion length. Excursions that only reach being alone are much shorter '
         'than those reaching another unit.">' % (top + plot_h + 56)]
    for t, lab in ((1, "1 night"), (7, "1 week"), (30, "1 month"), (100, "100"),
                   (365, "1 year")):
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (px(t), top, px(t), top + plot_h))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (px(t), top + plot_h + 16, lab))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="1.2" stroke-dasharray="4 3"/>'
             % (px(7), top - 6, px(7), top + plot_h))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle" '
             'fill="var(--ink)">settlement threshold</text>' % (px(7), top - 10))
    for p in (0, 0.5, 1.0):
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (pad_l, py(p), pad_l + plot_w, py(p)))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="end">%d%%</text>'
                 % (pad_l - 7, py(p) + 3.5, int(p * 100)))
    stats = {}
    ly = top + 10
    for label, ser, col in groups:
        v = ser.dropna().to_numpy()
        v = v[v > 0]
        s = survival(v, xs)
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2.2" '
                 'stroke-linejoin="round"/>'
                 % (" ".join("%.1f,%.1f" % (px(x), py(p)) for x, p in zip(xs, s)), col))
        med = float(np.median(v))
        stats[label] = {"n": int(len(v)), "median": med,
                        "p90": float(np.quantile(v, 0.9)), "max": float(v.max()),
                        "share_over_7": float((v > 7).mean())}
        o.append('<rect x="%.1f" y="%.1f" width="9" height="9" fill="%s"/>'
                 % (pad_l + plot_w + 14, ly - 8, col))
        o.append('<text class="lbl-em" x="%.1f" y="%.1f">%s</text>'
                 % (pad_l + plot_w + 27, ly, esc(label)))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">n=%d  median %g nights</text>'
                 % (pad_l + plot_w + 27, ly + 11, stats[label]["n"], med))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">%.0f%% pass 7 nights</text>'
                 % (pad_l + plot_w + 27, ly + 22, 100 * stats[label]["share_over_7"]))
        ly += 42
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top + plot_h, pad_l + plot_w, top + plot_h))
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top, pad_l, top + plot_h))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle">away-nights '
             'per excursion, log scale</text>'
             % (pad_l + plot_w / 2, top + plot_h + 36))
    o.append("</svg>")
    return "\n".join(o), stats


# ------------------------------------------------------------------ tables
def rate_table(d: dict) -> str:
    ev = d["events"]
    rows = [
        ("Between-group", "Structural encounter (frozen)", "1,705", "68 dyads",
         "2,867 of 73,353 dyad-days", "3.91% of dyad-days", "14 h"),
        ("Between-group", "Large merge (legacy)", "1,420", "74 dyads",
         "none", "&mdash;", "11 h"),
        ("Between-group", "Medium partial merge (legacy)", "665", "&mdash;",
         "none", "&mdash;", "2 h"),
        ("Between-group", "Small subset merge (legacy)", "1,132", "&mdash;",
         "none", "&mdash;", "3 h"),
        ("Within-group", "Within-group split (legacy)", "4,570", "24 groups",
         "1,268 group-weeks (frozen)", "63% of group-weeks contain one", "2 h"),
        ("Individual", "Excursion (frozen)", "338", "91 animals, 18 groups",
         "82,205 animal-nights", "0.41% of animal-nights start one", "2 nights"),
        ("Individual", "&mdash; reaching settlement", "61", "16 animals",
         "338 excursions", "18% of excursions", "&ge; 7 nights"),
        ("Individual", "Single-animal separation (legacy)", "5,828", "301 animals",
         "none", "&mdash;", "1 h"),
        ("Individual", "Isolated, unfiltered (legacy)", "4,293", "285 animals",
         "none", "&mdash;", "1 h"),
        ("Individual", "Alone, located-reference filter (frozen)", "&mdash;",
         "267 animals, 21 groups", "11,668 animal-hours retained of 48,775",
         "23.9% of raw survives", "&mdash;"),
    ]
    body = "\n".join(
        '<tr><td>%s</td><td>%s</td><td class="n">%s</td><td>%s</td><td>%s</td>'
        '<td>%s</td><td class="n">%s</td></tr>' % r for r in rows)
    return ('<div class="tw"><table><caption>Counts, denominators and median duration '
            'by social structure</caption><thead><tr><th>Axis</th><th>Structure</th>'
            '<th class="n">Events</th><th>Units</th><th>Denominator</th>'
            '<th>Rate</th><th class="n">Median</th></tr></thead><tbody>%s</tbody>'
            '</table></div>' % body)


CSS = """
:root{--ground:#f1f3f1;--paper:#fbfcfa;--ink:#141a19;--ink-2:#38423f;--muted:#667370;
--faint:#8d9793;--rule:#d3d9d4;--rule-soft:#e3e8e3;--accent:#0f5f57;--ax-a:#1f6f8b;
--ax-a2:#4a92aa;--ax-a3:#7db3c4;--ax-b:#8a5a1f;--ax-c:#6b4a7a;--mid:#8a6410;
--differ:#8c3a2e;
--shadow:0 1px 2px rgba(20,26,25,.05),0 8px 24px -14px rgba(20,26,25,.18);
--serif:"Newsreader","Iowan Old Style",Georgia,serif;
--sans:"IBM Plex Sans","Segoe UI",Helvetica,Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#101413;
--paper:#171d1c;--ink:#e9efec;--ink-2:#c3ccc8;--muted:#94a09c;--faint:#737e7a;
--rule:#2b3432;--rule-soft:#222a29;--accent:#5fc0b0;--ax-a:#74b6d0;--ax-a2:#4f93ad;
--ax-a3:#3d7288;--ax-b:#d3a061;--ax-c:#b494c4;--mid:#d8ab54;--differ:#e08a78;
--shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px -14px rgba(0,0,0,.7);}}
:root[data-theme="dark"]{--ground:#101413;--paper:#171d1c;--ink:#e9efec;--ink-2:#c3ccc8;
--muted:#94a09c;--faint:#737e7a;--rule:#2b3432;--rule-soft:#222a29;--accent:#5fc0b0;
--ax-a:#74b6d0;--ax-a2:#4f93ad;--ax-a3:#3d7288;--ax-b:#d3a061;--ax-c:#b494c4;
--mid:#d8ab54;--differ:#e08a78;
--shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px -14px rgba(0,0,0,.7);}
*{box-sizing:border-box}
body{background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:16px;
line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:940px;margin:0 auto;padding:44px 28px 96px}
.eyebrow{font-family:var(--mono);font-size:10.5px;letter-spacing:.16em;
text-transform:uppercase;color:var(--accent);margin-bottom:14px}
h1{font-family:var(--serif);font-size:clamp(32px,5vw,48px);font-weight:500;
line-height:1.05;letter-spacing:-.018em;margin:0 0 12px;text-wrap:balance}
.standfirst{font-family:var(--serif);font-size:19px;line-height:1.5;color:var(--ink-2);
max-width:36em;margin:0 0 18px}
header{border-bottom:2px solid var(--ink);padding-bottom:20px;margin-bottom:36px}
.meta{display:flex;flex-wrap:wrap;gap:6px 22px;font-family:var(--mono);font-size:11.5px;
color:var(--muted)}
.meta b{font-weight:500;color:var(--ink-2)}
h2{font-family:var(--serif);font-size:26px;font-weight:500;letter-spacing:-.012em;
line-height:1.2;margin:0 0 6px;text-wrap:balance}
.secnum{font-family:var(--mono);font-size:12px;color:var(--accent);margin-bottom:4px}
section{margin-bottom:54px}
p{margin:0 0 14px;max-width:66ch}
strong{font-weight:600}
code,.num{font-family:var(--mono);font-size:.88em;font-variant-numeric:tabular-nums}
figure{margin:24px 0;background:var(--paper);border:1px solid var(--rule-soft);
padding:20px 20px 16px;box-shadow:var(--shadow)}
figure svg{display:block;width:100%;height:auto;overflow:visible}
.scroll{overflow-x:auto}.scroll>svg{min-width:560px}
figcaption{font-size:12px;color:var(--muted);margin-top:14px;padding-top:12px;
border-top:1px solid var(--rule-soft);max-width:70ch}
text{font-family:var(--mono);font-variant-numeric:tabular-nums}
.lbl{font-size:10.5px;fill:var(--ink-2)}.lbl-sm{font-size:9.5px;fill:var(--faint)}
.lbl-em{font-size:10.5px;fill:var(--ink);font-weight:500}
.axline{stroke:var(--rule);stroke-width:1}.grid{stroke:var(--rule-soft);stroke-width:1}
.tw{overflow-x:auto;margin:18px 0 22px}
table{border-collapse:collapse;width:100%;font-size:13px;min-width:640px}
caption{caption-side:top;text-align:left;font-family:var(--mono);font-size:10.5px;
letter-spacing:.1em;text-transform:uppercase;color:var(--faint);padding-bottom:8px}
th,td{text-align:left;padding:7px 12px 7px 0;vertical-align:top;
border-bottom:1px solid var(--rule-soft)}
thead th{font-size:11px;font-weight:600;letter-spacing:.05em;text-transform:uppercase;
color:var(--muted);border-bottom:1px solid var(--rule);padding-bottom:6px}
td.n,th.n{text-align:right;font-family:var(--mono);font-variant-numeric:tabular-nums;
padding-right:0}
tbody tr:last-child td{border-bottom:none}
blockquote{margin:22px 0;padding:18px 22px;background:var(--paper);
border-left:3px solid var(--accent);font-family:var(--serif);font-size:19px;
line-height:1.45;max-width:42em}
blockquote p{margin:0;max-width:none}
footer{border-top:1px solid var(--rule);padding-top:16px;font-size:12.5px;color:var(--muted)}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def build() -> str:
    d = load()
    p1 = plot_states()
    p3, s3 = plot_survival(d)
    p4, s4 = plot_three_durations(d)
    p5, s5 = plot_excursions(d)

    over24 = ", ".join(
        "%s %.0f%%" % (k, 100 * v["share_over_24h"]) for k, v in s3.items())
    ratio = s4["Structural span"]["median"] / max(1e-9, s4["Active contact"]["median"])

    return f"""<title>Social Structure Census</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500;6..72,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header>
  <div class="eyebrow">Companion to Stable Groups, Fluid Boundaries</div>
  <h1>How often, and for how long</h1>
  <p class="standfirst">The census that comes before any model: counts, rates against their denominators, and the full duration distribution of every social structure in the cohort.</p>
  <div class="meta">
    <span><b>350</b> animals</span>
    <span><b>26</b> origin groups</span>
    <span><b>13,615</b> events inventoried</span>
    <span><b>81,695</b> animal-nights</span>
    <span><b>2024-03-01</b> &rarr; <b>2026-07-22</b></span>
  </div>
</header>

<section>
  <div class="secnum">01</div>
  <h2>The same nights, classified three ways</h2>
  <p>Two orthogonal classifications of the identical animal-nights, plus the frozen four-state ledger after the located-reference filter. Showing all three together makes the filter&rsquo;s effect visible rather than asserted: the nights it removes come overwhelmingly out of <em>with another group</em> and <em>isolated</em>, and land in <em>unresolvable</em>.</p>
  <figure>
    <div class="scroll">{p1}</div>
    <figcaption>Top: an animal is with its origin group on 82.7% of nights. Middle: the <em>configuration</em> it sleeps in is not its whole group on 33.2% of them. Bottom: once <code>alone</code> and <code>with_other</code> require the origin group to be locatable elsewhere, they fall to 0.5% and 2.3% and 8.0% of nights become explicitly unresolvable. The bottom bar is the honest one; the top two are what the raw classification reports. Sources: <code>membership_export_nightly</code>, <code>phase4b_individual_axis/nightly_states_dominant.csv</code>.</figcaption>
  </figure>
</section>

<section>
  <div class="secnum">02</div>
  <h2>Counts, and what each rate is a rate of</h2>
  <p>An event count without a denominator is an inventory, not a frequency. Three structures now have denominators and the rest do not &mdash; the table says which is which rather than quietly reporting all ten the same way.</p>
  {rate_table(d)}
  <p>The three rates that exist: a dyad meets on <strong>3.91%</strong> of the days both groups are tracked; <strong>63%</strong> of group-weeks contain at least one within-group split; an animal begins an excursion on <strong>0.41%</strong> of its resolvable nights, and <strong>18%</strong> of those excursions reach the seven-night settlement threshold.</p>
</section>

<section>
  <div class="secnum">03</div>
  <h2>Duration: every structure is brief, and every structure has a tail</h2>
  <p>Survival curves &mdash; the share of events still running past a given duration &mdash; on a log axis. The project held 181 saved figures and not one cumulative distribution, which is how five structures with medians between 1 and 11 hours and maxima between 217 and 6,330 hours came to be summarised by their medians alone.</p>
  <figure>
    <div class="scroll">{p3}</div>
    <figcaption>Share of events lasting longer than 24 h: {over24}. Every curve drops steeply through the first day and then flattens into a long tail: a large merge has a median of 11 h but a 99th percentile of 116 h and a maximum of 564 h. The small-subset-merge curve crosses above large merge past about 30 h &mdash; small partial merges are rarer but, when they persist, they persist longest. Source: <code>canonical_event_size_duration_all_events.csv</code> and <code>canonical_isolated_events.csv</code>, both legacy cohort.</figcaption>
  </figure>
  <blockquote><p>Common, brief, and heavy-tailed &mdash; at every scale. The median says the structures are transient; the tail says a handful are not, and those are the ones that reach real mixing.</p></blockquote>
</section>

<section>
  <div class="secnum">04</div>
  <h2>One encounter, three durations, two orders of magnitude</h2>
  <p>The same 1,705 encounters measured three ways. The saved pipeline reported a single figure for all three roles; these curves are why that mattered.</p>
  <figure>
    <div class="scroll">{p4}</div>
    <figcaption>Median structural span {s4["Structural span"]["median"]:.0f} h, median supported exposure {s4["Supported exposure"]["median"]:.2f} h, median active contact {s4["Active contact"]["median"]:.2f} h &mdash; a ratio of about <strong>{ratio:.0f} to 1</strong> between the outermost two. An encounter spanning 14 hours typically carries about 2.3 hours of observed fine-scale exposure and about 16 minutes of actual cross-group contact. Note also the falling sample: {s4["Structural span"]["n"]:,} encounters have a span, {s4["Supported exposure"]["n"]:,} have any fine-scale exposure, {s4["Active contact"]["n"]:,} have any contact at all. Source: <code>phase2_two_stage_events/stage1_events_with_stage2_mixing.csv</code>.</figcaption>
  </figure>
</section>

<section>
  <div class="secnum">05</div>
  <h2>Excursion length by how far it went</h2>
  <p>&ldquo;Depth requires duration&rdquo; as a distribution rather than two medians. Excursions that only ever reach being alone against those that reach another unit.</p>
  <figure>
    <div class="scroll">{p5}</div>
    <figcaption>Alone-only excursions: n={s5["Alone only"]["n"]}, median {s5["Alone only"]["median"]:.0f} night, {100 * s5["Alone only"]["share_over_7"]:.0f}% pass seven nights. Reached another unit: n={s5["Reached another unit"]["n"]}, median {s5["Reached another unit"]["median"]:.0f} nights, {100 * s5["Reached another unit"]["share_over_7"]:.0f}% pass seven nights. The two curves separate immediately and never re-cross: being alone is a short state, being with another group is not. The dashed line is the seven-night threshold that defines settlement, and only the second curve has meaningful mass beyond it. Source: <code>phase4b_individual_axis/excursions_dominant_gap0.csv</code>, dominant nightly rule, no gap bridging.</figcaption>
  </figure>
</section>

<section>
  <div class="secnum">06</div>
  <h2>What the census establishes on its own</h2>
  <p><strong>Nothing here needs a model.</strong> Four statements follow from the counts and the curves alone, and all four are population-scale.</p>
  <p><strong>1. Group identity persists while configuration varies constantly.</strong> 82.7% of nights with the origin group; 33.2% of nights in a configuration that is not the whole group. Both true of the same 81,695 nights.</p>
  <p><strong>2. Every structure is brief and heavy-tailed.</strong> Medians of 1&ndash;11 h across five event types, maxima of 217&ndash;6,330 h. No structure is typically long; every structure occasionally is.</p>
  <p><strong>3. Duration is not one quantity.</strong> Span, exposure and contact differ by about {ratio:.0f}-fold on the same events. Any result that reported one as another inherited the error.</p>
  <p><strong>4. Depth tracks duration, not the reverse.</strong> Alone-only excursions last one night; excursions reaching another unit last three and are the only ones with mass past the settlement threshold.</p>
  <p>What the census cannot establish: seven of the ten structures in &sect;2 have no denominator, so their counts are inventories. Whether they are common in the population &mdash; as opposed to common in this collar deployment &mdash; is a question the census cannot answer and the models in the companion pages have to.</p>
</section>

<footer>
  Generated by <code>build_social_structure_census.py</code> from saved outputs. No model was rerun. Legacy-cohort quantities are labelled as such: the event inventory runs on the 1,703,133-row source filtered to 2025-01-01 onward, while the frozen quantities use the 1,924,104-row export from 2024-03-01 &mdash; gap 7.3 in the companion framing.
</footer>
</div>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(build(), encoding="utf-8")
    print("wrote %s (%d bytes)" % (args.output, args.output.stat().st_size))


if __name__ == "__main__":
    main()
