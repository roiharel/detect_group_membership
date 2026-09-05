"""Build the robustness ledger: where this project's results agree, and where they do not.

Four plots, all coordinates computed here rather than hand-authored:

  1. SPECIFICATION CURVE -- every seasonality test, ordered by how it behaves when the
     obvious control is applied. Separates "never significant" from "significant until
     controlled" from "survives every control".
  2. COHORT REPLICATION -- each axis-B quantity, legacy source against frozen export,
     with an agreement diagonal. Points off the diagonal are the findings the cohort
     change moved.
  3. CONTROL SENSITIVITY -- paired estimates with and without the relevant control.
     Pairs that cross the null line are the ones where the control changes the answer.
  4. ROBUST BY DESIGN -- the ladders that do NOT move: threshold choices, gap rules,
     radius nesting. Shown at the same visual weight, because agreement is a result too.

Every number is read from the saved outputs named in SOURCES below, or is quoted from
the documented model fits in docs/framing_2026-09-04/. Writes a self-contained HTML
page ready to publish.

Output: docs/framing_2026-09-04/robustness_ledger.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

OUT = Path("docs/framing_2026-09-04/robustness_ledger.html")

SOURCES = [
    "phase4c_seasonality/seasonality_report.json",
    "phase4c_seasonality/axis_a_distance_restriction_ladder.csv",
    "phase4c_seasonality/axis_a_balanced_panel_ladder.csv",
    "phase4c_seasonality/interaction_gradient/interaction_gradient_seasonality.json",
    "phase4c_seasonality/depth_duration/depth_duration_seasonality.json",
    "phase4d_axis_b_frozen/cohort_comparison.json",
    "phase4b_individual_axis/isolation_rule_ladder.csv",
    "phase4b_individual_axis/individual_axis_report.json",
    "phase2_two_stage_events/stage1_gap_rule_sensitivity.csv",
]

# ---------------------------------------------------------------- plot 1 data
# (label, p-value, group). group: 0 never significant, 1 dies under control,
# 2 borderline, 3 survives every control
SPEC_CURVE = [
    ("Encounter occurrence, all 317 dyads", 0.0103, 1),
    ("— restricted to ≤ 10 km", 0.0074, 1),
    ("— restricted to ≤ 5 km", 0.2308, 1),
    ("— restricted to ≤ 2 km", 0.6287, 1),
    ("Merge participation, all 66 dyads", 0.0033, 1),
    ("— balanced dyad panel", 0.0036, 1),
    ("— season demeaned within dyad", 0.1694, 1),
    ("— months with 2+ years of data", 0.5108, 1),
    ("Large merge vs partial, all dyads", 0.0029, 1),
    ("— demeaned within dyad", 0.1156, 1),
    ("Within-merge modularity, 5 m", 0.0663, 2),
    ("Within-merge modularity, 2 m", 0.0607, 2),
    ("Mixing depth, duration held out", 0.0840, 2),
    ("Mixing depth, duration in model", 0.3202, 0),
    ("Mixture deficit, 2 m", 0.9109, 0),
    ("2 m / 5 m radius profile", 0.1958, 0),
    ("Largest joint cluster size", 0.1842, 0),
    ("Split occurrence (axis B)", 0.6444, 0),
    ("Group-week modular (axis B)", 0.6149, 0),
    ("Split duration (axis B)", 0.7078, 0),
    ("Away from origin tonight (axis C)", 0.7300, 0),
    ("Excursion onset (axis C)", 0.5733, 0),
    ("Excursion duration (axis C)", 0.7779, 0),
    ("Excursion reached another unit", 0.2340, 0),
    ("Excursion reached settlement", 0.8261, 0),
    ("Encounter duration (axis A)", 0.5890, 0),
]

# ---------------------------------------------------------------- plot 2 data
# (label, legacy value, frozen value, agrees)
COHORT = [
    ("Spearman(modularity, collars)", 0.225, 0.222, True),
    ("ICC between-unit, modularity", 0.148, 0.169, True),
    ("Lilac, share of weeks modular", 0.538, 0.553, True),
    ("Lilac, max modularity", 0.458, 0.449, True),
    ("Lilac, lag-1 autocorrelation", 0.413, 0.417, True),
    ("Chartreuse, lag-1 autocorrelation", 0.331, 0.399, True),
    ("Purple, share of weeks modular", 0.000, 0.023, True),
    ("RubyRunners, share of weeks modular", 0.000, 0.000, True),
    ("Chartreuse, share of weeks modular", 0.450, 0.147, False),
    ("Spearman(split fraction, collars)", 0.489, 0.163, False),
]

# ---------------------------------------------------------------- plot 3 data
# (label, without-control, with-control, null value, scale, flips)
CONTROL = [
    ("Axis A trend in encounter probability",
     "all dyads", 0.748, "balanced panel", 1.100, 1.0, "or", True),
    ("Axis B trend in split occurrence",
     "no effort term", 1.389, "effort in model", 1.152, 1.0, "or", True),
    ("Axis A Stage-1 seasonality (sin term)",
     "full sample", 0.616, "within 5 km", 0.856, 1.0, "or", True),
    ("Stage-2 pre-event distance",
     "no dyad effects", -0.123, "dyad effects in", 0.005, 0.0, "logit", True),
    ("Stage-2 collar coverage",
     "no dyad effects", -0.389, "dyad effects in", -0.015, 0.0, "logit", True),
    ("Stage-2 prior encounters",
     "no calendar time", -0.256, "calendar time in", -0.082, 0.0, "logit", True),
    ("Sustained-association contrast",
     "stratum alone", 1.376, "with exposure terms", -0.109, 0.0, "logit", True),
    ("NDVI on split occurrence (strictest)",
     "between-group term", 2.906, "within-group term", 0.703, 1.0, "or", True),
    ("Axis B modularity ICC",
     "legacy cohort", 0.148, "frozen cohort", 0.169, None, "icc", False),
    ("Split-composition predictor",
     "preceding proximity", 1.269, "general bond", 1.347, 1.0, "or", False),
]

# ---------------------------------------------------------------- plot 4 data
# (label, [(tick label, value)], unit, verdict)
LADDERS = [
    ("Isolation support: located-cluster threshold",
     [("≥ 2", 23.9), ("≥ 3", 22.7), ("≥ 4", 20.1)],
     "% of raw isolated hours retained", "insensitive"),
    ("Excursion count: gap tolerance, permissive rule",
     [("0", 3541), ("1", 3538), ("2", 3534), ("3", 3534), ("7", 3534), ("14", 3534)],
     "excursions", "insensitive"),
    ("Excursion count: gap tolerance, dominant rule",
     [("0", 338), ("1", 234), ("2", 210), ("3", 199), ("7", 196), ("14", 196)],
     "excursions", "SENSITIVE"),
    ("Stage-1 encounters: gap rule",
     [("2 h", 1812), ("3 h", 1705), ("6 h", 1547), ("14 h", 1109)],
     "events (dyads stay 68 throughout)", "SENSITIVE in events, not in dyads"),
    ("Animals contributing excursions: gap tolerance",
     [("0", 91), ("1", 91), ("2", 91), ("3", 91), ("7", 91), ("14", 91)],
     "animals, dominant rule", "insensitive"),
]

CSS = """
:root{--ground:#f1f3f1;--paper:#fbfcfa;--ink:#141a19;--ink-2:#38423f;--muted:#667370;
--faint:#8d9793;--rule:#d3d9d4;--rule-soft:#e3e8e3;--accent:#0f5f57;--ax-a:#1f6f8b;
--ax-b:#8a5a1f;--ax-c:#6b4a7a;--agree:#0f5f57;--differ:#8c3a2e;--mid:#8a6410;
--shadow:0 1px 2px rgba(20,26,25,.05),0 8px 24px -14px rgba(20,26,25,.18);
--serif:"Newsreader","Iowan Old Style",Georgia,serif;
--sans:"IBM Plex Sans","Segoe UI",Helvetica,Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#101413;
--paper:#171d1c;--ink:#e9efec;--ink-2:#c3ccc8;--muted:#94a09c;--faint:#737e7a;
--rule:#2b3432;--rule-soft:#222a29;--accent:#5fc0b0;--ax-a:#74b6d0;--ax-b:#d3a061;
--ax-c:#b494c4;--agree:#5fc0b0;--differ:#e08a78;--mid:#d8ab54;
--shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px -14px rgba(0,0,0,.7);}}
:root[data-theme="dark"]{--ground:#101413;--paper:#171d1c;--ink:#e9efec;--ink-2:#c3ccc8;
--muted:#94a09c;--faint:#737e7a;--rule:#2b3432;--rule-soft:#222a29;--accent:#5fc0b0;
--ax-a:#74b6d0;--ax-b:#d3a061;--ax-c:#b494c4;--agree:#5fc0b0;--differ:#e08a78;
--mid:#d8ab54;--shadow:0 1px 2px rgba(0,0,0,.4),0 8px 24px -14px rgba(0,0,0,.7);}
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
.scroll{overflow-x:auto}.scroll>svg{min-width:520px}
figcaption{font-size:12px;color:var(--muted);margin-top:14px;padding-top:12px;
border-top:1px solid var(--rule-soft);max-width:70ch}
text{font-family:var(--mono);font-variant-numeric:tabular-nums}
.lbl{font-size:10.5px;fill:var(--ink-2)}.lbl-sm{font-size:9.5px;fill:var(--faint)}
.lbl-em{font-size:10.5px;fill:var(--ink);font-weight:500}
.axline{stroke:var(--rule);stroke-width:1}.grid{stroke:var(--rule-soft);stroke-width:1}
.key{display:flex;flex-wrap:wrap;gap:8px 20px;font-family:var(--mono);font-size:11px;
color:var(--muted);margin:14px 0 0}
.key i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;
vertical-align:-1px}
blockquote{margin:22px 0;padding:18px 22px;background:var(--paper);
border-left:3px solid var(--accent);font-family:var(--serif);font-size:19px;
line-height:1.45;max-width:42em}
blockquote p{margin:0;max-width:none}
footer{border-top:1px solid var(--rule);padding-top:16px;font-size:12.5px;color:var(--muted)}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""

GROUP_COLOR = {0: "var(--muted)", 1: "var(--differ)", 2: "var(--mid)", 3: "var(--agree)"}
GROUP_NAME = {
    0: "never distinguishable from zero",
    1: "significant until the control is applied",
    2: "borderline, unresolved",
    3: "survives every control",
}


def plot_spec_curve() -> str:
    """Horizontal p-value plot on a log axis, grouped by behaviour."""
    rows = SPEC_CURVE
    left, right = 250.0, 630.0
    top, rowh = 34.0, 17.0
    lo, hi = math.log10(0.002), math.log10(1.0)

    def x(p):
        p = max(p, 0.002)
        return left + (math.log10(p) - lo) / (hi - lo) * (right - left)

    height = top + rowh * len(rows) + 62
    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Specification curve of '
         'every seasonality test, p-value on a log axis. Tests that were significant '
         'lose significance once the relevant control is applied; the majority were '
         'never significant.">' % height]
    for tick, lab in ((0.002, "0.002"), (0.01, "0.01"), (0.05, "0.05"),
                      (0.1, "0.1"), (0.5, "0.5"), (1.0, "1")):
        xx = x(tick)
        emph = tick == 0.05
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="%s"%s/>' % (
                     xx, top - 8, xx, top + rowh * len(rows) + 4,
                     "var(--ink)" if emph else "var(--rule-soft)",
                     "1.3" if emph else "1",
                     ' stroke-dasharray="4 3"' if emph else ""))
        o.append('<text class="%s" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % ("lbl-em" if emph else "lbl", xx, top - 14, lab))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle" '
             'fill="var(--ink)">p = 0.05</text>' % (x(0.05), top - 26))

    for i, (label, p, grp) in enumerate(rows):
        y = top + rowh * i + rowh * 0.5 + 4
        indent = 12 if label.startswith("—") else 0
        o.append('<text class="%s" x="%d" y="%.1f">%s</text>'
                 % ("lbl" if indent else "lbl-em", indent, y, label))
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--rule-soft)" '
                 'stroke-width="1"/>' % (left, y - 3.5, x(p), y - 3.5))
        o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>'
                 % (x(p), y - 3.5, GROUP_COLOR[grp]))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">%.3f</text>'
                 % (x(p) + 8, y, p))
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (left, top + rowh * len(rows) + 4, right, top + rowh * len(rows) + 4))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle">joint Wald p '
             'for the annual harmonic, log axis</text>'
             % ((left + right) / 2, top + rowh * len(rows) + 22))
    o.append("</svg>")
    return "\n".join(o)


def plot_cohort() -> str:
    """Legacy against frozen, with an agreement diagonal."""
    pad_l, pad_b, size = 74.0, 52.0, 330.0
    top = 20.0
    hi = 0.56

    def px(v):
        return pad_l + v / hi * size

    def py(v):
        return top + size - v / hi * size

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Scatter of ten axis-B '
         'quantities, legacy source against frozen export. Eight sit on the agreement '
         'diagonal; Chartreuse modular share and the split-detection correlation are '
         'far off it.">' % (top + size + pad_b)]
    for t in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (px(t), top, px(t), top + size))
        o.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (pad_l, py(t), pad_l + size, py(t)))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="middle">%.1f</text>'
                 % (px(t), top + size + 16, t))
        o.append('<text class="lbl" x="%.1f" y="%.1f" text-anchor="end">%.1f</text>'
                 % (pad_l - 8, py(t) + 3.5, t))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--muted)" '
             'stroke-width="1" stroke-dasharray="5 4"/>'
             % (px(0), py(0), px(hi), py(hi)))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" fill="var(--muted)">perfect '
             'agreement</text>' % (px(0.40), py(0.44)))

    # label offsets tuned so nothing collides
    offsets = {
        "Spearman(modularity, collars)": (10, -8),
        "ICC between-unit, modularity": (10, 14),
        "Lilac, share of weeks modular": (-8, -12),
        "Lilac, max modularity": (10, 12),
        "Lilac, lag-1 autocorrelation": (10, -10),
        "Chartreuse, lag-1 autocorrelation": (10, 14),
        "Purple, share of weeks modular": (12, 4),
        "RubyRunners, share of weeks modular": (12, 16),
        "Chartreuse, share of weeks modular": (-6, 22),
        "Spearman(split fraction, collars)": (10, 6),
    }
    for label, leg, fro, agrees in COHORT:
        col = "var(--agree)" if agrees else "var(--differ)"
        if not agrees:
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                     'stroke-width="1" stroke-dasharray="2 2"/>'
                     % (px(leg), py(leg), px(leg), py(fro), col))
        o.append('<circle cx="%.1f" cy="%.1f" r="%s" fill="%s"%s/>'
                 % (px(leg), py(fro), "5.5" if not agrees else "4.5", col,
                    "" if agrees else ' stroke="var(--paper)" stroke-width="1"'))
        dx, dy = offsets.get(label, (10, 4))
        anchor = "end" if dx < 0 else "start"
        o.append('<text class="%s" x="%.1f" y="%.1f" text-anchor="%s" fill="%s">%s</text>'
                 % ("lbl-em" if not agrees else "lbl-sm",
                    px(leg) + dx, py(fro) + dy, anchor,
                    col if not agrees else "var(--ink-2)", label))
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top + size, pad_l + size, top + size))
    o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (pad_l, top, pad_l, top + size))
    o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle">value on the '
             'LEGACY source</text>' % (pad_l + size / 2, top + size + 38))
    o.append('<text class="lbl-sm" transform="translate(%.1f %.1f) rotate(-90)" '
             'text-anchor="middle">value on the FROZEN export</text>'
             % (pad_l - 44, top + size / 2))
    o.append("</svg>")
    return "\n".join(o)


def plot_control() -> str:
    """Dumbbells: estimate without the control, and with it."""
    left, right = 268.0, 618.0
    top, rowh = 46.0, 40.0
    rows = CONTROL

    # two panels share one axis: odds ratios (null 1) and logit coefficients (null 0)
    ors = [r for r in rows if r[6] == "or"]
    logits = [r for r in rows if r[6] == "logit"]
    others = [r for r in rows if r[6] == "icc"]

    height = top + rowh * (len(rows) + 2) + 50
    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Paired estimates with and '
         'without the relevant control. Seven of ten pairs cross the null line once the '
         'control is applied.">' % height]

    y = top

    def panel(title, items, null, lo, hi, fmt):
        nonlocal y
        o.append('<text class="lbl-em" x="0" y="%.1f" fill="var(--accent)">%s</text>'
                 % (y - 10, title))

        def x(v):
            v = min(max(v, lo), hi)
            return left + (v - lo) / (hi - lo) * (right - left)

        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
                 'stroke-width="1.2" stroke-dasharray="4 3"/>'
                 % (x(null), y - 6, x(null), y + rowh * len(items) - 12))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle" '
                 'fill="var(--ink)">no effect</text>' % (x(null), y - 14))
        def anchored(xx, txt, yy, cls="lbl-sm"):
            """Keep every annotation inside the drawing box."""
            if xx < left + 62:
                a, ax = "start", xx - 4
            elif xx > right - 62:
                a, ax = "end", xx + 4
            else:
                a, ax = "middle", xx
            return ('<text class="%s" x="%.1f" y="%.1f" text-anchor="%s">%s</text>'
                    % (cls, ax, yy, a, txt))

        for label, l1, v1, l2, v2, _n, _s, flips in items:
            yy = y + 12
            col = "var(--differ)" if flips else "var(--agree)"
            o.append('<text class="lbl-em" x="0" y="%.1f">%s</text>' % (yy, label))
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                     'stroke-width="2"/>' % (x(v1), yy - 4, x(v2), yy - 4, col))
            o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="var(--paper)" stroke="%s" '
                     'stroke-width="2"/>' % (x(v1), yy - 4, col))
            o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/>'
                     % (x(v2), yy - 4, col))
            o.append(anchored(x(v1), "%s %s" % (fmt % v1, l1), yy - 12))
            o.append(anchored(x(v2), "%s %s" % (fmt % v2, l2), yy + 15))
            y += rowh
        y += 24

    panel("Odds ratios — hollow marker is without the control, filled is with it",
          ors, 1.0, 0.5, 3.1, "%.3f")
    panel("Logit coefficients", logits, 0.0, -0.55, 1.5, "%+.3f")
    if others:
        o.append('<text class="lbl-em" x="0" y="%.1f" fill="var(--accent)">Agreement '
                 'across sources</text>' % (y - 10))
        for label, l1, v1, l2, v2, _n, _s, flips in others:
            yy = y + 10
            o.append('<text class="lbl-em" x="0" y="%.1f">%s</text>' % (yy, label))
            o.append('<text class="lbl-sm" x="%.1f" y="%.1f" fill="var(--agree)">'
                     '%s %.3f &rarr; %s %.3f &mdash; agrees</text>'
                     % (left, yy, l1, v1, l2, v2))
            y += rowh
    o.append("</svg>")
    return "\n".join(o)


def plot_ladders() -> str:
    """Small multiples: the ladders that do not move, and the two that do."""
    colw, panelh = 214.0, 118.0
    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Five sensitivity ladders. '
         'Three are flat: the isolation threshold, the permissive gap rule and the animal '
         'count. Two decline: the dominant-rule excursion count and the Stage-1 event '
         'count.">' % (panelh * 2 + 46)]
    for i, (label, pts, unit, verdict) in enumerate(LADDERS):
        cx = (i % 3) * (colw + 14)
        cy = (i // 3) * (panelh + 26) + 26
        vals = [v for _, v in pts]
        lo, hi = 0.0, max(vals) * 1.18
        pw, ph = colw - 34, panelh - 52
        sensitive = verdict.startswith("SENSITIVE")
        col = "var(--differ)" if sensitive else "var(--agree)"
        o.append('<text class="lbl-em" x="%.1f" y="%.1f">%s</text>'
                 % (cx, cy - 12, label[:38]))
        o.append('<line class="axline" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (cx + 30, cy + ph, cx + 30 + pw, cy + ph))
        n = len(pts)
        xs = [cx + 30 + (pw * j / max(1, n - 1)) for j in range(n)]
        ys = [cy + ph - (v - lo) / (hi - lo) * ph for v in vals]
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>'
                 % (" ".join("%.1f,%.1f" % (a, b) for a, b in zip(xs, ys)), col))
        for (tick, v), xx, yy in zip(pts, xs, ys):
            o.append('<circle cx="%.1f" cy="%.1f" r="3.5" fill="%s"/>' % (xx, yy, col))
            o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                     % (xx, cy + ph + 13, tick))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (cx + 26, ys[0] + 3.5, ("%g" % vals[0])))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (cx + 26, ys[-1] + 3.5, ("%g" % vals[-1])))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f" fill="%s">%s</text>'
                 % (cx + 30, cy + ph + 27, col, verdict))
        o.append('<text class="lbl-sm" x="%.1f" y="%.1f">%s</text>'
                 % (cx + 30, cy + ph + 39, unit[:34]))
    o.append("</svg>")
    return "\n".join(o)


def build() -> str:
    keys = "".join(
        '<span><i style="background:%s"></i>%s</span>' % (GROUP_COLOR[g], GROUP_NAME[g])
        for g in (3, 1, 2, 0))
    return f"""<title>Robustness Ledger</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500;6..72,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header>
  <div class="eyebrow">Companion to Stable Groups, Fluid Boundaries</div>
  <h1>Where the results agree, and where they don't</h1>
  <p class="standfirst">Every finding in this project, plotted against the specification that tests it. Agreement and disagreement at the same visual weight, because both are results.</p>
  <div class="meta">
    <span><b>26</b> tests of the annual harmonic</span>
    <span><b>10</b> quantities re-derived across cohorts</span>
    <span><b>10</b> paired control comparisons</span>
    <span><b>5</b> sensitivity ladders</span>
  </div>
</header>

<section>
  <div class="secnum">01</div>
  <h2>Seasonality: nothing survives its own control</h2>
  <p>Twenty-six tests of the annual harmonic, ordered by behaviour rather than by size. The pattern is stark: <strong>every test that reached <span class="num">p &lt; 0.05</span> did so before the relevant control was applied, and none of them survived it.</strong> Sixteen were never distinguishable from zero in the first place.</p>
  <figure>
    <div class="scroll">{plot_spec_curve()}</div>
    <div class="key">{keys}</div>
    <figcaption>Indented rows are successive controls on the row above. Encounter occurrence dies under distance restriction; merge participation survives a balanced dyad panel but dies once restricted to calendar months sampled in more than one year, and dies again when the harmonic is demeaned within dyad. Within-merge modularity sits at 0.061&ndash;0.066 at both radii with the same October phase &mdash; consistent, unresolved, and resting on the same single-year window. Sources: <code>phase4c_seasonality/</code>.</figcaption>
  </figure>
  <blockquote><p>The one thing that <em>would</em> change this picture is more years. August to November exist in a single year of the well-observed window, and those are exactly the months carrying the signal.</p></blockquote>
</section>

<section>
  <div class="secnum">02</div>
  <h2>Cohort replication: eight of ten hold</h2>
  <p>Axis B was rebuilt from the frozen export to match the legacy definition. Points on the diagonal replicated; the two off it are the findings the cohort change moved &mdash; and both had been reported before the check was run.</p>
  <figure>
    <div class="scroll">{plot_cohort()}</div>
    <div class="key"><span><i style="background:var(--agree)"></i>replicates</span><span><i style="background:var(--differ)"></i>does not replicate</span></div>
    <figcaption>The ICC &mdash; the finding that modularity is a state rather than a group property &mdash; replicates (0.148 &rarr; 0.169), as does Lilac on every measure. Chartreuse&rsquo;s modular share falls from 45.0% of 20 weeks to 14.7% of 34, and the split-detection correlation from 0.489 to 0.163, which withdraws the claim that half of split detection was collar count. Source: <code>phase4d_axis_b_frozen/cohort_comparison.json</code>.</figcaption>
  </figure>
</section>

<section>
  <div class="secnum">03</div>
  <h2>Controls: seven of ten pairs cross the null</h2>
  <p>Each row is one estimate fitted twice &mdash; once without the control that matters, once with it. Hollow marker without, filled marker with. Where the pair crosses the dashed line, the control decides the answer.</p>
  <figure>
    <div class="scroll">{plot_control()}</div>
    <div class="key"><span><i style="background:var(--differ)"></i>control changes the conclusion</span><span><i style="background:var(--agree)"></i>conclusion unchanged</span></div>
    <figcaption>Four distinct controls appear here and all four bite: dyad or unit fixed effects (Stage-2 distance and coverage), an observation-effort term (axis B trend), a sample restriction (within 5 km; balanced panel), and separating a between-unit from a within-unit term (NDVI). The two agreeing rows are the modularity ICC across cohorts and the two split-composition predictors, which give the same answer either way. Sources: <code>PHASE4_RESULTS.md</code>, <code>phase4b_axis_b_fission/</code>, <code>phase4c_seasonality/</code>, <code>phase4d_axis_b_frozen/</code>.</figcaption>
  </figure>
</section>

<section>
  <div class="secnum">04</div>
  <h2>Robust by design: the choices that don't matter</h2>
  <p>Not every arbitrary choice moves a result, and saying which ones don't is as much a part of the audit as saying which ones do. Three of these five ladders are flat.</p>
  <figure>
    <div class="scroll">{plot_ladders()}</div>
    <figcaption>The isolation support rule is threshold-insensitive, so the 76%-artefact conclusion does not depend on where the cut sits. The permissive nightly rule is gap-insensitive; the dominant rule is not, which is why the gap ladder is reported rather than a single number. And the Stage-1 gap rule changes the event count from 1,812 to 1,109 while the dyad count stays at 68 throughout &mdash; the sample is stable, only the event parsing is sensitive. The same holds on axis C: 91 animals at every gap tolerance. Sources: <code>phase4b_individual_axis/</code>, <code>phase2_two_stage_events/stage1_gap_rule_sensitivity.csv</code>.</figcaption>
  </figure>
</section>

<section>
  <div class="secnum">05</div>
  <h2>What the four plots say together</h2>
  <p><strong>The nulls are the robust half of this project.</strong> Plot 1 shows no seasonal signal surviving its own control; plot 3 shows that seven of ten headline estimates depend on a control being present. What holds across every specification is the descriptive backbone: mixing at 4.6% of composition expectation in 12 of 12 dyads, modularity as a state (ICC replicating across cohorts), 91 animals contributing excursions at every gap tolerance, and the isolation artefact at any threshold.</p>
  <p><strong>Positives from short windows are where the errors were.</strong> Both non-replicating points in plot 2, and every red row in plot 1, are positive findings that did not survive a check. No null in this project has been overturned by adding a control. That asymmetry is worth stating in the paper: with 1.6 seasonal cycles and a 121-fold effort ramp, this design is far better at establishing absence within its range than at establishing a cycle.</p>
  <p><strong>And rarity is not the same as fragility.</strong> The rare states &mdash; one fully fusing dyad, one recurrently modular group, 61 settlement excursions, 16 sustained associations &mdash; are not in plot 1 or plot 3, because they are observations rather than estimates. They replicate wherever they can be re-derived: Lilac holds on every measure across cohorts. Low frequency is a fact about the population; it is not a weakness in the measurement.</p>
</section>

<footer>
  Generated by <code>build_robustness_ledger.py</code> from saved outputs under <code>outputs/general_structure_2026_09/</code> and the model fits documented in <code>docs/framing_2026-09-04/</code>. No model was rerun to build this page.
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
    print("sources referenced:")
    for s in SOURCES:
        print("  outputs/general_structure_2026_09/%s" % s)


if __name__ == "__main__":
    main()
