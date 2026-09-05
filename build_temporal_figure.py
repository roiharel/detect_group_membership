"""Manuscript figure: what in social structure is a group property, and what drifts.

  (a) A hierarchy of stability. For each measured quantity, the share of variance that
      sits between groups rather than within a group over time, with the null test that
      asks whether groups differ once quarter-to-quarter drift is allowed.
  (b) Merge-scale composition quarter by quarter, one line per group. The between-group
      ordering in Figure 3 does not hold: well-sampled groups traverse most of the range
      inside their own record.
  (c) Merge-partner fidelity. Jaccard similarity of a group's partner set between
      consecutive quarters, against a null that redraws from the partners actually
      available that quarter.

Output: docs/framing_2026-09-04/temporal_figure.html
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("outputs/general_structure_2026_09/phase4f_temporal_variation")
MODB = Path("outputs/general_structure_2026_09/phase4d_axis_b_frozen/cohort_comparison.json")
OUT = Path("docs/framing_2026-09-04/temporal_figure.html")

QUARTERS = ["2025Q1", "2025Q2", "2025Q3", "2025Q4", "2026Q1", "2026Q2"]


def load():
    rep = json.load(open(BASE / "temporal_variation_report.json", encoding="utf-8"))
    mrg = pd.read_csv(BASE / "composition_by_group_quarter_merge.csv")
    turn = pd.read_csv(BASE / "partner_turnover.csv")
    mod = json.load(open(MODB, encoding="utf-8"))
    return rep, mrg, turn, mod


# ------------------------------------------------------------------ panel a
def panel_a(rep, mod, y0):
    """Stability hierarchy: ICC per quantity, with the permutation verdict."""
    fa = rep["family"]["anova_permutation"]["by_proportion"]
    ma = rep["merge_scale"]["anova_permutation"]["by_proportion"]
    rows = [
        ("Event rate per group-week", rep["event_rate"]["icc_between_group"], None,
         "%d group-quarters, 23 groups" % rep["event_rate"]["group_quarters"], True),
        ("Merge-partner identity", None, None, "Jaccard %.2f vs null %.2f (%.1f×)"
         % (rep["partner_turnover"]["mean_jaccard"],
            rep["partner_turnover"]["null_mean_jaccard"],
            rep["partner_turnover"]["ratio_to_null"]), True),
        ("Family: individual separations", fa["individual"]["icc_anova"],
         fa["individual"]["permutation_p"], "76 group-quarters, 14 groups", True),
        ("Family: between-group merges", fa["group_merge"]["icc_anova"],
         fa["group_merge"]["permutation_p"], "", True),
        ("Family: within-group splits", fa["within_group"]["icc_anova"],
         fa["within_group"]["permutation_p"], "", True),
        ("Within-group modularity", mod["frozen"]["icc"]["icc_between_unit"], None,
         "368 group-weeks, 12 units", False),
        ("Merge scale: small subset", ma["small_subset_merge"]["icc_anova"],
         ma["small_subset_merge"]["permutation_p"], "54 group-quarters, 10 groups", False),
        ("Merge scale: large merge", ma["large_merge"]["icc_anova"],
         ma["large_merge"]["permutation_p"], "", False),
        ("Merge scale: medium partial", ma["medium_partial_merge"]["icc_anova"],
         ma["medium_partial_merge"]["permutation_p"], "", False),
    ]
    x0, w, rowh = 176.0, 300.0, 20.0
    o = []
    for t in (0, 0.25, 0.5, 0.75, 1.0):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0 + t * w, y0 - 6, x0 + t * w, y0 + rowh * len(rows)))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (x0 + t * w, y0 + rowh * len(rows) + 12,
                    "0" if t == 0 else ("1" if t == 1 else "%.2f" % t)))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="1" stroke-dasharray="3 2" opacity=".6"/>'
             % (x0 + 0.5 * w, y0 - 6, x0 + 0.5 * w, y0 + rowh * len(rows)))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">half</text>'
             % (x0 + 0.5 * w, y0 - 9))
    for r, (label, icc, p, note, stable) in enumerate(rows):
        y = y0 + rowh * r + rowh * 0.62
        col = "var(--c1)" if stable else "var(--c4)"
        o.append('<text class="gl" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 - 9, y, label))
        if icc is not None:
            o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="9" fill="%s" '
                     'opacity=".85"/>' % (x0, y - 8, icc * w, col))
            txt = "%.2f" % icc
            if p is not None:
                txt += "   p = %s" % ("< 0.001" if p < 0.001 else "%.3f" % p)
            o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>'
                     % (x0 + icc * w + 7, y, txt))
        else:
            o.append('<text class="ts" x="%.1f" y="%.1f" fill="%s">not a variance '
                     'ratio &mdash; see note</text>' % (x0 + 4, y, col))
        if note:
            o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>'
                     % (x0 + w + 58, y, note))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + rowh * len(rows), x0 + w, y0 + rowh * len(rows)))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">share of variance '
             'between groups rather than within a group over time</text>'
             % (x0 + w / 2, y0 + rowh * len(rows) + 27))
    return "\n".join(o), rowh * len(rows) + 34


# ------------------------------------------------------------------ panel b
def panel_b(mrg, y0):
    use = mrg[mrg["n"] >= 20].copy()
    use["qi"] = use["quarter"].map({q: i for i, q in enumerate(QUARTERS)})
    use = use[use["qi"].notna()]
    counts = use.groupby("group")["qi"].nunique()
    hi_n = counts[counts >= 5].index          # well-sampled groups get labelled
    x0, w, h = 176.0, 300.0, 168.0
    o = []
    for i, q in enumerate(QUARTERS):
        xx = x0 + i / (len(QUARTERS) - 1) * w
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (xx, y0, xx, y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (xx, y0 + h + 12, q.replace("20", "'")))
    for t in (0, 0.5, 1.0):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0, y0 + h - t * h, x0 + w, y0 + h - t * h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 - 6, y0 + h - t * h + 3, "0" if t == 0 else
                    ("1" if t == 1 else ".5")))
    for g, sub in use.groupby("group"):
        sub = sub.sort_values("qi")
        if len(sub) < 2:
            continue
        wide = g in hi_n
        col = "var(--c1)" if wide else "var(--muted)"
        pts = " ".join("%.1f,%.1f" % (x0 + r.qi / (len(QUARTERS) - 1) * w,
                                      y0 + h - r.p_large_merge * h)
                       for r in sub.itertuples())
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s" '
                 'opacity="%s" stroke-linejoin="round"/>'
                 % (pts, col, "1.6" if wide else "1", ".9" if wide else ".45"))
        for r in sub.itertuples():
            o.append('<circle cx="%.1f" cy="%.1f" r="%s" fill="%s" opacity="%s"/>'
                     % (x0 + r.qi / (len(QUARTERS) - 1) * w,
                        y0 + h - r.p_large_merge * h,
                        "2.4" if wide else "1.6", col, ".9" if wide else ".45"))
        if wide:
            last = sub.iloc[-1]
            o.append('<text class="ts" x="%.1f" y="%.1f" fill="var(--c1)">%s</text>'
                     % (x0 + last["qi"] / (len(QUARTERS) - 1) * w + 6,
                        y0 + h - last["p_large_merge"] * h + 3, g))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + h, x0 + w, y0 + h))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0, x0, y0 + h))
    o.append('<text class="ts" transform="translate(%.1f %.1f) rotate(-90)" '
             'text-anchor="middle">share of merges that are large</text>'
             % (x0 - 26, y0 + h / 2))
    return "\n".join(o), h + 26


# ------------------------------------------------------------------ panel c
def panel_c(turn, rep, y0):
    j = np.sort(turn["jaccard"].to_numpy())
    nullm = rep["partner_turnover"]["null_mean_jaccard"]
    obsm = rep["partner_turnover"]["mean_jaccard"]
    x0, w, h = 176.0, 300.0, 96.0
    o = []
    for t in (0, 0.25, 0.5, 0.75, 1.0):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0 + t * w, y0, x0 + t * w, y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (x0 + t * w, y0 + h + 12,
                    "0" if t == 0 else ("1" if t == 1 else "%.2f" % t)))
    # empirical CDF of observed Jaccard
    pts = " ".join("%.1f,%.1f" % (x0 + v * w, y0 + h - (i + 1) / len(j) * h)
                   for i, v in enumerate(j))
    o.append('<polyline points="%.1f,%.1f %s" fill="none" stroke="var(--c1)" '
             'stroke-width="1.8"/>' % (x0 + j[0] * w, y0 + h, pts))
    for val, col, lab in ((nullm, "var(--c4)", "null mean"),
                          (obsm, "var(--c1)", "observed mean")):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="1.2" stroke-dasharray="4 3"/>'
                 % (x0 + val * w, y0 - 5, x0 + val * w, y0 + h, col))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" fill="%s">'
                 '%s %.2f</text>' % (x0 + val * w, y0 - 8, col, lab, val))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + h, x0 + w, y0 + h))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0, x0, y0 + h))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">Jaccard '
             'similarity of merge-partner sets, consecutive quarters (n = %d)</text>'
             % (x0 + w / 2, y0 + h + 27, len(j)))
    o.append('<text class="ts" transform="translate(%.1f %.1f) rotate(-90)" '
             'text-anchor="middle">cumulative</text>' % (x0 - 26, y0 + h / 2))
    return "\n".join(o), h + 34


def build():
    rep, mrg, turn, mod = load()
    ya = 56.0
    pa, ha = panel_a(rep, mod, ya)
    yb = ya + ha + 60
    pb, hb = panel_b(mrg, yb)
    yc = yb + hb + 62
    pc, hc = panel_c(turn, rep, yc)
    H = yc + hc + 20
    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Three panels. Event rate '
         'and partner identity are stable group properties; merge-scale composition is '
         'not distinguishable from shuffled group labels and well-sampled groups '
         'traverse most of its range within their own record.">' % H]
    o.append('<text class="pl" x="0" y="%.1f">a</text>' % (ya - 30))
    o.append('<text class="ph" x="18" y="%.1f">A hierarchy of stability</text>'
             % (ya - 30))
    o.append(pa)
    o.append('<text class="pl" x="0" y="%.1f">b</text>' % (yb - 28))
    o.append('<text class="ph" x="18" y="%.1f">Merge scale, quarter by quarter</text>'
             % (yb - 28))
    o.append(pb)
    o.append('<text class="pl" x="0" y="%.1f">c</text>' % (yc - 28))
    o.append('<text class="ph" x="18" y="%.1f">Merge-partner fidelity</text>'
             % (yc - 28))
    o.append(pc)
    o.append("</svg>")
    return "\n".join(o), rep, mrg


CSS = """
:root{--ground:#f4f5f3;--paper:#fff;--ink:#15191a;--ink-2:#3a4244;--muted:#6b7476;
--rule:#d5dad7;--rule-soft:#e6eae7;--accent:#0f5f57;
--c1:#1f6f8b;--c2:#4a92aa;--c3:#84b6c6;--c4:#8a5a1f;--c5:#6b4a7a;
--sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;
--serif:"Newsreader",Georgia,serif;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0f1213;
--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;--muted:#8f9998;--rule:#2c3435;
--rule-soft:#232a2b;--accent:#5fc0b0;--c1:#74b6d0;--c2:#4f93ad;--c3:#356a7f;
--c4:#d3a061;--c5:#b494c4;}}
:root[data-theme="dark"]{--ground:#0f1213;--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;
--muted:#8f9998;--rule:#2c3435;--rule-soft:#232a2b;--accent:#5fc0b0;--c1:#74b6d0;
--c2:#4f93ad;--c3:#356a7f;--c4:#d3a061;--c5:#b494c4;}
*{box-sizing:border-box}
body{background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:15px;
line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:800px;margin:0 auto;padding:48px 28px 80px}
.kicker{font-family:var(--mono);font-size:10px;letter-spacing:.16em;
text-transform:uppercase;color:var(--muted);margin-bottom:10px}
h1{font-family:var(--serif);font-size:30px;font-weight:500;line-height:1.15;
letter-spacing:-.01em;margin:0 0 26px;text-wrap:balance}
.plate{background:var(--paper);border:1px solid var(--rule-soft);padding:22px 22px 8px}
.plate svg{display:block;width:100%;height:auto;overflow:visible}
.scroll{overflow-x:auto}.scroll>svg{min-width:640px}
.cap{font-size:12.5px;line-height:1.62;color:var(--ink-2);margin:16px 0 0;
padding-top:14px;border-top:1px solid var(--rule)}
.cap b{font-weight:600;color:var(--ink)}.cap .fignum{font-weight:600;color:var(--ink)}
.note{font-size:11.5px;line-height:1.6;color:var(--muted);margin-top:22px;
padding-top:14px;border-top:1px solid var(--rule-soft)}
.note code{font-family:var(--mono);font-size:.92em}
text{font-family:var(--sans)}
.gl{font-size:9.4px;fill:var(--ink-2)}
.ts{font-size:8.4px;fill:var(--muted)}
.ph{font-size:11px;fill:var(--ink);font-weight:600}
.pl{font-size:14px;font-weight:700;fill:var(--ink)}
.g{stroke:var(--rule-soft);stroke-width:.8}
.ax{stroke:var(--rule);stroke-width:1}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def page():
    fig, rep, mrg = build()
    fa = rep["family"]["anova_permutation"]["by_proportion"]
    ma = rep["merge_scale"]["anova_permutation"]["by_proportion"]
    pt = rep["partner_turnover"]
    er = rep["event_rate"]
    sw = pd.DataFrame(rep["within_group_swing_p_large_merge"])
    lap = sw[sw["group"].eq("Lapis")].iloc[0]
    lsp = sw[sw["group"].eq("LapisSplinter")].iloc[0]
    return f"""<title>Temporal Variation in Social Structure</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 4 &middot; draft</div>
  <h1>How often and with whom is a group property; how is not</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 4.</span> Temporal decomposition of social structure over six
      calendar quarters (2025Q1&ndash;2026Q2).
      <b>(a)</b> For each quantity, the share of variance sitting between groups rather than
      within a group over time, with a permutation test that shuffles group labels across
      group-quarters and so keeps quarter-to-quarter drift intact. Event rate per group-week is
      strongly a group property (ICC {er["icc_between_group"]:.2f}, and collar count does not
      predict it, p&nbsp;=&nbsp;{er["collar_pvalue"]:.2f}). Event <em>family</em> composition is
      about half a group property (ICC {fa["individual"]["icc_anova"]:.2f},
      {fa["group_merge"]["icc_anova"]:.2f}, {fa["within_group"]["icc_anova"]:.2f}; all
      p&nbsp;&le;&nbsp;0.001). Merge-<em>scale</em> composition is not
      (ICC {ma["large_merge"]["icc_anova"]:.2f}&ndash;{ma["small_subset_merge"]["icc_anova"]:.2f},
      permutation p&nbsp;=&nbsp;{ma["small_subset_merge"]["permutation_p"]:.2f}&ndash;{ma["large_merge"]["permutation_p"]:.2f}),
      nor is within-group modularity (ICC 0.17).
      <b>(b)</b> The same merge-scale proportion quarter by quarter, one line per group; groups
      with five or more usable quarters are drawn heavy and labelled. Lapis traverses
      {lap["min"]:.2f}&ndash;{lap["max"]:.2f} across its own six quarters on
      {int(lap["total_events"])} merges &mdash; a swing of {lap["swing"]:.2f}, wider than the
      whole between-group spread. LapisSplinter is the exception, holding
      {lsp["min"]:.2f}&ndash;{lsp["max"]:.2f} over six quarters.
      <b>(c)</b> Empirical distribution of merge-partner set overlap between consecutive
      quarters, {pt["consecutive_quarter_pairs"]} group-quarter transitions. Mean Jaccard
      {pt["mean_jaccard"]:.2f} against {pt["null_mean_jaccard"]:.2f} for a null that redraws the
      same number of partners from those actually available that quarter &mdash;
      {pt["ratio_to_null"]:.1f}&times;. Groups keep their partners.
    </p>
  </div>

  <p class="note">
    <b>This corrects Figure 3.</b> That figure tested between-group composition against a
    multinomial null in which each group drew from a fixed pooled distribution with no temporal
    structure. Panel (b) here shows why that null is too weak for merge scale: once
    quarter-to-quarter drift is retained, the between-group ordering is not distinguishable from
    shuffled labels. Figure 3's family-composition result stands; its merge-scale panel should be
    read as showing variation in the population, not differences between groups. Legacy hourly
    source, 2025-01-01 onward; group-quarters need 20 events and groups need 4 quarters to enter
    (a). Merges are counted for both participating groups. Modularity ICC is from the frozen
    cohort. Generated by <code>build_temporal_figure.py</code> from
    <code>analyze_temporal_variation.py</code>; seed 20260904.
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
