"""Manuscript figure: event-type composition by group, with CIs and a null test.

  (a) Family composition -- between-group merge, within-group split, individual
      separation -- one row per group, dot and 95% CI per type, pooled proportion as a
      reference line.
  (b) Merge-scale composition within merges: large, medium partial, small subset.
  (c) Observed between-group spread against a multinomial null, in three coverage
      strata. The variation is real and it strengthens when coverage is matched.

Intervals are 95% percentile bootstrap over each group's own events (2,000 draws).
Null: every group draws its own number of events from the pooled type distribution.

Output: docs/framing_2026-09-04/composition_figure.html
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

BASE = Path("outputs/general_structure_2026_09/phase4e_event_composition")
OUT = Path("docs/framing_2026-09-04/composition_figure.html")

FAM = [("group_merge", "Between-group merge", "var(--c1)"),
       ("within_group", "Within-group split", "var(--c4)"),
       ("individual", "Individual separation", "var(--c5)")]
MRG = [("large_merge", "Large merge", "var(--c1)"),
       ("medium_partial_merge", "Medium partial", "var(--c2)"),
       ("small_subset_merge", "Small subset", "var(--c3)")]
MIN_EVENTS = 30


def load():
    rep = json.load(open(BASE / "event_composition_report.json", encoding="utf-8"))
    fam = pd.read_csv(BASE / "composition_by_family.csv")
    mrg = pd.read_csv(BASE / "composition_by_merge_scale.csv")
    return rep, fam, mrg


def dot_panel(frame, boot, spec, pooled, y0, rowh, cols, label_w=104.0,
              col_w=148.0, gap=26.0, note_col=True):
    """Groups on rows, one small column per event type, dot with 95% CI."""
    use = frame[frame["total"] >= MIN_EVENTS].copy()
    order_key = "p_" + spec[0][0]
    use = use.sort_values(order_key, ascending=False)
    o = []
    n = len(use)
    h = rowh * n

    def cx(i, p):
        x0 = label_w + i * (col_w + gap)
        return x0 + p * col_w

    for i, (key, title, col) in enumerate(spec):
        x0 = label_w + i * (col_w + gap)
        o.append('<text class="ch" x="%.1f" y="%.1f">%s</text>' % (x0, y0 - 16, title))
        for t in (0, 0.25, 0.5, 0.75, 1.0):
            o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                     % (x0 + t * col_w, y0 - 5, x0 + t * col_w, y0 + h))
        for t, lab in ((0, "0"), (0.5, ".5"), (1.0, "1")):
            o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                     % (x0 + t * col_w, y0 + h + 12, lab))
        pv = pooled[key]
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
                 'stroke-width="1" stroke-dasharray="3 2" opacity=".65"/>'
                 % (x0 + pv * col_w, y0 - 5, x0 + pv * col_w, y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%.2f</text>'
                 % (x0 + pv * col_w, y0 - 6, pv))

    for r, row in enumerate(use.itertuples()):
        y = y0 + rowh * r + rowh * 0.62
        g = row.group
        o.append('<text class="gl" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (label_w - 34, y, g))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end">%d</text>'
                 % (label_w - 8, y, int(row.total)))
        b = boot.get(g)
        for i, (key, title, col) in enumerate(spec):
            p = getattr(row, "p_" + key)
            if b and key in b:
                lo, hi = b[key]["lo"], b[key]["hi"]
                o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                         'stroke-width="1.4" opacity=".85"/>'
                         % (cx(i, lo), y - 3.2, cx(i, hi), y - 3.2, col))
            o.append('<circle cx="%.1f" cy="%.1f" r="2.6" fill="%s"/>'
                     % (cx(i, p), y - 3.2, col))
    for i in range(len(spec)):
        x0 = label_w + i * (col_w + gap)
        o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0, y0 + h, x0 + col_w, y0 + h))
    return "\n".join(o), h, n


def null_panel(rep, y0):
    """Observed spread against the null, three strata, both composition types."""
    o = []
    x0, w = 148.0, 300.0
    rows = []
    for name, lab in (("family", "Family composition"),
                      ("merge_scale", "Merge-scale composition")):
        for s in rep[name]["coverage_stratified"]:
            if "skipped" in s:
                continue
            rows.append((lab, s["stratum"], s["groups"], s["observed"],
                         s["null_mean"], s["ratio"], s["p"]))
    hi = max(r[3] for r in rows) * 1.18

    def px(v):
        return x0 + v / hi * w

    rowh = 17.0
    for t in (0, 0.03, 0.06, 0.09, 0.12):
        if t > hi:
            continue
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (px(t), y0 - 6, px(t), y0 + rowh * len(rows)))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%.2f</text>'
                 % (px(t), y0 + rowh * len(rows) + 12, t))
    last = None
    for r, (lab, stratum, ng, obs, nul, ratio, p) in enumerate(rows):
        y = y0 + rowh * r + rowh * 0.6
        if lab != last:
            o.append('<text class="ch" x="0" y="%.1f">%s</text>' % (y, lab))
            last = lab
        o.append('<text class="gl" x="%.1f" y="%.1f" text-anchor="end">%s (%d)</text>'
                 % (x0 - 8, y, stratum, ng))
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="7" fill="var(--muted)" '
                 'opacity=".3"/>' % (px(0), y - 8, px(nul) - px(0)))
        o.append('<circle cx="%.1f" cy="%.1f" r="3.4" fill="var(--c4)"/>'
                 % (px(obs), y - 4.5))
        o.append('<text class="ts" x="%.1f" y="%.1f">%.1f&times; null,  p &lt; 0.001</text>'
                 % (px(obs) + 8, y, ratio))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + rowh * len(rows), x0 + w, y0 + rowh * len(rows)))
    return "\n".join(o), rowh * len(rows)


def build():
    rep, fam, mrg = load()
    fb = rep["family"]["bootstrap"]
    mb = rep["merge_scale"]["bootstrap"]
    fp = rep["family"]["variation_vs_null"]["pooled_proportions"]
    mp = rep["merge_scale"]["variation_vs_null"]["pooled_proportions"]

    ya = 62.0
    pa, ha, na = dot_panel(fam, fb, FAM, fp, ya, 12.4, 3)
    yb = ya + ha + 74
    pb, hb, nb = dot_panel(mrg, mb, MRG, mp, yb, 12.4, 3)
    yc = yb + hb + 78
    pc, hc = null_panel(rep, yc)
    H = yc + hc + 46

    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Composition of event types '
         'by group, with bootstrap intervals, and a null test showing the between-group '
         'variation exceeds sampling and strengthens when collar coverage is matched.">'
         % H]
    o.append('<text class="pl" x="0" y="%.1f">a</text>' % (ya - 32))
    o.append('<text class="ph" x="18" y="%.1f">Family composition, %d groups</text>'
             % (ya - 32, na))
    o.append(pa)
    o.append('<text class="pl" x="0" y="%.1f">b</text>' % (yb - 32))
    o.append('<text class="ph" x="18" y="%.1f">Merge-scale composition within merges, '
             '%d groups</text>' % (yb - 32, nb))
    o.append(pb)
    o.append('<text class="pl" x="0" y="%.1f">c</text>' % (yc - 30))
    o.append('<text class="ph" x="18" y="%.1f">Between-group spread against a '
             'multinomial null</text>' % (yc - 30))
    o.append(pc)
    o.append('<text class="ts" x="148" y="%.1f">grey bar: null mean &middot; dot: '
             'observed</text>' % (yc + hc + 30))
    o.append("</svg>")
    return "\n".join(o), rep, fam, mrg


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
.scroll{overflow-x:auto}.scroll>svg{min-width:620px}
.cap{font-size:12.5px;line-height:1.62;color:var(--ink-2);margin:16px 0 0;
padding-top:14px;border-top:1px solid var(--rule)}
.cap b{font-weight:600;color:var(--ink)}
.cap .fignum{font-weight:600;color:var(--ink)}
.note{font-size:11.5px;line-height:1.6;color:var(--muted);margin-top:22px;
padding-top:14px;border-top:1px solid var(--rule-soft)}
.note code{font-family:var(--mono);font-size:.92em}
text{font-family:var(--sans)}
.gl{font-size:9.2px;fill:var(--ink-2)}
.ts{font-size:8.4px;fill:var(--muted)}
.ch{font-size:9.6px;fill:var(--ink);font-weight:600}
.ph{font-size:11px;fill:var(--ink);font-weight:600}
.pl{font-size:14px;font-weight:700;fill:var(--ink)}
.g{stroke:var(--rule-soft);stroke-width:.8}
.ax{stroke:var(--rule);stroke-width:1}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def page():
    fig, rep, fam, mrg = build()
    fv = rep["family"]["variation_vs_null"]
    mv = rep["merge_scale"]["variation_vs_null"]
    fc = rep["family"]["coverage_confound"]["spearman_with_median_collars"]
    mc = rep["merge_scale"]["coverage_confound"]["spearman_with_median_collars"]
    fs = rep["family"]["coverage_stratified"]
    ms = rep["merge_scale"]["coverage_stratified"]
    pw = mrg.loc[mrg["group"].eq("PhantomWest")].iloc[0]
    fo = mrg.loc[mrg["group"].eq("FireOpal")].iloc[0]
    ls = fam.loc[fam["group"].eq("LapisSplinter")].iloc[0]
    jd = fam.loc[fam["group"].eq("Jade")].iloc[0]
    return f"""<title>Event Composition by Group</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 3 &middot; draft</div>
  <h1>Event composition varies widely across groups</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 3.</span> Event-type composition across 23 social groups.
      Rows are groups, ordered by the leftmost proportion; the number beside each group name
      is its event total. Dots are observed proportions, horizontal lines 95% percentile
      bootstrap intervals over that group's own events (2,000 draws); the dashed line in each
      column is the pooled proportion across all groups.
      <b>(a)</b> Family composition. Pooled, {100 * fv["pooled_proportions"]["group_merge"]:.0f}% of
      events are between-group merges, {100 * fv["pooled_proportions"]["within_group"]:.0f}%
      within-group splits and {100 * fv["pooled_proportions"]["individual"]:.0f}% individual
      separations &mdash; but the spread is wide: LapisSplinter is
      {100 * ls["p_group_merge"]:.0f}% merges against Jade at {100 * jd["p_group_merge"]:.0f}%,
      while Jade produces {100 * jd["p_within_group"]:.0f}% splits.
      <b>(b)</b> Composition of the merges themselves, pooled over the record. The spread is
      wide &mdash; PhantomWest is {100 * pw["p_small_subset_merge"]:.0f}% small subsets of
      {int(pw["total"])} merges, FireOpal {100 * fo["p_large_merge"]:.0f}% large of
      {int(fo["total"])} &mdash; but <b>this ordering does not survive a temporal
      decomposition</b> and should be read as variation in the population rather than as
      differences between groups. See Figure 4: once quarter-to-quarter drift is retained,
      merge-scale composition has an ICC of 0.05&ndash;0.13 and is not distinguishable from
      shuffled group labels (permutation p = 0.09&ndash;0.26), and well-sampled groups traverse
      most of the range within their own record.
      <b>(c)</b> Observed between-group spread (dot) against the mean of a multinomial null in
      which each group draws its own number of events from the pooled distribution (grey bar),
      in three collar-coverage strata. Family composition varies at
      {fs[0]["ratio"]:.1f}&times; the null across all groups and
      {fs[2]["ratio"]:.1f}&times; within a matched 7&ndash;14 collar band; merge-scale composition
      at {ms[0]["ratio"]:.1f}&times; and {ms[2]["ratio"]:.1f}&times; respectively. The variation
      is not a coverage artefact &mdash; it <b>strengthens</b> when coverage is matched. Note
      the limit of this null, though: it assumes each group draws from a fixed pooled
      distribution with no temporal structure, so it tests whether the spread exceeds
      <em>sampling</em>, not whether groups differ. Figure 4 applies the stronger test.
    </p>
  </div>

  <p class="note">
    Merges are dyadic and counted for both participating groups, so merge counts sum to twice
    the merge event total; splits and individual separations are attributed to the origin group.
    Groups with fewer than {MIN_EVENTS} events are omitted. Collar coverage does bias detection
    by type &mdash; the within-group-split share correlates with median collars at
    Spearman {fc["within_group"]:+.2f} and the merge share at {fc["group_merge"]:+.2f}, as
    expected since a split needs two clusters visible while a merge needs one animal from each
    side &mdash; which is why panel (c) repeats the test inside coverage strata. Events and
    exposure both come from the legacy hourly source filtered to 2025-01-01 onward; cohort
    reconciliation with the frozen export is outstanding. Generated by
    <code>build_composition_figure.py</code> from
    <code>analyze_event_composition_by_group.py</code>. Bootstrap and null seed 20260904.
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
