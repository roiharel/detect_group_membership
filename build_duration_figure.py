"""Manuscript figure: duration distributions of every social structure, with CIs.

Three panels, one point each, minimal on-figure text:

  (a) Survival curves for six social structures. Every structure is brief at the
      median and heavy-tailed.
  (b) The same 1,705 encounters measured three ways -- structural span, supported
      exposure, active contact. Duration is not one quantity.
  (c) Excursion length by the depth it reached. Depth tracks duration.

CONFIDENCE INTERVALS. Ribbons are 95% percentile intervals from a cluster bootstrap,
1,000 draws, resampling the natural unit of replication for each structure rather than
the events themselves:

  merges                       -> dyad (`pair`)
  within-group split           -> origin group
  single-animal separation     -> animal
  isolated                     -> animal
  encounter durations          -> dyad (`pair_key`)
  excursions                   -> animal

Resampling events would understate the interval wherever one dyad or one animal
contributes many events, which is the case throughout. Curves are also distinguished by
dash pattern so the figure survives greyscale reproduction.

Output: docs/framing_2026-09-04/duration_figure.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("outputs/general_structure_2026_09")
LEGACY = Path("outputs/canonical_group_merge_scale_log_scatter")
OUT = Path("docs/framing_2026-09-04/duration_figure.html")

N_BOOT = 1000
SEED = 20260904


# ------------------------------------------------------------------ statistics
def survival_at(v: np.ndarray, xs: np.ndarray) -> np.ndarray:
    """P(T > x), vectorised via the sorted sample."""
    s = np.sort(v)
    return 1.0 - np.searchsorted(s, xs, side="right") / len(s)


def cluster_bootstrap(vals: np.ndarray, clusters: np.ndarray, xs: np.ndarray,
                      rng: np.random.Generator, n_boot: int = N_BOOT
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Point estimate and 95% percentile band, resampling whole clusters."""
    vals = np.asarray(vals, dtype=float)
    clusters = np.asarray(clusters)
    keys, inv = np.unique(clusters, return_inverse=True)
    buckets = [vals[inv == i] for i in range(len(keys))]
    point = survival_at(vals, xs)
    draws = np.empty((n_boot, len(xs)), dtype=float)
    n_k = len(keys)
    for b in range(n_boot):
        pick = rng.integers(0, n_k, n_k)
        sample = np.concatenate([buckets[i] for i in pick])
        draws[b] = survival_at(sample, xs)
    lo = np.percentile(draws, 2.5, axis=0)
    hi = np.percentile(draws, 97.5, axis=0)
    return point, lo, hi, n_k


# ------------------------------------------------------------------ data
def load(rng) -> dict:
    ev = pd.read_csv(LEGACY / "canonical_event_size_duration_all_events.csv")
    iso = pd.read_csv(LEGACY / "canonical_isolated_events.csv")
    st = pd.read_csv(BASE / "phase2_two_stage_events/stage1_events_with_stage2_mixing.csv")
    ex = pd.read_csv(BASE / "phase4b_individual_axis/excursions_dominant_gap0.csv")

    def clean(frame, valcol, clustercol):
        f = frame[[valcol, clustercol]].copy()
        f[valcol] = pd.to_numeric(f[valcol], errors="coerce")
        f = f.dropna()
        f = f[f[valcol] > 0]
        return f[valcol].to_numpy(), f[clustercol].astype(str).to_numpy()

    panel_a = [
        ("Large merge", "solid", "var(--c1)",
         *clean(ev[ev.event_type.eq("large_merge")], "duration_hours", "pair")),
        ("Medium partial merge", "6 3", "var(--c2)",
         *clean(ev[ev.event_type.eq("medium_partial_merge")], "duration_hours", "pair")),
        ("Small subset merge", "2 3", "var(--c3)",
         *clean(ev[ev.event_type.eq("small_subset_merge")], "duration_hours", "pair")),
        ("Within-group split", "solid", "var(--c4)",
         *clean(ev[ev.event_type.eq("within_group_split")], "duration_hours",
                "origin_group")),
        ("Single-animal separation", "6 3", "var(--c5)",
         *clean(ev[ev.event_type.eq("single_animal_separation")], "duration_hours",
                "animal_id")),
        ("Isolated", "2 3", "var(--c6)",
         *clean(iso, "duration_hours", "animal_id")),
    ]
    panel_b = [
        ("Structural span", "solid", "var(--c1)",
         *clean(st, "structural_span_hours", "pair_key")),
        ("Supported exposure", "6 3", "var(--c2)",
         *clean(st, "5m_supported_exposure_hours", "pair_key")),
        ("Active contact", "2 3", "var(--c5)",
         *clean(st, "5m_active_contact_hours", "pair_key")),
    ]
    ex_alone = ex[ex.depth.eq("alone_only")]
    ex_join = ex[ex.depth.ne("alone_only")]
    panel_c = [
        ("Alone only", "2 3", "var(--c6)",
         *clean(ex_alone, "away_nights", "animal_id")),
        ("Reached another unit", "solid", "var(--c5)",
         *clean(ex_join, "away_nights", "animal_id")),
    ]
    return {"a": panel_a, "b": panel_b, "c": panel_c}


# ------------------------------------------------------------------ drawing
def draw_panel(series, box, xlim, xticks, rng, legend_xy=(None, 0.06),
               ylab="P(duration > x)", xlab="", refs=()):
    """One log-x survival panel with bootstrap ribbons. Returns svg, stats."""
    x0, y0, w, h = box
    lo, hi = xlim
    xs = np.geomspace(lo, hi, 90)

    def px(v):
        v = min(max(v, lo), hi)
        return x0 + (math.log10(v) - math.log10(lo)) / (
            math.log10(hi) - math.log10(lo)) * w

    def py(p):
        return y0 + h - p * h

    o, stats = [], {}
    for tick, lab in xticks:
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (px(tick), y0, px(tick), y0 + h))
        o.append('<text class="t" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (px(tick), y0 + h + 13, lab))
    for p in (0, 0.5, 1.0):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0, py(p), x0 + w, py(p)))
        o.append('<text class="t" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (x0 - 6, py(p) + 3.2, ("0", "0.5", "1")[(0, 0.5, 1.0).index(p)]))
    for rv, rlab in refs:
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
                 'stroke-width=".8" stroke-dasharray="2 3" opacity=".5"/>'
                 % (px(rv), y0 - 4, px(rv), y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (px(rv), y0 - 7, rlab))

    # ribbons first, then lines on top
    computed = []
    for label, dash, col, vals, clus in series:
        pt, blo, bhi, nk = cluster_bootstrap(vals, clus, xs, rng)
        computed.append((label, dash, col, pt, blo, bhi, vals, nk))
        fwd = " ".join("%.1f,%.1f" % (px(x), py(p)) for x, p in zip(xs, bhi))
        rev = " ".join("%.1f,%.1f" % (px(x), py(p)) for x, p in
                       zip(xs[::-1], blo[::-1]))
        o.append('<polygon points="%s %s" fill="%s" opacity=".16"/>' % (fwd, rev, col))
    for label, dash, col, pt, blo, bhi, vals, nk in computed:
        pts = " ".join("%.1f,%.1f" % (px(x), py(p)) for x, p in zip(xs, pt))
        da = "" if dash == "solid" else ' stroke-dasharray="%s"' % dash
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.7"%s '
                 'stroke-linejoin="round"/>' % (pts, col, da))
        med = float(np.median(vals))
        stats[label] = {"n": int(len(vals)), "clusters": int(nk), "median": med,
                        "p90": float(np.quantile(vals, 0.9)),
                        "max": float(vals.max())}

    # Legend is placed by hand into whichever corner the curves leave empty: the
    # upper right where they have all decayed, or the lower left where they are all
    # still high. `legend_xy` is (data x, fraction of panel height from the top).
    leg_x, leg_yf = legend_xy
    lx = px(leg_x) if leg_x is not None else x0 + 8
    ly = y0 + leg_yf * h + 7
    for label, dash, col, pt, blo, bhi, vals, nk in computed:
        da = "" if dash == "solid" else ' stroke-dasharray="%s"' % dash
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="1.7"%s/>' % (lx, ly - 3, lx + 17, ly - 3, col, da))
        o.append('<text class="t" x="%.1f" y="%.1f">%s</text>'
                 % (lx + 22, ly, label))
        ly += 11

    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + h, x0 + w, y0 + h))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0, x0, y0 + h))
    if xlab:
        o.append('<text class="t" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (x0 + w / 2, y0 + h + 30, xlab))
    o.append('<text class="t" transform="translate(%.1f %.1f) rotate(-90)" '
             'text-anchor="middle">%s</text>' % (x0 - 24, y0 + h / 2, ylab))
    return "\n".join(o), stats


def build() -> tuple[str, dict]:
    rng = np.random.default_rng(SEED)
    data = load(rng)

    H = 560.0
    o = ['<svg viewBox="0 0 700 %.0f" role="img" aria-label="Three-panel figure of '
         'duration distributions with bootstrap confidence bands. Panel a: six social '
         'structures, all brief at the median with long tails. Panel b: the same '
         'encounters measured as structural span, supported exposure and active contact, '
         'separated by about fifty-fold. Panel c: excursions that only reach being alone '
         'are shorter than those reaching another group.">' % H]

    svg_a, st_a = draw_panel(
        data["a"], (60.0, 34.0, 596.0, 196.0), (1.0, 8000.0),
        [(1, "1"), (10, "10"), (100, "100"), (1000, "1,000")],
        rng, legend_xy=(105, 0.05), xlab="duration (h)",
        refs=[(24, "1 d"), (168, "1 wk"), (730, "1 mo")])
    o.append('<text class="pl" x="18" y="30">a</text>')
    o.append(svg_a)

    svg_b, st_b = draw_panel(
        data["b"], (60.0, 330.0, 260.0, 158.0), (0.03, 8000.0),
        [(0.03, "2 min"), (1, "1 h"), (100, "100"), (8000, "8,000")],
        rng, legend_xy=(0.045, 0.60), xlab="duration (h)")
    o.append('<text class="pl" x="18" y="326">b</text>')
    o.append(svg_b)

    svg_c, st_c = draw_panel(
        data["c"], (420.0, 330.0, 236.0, 158.0), (1.0, 520.0),
        [(1, "1"), (7, "7"), (30, "30"), (365, "365")],
        rng, legend_xy=(11, 0.05), ylab="P(nights > x)", xlab="away-nights",
        refs=[(7, "settle")])
    o.append('<text class="pl" x="378" y="326">c</text>')
    o.append(svg_c)
    o.append("</svg>")
    return "\n".join(o), {"a": st_a, "b": st_b, "c": st_c}


CSS = """
:root{--ground:#f4f5f3;--paper:#fff;--ink:#15191a;--ink-2:#3a4244;--muted:#6b7476;
--rule:#d5dad7;--rule-soft:#e6eae7;--accent:#0f5f57;
--c1:#1f6f8b;--c2:#4a92aa;--c3:#84b6c6;--c4:#8a5a1f;--c5:#6b4a7a;--c6:#0f5f57;
--sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;
--serif:"Newsreader",Georgia,serif;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0f1213;
--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;--muted:#8f9998;--rule:#2c3435;
--rule-soft:#232a2b;--accent:#5fc0b0;--c1:#74b6d0;--c2:#4f93ad;--c3:#356a7f;
--c4:#d3a061;--c5:#b494c4;--c6:#5fc0b0;}}
:root[data-theme="dark"]{--ground:#0f1213;--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;
--muted:#8f9998;--rule:#2c3435;--rule-soft:#232a2b;--accent:#5fc0b0;--c1:#74b6d0;
--c2:#4f93ad;--c3:#356a7f;--c4:#d3a061;--c5:#b494c4;--c6:#5fc0b0;}
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
.scroll{overflow-x:auto}.scroll>svg{min-width:600px}
.cap{font-size:12.5px;line-height:1.62;color:var(--ink-2);margin:16px 0 0;
padding-top:14px;border-top:1px solid var(--rule);max-width:none}
.cap b{font-weight:600;color:var(--ink)}
.cap .fignum{font-weight:600;color:var(--ink)}
.note{font-size:11.5px;line-height:1.6;color:var(--muted);margin-top:22px;
padding-top:14px;border-top:1px solid var(--rule-soft)}
.note code{font-family:var(--mono);font-size:.92em}
text{font-family:var(--sans)}
.t{font-size:9.5px;fill:var(--ink-2)}
.ts{font-size:8.5px;fill:var(--muted)}
.pl{font-size:14px;font-weight:700;fill:var(--ink);font-family:var(--sans)}
.g{stroke:var(--rule-soft);stroke-width:.8}
.ax{stroke:var(--rule);stroke-width:1}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def page() -> str:
    fig, st = build()
    a, b, c = st["a"], st["b"], st["c"]
    ratio = b["Structural span"]["median"] / b["Active contact"]["median"]
    return f"""<title>Duration of Social Structures</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 2 &middot; draft</div>
  <h1>Every social structure is brief, and duration is not one quantity</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 2.</span> Duration distributions of social structures in a
      350-animal, 26-group baboon cohort. Curves are complementary cumulative distributions,
      P(duration &gt; x); shaded bands are 95% percentile intervals from a cluster bootstrap
      (1,000 draws) resampling the natural unit of replication &mdash; dyads for merges, groups for
      within-group splits, animals for individual structures.
      <b>(a)</b> Six structures, log duration axis. Every structure has a median at or below
      11 h and a tail beyond one month: large merge
      median {a["Large merge"]["median"]:.0f} h (max {a["Large merge"]["max"]:.0f} h,
      {a["Large merge"]["clusters"]} dyads), within-group split
      {a["Within-group split"]["median"]:.0f} h (max {a["Within-group split"]["max"]:.0f} h,
      {a["Within-group split"]["clusters"]} groups), single-animal separation
      {a["Single-animal separation"]["median"]:.0f} h (max {a["Single-animal separation"]["max"]:.0f} h,
      {a["Single-animal separation"]["clusters"]} animals). Small subset merges are the least
      frequent but the most persistent, crossing above large merges beyond roughly 30 h.
      <b>(b)</b> The same {b["Structural span"]["n"]:,} between-group encounters measured three
      ways: structural span (median {b["Structural span"]["median"]:.0f} h,
      n&nbsp;=&nbsp;{b["Structural span"]["n"]:,}), fine-scale supported exposure
      ({b["Supported exposure"]["median"]:.2f} h, n&nbsp;=&nbsp;{b["Supported exposure"]["n"]:,}) and
      active cross-group contact ({b["Active contact"]["median"]:.2f} h,
      n&nbsp;=&nbsp;{b["Active contact"]["n"]:,}) &mdash; a {ratio:.0f}-fold separation between the
      outermost two, on identical events. The falling sample reflects fine-scale coverage, not
      shorter encounters.
      <b>(c)</b> Individual excursions by the depth reached. Excursions that only reach being
      alone (median {c["Alone only"]["median"]:.0f} night, n&nbsp;=&nbsp;{c["Alone only"]["n"]}) are
      shorter than those reaching another social unit
      ({c["Reached another unit"]["median"]:.0f} nights, n&nbsp;=&nbsp;{c["Reached another unit"]["n"]});
      the curves separate immediately and do not re-cross, and only the latter carries appreciable
      mass beyond the 7-night threshold that defines settlement.
    </p>
  </div>

  <p class="note">
    Panel (a) derives from the legacy hourly inventory filtered to 2025-01-01 onward; panels (b)
    and (c) from the frozen 1,924,104-row export (2024-03-01 to 2026-07-22). Cohort reconciliation
    is outstanding for (a). Curves are distinguished by dash pattern as well as hue so the figure
    reproduces in greyscale. Generated by <code>build_duration_figure.py</code>; sources
    <code>canonical_event_size_duration_all_events.csv</code>,
    <code>canonical_isolated_events.csv</code>,
    <code>stage1_events_with_stage2_mixing.csv</code>,
    <code>excursions_dominant_gap0.csv</code>. Bootstrap seed {SEED}.
  </p>
</div>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    html = page()
    args.output.write_text(html, encoding="utf-8")
    print("wrote %s (%d bytes)" % (args.output, args.output.stat().st_size))


if __name__ == "__main__":
    main()
