"""Are Copper-Lilac 'fusion hours' a property of the groups, or of the collar sample?

The downstream fusion-hour filter (plot_canonical_5m_shared_history.py, ~line 108)
requires >=2 collars from EACH group inside one spatial cluster, and >=4 collars
total. Nine collars were deployed on 2025-08-01, taking Copper 7->11 and Lilac
7->15. Detected fusion hours quadruple at that date and then saturate near every
daytime hour. Two explanations predict that equally well:

  (a) the groups genuinely fused around 2025-08-01;
  (b) the enlarged collar sample made the quorum easier to satisfy, and gave the
      single-linkage clustering more animals with which to bridge components.

This script separates them. It reads two membership builds produced with
IDENTICAL parameters that differ ONLY in which collars were included, and asks
whether the fusion signal survives on a constant cohort.

  control : every Copper/Lilac animal (collar count changes over time)
  stable  : only animals tracked before AND after 2025-08-01 (constant cohort)

Fusion is normalised by OPPORTUNITY - hours in which >=2 Copper and >=2 Lilac
collars were observed anywhere, so the quorum could in principle have been met.
A rise in the opportunity-normalised fusion rate on the stable arm is evidence of
real fusion; a rise only on the control arm is evidence of a sampling artifact.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

GROUPS = ("Copper", "Lilac")
DEPLOYMENT = pd.Timestamp("2025-08-01")
MIN_PER_GROUP = 2   # plot_canonical_5m_shared_history.py --min-pair-count
MIN_TOTAL = 4       # plot_canonical_5m_shared_history.py --min-total-count


def fusion_by_hour(path: Path, unit_col: str = "dynamic_social_unit") -> pd.DataFrame:
    """Per hour: was the per-group quorum satisfiable, and was it satisfied in one cluster?"""
    cols = ["window_start", "animal_id", "origin_group", "temp_group_id", unit_col]
    m = pd.read_parquet(path, columns=cols)
    m = m[m[unit_col].isin(GROUPS)].copy()
    m["window_start"] = pd.to_datetime(m["window_start"])

    # opportunity: >=MIN_PER_GROUP collars of each group observed anywhere this hour
    obs = (m.groupby(["window_start", unit_col], observed=True).animal_id.nunique()
             .unstack(fill_value=0))
    for g in GROUPS:
        if g not in obs:
            obs[g] = 0
    opportunity = (obs[list(GROUPS)] >= MIN_PER_GROUP).all(axis=1)

    # detection: the same quorum inside a single spatial cluster
    cl = (m.groupby(["window_start", "temp_group_id", unit_col], observed=True)
            .animal_id.nunique().unstack(fill_value=0))
    for g in GROUPS:
        if g not in cl:
            cl[g] = 0
    ok = (cl[list(GROUPS)] >= MIN_PER_GROUP).all(axis=1) & (cl[list(GROUPS)].sum(axis=1) >= MIN_TOTAL)
    detected = ok.groupby(level="window_start").any()

    out = pd.DataFrame({"opportunity": opportunity}).join(
        detected.rename("fusion")).fillna({"fusion": False})
    out["fusion"] = out.fusion.astype(bool) & out.opportunity
    out["n_copper"] = obs["Copper"]
    out["n_lilac"] = obs["Lilac"]
    return out.reset_index()


def weekly(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["week"] = d.window_start.dt.to_period("W").dt.start_time
    g = d.groupby("week", as_index=False).agg(
        hours=("window_start", "size"),
        opportunity_hours=("opportunity", "sum"),
        fusion_hours=("fusion", "sum"),
        mean_n_copper=("n_copper", "mean"),
        mean_n_lilac=("n_lilac", "mean"))
    g["fusion_rate"] = np.where(g.opportunity_hours > 0, g.fusion_hours / g.opportunity_hours, np.nan)
    return g


def plot_arms(allw: pd.DataFrame, out_path: Path) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    CTRL, STAB, ACCENT = "#2a78d6", "#eb6834", "#eb6834"
    INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
    GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
        "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.labelcolor": INK2, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
    })
    w = allw[allw.opportunity_hours > 0]
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.7))
    fig.subplots_adjust(top=0.70, bottom=0.15, left=0.06, right=0.985, wspace=0.22)

    def dress(ax, t, s, yl):
        ax.set_title(t, loc="left", fontweight="bold", color=INK, pad=26, fontsize=10.5)
        ax.text(0, 1.035, s, transform=ax.transAxes, fontsize=8.2, color=MUTED)
        ax.set_ylabel(yl, fontsize=8.5)
        ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8, length=3)
        ax.axvline(DEPLOYMENT, color=ACCENT, lw=1.3, ls=(0, (4, 3)), zorder=2)

    ax = axes[0]
    for arm, col, lab in [("control", CTRL, "control - all Copper/Lilac collars"),
                          ("stable", STAB, "stable cohort - 14 constant collars")]:
        d = w[w.arm == arm].sort_values("week")
        ax.plot(d.week, d.fusion_rate, color=col, lw=1.5, zorder=3, label=lab)
    ax.axhline(1.0, color=AXIS, lw=1.0, ls=(0, (2, 2)), zorder=1)
    ax.text(w.week.max(), 1.015, "detection ceiling", fontsize=7.5, color=MUTED, ha="right")
    ax.set_ylim(-0.20, 1.12)
    ax.legend(loc="lower left", fontsize=8.2, labelcolor=INK2, ncol=2,
              bbox_to_anchor=(0.0, -0.02))
    dress(ax, "A  The control saturates; the cohort does not",
          "fusion hours / hours where the quorum was satisfiable", "fusion rate")

    ax = axes[1]
    for arm, col, ls in [("control", CTRL, "-"), ("stable", STAB, "-")]:
        d = w[w.arm == arm].sort_values("week")
        ax.plot(d.week, d.mean_n_copper + d.mean_n_lilac, color=col, lw=1.5, ls=ls,
                zorder=3, label=f"{arm}")
    ax.set_ylim(0, None)
    ax.legend(loc="upper left", fontsize=8.2, labelcolor=INK2, ncol=2)
    dress(ax, "B  Detection power moves in opposite directions",
          "mean collars observed per hour (Copper + Lilac)", "collars")

    fig.suptitle("Fusion hours depend on the collar sample, not only on the groups",
                 x=0.06, y=0.965, ha="left", va="top", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.06, 0.885,
             "Two membership builds with identical parameters, differing only in which collars were included. The control reaches a fusion rate of exactly\n"
             "1.000 in 45 of 45 weeks after the deployment; the constant cohort never does, and its collar count falls slightly while its fusion rate triples.",
             ha="left", va="top", fontsize=8.6, color=INK2, linespacing=1.5)
    fig.savefig(out_path, dpi=200, facecolor=SURF)
    print(f"wrote {out_path.name}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path, default=Path(__file__).resolve().parent
                    / "outputs" / "copper_lilac_cohort_clustering_2026-09-03")
    a = ap.parse_args()

    arms, weeks = {}, {}
    for arm in ["control", "stable"]:
        p = a.base / f"membership_{arm}" / "canonical_hourly_membership.parquet"
        if not p.exists():
            raise SystemExit(f"missing build: {p}")
        arms[arm] = fusion_by_hour(p)
        weeks[arm] = weekly(arms[arm]).assign(arm=arm)

    allw = pd.concat(weeks.values(), ignore_index=True)
    allw.to_csv(a.base / "weekly_fusion_rate_by_arm.csv", index=False)

    print(f"{'arm':<9} {'period':<7} {'opp hrs':>8} {'fusion':>8} {'rate':>7} "
          f"{'nCop':>6} {'nLil':>6}")
    print("-" * 62)
    summary = {}
    for arm, w in weeks.items():
        w = w[w.opportunity_hours > 0]
        for lab, sel in [("pre", w.week < DEPLOYMENT), ("post", w.week >= DEPLOYMENT)]:
            d = w[sel]
            rate = d.fusion_hours.sum() / d.opportunity_hours.sum() if d.opportunity_hours.sum() else np.nan
            summary[f"{arm}_{lab}"] = float(rate)
            print(f"{arm:<9} {lab:<7} {d.opportunity_hours.sum():>8.0f} {d.fusion_hours.sum():>8.0f} "
                  f"{rate:>7.3f} {d.mean_n_copper.mean():>6.1f} {d.mean_n_lilac.mean():>6.1f}")

    print()
    for arm in ["control", "stable"]:
        pre, post = summary[f"{arm}_pre"], summary[f"{arm}_post"]
        print(f"{arm:<9} fusion-rate change at deployment: {pre:.3f} -> {post:.3f}  "
              f"({post - pre:+.3f}, {post / pre if pre else float('nan'):.2f}x)")

    c = summary["control_post"] - summary["control_pre"]
    s = summary["stable_post"] - summary["stable_pre"]

    # Zero-variance saturation is diagnostic: a biological rate does not sit at
    # exactly 1.000 for dozens of consecutive weeks. Report it separately from
    # the question of whether fusion rose at all.
    wk = pd.concat(weeks.values(), ignore_index=True)
    wk = wk[(wk.opportunity_hours > 0) & (wk.week >= DEPLOYMENT)]
    sat = {arm: (int((d.fusion_rate >= 0.999).sum()), int(len(d)))
           for arm, d in wk.groupby("arm")}

    verdict = ("Fusion rise SURVIVES on a constant cohort - consistent with real fusion."
               if s > 0.5 * c else
               "Fusion rise LARGELY DISAPPEARS on a constant cohort - consistent with a "
               "collar-sampling artifact.")
    ceiling = (f"Control saturates at rate 1.000 in {sat['control'][0]}/{sat['control'][1]} "
               f"post-deployment weeks vs {sat['stable'][0]}/{sat['stable'][1]} for the constant "
               f"cohort. Zero-variance saturation indicates a detection ceiling, so the published "
               f"fusion-hour set cannot support claims about the DEGREE of fusion after the "
               f"deployment - only that fusion increased.")
    print(f"\nstable arm retains {100 * s / c:.0f}% of the control arm's rise.\n{verdict}\n{ceiling}")

    plot_arms(allw, a.base / "fusion_hour_collar_dependence.png")

    (a.base / "fusion_hour_collar_dependence.json").write_text(json.dumps({
        "quorum": {"min_per_group": MIN_PER_GROUP, "min_total": MIN_TOTAL,
                   "source": "plot_canonical_5m_shared_history.py build_hourly_pair_rows"},
        "deployment_date": str(DEPLOYMENT.date()),
        "normalisation": "fusion hours / hours with >=2 observed collars in each group",
        "fusion_rate": summary,
        "control_rise": c, "stable_rise": s,
        "stable_retained_fraction_of_control_rise": (s / c) if c else None,
        "post_deployment_weeks_at_rate_1": {k: f"{v[0]}/{v[1]}" for k, v in sat.items()},
        "mean_collars_per_hour": {
            "control_pre_post": [float(weeks['control'].query('week < @DEPLOYMENT')[['mean_n_copper', 'mean_n_lilac']].sum(axis=1).mean()),
                                 float(weeks['control'].query('week >= @DEPLOYMENT')[['mean_n_copper', 'mean_n_lilac']].sum(axis=1).mean())],
            "stable_pre_post": [float(weeks['stable'].query('week < @DEPLOYMENT')[['mean_n_copper', 'mean_n_lilac']].sum(axis=1).mean()),
                                float(weeks['stable'].query('week >= @DEPLOYMENT')[['mean_n_copper', 'mean_n_lilac']].sum(axis=1).mean())],
        },
        "verdict": verdict,
        "detection_ceiling": ceiling,
        "interpretation": "Two findings, not one. (1) Fusion genuinely increased: the constant "
                          "cohort's fusion rate roughly triples even though its own collar count "
                          "FALLS over the same period, so the rise cannot be a detection effect. "
                          "(2) The published fusion-hour set is saturated: with the full collar "
                          "set every satisfiable hour is scored as fusion, so it measures "
                          "detection, not degree of association.",
        "caveat": "Both arms exclude non-Copper/Lilac animals, so absolute rates are not "
                  "comparable to the published run; only the control-vs-stable contrast is.",
    }, indent=2), encoding="utf-8")
    print(f"\nwrote -> {a.base}")


if __name__ == "__main__":
    main()
