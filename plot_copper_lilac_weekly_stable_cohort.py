"""Figures for the weekly, constant-cohort Copper-Lilac re-analysis.

Reads the dated output directory written by
`analyze_copper_lilac_weekly_stable_cohort.py` and renders:

  copper_lilac_weekly_stable_main.png    3-panel candidate main figure
  copper_lilac_weekly_components.png     component decomposition

Colours follow a CVD-validated categorical palette (blue/orange adjacent slots)
and a single-hue blue ordinal ramp for radius.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
DEPLOYMENT = pd.Timestamp("2025-08-01")
PANEL_RADII = [2.0, 5.0, 20.0, 100.0]

COPPER, LILAC = "#2a78d6", "#eb6834"
RAMP = ["#86b6ef", "#3987e5", "#256abf", "#0d366b"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
ACCENT = "#eb6834"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def dress(ax, title, sub, ylab, pad=26):
    ax.set_title(title, loc="left", fontweight="bold", color=INK, pad=pad, fontsize=10.5)
    ax.text(0, 1.035, sub, transform=ax.transAxes, fontsize=8.2, color=MUTED)
    ax.set_ylabel(ylab, fontsize=8.5)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, length=3)


def main_figure(weekly, tests, out_path):
    ok = weekly[weekly.complete_cell & weekly.integration_ratio.notna()]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.1))
    fig.subplots_adjust(top=0.755, bottom=0.145, left=0.052, right=0.985, wspace=0.28)

    # A - scale gradient before vs after, with bootstrap CIs
    ax = axes[0]
    era = {"before 1 Aug 2025": ok[ok.week < DEPLOYMENT], "from 1 Aug 2025": ok[ok.week >= DEPLOYMENT]}
    for (lab, d), col, mk in zip(era.items(), [RAMP[0], RAMP[3]], ["o", "s"]):
        g = d.groupby("radius_m").agg(r=("integration_ratio", "mean"),
                                      lo=("ci_95_low", "mean"), hi=("ci_95_high", "mean")).reset_index()
        ax.fill_between(g.radius_m, g.lo, g.hi, color=col, alpha=0.16, lw=0, zorder=2)
        ax.plot(g.radius_m, g.r, color=col, lw=2.2, marker=mk, ms=6, mec=SURF, mew=1.0,
                zorder=3, label=lab)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 5, 10, 20, 50, 100, 200, 400])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.axhline(1.0, color=AXIS, lw=1.0, ls=(0, (2, 2)), zorder=1)
    ax.set_ylim(0, 1.12)
    ax.set_xlabel("proximity radius (m, log scale)", fontsize=8.5)
    ax.legend(loc="upper left", fontsize=8.2, labelcolor=INK2)
    dress(ax, "A  Mixing is scale-dependent", "integration ratio by radius, mean of weekly values",
          "integration ratio")

    # B - weekly series at 5 m with CI ribbon and trend
    ax = axes[1]
    d = ok[ok.radius_m == 5.0].sort_values("week")
    ax.axvline(DEPLOYMENT, color=ACCENT, lw=1.4, ls=(0, (4, 3)), zorder=2)
    ax.fill_between(d.week, d.ci_95_low, d.ci_95_high, color=RAMP[1], alpha=0.20, lw=0, zorder=2)
    ax.plot(d.week, d.integration_ratio, color=RAMP[1], lw=1.3, marker="o", ms=3.4,
            mec=SURF, mew=0.6, zorder=3)
    t = (d.week - d.week.min()).dt.days / 365.25
    b1, b0 = np.polyfit(t, d.integration_ratio, 1)
    xs = np.linspace(t.min(), t.max(), 100)
    ax.plot(d.week.min() + pd.to_timedelta(xs * 365.25, "D"), b0 + b1 * xs,
            color=RAMP[3], lw=2.4, zorder=4)
    row = tests[(tests.radius_m == 5.0) & (tests.response == "integration_ratio")].iloc[0]
    ax.set_ylim(0, 1.12)
    ax.text(DEPLOYMENT, 1.07, " collar deployment\n 1 Aug 2025", color=ACCENT, fontsize=8, va="top")
    ax.text(0.03, 0.955,
            f"trend {row.trend_per_year:+.2f}/yr (p<0.001)\n"
            f"step at deployment  p = {row.step_p:.2f}",
            transform=ax.transAxes, fontsize=8.4, color=INK2, va="top")
    dress(ax, "B  A smooth trend, not a step",
          "weekly integration ratio at 5 m, 95% bootstrap CI", "integration ratio")

    # C - identification: trend vs step by radius
    ax = axes[2]
    t = tests[tests.response == "integration_ratio"].sort_values("radius_m")
    y = np.arange(len(t))
    ax.axvline(0, color=AXIS, lw=1.2, zorder=1)
    ax.errorbar(t.trend_per_year, y + 0.17,
                xerr=[t.trend_per_year - t.trend_ci_low, t.trend_ci_high - t.trend_per_year],
                fmt="o", color=RAMP[3], ms=6, lw=1.8, capsize=2.5, zorder=3, label="trend per year")
    ax.errorbar(t.step_at_deployment, y - 0.17,
                xerr=[t.step_at_deployment - t.step_ci_low, t.step_ci_high - t.step_at_deployment],
                fmt="o", color=ACCENT, ms=6, lw=1.8, capsize=2.5, zorder=3, label="step at deployment")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r:.0f} m" for r in t.radius_m])
    ax.set_ylim(-0.7, len(t) + 0.15)
    ax.invert_yaxis()
    ax.set_xlabel("change in integration ratio", fontsize=8.5)
    ax.legend(loc="lower center", fontsize=8.2, labelcolor=INK2, ncol=2, bbox_to_anchor=(0.5, -0.02))
    dress(ax, "C  Every step estimate straddles zero",
          f"OLS, HAC errors; {int(row.n_pre)} pre / {int(row.n_post)} post weeks", "")
    ax.grid(axis="y", lw=0)

    fig.suptitle("Copper-Lilac integration among consistently tracked animals",
                 x=0.052, y=0.978, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.052, 0.905,
             "Stable cohort: 14 animals (7 Copper, 7 Lilac) collared before and after 1 Aug 2025, when nine new collars were "
             "deployed on a single day.\nWeekly resolution identifies the deployment separately from the trend; monthly bins cannot.",
             ha="left", va="top", fontsize=8.8, color=INK2, linespacing=1.5)
    fig.savefig(out_path, dpi=200, facecolor=SURF)
    print("wrote", out_path.name)


def component_figure(weekly, out_path):
    ok = weekly[weekly.complete_cell]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.7))
    fig.subplots_adjust(top=0.72, bottom=0.155, left=0.052, right=0.985, wspace=0.27)

    ax = axes[0]
    ax.axvline(DEPLOYMENT, color=ACCENT, lw=1.3, ls=(0, (4, 3)), zorder=2)
    for c, R in zip(RAMP, PANEL_RADII):
        d = ok[ok.radius_m == R].sort_values("week")
        ax.plot(d.week, d.cross_balanced.replace(0, np.nan), color=c, lw=1.4, zorder=3)
        ax.annotate(f"{R:.0f} m", (d.week.iloc[-1], d.cross_balanced.iloc[-1]), color=c,
                    fontsize=8, fontweight="bold", xytext=(5, -2), textcoords="offset points")
    ax.set_yscale("log")
    dress(ax, "A  Intergroup proximity", "weekly cross-origin contact rate (log)", "contact rate")

    ax = axes[1]
    d = ok[ok.radius_m == 5.0].sort_values("week")
    ax.axvline(DEPLOYMENT, color=ACCENT, lw=1.3, ls=(0, (4, 3)), zorder=2)
    for col, c, lab in [("within_Copper", COPPER, "within Copper"),
                        ("within_Lilac", LILAC, "within Lilac")]:
        ax.plot(d.week, d[col], color=c, lw=1.5, zorder=3, label=lab)
    ax.plot(d.week, d.cross_balanced, color=INK2, lw=1.4, ls=(0, (3, 2)), zorder=3, label="cross-origin")
    ax.legend(loc="upper right", fontsize=8, labelcolor=INK2)
    dress(ax, "B  Inner-group cohesion holds", "within- and cross-origin rates at 5 m", "contact rate")

    ax = axes[2]
    ax.axvline(DEPLOYMENT, color=ACCENT, lw=1.3, ls=(0, (4, 3)), zorder=2)
    for col, c, lab in [("n_within_Copper", COPPER, "Copper"), ("n_within_Lilac", LILAC, "Lilac")]:
        ax.plot(d.week, d[col], color=c, lw=1.5, zorder=3, label=lab)
    ax.axhline(2, color=MUTED, lw=1.0, ls=(0, (2, 2)), zorder=2)
    ax.legend(loc="upper right", fontsize=8, labelcolor=INK2, ncol=2)
    ax.set_ylim(0, 8)
    dress(ax, "C  Support thins late", "animals contributing a within-group estimate", "animals")

    fig.suptitle("Components of the integration ratio, stable cohort",
                 x=0.052, y=0.972, ha="left", va="top", fontsize=13,
                 fontweight="bold", color=INK)
    fig.savefig(out_path, dpi=200, facecolor=SURF)
    print("wrote", out_path.name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", type=Path, required=True)
    a = ap.parse_args()
    weekly = pd.read_csv(a.output_dir / "weekly_integration_stable_cohort.csv", parse_dates=["week"])
    tests = pd.read_csv(a.output_dir / "trend_vs_deployment_step_tests.csv")
    main_figure(weekly, tests, a.output_dir / "copper_lilac_weekly_stable_main.png")
    component_figure(weekly, a.output_dir / "copper_lilac_weekly_components.png")


if __name__ == "__main__":
    main()
