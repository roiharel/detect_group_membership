"""Figure: the fusion-hour detection ceiling is specific to Copper-Lilac."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
DEPLOYMENT = pd.Timestamp("2025-08-01")
FOCUS, OTHER = "#eb6834", "#2a78d6"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.labelcolor": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})


def dress(ax, t, s, yl):
    ax.set_title(t, loc="left", fontweight="bold", color=INK, pad=26, fontsize=10.5)
    ax.text(0, 1.035, s, transform=ax.transAxes, fontsize=8.2, color=MUTED)
    ax.set_ylabel(yl, fontsize=8.5)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, length=3)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=PROJECT / "outputs" / "fusion_quorum_all_dyads_2026-09-03")
    a = ap.parse_args()

    d = pd.read_csv(a.dir / "fusion_rate_by_dyad.csv")
    m = pd.read_csv(a.dir / "monthly_fusion_rate_by_dyad.csv", parse_dates=["month"])
    top = d[d.fusion_hours > 0].sort_values("fusion_rate", ascending=False).head(14)

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0))
    fig.subplots_adjust(top=0.72, bottom=0.13, left=0.155, right=0.985, wspace=0.42)

    ax = axes[0]
    y = np.arange(len(top))[::-1]
    cols = [FOCUS if n == "Copper - Lilac" else OTHER for n in top.dyad]
    ax.barh(y, top.fusion_rate, color=cols, height=0.72, zorder=3, linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels(top.dyad, fontsize=8)
    for yi, (rate, sm, sc) in enumerate(zip(top.fusion_rate, top.months_saturated,
                                            top.months_scored)):
        lbl = f"{rate:.2f}" + (f"   {int(sm)}/{int(sc)} mo at ceiling" if sm else "")
        ax.text(rate + 0.012, y[yi], lbl, va="center", fontsize=7.8,
                color=INK if sm else MUTED, fontweight="bold" if sm else "normal")
    ax.set_xlim(0, 0.92)
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.set_xlabel("fusion hours / opportunity hours", fontsize=8.5)
    dress(ax, "A  One dyad stands apart", "all dyads with any detected fusion, ranked", "")

    ax = axes[1]
    names = list(top.dyad.head(6))
    for n in names:
        s = m[(m.dyad == n) & (m.opportunity_hours >= 20)].sort_values("month")
        focus = n == "Copper - Lilac"
        ax.plot(s.month, s.fusion_rate, color=FOCUS if focus else OTHER,
                lw=2.2 if focus else 1.1, alpha=1.0 if focus else 0.45,
                zorder=4 if focus else 3, label=n if focus else None)
    ax.axhline(1.0, color=AXIS, lw=1.0, ls=(0, (2, 2)), zorder=1)
    ax.axvline(DEPLOYMENT, color=MUTED, lw=1.2, ls=(0, (4, 3)), zorder=2)
    ax.text(DEPLOYMENT, 1.10, " collar deployment", color=MUTED, fontsize=7.8, va="top")
    ax.text(m.month.min(), 1.02, "detection ceiling", fontsize=7.5, color=MUTED)
    ax.set_ylim(-0.03, 1.16)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.legend(loc="center left", fontsize=8.2, labelcolor=INK2)
    ax.text(0.02, 0.40, "other five dyads", transform=ax.transAxes,
            fontsize=8.2, color=OTHER)
    dress(ax, "B  Only Copper-Lilac reaches the ceiling",
          "monthly fusion rate, six highest-rate dyads", "fusion rate")

    fig.suptitle("The fusion-hour ceiling is specific to Copper-Lilac",
                 x=0.02, y=0.965, ha="left", va="top", fontsize=13,
                 fontweight="bold", color=INK)
    fig.text(0.02, 0.885,
             "Quorum: >=2 collars of each unit inside one spatial cluster, >=4 total. Across 248 dyads there are 10 months at a fusion rate of exactly\n"
             "1.000 in the whole project, and every one is Copper-Lilac - beginning the month nine collars were deployed.",
             ha="left", va="top", fontsize=8.6, color=INK2, linespacing=1.5)
    out = a.dir / "fusion_quorum_all_dyads.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
