"""Figure: who brokers, how stable the role is, and how concentrated."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

PROJECT = Path(__file__).resolve().parent
DEPLOYMENT = pd.Timestamp("2025-08-01")
COPPER, LILAC = "#2a78d6", "#eb6834"
COPPER_DK, LILAC_DK = "#0d366b", "#a33c12"
RAMP = ["#cde2fb", "#9ec5f4", "#3987e5", "#0d366b"]
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=PROJECT / "outputs" / "brokerage_individuals_2026-09-03")
    ap.add_argument("--radii", type=float, nargs="+", default=[1.0, 2.0, 5.0, 20.0])
    a = ap.parse_args()

    per = pd.read_csv(a.dir / "brokerage_by_animal_radius.csv")
    info = pd.read_csv(a.dir / "broker_concentration_weekly.csv", parse_dates=["week"])

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.4))
    fig.subplots_adjust(top=0.71, bottom=0.20, left=0.055, right=0.985, wspace=0.30)

    # A - per-animal brokerage at 5 m
    ax = axes[0]
    s = per[(per.radius_m == 5.0) & (per.weeks >= 10)].sort_values(
        "mean_brokerage", ascending=False).head(20)
    y = np.arange(len(s))[::-1]
    cols = [COPPER if o == "Copper" else LILAC for o in s.origin_group]
    edges = ["#0b0b0b" if c == "new_Aug2025" else "none" for c in s.cohort]
    ax.barh(y, s.mean_brokerage, color=cols, edgecolor=edges, linewidth=1.6,
            height=0.74, zorder=3)
    labs = [f"{aid.split('_')[0][2:]}  {str(sx)[:1]}/{str(ag)[:3]}"
            for aid, sx, ag in zip(s.animal_id, s.sex, s.age_class)]
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=7.4)
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.set_xlabel("mean brokerage at 5 m", fontsize=8.5)
    ax.legend(handles=[
        Line2D([], [], marker="s", ls="", ms=9, color=COPPER, label="Copper"),
        Line2D([], [], marker="s", ls="", ms=9, color=LILAC, label="Lilac"),
        Line2D([], [], marker="s", ls="", ms=9, mfc="#dddcd6", mec="#0b0b0b", mew=1.6,
               label="collared 2025-07-31")],
        loc="lower right", fontsize=8, labelcolor=INK2)
    dress(ax, "A  A few animals carry the connection",
          "top 20 by mean brokerage; label = ID, sex/age", "")

    # B - rank agreement between radii
    ax = axes[1]
    wide = per.pivot_table(index="animal_id", columns="radius_m",
                           values="mean_brokerage").dropna()
    from scipy.stats import spearmanr
    n = len(a.radii)
    M = np.ones((n, n))
    for i, r1 in enumerate(a.radii):
        for j, r2 in enumerate(a.radii):
            M[i, j] = spearmanr(wide[r1], wide[r2]).statistic
    im = ax.imshow(M, cmap=mpl.colors.LinearSegmentedColormap.from_list("b", RAMP),
                   vmin=0.3, vmax=1.0, zorder=2)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=9,
                    color="#ffffff" if M[i, j] > 0.72 else INK, fontweight="bold", zorder=3)
    ax.set_xticks(range(n)); ax.set_xticklabels([f"{r:g} m" for r in a.radii], fontsize=8.5)
    ax.set_yticks(range(n)); ax.set_yticklabels([f"{r:g} m" for r in a.radii], fontsize=8.5)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label(
        "Spearman rho", fontsize=8, color=INK2)
    dress(ax, "B  Same brokers at fine scales", f"rank agreement, n={len(wide)} animals", "")

    # C - concentration over time
    ax = axes[2]
    for c, r in zip(RAMP, a.radii):
        d = info[info.radius_m == r].sort_values("week")
        ax.plot(d.week, d.fraction_of_network.rolling(5, center=True, min_periods=2).mean(),
                color=c, lw=2.2, zorder=3, label=f"{r:g} m")
    ax.axvline(DEPLOYMENT, color=LILAC, lw=1.3, ls=(0, (4, 3)), zorder=2)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.legend(loc="upper left", fontsize=8.2, labelcolor=INK2, ncol=2)
    ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=5))
    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y-%m"))
    ax.set_ylim(0, None)
    dress(ax, "C  The role spreads at close range",
          "effective brokers / tracked animals, 5-week mean", "fraction of network")

    fig.suptitle("Who brokers between Copper and Lilac, and does the role stay with them?",
                 x=0.055, y=0.965, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.055, 0.885,
             "Brokerage = betweenness on Copper->Lilac shortest paths, edge distance 1/association. Effective brokers = exp(Shannon entropy) of the brokerage share:\n"
             "the equivalent number of equally-important brokers. Nothing was measured travelling these paths - this is network position, not a flow.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    out = a.dir / "brokerage_individuals.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
