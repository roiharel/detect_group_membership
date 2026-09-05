"""Figure: sleep-site membership from the EAS night-location table."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
COPPER, LILAC = "#2a78d6", "#eb6834"
COPPER_DK, LILAC_DK = "#0d366b", "#a33c12"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
ORDER = ["Copper_original", "Copper_new", "Lilac_original", "Lilac_new"]
NICE = {"Copper_original": "Copper\n(original)", "Copper_new": "Copper\n(Jul-2025)",
        "Lilac_original": "Lilac\n(original)", "Lilac_new": "'Lilac'\n(Jul-2025)"}

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
                    default=PROJECT / "outputs" / "sleep_sites_eas_2026-09-03")
    ap.add_argument("--focal", default="24AE04_6L7M")
    a = ap.parse_args()
    aff = pd.read_csv(a.dir / "sleep_site_affinity.csv")
    case = pd.read_csv(a.dir / f"case_{a.focal}_monthly.csv", parse_dates=["month"])
    aff = aff[aff.Copper_original > 0.4]   # 25AB44/26AB52: too few nights to place

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2))
    fig.subplots_adjust(top=0.72, bottom=0.16, left=0.065, right=0.985, wspace=0.26)

    ax = axes[0]
    for i, c in enumerate(ORDER):
        s = aff[aff.focal_cohort == c].sort_values("lilac_minus_copper")
        if not len(s):
            continue
        x = np.full(len(s), i) + np.linspace(-0.23, 0.23, len(s))
        col = [LILAC if v > 0 else COPPER for v in s.lilac_minus_copper]
        ax.scatter(x, s.lilac_minus_copper, s=85, c=col, edgecolor="#ffffff",
                   linewidth=1.0, zorder=4)
    ax.axhline(0, color=INK2, lw=1.4, zorder=2)
    ax.annotate("24AE04\n(transferred)", xy=(2 - 0.23, -0.2844), xytext=(1.30, -0.245),
                fontsize=7.8, color=INK, ha="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    ax.text(3.45, 0.028, "sleeps at LILAC sites", fontsize=8, color=LILAC_DK,
            ha="right", fontweight="bold")
    ax.text(3.45, -0.045, "sleeps at COPPER sites", fontsize=8, color=COPPER_DK,
            ha="right", fontweight="bold")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([NICE[c] for c in ORDER], fontsize=8.2)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_xlim(-0.5, 3.5)
    dress(ax, "A  All nine sleep at Copper sites",
          "same-site rate with Lilac minus with Copper", "difference in same-site rate")

    ax = axes[1]
    ax.plot(case.month, case.with_lilac, color=LILAC, lw=2.4, marker="o", ms=5,
            mec=SURF, mew=0.9, zorder=4, label="same site as Lilac")
    ax.plot(case.month, case.with_copper, color=COPPER, lw=2.4, marker="s", ms=5,
            mec=SURF, mew=0.9, zorder=4, label="same site as Copper")
    cross = pd.Timestamp("2025-06-01")
    ax.axvline(cross, color=INK2, lw=1.3, ls=(0, (4, 3)), zorder=2)
    ax.text(cross, 1.06, " switches", fontsize=8, color=INK2, va="top", fontweight="bold")
    ax.axvline(pd.Timestamp("2025-07-31"), color=MUTED, lw=1.1, ls=(0, (2, 2)), zorder=2)
    ax.text(pd.Timestamp("2025-07-31"), 0.06, " collar deployment", fontsize=7.4,
            color=MUTED, va="bottom")
    ax.set_ylim(-0.04, 1.12)
    ax.legend(loc="lower left", fontsize=8.2, labelcolor=INK2)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y-%m"))
    dress(ax, "B  24AE04_6L7M transferred in June 2025",
          "adult male, Lilac-born; fraction of nights at the same site", "fraction of nights")

    fig.suptitle("Sleep-site membership from the EAS night-location table",
                 x=0.065, y=0.965, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.065, 0.885,
             "Same site = identical cluster_united in processed/2025/gps/individual_night_locations - the project's own sleep-site assignment, not a distance threshold.\n"
             "Original Lilac animals still prefer Lilac sites (6 of 9), so the measure detects group membership where it exists. The July-2025 collars show none of it.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    out = a.dir / "sleep_sites_eas.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
