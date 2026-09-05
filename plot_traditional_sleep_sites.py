"""Figure: sleep-site fusion, and whose traditional ground animals sleep on."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
COPPER, LILAC, SHARED = "#2a78d6", "#eb6834", "#c9c8c2"
COPPER_DK, LILAC_DK = "#0d366b", "#a33c12"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
FUSION = pd.Timestamp("2025-08-01")
# First week in which 24AE04 favoured Copper sites (weekly night-location series).
TRANSFER = pd.Timestamp("2025-06-02")
# Its collar was swapped 2025-08-01 (tag 15508 off 11:00, tag 16463 on 13:00), the
# same field operation that deployed the nine new 'Lilac' collars. The transfer
# precedes that swap by 60 days, so it cannot be a collar-change artefact.
COLLAR_SWAP = pd.Timestamp("2025-08-01")

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.labelcolor": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})


def dress(ax, t, s, yl, dates=True):
    ax.set_title(t, loc="left", fontweight="bold", color=INK, pad=26, fontsize=10.5)
    ax.text(0, 1.035, s, transform=ax.transAxes, fontsize=8.2, color=MUTED)
    ax.set_ylabel(yl, fontsize=8.5)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, length=3)
    if dates:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=4))
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y-%m"))


def shade_eras(ax, lo, hi):
    ax.axvspan(lo, FUSION, color="#f2f1ec", zorder=0)
    ax.axvline(FUSION, color=INK2, lw=1.5, zorder=5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=PROJECT / "outputs" / "traditional_sleep_sites_2026-09-03")
    ap.add_argument("--focal", default="24AE04_6L7M")
    a = ap.parse_args()
    fus = pd.read_csv(a.dir / "group_site_sharing_monthly.csv", parse_dates=["month"])
    foc = pd.read_csv(a.dir / f"case_{a.focal}_monthly_site_type.csv", parse_dates=["month"])
    coh = pd.read_csv(a.dir / "cohort_monthly_site_type.csv", parse_dates=["month"])

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.3))
    fig.subplots_adjust(top=0.66, bottom=0.17, left=0.055, right=0.985, wspace=0.28)

    # A - when the groups fused, from night locations alone
    ax = axes[0]
    lo, hi = fus.month.min(), fus.month.max()
    shade_eras(ax, lo, hi)
    ax.plot(fus.month, fus.same_site, color=INK, lw=2.4, marker="o", ms=4.5,
            mec=SURF, mew=0.9, zorder=4)
    ax.set_ylim(-0.04, 1.12)
    ax.text(lo, 1.06, "  separate groups", fontsize=8.4, color=INK2, fontweight="bold")
    ax.text(FUSION, 1.06, "  fused", fontsize=8.4, color=INK2, fontweight="bold")
    dress(ax, "A  The groups fuse in August 2025",
          "nights on which both groups' modal site was the same", "fraction of nights")

    # B - focal animal on traditional ground
    ax = axes[1]
    lo, hi = foc.month.min(), foc.month.max()
    shade_eras(ax, lo, hi)
    ax.stackplot(foc.month, foc.Lilac_traditional, foc.shared, foc.Copper_traditional,
                 colors=[LILAC, SHARED, COPPER], zorder=3,
                 labels=["Lilac-traditional sites", "shared sites", "Copper-traditional sites"])
    ax.annotate("", xy=(TRANSFER, 1.16), xytext=(TRANSFER, 1.02),
                annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.5))
    ax.text(TRANSFER, 1.19, "transfer\nJun 2025", fontsize=8, color=INK,
            fontweight="bold", ha="center", va="bottom", clip_on=False)
    ax.axvline(TRANSFER, color=INK, lw=1.4, ls=(0, (4, 3)), zorder=6)
    ax.text(COLLAR_SWAP, 1.19, "collar swap\n+ fusion", fontsize=8, color=INK2,
            fontweight="bold", ha="center", va="bottom", clip_on=False)
    ax.annotate("", xy=(COLLAR_SWAP, 1.16), xytext=(COLLAR_SWAP, 1.02),
                annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.5))
    ax.annotate("", xy=(COLLAR_SWAP, 1.09), xytext=(TRANSFER, 1.09),
                annotation_clip=False,
                arrowprops=dict(arrowstyle="<|-|>", color=MUTED, lw=1.0,
                                shrinkA=0, shrinkB=0))
    ax.text(TRANSFER + (COLLAR_SWAP - TRANSFER) / 2, 1.105, "60 days", fontsize=7.4,
            color=MUTED, ha="center", va="bottom", clip_on=False)
    ax.set_ylim(0, 1.02)
    ax.set_xlim(lo, hi)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1], loc="lower left", fontsize=7.8, labelcolor=INK2)
    dress(ax, f"B  {a.focal} leaves Lilac ground before its collar changed",
          "adult male, Lilac-born; nights per month by site type", "proportion of nights")

    # C - cohorts on each group's traditional ground
    ax = axes[2]
    coh["era"] = np.where(coh.month < FUSION, "before fusion", "after fusion")
    order = ["Copper_original", "Copper_new", "Lilac_original", "Lilac_new"]
    nice = ["Copper\n(original)", "Copper\n(Jul-25)", "Lilac\n(original)", "'Lilac'\n(Jul-25)"]
    g = coh.groupby(["cohort", "era"])[["Lilac_traditional", "Copper_traditional"]].mean()
    x = np.arange(len(order))
    w = 0.36
    for k, (era, off) in enumerate([("before fusion", -w / 2), ("after fusion", w / 2)]):
        lil = [g.loc[(c, era), "Lilac_traditional"] if (c, era) in g.index else np.nan
               for c in order]
        cop = [g.loc[(c, era), "Copper_traditional"] if (c, era) in g.index else np.nan
               for c in order]
        ax.bar(x + off, lil, w, color=LILAC, alpha=1.0 if k else 0.45, zorder=3,
               edgecolor=SURF, linewidth=0.8)
        ax.bar(x + off, cop, w, bottom=np.nan_to_num(lil), color=COPPER,
               alpha=1.0 if k else 0.45, zorder=3, edgecolor=SURF, linewidth=0.8)
        for xi, (lv, cv) in enumerate(zip(lil, cop)):
            if np.isfinite(lv) or np.isfinite(cv):
                ax.text(x[xi] + off, np.nan_to_num(lv) + np.nan_to_num(cv) + 0.012,
                        "before" if k == 0 else "after", ha="center", fontsize=6.4,
                        color=MUTED, rotation=90, va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(nice, fontsize=8.2)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_ylim(0, 0.72)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=LILAC, label="Lilac-traditional"),
                       plt.Rectangle((0, 0), 1, 1, color=COPPER, label="Copper-traditional")],
              loc="upper left", fontsize=8, labelcolor=INK2)
    dress(ax, "C  The nine never use Lilac ground",
          "mean proportion of nights, pale = before fusion", "proportion of nights", dates=False)

    fig.suptitle("Sleep sites: when the groups fused, and whose traditional ground they use",
                 x=0.055, y=0.965, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.055, 0.885,
             "Sites are classified once, from nights before 2025-05: a site is group-traditional when that group's usage share is 3x the other's. Usage is a within-group\n"
             "share, so a group with more collars cannot appear to own more sites. Fixing site identity keeps the measure meaningful after the groups merge.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    out = a.dir / "traditional_sleep_sites.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)


if __name__ == "__main__":
    main()

