"""Figure: do the July-2025 'Lilac' collars roost with Lilac or with Copper?"""
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
                    default=PROJECT / "outputs" / "lilac_range_sleep_sites_2026-09-03")
    a = ap.parse_args()
    sleep = pd.read_csv(a.dir / "affinity_frac_same_site.csv")
    rng = pd.read_csv(a.dir / "affinity_range_overlap.csv")
    # 25AB44 was tracked 12 days; its estimates are unusable
    sleep = sleep[sleep.Copper_original > 0.3]

    fig, axes = plt.subplots(1, 3, figsize=(14.6, 5.2))
    fig.subplots_adjust(top=0.71, bottom=0.19, left=0.06, right=0.985, wspace=0.30)

    # A - roost affinity per animal
    ax = axes[0]
    for i, c in enumerate(ORDER):
        s = sleep[sleep.focal_cohort == c].sort_values("lilac_minus_copper")
        if not len(s):
            continue
        x = np.full(len(s), i) + np.linspace(-0.22, 0.22, len(s))
        col = [LILAC if v > 0 else COPPER for v in s.lilac_minus_copper]
        ax.scatter(x, s.lilac_minus_copper, s=70, c=col, edgecolor="#ffffff",
                   linewidth=0.9, zorder=4)
    ax.axhline(0, color=INK2, lw=1.4, zorder=2)
    ax.text(3.42, 0.035, "roosts with LILAC", fontsize=7.8, color=LILAC_DK,
            ha="right", fontweight="bold")
    ax.text(3.42, -0.06, "roosts with COPPER", fontsize=7.8, color=COPPER_DK,
            ha="right", fontweight="bold")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([NICE[c] for c in ORDER], fontsize=8)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_xlim(-0.5, 3.5)
    dress(ax, "A  The July collars roost with Copper",
          "roost sharing with Lilac minus with Copper", "difference in shared-roost fraction")

    # B - class means
    ax = axes[1]
    pairs = pd.read_csv(a.dir / "sleep_site_sharing_pairs.csv")
    meta = pd.read_csv(a.dir / "affinity_frac_same_site.csv")[["focal", "focal_cohort"]]
    m = dict(zip(meta.focal, meta.focal_cohort))
    pairs["cls"] = [" x ".join(sorted([m.get(x, "?"), m.get(y, "?")]))
                    for x, y in zip(pairs.animal_a, pairs.animal_b)]
    key = ["Lilac_new x Lilac_original", "Lilac_new x Copper_original",
           "Lilac_original x Lilac_original", "Copper_original x Copper_original"]
    labels = {"Lilac_new x Lilac_original": "'Lilac' Jul-2025\nvs Lilac original",
              "Copper_original x Lilac_new": "'Lilac' Jul-2025\nvs Copper original",
              "Lilac_original x Lilac_original": "Lilac original\nvs itself",
              "Copper_original x Copper_original": "Copper original\nvs itself"}
    g = (pairs[pairs.cls.isin(labels)].groupby("cls").frac_same_site.agg(["mean", "size"])
         .reindex(list(labels)).dropna())
    y = np.arange(len(g))[::-1]
    cols = [LILAC, COPPER, LILAC_DK, COPPER_DK][:len(g)]
    ax.barh(y, g["mean"], color=cols, height=0.68, zorder=3)
    for yi, (v, n) in zip(y, zip(g["mean"], g["size"])):
        ax.text(v + 0.008, yi, f"{v:.3f}  (n={int(n)})", va="center", fontsize=8, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[i] for i in g.index], fontsize=7.8)
    ax.set_xlim(0, 1.02)
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.set_xlabel("fraction of nights sharing a roost", fontsize=8.5)
    dress(ax, "B  Lowest of all: 'Lilac' with Lilac",
          "pairwise roost sharing by class", "")

    # C - range overlap is uninformative
    ax = axes[2]
    for i, c in enumerate(ORDER):
        s = rng[rng.focal_cohort == c]
        if not len(s):
            continue
        x = np.full(len(s), i) + np.linspace(-0.22, 0.22, len(s))
        col = [LILAC if v > 0 else COPPER for v in s.lilac_minus_copper]
        ax.scatter(x, s.lilac_minus_copper, s=70, c=col, edgecolor="#ffffff",
                   linewidth=0.9, zorder=4)
    ax.axhline(0, color=INK2, lw=1.4, zorder=2)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([NICE[c] for c in ORDER], fontsize=8)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(-0.09, 0.09)
    dress(ax, "C  Home range says nothing",
          "range overlap with Lilac minus with Copper", "difference in Bhattacharyya overlap")

    fig.suptitle("Where do the July-2025 'Lilac' collars sleep, and whose range do they share?",
                 x=0.06, y=0.965, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.06, 0.885,
             "Full GPS record, post-deployment, no fusion-hour conditioning. Collars record 03:00-16:00 UTC only, so a roost is scored when a pair is within 100 m at BOTH\n"
             "dusk and the following dawn. Original Lilac animals do show Lilac roost affinity, so the measure detects group membership where it exists - the July collars show none.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    out = a.dir / "lilac_range_sleep_sites.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
