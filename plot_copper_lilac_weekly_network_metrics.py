"""Figures for weekly betweenness and origin-blind community structure.

  network_betweenness_weekly.png   betweenness at three levels, plus Q_origin/ARI
  network_origin_blind_5m_100m.png communities detected WITHOUT origin, at both
                                   scales. Node fill = detected community;
                                   node ring = true origin. When fill stops
                                   tracking ring, origin has stopped organising
                                   the network.
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

PROJECT = Path(__file__).resolve().parent
DEPLOYMENT = pd.Timestamp("2025-08-01")
COPPER, LILAC = "#2a78d6", "#eb6834"
COPPER_DK, LILAC_DK = "#0d366b", "#a33c12"
# community fills deliberately avoid blue/orange, which are reserved for origin rings
COMM = ["#4a3aa7", "#008300", "#1baf7a", "#eda100", "#e87ba4", "#5b5b58", "#8c5a2b"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"

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
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, length=3)
    ax.axvline(DEPLOYMENT, color=LILAC, lw=1.3, ls=(0, (4, 3)), zorder=2)
    if dates:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=5))
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y-%m"))


def smooth(s, k=5):
    return s.rolling(k, center=True, min_periods=2).mean()


def fig_metrics(d, radii, out_path):
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 8.0))
    fig.subplots_adjust(top=0.80, bottom=0.07, left=0.055, right=0.985,
                        wspace=0.26, hspace=0.42)
    for row, r in enumerate(radii):
        s = d[d.radius_m == r].sort_values("week")
        ax = axes[row, 0]
        for col, c, lab in [("between_within_copper_mean", COPPER, "within Copper"),
                            ("between_within_lilac_mean", LILAC, "within Lilac"),
                            ("between_cross_mean", "#4a3aa7", "Copper<->Lilac brokerage")]:
            ax.plot(s.week, s[col], color=c, lw=0.8, alpha=0.30, zorder=2)
            ax.plot(s.week, smooth(s[col]), color=c, lw=2.2, zorder=3, label=lab)
        ax.legend(loc="upper right", fontsize=7.8, labelcolor=INK2)
        dress(ax, f"{'ABC'[0]}{row + 1}  Betweenness, {r:g} m",
              "weekly, 5-week rolling mean", "normalised betweenness")

        ax = axes[row, 1]
        ax.plot(s.week, s.Q_origin, color="#256abf", lw=0.8, alpha=0.30, zorder=2)
        ax.plot(s.week, smooth(s.Q_origin), color="#0d366b", lw=2.2, zorder=3)
        ax.axhline(0, color=AXIS, lw=1.0, zorder=1)
        dress(ax, f"B{row + 1}  Origin modularity, {r:g} m",
              "Q of the Copper | Lilac partition", "Q_origin")

        ax = axes[row, 2]
        ax.plot(s.week, s.ARI_vs_origin, color="#008300", lw=0.8, alpha=0.30, zorder=2)
        ax.plot(s.week, smooth(s.ARI_vs_origin), color="#005f00", lw=2.2, zorder=3)
        ax.axhline(0, color=AXIS, lw=1.0, zorder=1)
        ax.set_ylim(-0.15, 1.05)
        dress(ax, f"C{row + 1}  Communities vs origin, {r:g} m",
              "ARI of origin-blind communities against origin", "adjusted Rand index")

    fig.suptitle("Weekly network structure: brokerage, modularity, and whether origin still explains it",
                 x=0.055, y=0.975, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.055, 0.905,
             "Top row 5 m (close-range mixing), bottom row 100 m (shared space). Betweenness uses 1/association as edge distance and is normalised, so values compare\n"
             "across weeks with different numbers of tracked animals. ARI = 1 means Louvain communities, detected without ever seeing origin, still reproduce it exactly.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    fig.savefig(out_path, dpi=200, facecolor=SURF)
    print("wrote", out_path.name)


def fig_origin_blind(pairs, origin, radii, out_path, bins=6, seed=20260903):
    fig, axes = plt.subplots(len(radii), bins, figsize=(16.2, 3.1 * len(radii) + 2.0))
    fig.subplots_adjust(top=0.76, bottom=0.03, left=0.035, right=0.99,
                        wspace=0.04, hspace=0.16)
    axes = np.atleast_2d(axes)
    weeks = np.sort(pairs.week.unique())
    chunks = np.array_split(weeks, bins)

    for row, r in enumerate(radii):
        for col, chunk in enumerate(chunks):
            ax = axes[row, col]
            f = pairs[(pairs.radius_m == r) & (pairs.week.isin(chunk))]
            g = nx.Graph()
            agg = f.groupby(["animal_a", "animal_b"], as_index=False).agg(
                opportunity=("opportunity", "sum"), contact=("contact", "sum"))
            agg["association"] = agg.contact / agg.opportunity
            for x in agg.itertuples(index=False):
                if x.association > 0:
                    g.add_edge(x.animal_a, x.animal_b, weight=float(x.association))
            if g.number_of_edges() < 3:
                ax.axis("off")
                continue
            try:
                comms = nx.community.louvain_communities(g, weight="weight", seed=seed)
            except Exception:
                comms = nx.community.greedy_modularity_communities(g, weight="weight")
            comms = [set(c) for c in comms if c]
            cmap = {n: i for i, c in enumerate(comms) for n in c}
            pos = nx.spring_layout(g, weight="weight", seed=seed, k=0.55, iterations=400)
            xy = np.array(list(pos.values()))
            c0, sc = xy.mean(axis=0), np.abs(xy - xy.mean(axis=0)).max() or 1.0
            pos = {k: ((v[0] - c0[0]) / sc * 0.85, (v[1] - c0[1]) / sc * 0.85)
                   for k, v in pos.items()}
            ref = agg.association.quantile(0.95) or 1.0
            for x in agg.itertuples(index=False):
                if x.animal_a in pos and x.animal_b in pos and x.association > 0:
                    p1, p2 = pos[x.animal_a], pos[x.animal_b]
                    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#b9c9dd",
                            lw=0.25 + 2.6 * min(x.association / ref, 1.0),
                            alpha=0.55, zorder=2, solid_capstyle="round")
            for n, (x, y) in pos.items():
                ax.add_patch(plt.Circle((x, y), 0.062,
                                        facecolor=COMM[cmap.get(n, 0) % len(COMM)],
                                        edgecolor=COPPER_DK if origin.get(n) == "Copper" else LILAC_DK,
                                        lw=2.0, zorder=4))
            lo = pd.Timestamp(chunk[0]).strftime("%b %Y")
            hi = pd.Timestamp(chunk[-1]).strftime("%b %Y")
            ax.set_title(f"{lo} - {hi}" if row == 0 else "", loc="left",
                         fontsize=8.6, fontweight="bold", color=INK2, pad=5)
            ax.text(0, -1.18, f"{len(comms)} communities", ha="center", fontsize=7.6,
                    color=INK2, fontweight="bold")
            if col == 0:
                ax.text(-1.16, 0, f"{r:g} m", rotation=90, va="center", ha="center",
                        fontsize=12, fontweight="bold", color=INK)
            ax.set_xlim(-1.25, 1.15)
            ax.set_ylim(-1.30, 1.12)
            ax.set_aspect("equal")
            ax.axis("off")

    handles = [Line2D([], [], marker="o", ls="", ms=9, mfc=COMM[i], mec="#7a7975",
                      label=f"community {i + 1}") for i in range(4)]
    handles += [Line2D([], [], marker="o", ls="", ms=9, mfc="#dddcd6", mec=COPPER_DK,
                       mew=2.2, label="ring: Copper origin"),
                Line2D([], [], marker="o", ls="", ms=9, mfc="#dddcd6", mec=LILAC_DK,
                       mew=2.2, label="ring: Lilac origin")]
    fig.legend(handles=handles, loc="upper right", ncol=3, fontsize=8.4,
               labelcolor=INK2, bbox_to_anchor=(0.99, 0.90))
    fig.suptitle("Community structure detected without using origin, at two spatial scales",
                 x=0.035, y=0.975, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.035, 0.905,
             "Louvain communities on the effort-corrected association network, with no knowledge of group identity. Node FILL is the detected community;\n"
             "node RING is the true origin. Early panels: fill tracks ring, so the network organises by origin. Later: fill and ring come apart.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    fig.savefig(out_path, dpi=200, facecolor=SURF)
    print("wrote", out_path.name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=PROJECT / "outputs" / "copper_lilac_weekly_network_metrics_2026-09-03")
    ap.add_argument("--radii", type=float, nargs="+", default=[5.0, 100.0])
    a = ap.parse_args()
    d = pd.read_csv(a.dir / "weekly_network_metrics.csv", parse_dates=["week"])
    pairs = pd.read_csv(a.dir / "weekly_pair_rates.csv", parse_dates=["week"])
    origin = {}
    for c in ("a", "b"):
        origin.update(dict(zip(pairs[f"animal_{c}"], pairs[f"origin_{c}"])))
    fig_metrics(d, a.radii, a.dir / "network_betweenness_weekly.png")
    fig_origin_blind(pairs, origin, a.radii, a.dir / "network_origin_blind_5m_100m.png")


if __name__ == "__main__":
    main()
