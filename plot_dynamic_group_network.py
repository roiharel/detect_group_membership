"""Draw the dynamic, weighted, directional, effort-corrected group network.

A circular (chord) layout is used rather than a force layout: with 26 groups and
a handful of strong dyads, spring layouts collapse the weakly-connected majority
into an unreadable clump. Nodes are ordered by hierarchical clustering of the
association matrix so related groups sit adjacent and chords stay short.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

PROJECT = Path(__file__).resolve().parent
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
AXIS, SURF = "#c3c2b7", "#fcfcfb"
RAMP = ["#cde2fb", "#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
NODE, NODE_EDGE, HILITE = "#9ec5f4", "#1c5cab", "#eb6834"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "legend.frameon": False,
})


def ring_positions(order: list[str], radius: float = 1.0) -> dict[str, tuple[float, float]]:
    n = len(order)
    return {u: (radius * np.cos(np.pi / 2 - 2 * np.pi * i / n),
                radius * np.sin(np.pi / 2 - 2 * np.pi * i / n))
            for i, u in enumerate(order)}


def cluster_order(units: list[str], agg: pd.DataFrame) -> list[str]:
    idx = {u: i for i, u in enumerate(units)}
    M = np.zeros((len(units), len(units)))
    for r in agg.itertuples(index=False):
        i, j = idx[r.unit_a], idx[r.unit_b]
        M[i, j] = M[j, i] = r.association
    mx = M.max() or 1.0
    D = 1.0 - M / mx
    np.fill_diagonal(D, 0.0)
    try:
        z = linkage(squareform((D + D.T) / 2, checks=False), method="average")
        return [units[i] for i in leaves_list(z)]
    except Exception:
        return units


def edge_color(w: float, vmax: float) -> str:
    if vmax <= 0:
        return RAMP[0]
    return RAMP[int(np.clip(np.searchsorted([0.04, 0.12, 0.30, 0.60], w / vmax), 0, 4))]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=PROJECT / "outputs" / "dynamic_group_network_2026-09-03")
    ap.add_argument("--freq", default="monthly")
    ap.add_argument("--min-assoc", type=float, default=0.015)
    ap.add_argument("--bins", type=int, default=6)
    a = ap.parse_args()

    und = pd.read_csv(a.dir / f"group_network_{a.freq}_undirected.csv", parse_dates=["period"])
    dire = pd.read_csv(a.dir / f"group_network_{a.freq}_directed.csv", parse_dates=["period"])
    nodes = pd.read_csv(a.dir / "node_codes.csv")
    code = dict(zip(nodes.unit, nodes.code))
    sizes = dict(zip(nodes.unit, nodes["group_size"])) if "group_size" in nodes else {}
    med = float(np.nanmedian([v for v in sizes.values() if np.isfinite(v)] or [50.0]))

    agg = (und.groupby(["unit_a", "unit_b"], as_index=False)
              .agg(opportunity=("opportunity", "sum"), contact=("contact", "sum")))
    agg["association"] = agg.contact / agg.opportunity
    units = sorted(set(und.unit_a) | set(und.unit_b))
    pos = ring_positions(cluster_order(units, agg))
    vmax = agg[agg.association >= a.min_assoc].association.max()

    dagg = (dire.groupby(["source", "target"], as_index=False)
                .agg(opportunity=("opportunity", "sum"), contact=("contact", "sum")))
    dagg["association"] = dagg.contact / dagg.opportunity
    tot = dagg.groupby("source").association.transform("sum")
    dagg["out_share"] = np.where(tot > 0, dagg.association / tot, np.nan)

    def node_radius(u, scale=1.0):
        s = sizes.get(u, med)
        s = med if not np.isfinite(s) else s
        return (0.035 + 0.030 * np.sqrt(s / 110.0)) * scale

    def draw_nodes(ax, scale=1.0, label=True, fs=7.4, pad=1.13):
        for u, (x, y) in pos.items():
            hot = u in ("Copper", "Lilac")
            ax.add_patch(plt.Circle((x, y), node_radius(u, scale),
                                    facecolor=HILITE if hot else NODE,
                                    edgecolor=NODE_EDGE, lw=0.8, zorder=4))
            if label:
                ang = np.degrees(np.arctan2(y, x))
                ha = "left" if -90 <= ang <= 90 else "right"
                rot = ang if -90 <= ang <= 90 else ang + 180
                ax.text(x * pad, y * pad, code[u], ha=ha, va="center", fontsize=fs,
                        rotation=rot, rotation_mode="anchor",
                        fontweight="bold" if hot else "normal",
                        color=HILITE if hot else INK2, zorder=5)

    fig = plt.figure(figsize=(15.2, 8.0))
    gs = fig.add_gridspec(2, 5, width_ratios=[1.75, 1, 1, 1, 1],
                          left=0.015, right=0.99, top=0.845, bottom=0.02,
                          wspace=0.05, hspace=0.14)
    axm = fig.add_subplot(gs[:, 0])

    for r in dagg[dagg.association >= a.min_assoc].itertuples(index=False):
        x1, y1 = pos[r.source]
        x2, y2 = pos[r.target]
        axm.add_patch(FancyArrowPatch(
            (x1, y1), (x2, y2), connectionstyle="arc3,rad=0.20", arrowstyle="-|>",
            mutation_scale=7 + 20 * float(r.out_share),
            linewidth=0.6 + 6.0 * float(r.out_share),
            color=edge_color(r.association, vmax), alpha=0.9,
            shrinkA=7, shrinkB=9, zorder=2))
    draw_nodes(axm, scale=1.0, fs=8.2)
    axm.set_title("Pooled, directed", loc="left", fontweight="bold", color=INK,
                  fontsize=11.5, pad=16)
    axm.text(0, 1.017, "arrow width = share of outward association",
             transform=axm.transAxes, fontsize=8.2, color=MUTED)

    periods = np.sort(und.period.unique())
    for slot, chunk in zip([gs[i, j] for i in range(2) for j in range(1, 5)],
                           np.array_split(periods, a.bins)):
        ax = fig.add_subplot(slot)
        d = und[und.period.isin(chunk)]
        g = d.groupby(["unit_a", "unit_b"], as_index=False).agg(
            opportunity=("opportunity", "sum"), contact=("contact", "sum"))
        g["association"] = g.contact / g.opportunity
        for r in g[g.association >= a.min_assoc].itertuples(index=False):
            x1, y1 = pos[r.unit_a]
            x2, y2 = pos[r.unit_b]
            ax.add_patch(FancyArrowPatch(
                (x1, y1), (x2, y2), connectionstyle="arc3,rad=0.20", arrowstyle="-",
                linewidth=0.5 + 5.5 * float(r.association),
                color=edge_color(r.association, vmax), alpha=0.9,
                shrinkA=4, shrinkB=4, zorder=2))
        draw_nodes(ax, scale=0.72, fs=5.0, pad=1.16)
        lo = pd.Timestamp(chunk[0]).strftime("%b %Y")
        hi = pd.Timestamp(chunk[-1]).strftime("%b %Y")
        ax.set_title(f"{lo} – {hi}", loc="left", fontsize=8.8, fontweight="bold",
                     color=INK2, pad=5)

    for ax in fig.axes:
        ax.set_xlim(-1.34, 1.34)
        ax.set_ylim(-1.30, 1.30)
        ax.set_aspect("equal")
        ax.axis("off")

    fig.suptitle("Between-group association, corrected for shared collaring effort",
                 x=0.015, y=0.985, ha="left", va="top", fontsize=14.5,
                 fontweight="bold", color=INK)
    fig.text(0.015, 0.928,
             "Edge weight = co-observed animal-pairs sharing an hourly spatial cluster / all co-observed animal-pairs, so a dyad with 2 collars each and one with 20 each are\n"
             "comparable (coverage across groups runs 3%–42%). Node area = demographic group size from EAS, never collar count. Copper and Lilac in orange.",
             ha="left", va="top", fontsize=8.4, color=INK2, linespacing=1.55)

    out = a.dir / f"dynamic_group_network_{a.freq}.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
