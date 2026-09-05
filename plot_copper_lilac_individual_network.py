"""Draw the individual-level Copper-Lilac network.

Two layouts:

  --layout arc        Copper on the left arc, Lilac on the right, positions fixed
                      across panels. Cross-origin association crosses the middle,
                      so integration is visible as the centre filling in. Best for
                      comparing periods, because nodes never move.

  --layout proximity  Force-directed per period, using association as attraction.
                      Position carries meaning: animals that associate sit close.
                      Origin is colour only, never position, so the two colours
                      separating or mixing is an emergent result rather than
                      something the layout imposed. Positions are NOT comparable
                      between panels.

Animals with no co-observation in a period are drawn grey and small, so a thinning
network is distinguishable from animals simply not being tracked.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

PROJECT = Path(__file__).resolve().parent
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
SURF, OFF = "#fcfcfb", "#c9c8c2"
COPPER, LILAC = "#2a78d6", "#eb6834"
COPPER_DK, LILAC_DK = "#0d366b", "#a33c12"
CROSS = ["#cde2fb", "#9ec5f4", "#5598e7", "#256abf", "#0d366b"]

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "legend.frameon": False,
})


def arc_positions(copper, lilac, gap_deg=14.0):
    pos = {}
    for ids, lo, hi in ((lilac, -90 + gap_deg / 2, 90 - gap_deg / 2),
                        (copper, 90 + gap_deg / 2, 270 - gap_deg / 2)):
        n = max(len(ids), 1)
        for i, u in enumerate(ids):
            ang = np.radians(lo + (hi - lo) * (i / max(n - 1, 1)))
            pos[u] = (np.cos(ang), np.sin(ang))
    return pos


def proximity_positions(frame, all_nodes, seed=7):
    """Force-directed on association; isolated animals pushed to the rim."""
    g = nx.Graph()
    g.add_nodes_from(all_nodes)
    for r in frame.itertuples(index=False):
        if r.association > 0:
            g.add_edge(r.animal_a, r.animal_b, weight=float(r.association) ** 0.5)
    active = [n for n in g.nodes if g.degree(n) > 0]
    sub = g.subgraph(active)
    p = nx.spring_layout(sub, weight="weight", seed=seed, k=0.55, iterations=500) if active else {}
    if p:
        xy = np.array(list(p.values()))
        c, s = xy.mean(axis=0), np.abs(xy - xy.mean(axis=0)).max() or 1.0
        p = {k: ((v[0] - c[0]) / s * 0.82, (v[1] - c[1]) / s * 0.82) for k, v in p.items()}
    idle = [n for n in all_nodes if n not in p]
    for i, n in enumerate(idle):
        ang = 2 * np.pi * i / max(len(idle), 1)
        p[n] = (1.12 * np.cos(ang), 1.12 * np.sin(ang))
    return p


def ecolor(w, vmax):
    if vmax <= 0:
        return CROSS[0]
    return CROSS[int(np.clip(np.searchsorted([0.05, 0.15, 0.35, 0.65], w / vmax), 0, 4))]


# Association values span roughly 0.01-0.19, so a linear map like `0.25 + 3*w`
# renders every edge at well under 1 px with no visible variation. Normalise by a
# reference (the 95th percentile of all pooled edges) and expand the low end with
# a gamma < 1, so differences among the many weak edges stay legible while the
# strongest still dominate. One shared reference across all pair types, so a
# thicker Lilac edge really does mean stronger association than a Copper one.
W_MIN, W_MAX, GAMMA = 0.35, 7.0, 0.65


def width(assoc: float, ref: float, scale: float = 1.0) -> float:
    if ref <= 0:
        return W_MIN * scale
    x = float(np.clip(assoc / ref, 0.0, 1.0)) ** GAMMA
    return (W_MIN + (W_MAX - W_MIN) * x) * scale


def top_by_type(frame: pd.DataFrame, frac: float) -> pd.DataFrame:
    """Keep the strongest `frac` WITHIN each pair type.

    A global quantile silently deletes whole pair types when one is much more
    cohesive than another - with a global 80th percentile on this data every
    within-Copper and every cross-origin edge disappeared while within-Lilac
    survived, making Copper look disconnected when it is merely less cohesive.
    """
    if frac >= 1.0:
        return frame
    out = []
    for _, sub in frame.groupby("pair_type"):
        out.append(sub[sub.association >= sub.association.quantile(1.0 - frac)])
    return pd.concat(out, ignore_index=True) if out else frame


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=PROJECT / "outputs" / "copper_lilac_individual_network_2026-09-03")
    ap.add_argument("--freq", default="monthly")
    ap.add_argument("--layout", choices=["arc", "proximity"], default="arc")
    ap.add_argument("--min-assoc", type=float, default=0.005)
    ap.add_argument("--top-frac", type=float, default=0.6,
                    help="pooled panel: keep the strongest this fraction WITHIN each pair "
                         "type. Applied per type because a global cut deletes whole types "
                         "when one group is much more cohesive than the other. 1.0 = show all.")
    ap.add_argument("--bins", type=int, default=6)
    a = ap.parse_args()

    und = pd.read_csv(a.dir / f"individual_network_{a.freq}_undirected.csv", parse_dates=["period"])
    dire = pd.read_csv(a.dir / f"individual_network_{a.freq}_directed.csv", parse_dates=["period"])
    nodes = pd.read_csv(a.dir / "node_codes.csv")
    code = dict(zip(nodes.animal_id, nodes.code))
    cohort = dict(zip(nodes.animal_id, nodes.cohort))
    origin = dict(zip(nodes.animal_id, nodes.origin_group))
    all_nodes = list(nodes.animal_id)

    meta_path = a.dir / "metadata.json"
    radius_m = 5.0
    if meta_path.exists():
        import json
        radius_m = float(json.loads(meta_path.read_text(encoding="utf-8")).get("radius_m", 5.0))

    cop = list(nodes[nodes.origin_group == "Copper"].sort_values(["cohort", "animal_id"]).animal_id)
    lil = list(nodes[nodes.origin_group == "Lilac"]
               .assign(k=lambda d: (d.cohort == "new_Aug2025").astype(int))
               .sort_values(["k", "animal_id"]).animal_id)
    arc = arc_positions(cop, lil)

    agg = (und.groupby(["animal_a", "animal_b", "pair_type"], as_index=False)
              .agg(opportunity=("opportunity", "sum"), contact=("contact", "sum")))
    agg["association"] = agg.contact / agg.opportunity
    vmax = agg.association.quantile(0.97)

    dagg = (dire.groupby(["source", "target", "pair_type"], as_index=False)
                .agg(opportunity=("opportunity", "sum"), contact=("contact", "sum")))
    dagg["association"] = dagg.contact / dagg.opportunity
    tot = dagg.groupby("source").association.transform("sum")
    dagg["out_share"] = np.where(tot > 0, dagg.association / tot, np.nan)

    def draw_nodes(ax, pos, active=None, r=0.036, fs=6.4, pad=1.10, label=True):
        for u, (x, y) in pos.items():
            on = active is None or u in active
            is_cop = origin[u] == "Copper"
            new = cohort[u] == "new_Aug2025"
            fc = (COPPER if is_cop else LILAC) if on else OFF
            ec = (COPPER_DK if is_cop else LILAC_DK) if on else "#a9a8a2"
            ax.add_patch(plt.Circle((x, y), r if on else r * 0.62, facecolor=fc, edgecolor=ec,
                                    lw=(2.0 if new else 0.8) if on else 0.6,
                                    alpha=1.0 if on else 0.55, zorder=5 if on else 4))
            if label:
                ang = np.degrees(np.arctan2(y, x))
                ha = "left" if -90 <= ang <= 90 else "right"
                rot = ang if -90 <= ang <= 90 else ang + 180
                ax.text(x * pad, y * pad, code[u], ha=ha, va="center", fontsize=fs,
                        rotation=rot, rotation_mode="anchor",
                        color=(COPPER_DK if is_cop else LILAC_DK) if on else "#b3b2ac",
                        zorder=6)

    ref = float(agg.association.quantile(0.95))

    def draw_edges(ax, pos, frame, directed=False, wscale=1.0, rad_within=0.38, rad_cross=0.13):
        for r in frame.itertuples(index=False):
            s, t = (r.source, r.target) if directed else (r.animal_a, r.animal_b)
            if s not in pos or t not in pos:
                continue
            x1, y1 = pos[s]
            x2, y2 = pos[t]
            cross = r.pair_type == "cross_origin"
            lw = width(r.association, ref, wscale)
            if cross and directed:
                ax.add_patch(FancyArrowPatch(
                    (x1, y1), (x2, y2), connectionstyle=f"arc3,rad={rad_cross}",
                    arrowstyle="-|>", mutation_scale=6 + 2.2 * lw / max(wscale, 0.01),
                    linewidth=lw, color=ecolor(r.association, vmax), alpha=0.9,
                    zorder=3, shrinkA=5, shrinkB=7))
            elif cross:
                ax.add_patch(FancyArrowPatch(
                    (x1, y1), (x2, y2), connectionstyle=f"arc3,rad={rad_cross}",
                    arrowstyle="-", linewidth=lw,
                    color=ecolor(r.association, vmax), alpha=0.85, zorder=3,
                    shrinkA=3, shrinkB=3))
            else:
                is_cop = origin[s] == "Copper"
                ax.add_patch(FancyArrowPatch(
                    (x1, y1), (x2, y2), connectionstyle=f"arc3,rad={rad_within}",
                    arrowstyle="-|>" if directed else "-",
                    mutation_scale=6 + 2.2 * lw / max(wscale, 0.01) if directed else 1,
                    linewidth=lw, color=COPPER_DK if is_cop else LILAC_DK,
                    alpha=0.55, zorder=2, shrinkA=3, shrinkB=5 if directed else 3))

    prox = a.layout == "proximity"
    fig = plt.figure(figsize=(15.2, 8.2))
    gs = fig.add_gridspec(2, 5, width_ratios=[1.75, 1, 1, 1, 1],
                          left=0.015, right=0.99, top=0.835, bottom=0.02,
                          wspace=0.05, hspace=0.14)

    axm = fig.add_subplot(gs[:, 0])
    pooled_pos = proximity_positions(agg, all_nodes) if prox else arc
    draw_edges(axm, pooled_pos, top_by_type(dagg, a.top_frac), directed=True,
               rad_within=0.38 if not prox else 0.10, rad_cross=0.13 if not prox else 0.08)
    draw_nodes(axm, pooled_pos, r=0.040, fs=7.4)
    axm.set_title("Pooled, directed", loc="left", fontweight="bold", color=INK,
                  fontsize=11.5, pad=16)
    pct = f"strongest {a.top_frac:.0%} within each pair type" if a.top_frac < 1 else "all pairs"
    axm.text(0, 1.017, f"{pct}; edge width = association strength",
             transform=axm.transAxes, fontsize=8.2, color=MUTED)
    axm.text(0.34, -1.22, f"association within {radius_m:g} m", fontsize=7.6, color=MUTED)
    for i, (f, lab) in enumerate([(0.25, "weak"), (0.60, "medium"), (1.00, "strong")]):
        y = -1.32 - 0.085 * i
        axm.plot([0.36, 0.70], [y, y], color=INK2, lw=width(f * ref, ref),
                 solid_capstyle="round")
        axm.text(0.75, y, f"{lab}  {f * ref:.3f}", fontsize=7.2, color=INK2, va="center")
    axm.legend(handles=[
        Line2D([], [], marker="o", ls="", ms=8, mfc=COPPER, mec=COPPER_DK, label="Copper"),
        Line2D([], [], marker="o", ls="", ms=8, mfc=LILAC, mec=LILAC_DK, label="Lilac (original)"),
        Line2D([], [], marker="o", ls="", ms=8, mfc=LILAC, mec=LILAC_DK, mew=2.2,
               label="Lilac (collared 2025-07-31)"),
        Line2D([], [], marker="o", ls="", ms=6.5, mfc=OFF, mec="#a9a8a2",
               label="not co-observed in period"),
    ], loc="lower left", fontsize=8, labelcolor=INK2, bbox_to_anchor=(-0.02, -0.01))

    periods = np.sort(und.period.unique())
    for slot, chunk in zip([gs[i, j] for i in range(2) for j in range(1, 5)],
                           np.array_split(periods, a.bins)):
        ax = fig.add_subplot(slot)
        d = und[und.period.isin(chunk)]
        g = d.groupby(["animal_a", "animal_b", "pair_type"], as_index=False).agg(
            opportunity=("opportunity", "sum"), contact=("contact", "sum"))
        g["association"] = g.contact / g.opportunity
        active = set(g.animal_a) | set(g.animal_b)
        pos = proximity_positions(g, all_nodes) if prox else arc
        draw_edges(ax, pos, g[g.association >= a.min_assoc], directed=False, wscale=0.72,
                   rad_within=0.38 if not prox else 0.10, rad_cross=0.13 if not prox else 0.08)
        draw_nodes(ax, pos, active=active, r=0.026, fs=4.2, pad=1.11)
        lo = pd.Timestamp(chunk[0]).strftime("%b %Y")
        hi = pd.Timestamp(chunk[-1]).strftime("%b %Y")
        cross = g[g.pair_type == "cross_origin"]
        rate = cross.contact.sum() / cross.opportunity.sum() if len(cross) else np.nan
        ax.set_title(f"{lo} - {hi}", loc="left", fontsize=8.8, fontweight="bold",
                     color=INK2, pad=5)
        ax.text(0.0, -1.46, f"cross-origin  {rate:.3f}   --   {len(active)} tracked",
                ha="center", fontsize=7.8, color=INK2, fontweight="bold")

    for ax in fig.axes:
        ax.set_xlim(-1.34, 1.34)
        ax.set_ylim(-1.56, 1.26)
        ax.set_aspect("equal")
        ax.axis("off")

    common = (f"Edge weight = 2-minute bins in which the pair was within {radius_m:g} m, divided by all bins in which BOTH were observed. "
              f"Grey = not co-observed in that period.")
    sub = (f"Force-directed per panel: animals that associate are drawn close together, and origin is colour only - never position, so the two colours mixing is a result\n"
           f"rather than something the layout imposed. Positions are not comparable between panels. {common}"
           if prox else
           f"Copper on the left arc, Lilac on the right, so cross-origin association crosses the middle and within-origin stays on its own side.\n{common}")
    fig.suptitle(f"Copper-Lilac at the individual level - association within {radius_m:g} m, "
                 f"corrected for shared observation effort",
                 x=0.015, y=0.985, ha="left", va="top", fontsize=14.5,
                 fontweight="bold", color=INK)
    fig.text(0.015, 0.925, sub, ha="left", va="top", fontsize=8.4, color=INK2, linespacing=1.55)

    out = a.dir / f"copper_lilac_individual_network_{a.freq}_{a.layout}.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)


if __name__ == "__main__":
    main()
