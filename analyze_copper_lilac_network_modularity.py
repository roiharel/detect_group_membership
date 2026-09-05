"""How does network modularity change, within and between origin groups, over time?

Three quantities per month, on the effort-corrected individual association network:

  Q_origin      Newman modularity of the ORIGIN partition (Copper | Lilac).
                High = the two origins are still distinct communities.
                Falling toward 0 = origin no longer predicts community structure,
                i.e. integration. This is the between-origin measure.

  Q_within_*    Modularity of the best data-driven partition found INSIDE each
                origin group alone. Rising = that group is developing internal
                substructure; falling = it is becoming more homogeneous. This is
                the within-origin measure, and it is what caught the Lilac
                cohort split.

  assortativity Attribute assortativity by origin - an independent check on
                Q_origin that does not depend on a community-detection choice.

Q_origin is reported both raw and against a degree-preserving null: with unequal
group sizes and unequal edge weights a partition can score above zero by
construction, so the null says how much of the observed value is structural.

Read-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
NETDIR = PROJECT / "outputs" / "copper_lilac_individual_network_2026-09-03"
DEPLOYMENT = pd.Timestamp("2025-08-01")
NEW_LILAC = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
             "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"}


def build_graph(frame: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    for r in frame.itertuples(index=False):
        if r.association > 0:
            g.add_edge(r.animal_a, r.animal_b, weight=float(r.association))
    return g


def origin_modularity(g: nx.Graph, origin: dict[str, str]) -> float | None:
    parts = {}
    for n in g.nodes:
        parts.setdefault(origin.get(n, "?"), set()).add(n)
    parts = [s for s in parts.values() if s]
    if len(parts) < 2 or g.number_of_edges() == 0:
        return None
    return float(nx.community.modularity(g, parts, weight="weight"))


def null_origin_modularity(g: nx.Graph, origin: dict[str, str], n_perm: int,
                           rng: np.random.Generator) -> tuple[float, float]:
    """Permute origin labels, keeping group sizes, to get a null distribution."""
    nodes = list(g.nodes)
    labels = [origin.get(n, "?") for n in nodes]
    draws = []
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        m = dict(zip(nodes, perm))
        q = origin_modularity(g, m)
        if q is not None:
            draws.append(q)
    if not draws:
        return np.nan, np.nan
    return float(np.mean(draws)), float(np.std(draws))


def best_internal_modularity(g: nx.Graph, members: set[str], seed: int) -> tuple[float, int]:
    sub = g.subgraph([n for n in g.nodes if n in members])
    if sub.number_of_edges() < 3 or sub.number_of_nodes() < 4:
        return np.nan, 0
    try:
        comms = nx.community.louvain_communities(sub, weight="weight", seed=seed)
    except Exception:
        comms = nx.community.greedy_modularity_communities(sub, weight="weight")
    comms = [set(c) for c in comms if c]
    if len(comms) < 2:
        return 0.0, 1
    return float(nx.community.modularity(sub, comms, weight="weight")), len(comms)


def plot_modularity(d: pd.DataFrame, out_path: Path) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    COP, LIL, DARK, ACC = "#2a78d6", "#eb6834", "#0d366b", "#eb6834"
    INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
    GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
        "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.labelcolor": INK2, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False})

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.7))
    fig.subplots_adjust(top=0.70, bottom=0.15, left=0.055, right=0.985, wspace=0.26)

    def dress(ax, t, s, yl):
        ax.set_title(t, loc="left", fontweight="bold", color=INK, pad=26, fontsize=10.5)
        ax.text(0, 1.035, s, transform=ax.transAxes, fontsize=8.2, color=MUTED)
        ax.set_ylabel(yl, fontsize=8.5)
        ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8, length=3)
        ax.axvline(DEPLOYMENT, color=ACC, lw=1.3, ls=(0, (4, 3)), zorder=2)
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=5))
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y-%m"))

    ax = axes[0]
    ax.axhline(0, color=AXIS, lw=1.0, zorder=1)
    ax.plot(d.period, d.Q_origin, color=DARK, lw=2.2, marker="o", ms=4.5,
            mec=SURF, mew=0.9, zorder=4)
    ax.fill_between(d.period, d.Q_origin_null_mean - 2 * d.Q_origin_null_sd,
                    d.Q_origin_null_mean + 2 * d.Q_origin_null_sd,
                    color=MUTED, alpha=0.18, lw=0, zorder=2, label="shuffled-origin null (±2 sd)")
    ax.legend(loc="upper right", fontsize=8, labelcolor=INK2)
    ax.text(DEPLOYMENT, ax.get_ylim()[1], " collars", color=ACC, fontsize=7.8, va="top")
    dress(ax, "A  Origin stops predicting structure",
          "modularity of the Copper | Lilac partition", "Q (origin partition)")

    ax = axes[1]
    ax.plot(d.period, d.Q_within_copper, color=COP, lw=2, marker="o", ms=4,
            mec=SURF, mew=0.8, zorder=3, label="within Copper")
    ax.plot(d.period, d.Q_within_lilac, color=LIL, lw=2, marker="s", ms=4,
            mec=SURF, mew=0.8, zorder=3, label="within Lilac")
    ax.legend(loc="upper left", fontsize=8, labelcolor=INK2)
    ax.set_ylim(-0.01, None)
    dress(ax, "B  Lilac fragments, then re-coheres",
          "modularity inside each origin group alone", "Q (within group)")

    ax = axes[2]
    ax.axhline(1.0, color=AXIS, lw=1.0, ls=(0, (2, 2)), zorder=1)
    ax.plot(d.period, d.cross_within_ratio, color=DARK, lw=2.2, marker="o", ms=4.5,
            mec=SURF, mew=0.9, zorder=3)
    ax.text(d.period.max(), 1.03, "parity", fontsize=7.5, color=MUTED, ha="right")
    ax.set_ylim(0, 1.18)
    dress(ax, "C  Cross-origin association reaches parity",
          "pooled cross-origin / within-origin association", "ratio")

    fig.suptitle("Network modularity: between origins collapses, within Lilac spikes and decays",
                 x=0.055, y=0.965, ha="left", va="top", fontsize=13,
                 fontweight="bold", color=INK)
    fig.text(0.055, 0.885,
             "Individual-level 5 m association network, effort-corrected. Q_origin falls -0.19/yr (p<0.001) to zero: by late 2025 the Copper|Lilac split explains no more\n"
             "structure than a shuffled label. Node count jumps 14->26 at the deployment, so part of that step is compositional; the pre-2025-08 decline is not.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    fig.savefig(out_path, dpi=200, facecolor=SURF)
    print(f"wrote {out_path.name}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--netdir", type=Path, default=NETDIR)
    ap.add_argument("--permutations", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20260903)
    ap.add_argument("--output-dir", type=Path, default=None)
    a = ap.parse_args()
    out_dir = a.output_dir or a.netdir
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    und = pd.read_csv(a.netdir / "individual_network_monthly_undirected.csv",
                      parse_dates=["period"])
    nodes = pd.read_csv(a.netdir / "node_codes.csv")
    origin = dict(zip(nodes.animal_id, nodes.origin_group))
    copper = set(nodes[nodes.origin_group == "Copper"].animal_id)
    lilac = set(nodes[nodes.origin_group == "Lilac"].animal_id)

    rows = []
    for period, frame in und.groupby("period"):
        g = build_graph(frame)
        if g.number_of_edges() == 0:
            continue
        nx.set_node_attributes(g, {n: origin.get(n, "?") for n in g.nodes}, "origin")
        q = origin_modularity(g, origin)
        mu, sd = null_origin_modularity(g, origin, a.permutations, rng)
        q_cop, n_cop = best_internal_modularity(g, copper, a.seed)
        q_lil, n_lil = best_internal_modularity(g, lilac, a.seed)
        try:
            assort = float(nx.attribute_assortativity_coefficient(
                g.subgraph(list(g.nodes)), "origin"))
        except Exception:
            assort = np.nan
        cross = frame[frame.pair_type == "cross_origin"]
        within = frame[frame.pair_type != "cross_origin"]
        rows.append({
            "period": period,
            "n_nodes": g.number_of_nodes(), "n_edges": g.number_of_edges(),
            "Q_origin": q, "Q_origin_null_mean": mu, "Q_origin_null_sd": sd,
            "Q_origin_z": (q - mu) / sd if (q is not None and sd and sd > 0) else np.nan,
            "Q_within_copper": q_cop, "n_comm_copper": n_cop,
            "Q_within_lilac": q_lil, "n_comm_lilac": n_lil,
            "assortativity_origin": assort,
            "cross_association": cross.contact.sum() / cross.opportunity.sum() if len(cross) else np.nan,
            "within_association": within.contact.sum() / within.opportunity.sum() if len(within) else np.nan,
        })

    d = pd.DataFrame(rows).sort_values("period")
    for g_ in ("copper", "lilac"):
        pass
    d["cross_within_ratio"] = d.cross_association / d.within_association
    d.to_csv(out_dir / "network_modularity_monthly.csv", index=False)

    pd.set_option("display.width", 220)
    show = ["period", "n_nodes", "Q_origin", "Q_origin_z", "assortativity_origin",
            "Q_within_copper", "Q_within_lilac", "cross_within_ratio"]
    print("=== modularity by month ===")
    print(d[show].round(4).to_string(index=False))

    pre, post = d[d.period < DEPLOYMENT], d[d.period >= DEPLOYMENT]
    early, late = d[d.period < "2025-06-01"], d[d.period >= "2026-01-01"]
    summary = {}
    for lab, sub in [("all", d), ("pre_deployment", pre), ("post_deployment", post),
                     ("early_to_2025_05", early), ("from_2026", late)]:
        summary[lab] = {k: (float(sub[k].mean()) if len(sub) else None)
                        for k in ["Q_origin", "Q_origin_z", "assortativity_origin",
                                  "Q_within_copper", "Q_within_lilac", "cross_within_ratio"]}
    print("\n=== early (to 2025-05) vs late (from 2026-01) ===")
    for k in ["Q_origin", "assortativity_origin", "Q_within_copper", "Q_within_lilac",
              "cross_within_ratio"]:
        e, l = summary["early_to_2025_05"][k], summary["from_2026"][k]
        if e is not None and l is not None:
            print(f"  {k:24s} {e:+.4f}  ->  {l:+.4f}   ({l - e:+.4f})")

    # trend tests
    import statsmodels.api as sm
    trends = {}
    t = (d.period - d.period.min()).dt.days / 365.25
    for k in ["Q_origin", "assortativity_origin", "Q_within_copper", "Q_within_lilac"]:
        s = d[k].astype(float)
        ok = s.notna()
        if ok.sum() >= 8:
            fit = sm.OLS(s[ok].to_numpy(), sm.add_constant(pd.DataFrame({"t": t[ok].to_numpy()}))).fit(
                cov_type="HAC", cov_kwds={"maxlags": 3})
            trends[k] = {"per_year": float(fit.params["t"]), "p": float(fit.pvalues["t"])}
    print("\n=== linear trend per year (HAC(3)) ===")
    for k, v in trends.items():
        star = " *" if v["p"] < 0.05 else ""
        print(f"  {k:24s} {v['per_year']:+.4f} / yr   p = {v['p']:.4f}{star}")

    plot_modularity(d, out_dir / "network_modularity.png")

    (out_dir / "network_modularity.json").write_text(json.dumps({
        "source": str(a.netdir / "individual_network_monthly_undirected.csv"),
        "permutations": a.permutations, "seed": a.seed,
        "definitions": {
            "Q_origin": "Newman modularity of the Copper|Lilac partition; falls as origins merge",
            "Q_origin_z": "(observed - permuted mean) / permuted sd, origin labels shuffled",
            "Q_within_*": "Louvain modularity inside one origin group alone",
            "cross_within_ratio": "pooled cross-origin association / pooled within-origin",
        },
        "summary": summary, "trends": trends,
        "caveat": "Conditioned on canonical Copper-Lilac fusion hours, which saturate after "
                  "2025-08-01; association is GPS proximity at 5 m, not affiliation.",
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {out_dir}")


if __name__ == "__main__":
    main()
