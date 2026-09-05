"""Weekly betweenness and origin-blind community structure for Copper-Lilac.

BETWEENNESS - three levels
--------------------------
  within Copper   betweenness inside the Copper-only subgraph
  within Lilac    betweenness inside the Lilac-only subgraph
  Copper<->Lilac  betweenness_subset with sources=Copper, targets=Lilac: how much
                  of the shortest-path flow BETWEEN the two groups runs through
                  each animal. This is brokerage - it names the individuals the
                  two groups connect through, which a within-group measure cannot.

Edge distance for shortest paths is 1/association: strongly associated animals
are "close". Betweenness is normalised, so values are comparable across periods
with different numbers of tracked animals.

ORIGIN-BLIND COMMUNITIES
------------------------
Communities are detected by Louvain on the association network with NO knowledge
of origin. Their agreement with the true origin split is then measured by
adjusted Rand index and normalised mutual information. ARI ~1 means the network
still organises by origin; ARI ~0 means origin has stopped being the organising
principle - integration, established without ever using the labels.

Run at two scales: 5 m (close-range mixing) and 100 m (shared space). The
contrast is the point - groups can share space long before they mix.

Pair rates are rebuilt weekly from the 2-minute positions, because the saved
pair table is monthly.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

PROJECT = Path(__file__).resolve().parent
SRC = PROJECT / "outputs" / "copper_lilac_effort_corrected_integration"
POSITIONS = SRC / "copper_lilac_fusion_2min_positions.parquet"
R_EARTH_M = 6371000.0
DEPLOYMENT = pd.Timestamp("2025-08-01")


def haversine_m(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2.0) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2.0) ** 2)
    return 2.0 * R_EARTH_M * np.arcsin(np.sqrt(a))


def weekly_pair_rates(pos: pd.DataFrame, radii: list[float], min_bins: int) -> pd.DataFrame:
    """Rebuild effort-corrected pair rates at weekly resolution, for every radius."""
    animals = sorted(pos.animal_id.unique())
    lat = pos.pivot_table(index="bin_2min", columns="animal_id", values="latitude", aggfunc="mean")
    lon = pos.pivot_table(index="bin_2min", columns="animal_id", values="longitude", aggfunc="mean")
    weeks = lat.index.to_period("W")
    codes, uniq = pd.factorize(weeks, sort=True)
    n_weeks = len(uniq)
    origin = pos.drop_duplicates("animal_id").set_index("animal_id").origin_group

    rows = []
    for a, b in combinations(animals, 2):
        if a not in lat.columns or b not in lat.columns:
            continue
        d = haversine_m(lat[a].values, lon[a].values, lat[b].values, lon[b].values)
        ok = ~np.isnan(d)
        if not ok.any():
            continue
        wc, dv = codes[ok], d[ok]
        opp = np.bincount(wc, minlength=n_weeks)
        for r in radii:
            con = np.bincount(wc, weights=(dv <= r).astype(float), minlength=n_weeks)
            keep = opp >= min_bins
            if not keep.any():
                continue
            idx = np.flatnonzero(keep)
            rows.append(pd.DataFrame({
                "week": uniq[idx].to_timestamp(), "radius_m": r,
                "animal_a": a, "animal_b": b,
                "origin_a": origin[a], "origin_b": origin[b],
                "opportunity": opp[idx], "contact": con[idx],
            }))
    out = pd.concat(rows, ignore_index=True)
    out["association"] = out.contact / out.opportunity
    out["pair_type"] = np.where(out.origin_a == out.origin_b, "within_" + out.origin_a.str.lower(),
                                "cross_origin")
    return out


def graph_for(frame: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    for r in frame.itertuples(index=False):
        if r.association > 0:
            g.add_edge(r.animal_a, r.animal_b, weight=float(r.association),
                       distance=1.0 / float(r.association))
    return g


def rarefied_betweenness(g: nx.Graph, cop: list, lil: list, k: int, draws: int,
                         rng: np.random.Generator) -> dict:
    """Betweenness on repeated subsamples of exactly k collars per origin group.

    Normalised betweenness still falls as a network grows: more nodes means more
    alternative routes, so no single animal is essential. Collar count here
    roughly doubles at the 2025-07-31 deployment, so an uncorrected decline in
    betweenness is confounded with sample size. Holding the node count fixed at
    k per group removes that entirely - any remaining change is topological.
    """
    if len(cop) < k or len(lil) < k:
        return {"rar_within_copper": np.nan, "rar_within_lilac": np.nan,
                "rar_cross": np.nan, "rar_draws": 0}
    wc, wl, xr = [], [], []
    for _ in range(draws):
        cs = list(rng.choice(cop, k, replace=False))
        ls = list(rng.choice(lil, k, replace=False))
        sub = g.subgraph(cs + ls)
        if sub.number_of_edges() < 3:
            continue
        bs = nx.betweenness_centrality_subset(sub, sources=cs, targets=ls,
                                              normalized=True, weight="distance")
        xr.append(float(np.mean(list(bs.values()))))
        for members, acc in ((cs, wc), (ls, wl)):
            s2 = g.subgraph(members)
            if s2.number_of_edges() >= 3:
                acc.append(float(np.mean(list(nx.betweenness_centrality(
                    s2, normalized=True, weight="distance").values()))))
    return {"rar_within_copper": float(np.mean(wc)) if wc else np.nan,
            "rar_within_lilac": float(np.mean(wl)) if wl else np.nan,
            "rar_cross": float(np.mean(xr)) if xr else np.nan,
            "rar_draws": len(xr)}


def metrics_for_period(frame: pd.DataFrame, origin: dict, seed: int,
                       rarefy_k: int = 0, rarefy_draws: int = 0,
                       rng: np.random.Generator | None = None) -> dict:
    g = graph_for(frame)
    if g.number_of_edges() < 3:
        return {}
    cop = [n for n in g.nodes if origin.get(n) == "Copper"]
    lil = [n for n in g.nodes if origin.get(n) == "Lilac"]
    out = {"n_nodes": g.number_of_nodes(), "n_edges": g.number_of_edges(),
           "n_copper": len(cop), "n_lilac": len(lil)}

    # brokerage between the two groups
    if cop and lil:
        bs = nx.betweenness_centrality_subset(g, sources=cop, targets=lil,
                                              normalized=True, weight="distance")
        out["between_cross_mean"] = float(np.mean(list(bs.values())))
        out["between_cross_max"] = float(np.max(list(bs.values())))
        out["between_cross_copper"] = float(np.mean([bs[n] for n in cop]))
        out["between_cross_lilac"] = float(np.mean([bs[n] for n in lil]))
        top = max(bs, key=bs.get)
        out["top_broker"] = top
        out["top_broker_origin"] = origin.get(top)
        out["top_broker_value"] = float(bs[top])

    if rarefy_k and rng is not None:
        out.update(rarefied_betweenness(g, cop, lil, rarefy_k, rarefy_draws, rng))

    # within-group betweenness, each subgraph on its own
    for lab, members in (("copper", cop), ("lilac", lil)):
        sub = g.subgraph(members)
        if sub.number_of_edges() >= 3 and sub.number_of_nodes() >= 3:
            bw = nx.betweenness_centrality(sub, normalized=True, weight="distance")
            vals = list(bw.values())
            out[f"between_within_{lab}_mean"] = float(np.mean(vals))
            out[f"between_within_{lab}_max"] = float(np.max(vals))
            out[f"between_within_{lab}_sd"] = float(np.std(vals))
        else:
            out[f"between_within_{lab}_mean"] = np.nan
            out[f"between_within_{lab}_max"] = np.nan
            out[f"between_within_{lab}_sd"] = np.nan

    # origin-blind community structure
    try:
        comms = nx.community.louvain_communities(g, weight="weight", seed=seed)
    except Exception:
        comms = nx.community.greedy_modularity_communities(g, weight="weight")
    comms = [set(c) for c in comms if c]
    nodes = list(g.nodes)
    lab_detected = np.array([next(i for i, c in enumerate(comms) if n in c) for n in nodes])
    lab_origin = np.array([0 if origin.get(n) == "Copper" else 1 for n in nodes])
    out["n_communities"] = len(comms)
    out["Q_detected"] = float(nx.community.modularity(g, comms, weight="weight"))
    out["Q_origin"] = float(nx.community.modularity(g, [set(cop), set(lil)], weight="weight")) \
        if cop and lil else np.nan
    out["ARI_vs_origin"] = float(adjusted_rand_score(lab_origin, lab_detected))
    out["NMI_vs_origin"] = float(normalized_mutual_info_score(lab_origin, lab_detected))

    cross = frame[frame.pair_type == "cross_origin"]
    within = frame[frame.pair_type != "cross_origin"]
    out["cross_association"] = (cross.contact.sum() / cross.opportunity.sum()
                                if cross.opportunity.sum() else np.nan)
    out["within_association"] = (within.contact.sum() / within.opportunity.sum()
                                 if within.opportunity.sum() else np.nan)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--radii", type=float, nargs="+", default=[2.0, 5.0, 100.0, 200.0, 400.0])
    ap.add_argument("--min-bins", type=int, default=15,
                    help="minimum co-observed 2-min bins per pair-week (15 = 30 min)")
    ap.add_argument("--rarefy-k", type=int, default=5,
                    help="collars per origin group in the rarefied betweenness; 0 disables")
    ap.add_argument("--rarefy-draws", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20260903)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "copper_lilac_weekly_network_metrics_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    pos = pd.read_parquet(POSITIONS, columns=["bin_2min", "animal_id", "longitude",
                                              "latitude", "origin_group"])
    origin = pos.drop_duplicates("animal_id").set_index("animal_id").origin_group.to_dict()
    print(f"animals {len(origin)}  bins {pos.bin_2min.nunique():,}")

    pairs = weekly_pair_rates(pos, a.radii, a.min_bins)
    pairs.to_csv(a.output_dir / "weekly_pair_rates.csv", index=False)
    print(f"weekly pair-rows {len(pairs):,}  weeks {pairs.week.nunique()}")

    rng = np.random.default_rng(a.seed)
    rows = []
    for (radius, week), frame in pairs.groupby(["radius_m", "week"]):
        m = metrics_for_period(frame, origin, a.seed, a.rarefy_k, a.rarefy_draws, rng)
        if m:
            rows.append({"radius_m": radius, "week": week, **m})
    d = pd.DataFrame(rows).sort_values(["radius_m", "week"])
    d.to_csv(a.output_dir / "weekly_network_metrics.csv", index=False)

    pd.set_option("display.width", 240)
    for r in a.radii:
        s = d[d.radius_m == r]
        print(f"\n{'=' * 78}\n{r:g} m   ({len(s)} weeks)\n{'=' * 78}")
        cols = ["n_copper", "n_lilac",
                "between_within_copper_mean", "between_within_lilac_mean",
                "between_cross_mean",
                "rar_within_copper", "rar_within_lilac", "rar_cross",
                "Q_origin", "ARI_vs_origin", "n_communities"]
        pre, post = s[s.week < DEPLOYMENT], s[s.week >= DEPLOYMENT]
        print(f"{'metric':<32}{'pre-Aug2025':>14}{'post':>12}{'change':>12}")
        for c in cols:
            if c in s and s[c].notna().any():
                p, q = pre[c].mean(), post[c].mean()
                print(f"{c:<32}{p:>14.4f}{q:>12.4f}{q - p:>+12.4f}")
        brokers = s.dropna(subset=["top_broker"]).top_broker.value_counts().head(5)
        if len(brokers):
            print("\nmost frequent top broker between the groups:")
            for k, v in brokers.items():
                print(f"  {k:<16} {origin.get(k):<8} {v:>3} weeks")

    import statsmodels.api as sm
    trends = {}
    print(f"\n{'=' * 78}\nlinear trends per year (HAC(4))\n{'=' * 78}")
    for r in a.radii:
        s = d[d.radius_m == r].dropna(subset=["Q_origin"])
        t = (s.week - s.week.min()).dt.days / 365.25
        print(f"\n-- {r:g} m --")
        for c in ["between_within_copper_mean", "between_within_lilac_mean",
                  "between_cross_mean", "rar_within_copper", "rar_within_lilac",
                  "rar_cross", "Q_origin", "ARI_vs_origin"]:
            y = s[c].astype(float)
            ok = y.notna()
            if ok.sum() < 12:
                continue
            fit = sm.OLS(y[ok].to_numpy(),
                         sm.add_constant(pd.DataFrame({"t": t[ok].to_numpy()}))).fit(
                cov_type="HAC", cov_kwds={"maxlags": 4})
            trends[f"{r:g}m_{c}"] = {"per_year": float(fit.params["t"]),
                                     "p": float(fit.pvalues["t"])}
            star = " *" if fit.pvalues["t"] < 0.05 else ""
            print(f"  {c:<32}{fit.params['t']:+.4f}/yr   p={fit.pvalues['t']:.4f}{star}")

    (a.output_dir / "weekly_network_metrics.json").write_text(json.dumps({
        "radii_m": a.radii, "min_bins_per_pair_week": a.min_bins, "seed": a.seed,
        "weeks": int(d.week.nunique()), "animals": len(origin),
        "definitions": {
            "between_within_*": "normalised weighted betweenness inside that origin group's subgraph",
            "between_cross_*": "betweenness_centrality_subset, sources=Copper targets=Lilac - brokerage",
            "ARI_vs_origin": "adjusted Rand index of Louvain communities vs the origin split; "
                             "communities are detected without using origin",
            "edge_distance": "1/association",
        },
        "trends": trends,
        "caveat": "Conditioned on canonical Copper-Lilac fusion hours, which saturate after "
                  "2025-08-01. Association is GPS proximity, not affiliation.",
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
