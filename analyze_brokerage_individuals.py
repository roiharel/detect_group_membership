"""Who brokers between Copper and Lilac, and is it the same animals over time?

Brokerage here is betweenness_centrality_subset with sources = Copper collars and
targets = Lilac collars, on the effort-corrected association network with edge
distance 1/association. It answers: of all the shortest network paths running
from the Copper side to the Lilac side, what fraction pass through this animal?

Three questions:

1. WHO.       Per-animal brokerage, broken down by origin group, sex and age class.
2. STABLE?    Is there a consistent rank order across weeks, or does the role
              rotate? Measured by Kendall's W (coefficient of concordance) over
              weeks, and by the week-to-week Spearman correlation of ranks.
3. SCALE.     Does the same rank order hold at 1, 2, 5 and 20 m? Spearman between
              radii on per-animal mean brokerage.

IN INFORMATION TERMS
--------------------
Betweenness is the standard structural proxy for control over flow: a high-
brokerage animal sits on many of the shortest routes connecting the two groups,
so anything moving along association ties would tend to pass through it. To
summarise how concentrated that role is, brokerage is normalised to a
distribution over animals and its Shannon entropy taken. The reported
`effective_brokers` is exp(H) - a Hill number giving the equivalent number of
equally-important brokers. effective_brokers close to the number of tracked
animals means the role is spread evenly and the network has no bottleneck; a
value near 1 means a single animal carries the connection.

IMPORTANT: no information, disease or behaviour was measured travelling these
paths. Betweenness is a structural descriptor of a proximity network, and GPS
proximity is not affiliation. Treat "information" as a metaphor for network
position, not a measured flow.

Read-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

from analyze_copper_lilac_weekly_network_metrics import weekly_pair_rates

PROJECT = Path(__file__).resolve().parent
SRC = PROJECT / "outputs" / "copper_lilac_effort_corrected_integration"
POSITIONS = SRC / "copper_lilac_fusion_2min_positions.parquet"
EAS_META = Path(r"Z:\baboon\working\data\processed\2025\acc\movebank_metadata.csv")
LOCAL_META = Path(r"C:\Users\rharel\Documents\Github\MBRP_basic\data\metadata_animal_name.csv")
DEPLOYMENT = pd.Timestamp("2025-08-01")
NEW_LILAC = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
             "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"}


def load_demography() -> tuple[pd.DataFrame, str]:
    path, src = (EAS_META, "EAS") if EAS_META.exists() else (LOCAL_META, "local copy")
    m = pd.read_csv(path, dtype=str, low_memory=False)
    keep = m[["animal-id", "animal-sex", "animal-comments"]].rename(columns={
        "animal-id": "animal_id", "animal-sex": "sex", "animal-comments": "age_class"})
    keep = keep.drop_duplicates("animal_id")
    keep["sex"] = keep.sex.str.strip().str.lower().map({"m": "male", "f": "female"})
    keep["age_class"] = keep.age_class.str.strip().str.lower()
    return keep, f"{src}: {path}"


def brokerage_table(pairs: pd.DataFrame, origin: dict) -> pd.DataFrame:
    rows = []
    for (radius, week), frame in pairs.groupby(["radius_m", "week"]):
        g = nx.Graph()
        for r in frame.itertuples(index=False):
            if r.association > 0:
                g.add_edge(r.animal_a, r.animal_b, distance=1.0 / float(r.association))
        cop = [n for n in g.nodes if origin.get(n) == "Copper"]
        lil = [n for n in g.nodes if origin.get(n) == "Lilac"]
        if not cop or not lil or g.number_of_edges() < 3:
            continue
        bs = nx.betweenness_centrality_subset(g, sources=cop, targets=lil,
                                              normalized=True, weight="distance")
        tot = sum(bs.values())
        for n, v in bs.items():
            rows.append({"radius_m": radius, "week": week, "animal_id": n,
                         "origin_group": origin.get(n), "brokerage": float(v),
                         "share": float(v / tot) if tot > 0 else np.nan,
                         "n_nodes": g.number_of_nodes()})
    return pd.DataFrame(rows)


def effective_brokers(shares: np.ndarray) -> float:
    """exp(Shannon entropy): equivalent number of equally-important brokers."""
    p = shares[np.isfinite(shares) & (shares > 0)]
    if p.size == 0:
        return np.nan
    p = p / p.sum()
    return float(np.exp(-np.sum(p * np.log(p))))


def core_matrix(mat: pd.DataFrame, min_coverage: float = 0.8) -> pd.DataFrame:
    """Largest near-complete submatrix.

    Collars come and go, so no single period contains every animal; requiring
    complete rows outright leaves nothing. Keep animals present in at least
    `min_coverage` of periods, then drop the periods still missing any of them.
    """
    if mat.empty:
        return mat
    cov = mat.notna().mean(axis=0)
    keep = cov[cov >= min_coverage].index
    if len(keep) < 3:
        keep = cov.nlargest(min(8, len(cov))).index
    return mat[keep].dropna(axis=0, how="any")


def kendalls_w(mat: pd.DataFrame) -> tuple[float, int, int]:
    """Concordance of rank orders across periods. mat: rows=period, cols=animal."""
    m = core_matrix(mat)
    n_periods, n_items = m.shape
    if n_periods < 3 or n_items < 3:
        return np.nan, n_periods, n_items
    ranks = m.rank(axis=1)
    s = ((ranks.sum(axis=0) - ranks.sum(axis=0).mean()) ** 2).sum()
    w = 12 * s / (n_periods ** 2 * (n_items ** 3 - n_items))
    return float(w), n_periods, n_items


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--radii", type=float, nargs="+", default=[1.0, 2.0, 5.0, 20.0])
    ap.add_argument("--min-bins", type=int, default=15)
    ap.add_argument("--min-weeks", type=int, default=25,
                    help="animals present in fewer weeks are excluded from rank stability")
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "brokerage_individuals_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    demo, demo_src = load_demography()
    print(f"demography from {demo_src}")

    pos = pd.read_parquet(POSITIONS, columns=["bin_2min", "animal_id", "longitude",
                                              "latitude", "origin_group"])
    origin = pos.drop_duplicates("animal_id").set_index("animal_id").origin_group.to_dict()
    pairs = weekly_pair_rates(pos, a.radii, a.min_bins)
    print(f"weekly pair-rows {len(pairs):,}  weeks {pairs.week.nunique()}  radii {a.radii}")

    bk = brokerage_table(pairs, origin)
    bk = bk.merge(demo, on="animal_id", how="left")
    bk["cohort"] = np.where(bk.animal_id.isin(NEW_LILAC), "new_Aug2025", "original")
    bk["era"] = np.where(bk.week < DEPLOYMENT, "pre", "post")
    bk.to_csv(a.output_dir / "weekly_brokerage_by_animal.csv", index=False)

    # ---- 1. who ----
    print(f"\n{'=' * 74}\nMEAN BROKERAGE BY GROUP (all radii pooled)\n{'=' * 74}")
    for key in ["origin_group", "sex", "age_class", "cohort"]:
        g = (bk.groupby([key], dropna=False)
               .agg(animals=("animal_id", "nunique"), obs=("brokerage", "size"),
                    mean=("brokerage", "mean"), median=("brokerage", "median")))
        print(f"\n-- by {key} --")
        print(g.round(5).to_string())

    per = (bk.groupby(["radius_m", "animal_id", "origin_group", "sex", "age_class", "cohort"],
                      dropna=False)
             .agg(weeks=("week", "nunique"), mean_brokerage=("brokerage", "mean"),
                  mean_share=("share", "mean")).reset_index())
    per["rank"] = per.groupby("radius_m").mean_brokerage.rank(ascending=False)
    per.to_csv(a.output_dir / "brokerage_by_animal_radius.csv", index=False)

    print(f"\n{'=' * 74}\nTOP BROKERS at 5 m (weeks >= {a.min_weeks})\n{'=' * 74}")
    t = per[(per.radius_m == 5.0) & (per.weeks >= a.min_weeks)].nsmallest(12, "rank")
    print(t[["rank", "animal_id", "origin_group", "sex", "age_class", "cohort",
             "weeks", "mean_brokerage"]].round(5).to_string(index=False))

    # ---- 2. stability over time ----
    print(f"\n{'=' * 74}\nRANK STABILITY OVER WEEKS\n{'=' * 74}")
    stab = {}
    print("(monthly aggregation: weekly matrices are too gappy for a concordance test)")
    bk["month"] = bk.week.values.astype("datetime64[M]")
    for r in a.radii:
        s = bk[bk.radius_m == r]
        keep = s.groupby("animal_id").week.nunique()
        keep = keep[keep >= a.min_weeks].index
        mat = (s[s.animal_id.isin(keep)]
               .pivot_table(index="month", columns="animal_id", values="brokerage"))
        w, npd, nit = kendalls_w(mat)
        core = core_matrix(mat)
        rk = core.rank(axis=1)
        cons = [spearmanr(rk.iloc[i], rk.iloc[i + 1]).statistic for i in range(len(rk) - 1)]
        cons = [c for c in cons if np.isfinite(c)]
        # how much of the ranking is degenerate ties (many animals at exactly zero)?
        zero_frac = float((s.brokerage == 0).mean())
        stab[r] = {"kendall_w": w, "periods": npd, "animals": nit,
                   "median_period_to_period_spearman": float(np.median(cons)) if cons else np.nan,
                   "fraction_zero_brokerage": zero_frac}
        wtxt = f"{w:.3f}" if np.isfinite(w) else "  n/a"
        rtxt = (f"{stab[r]['median_period_to_period_spearman']:+.3f}"
                if cons else "  n/a")
        print(f"  {r:>5g} m   Kendall W = {wtxt}  ({npd} months x {nit} animals)   "
              f"month-to-month rho = {rtxt}   zero-brokerage rows {zero_frac:.0%}")

    # ---- 3. stability across radii ----
    print(f"\n{'=' * 74}\nRANK AGREEMENT BETWEEN RADII (per-animal mean brokerage)\n{'=' * 74}")
    wide = per.pivot_table(index="animal_id", columns="radius_m", values="mean_brokerage")
    wide = wide.dropna()
    print(f"{'':>8}" + "".join(f"{r:>9g} m" for r in a.radii))
    cross = {}
    for r1 in a.radii:
        line = f"{r1:>6g} m"
        for r2 in a.radii:
            rho = spearmanr(wide[r1], wide[r2]).statistic
            cross[f"{r1:g}-{r2:g}"] = float(rho)
            line += f"{rho:>11.3f}"
        print(line)
    print(f"\n(n = {len(wide)} animals present at every radius)")

    # ---- information view ----
    print(f"\n{'=' * 74}\nCONCENTRATION OF THE BROKER ROLE\n{'=' * 74}")
    info_rows = []
    for (r, week), s in bk.groupby(["radius_m", "week"]):
        info_rows.append({"radius_m": r, "week": week,
                          "effective_brokers": effective_brokers(s.share.to_numpy()),
                          "n_nodes": s.n_nodes.iloc[0]})
    info = pd.DataFrame(info_rows)
    info["fraction_of_network"] = info.effective_brokers / info.n_nodes
    info.to_csv(a.output_dir / "broker_concentration_weekly.csv", index=False)
    print(f"{'radius':>8}{'pre eff.':>11}{'post eff.':>11}{'pre frac':>11}{'post frac':>11}")
    for r in a.radii:
        s = info[info.radius_m == r]
        pre, post = s[s.week < DEPLOYMENT], s[s.week >= DEPLOYMENT]
        print(f"{r:>6g} m{pre.effective_brokers.mean():>11.2f}{post.effective_brokers.mean():>11.2f}"
              f"{pre.fraction_of_network.mean():>11.3f}{post.fraction_of_network.mean():>11.3f}")

    (a.output_dir / "brokerage_individuals.json").write_text(json.dumps({
        "demography_source": demo_src, "radii_m": a.radii,
        "min_weeks_for_stability": a.min_weeks,
        "rank_stability_over_weeks": stab,
        "rank_agreement_between_radii_spearman": cross,
        "definitions": {
            "brokerage": "betweenness_centrality_subset, sources=Copper targets=Lilac, "
                         "normalised, edge distance 1/association",
            "effective_brokers": "exp(Shannon entropy of the brokerage share distribution); "
                                 "equivalent number of equally-important brokers",
            "kendall_w": "concordance of per-week rank orders; 1 = identical ranking every week",
        },
        "caveats": [
            "GPS proximity is not affiliation, and nothing was measured travelling these paths; "
            "betweenness is a structural descriptor, not a flow measurement.",
            "Sex/age come from collar deployment metadata: 33 adult vs 7 subadult across 38 "
            "animals, so age-class contrasts rest on very few individuals.",
            "Conditioned on canonical Copper-Lilac fusion hours, which saturate after 2025-08-01.",
        ],
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
