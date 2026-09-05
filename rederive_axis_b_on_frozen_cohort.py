"""Gap 7.3: re-derive axis B's weekly network metrics on the FROZEN cohort.

THE PROBLEM. Axis B (within-group splits and modularity) runs on the legacy hourly
source -- `canonical_robust_hourly_membership_local_2h_support`, 1,703,133 rows, 347
animals, 21 units -- filtered to 2025-01-01 onward. Axes A and C run on the frozen
export, 1,924,104 rows, 350 animals, 26 units, from 2024-03-01. Result 2 puts all three
axes in one table, so as it stands it mixes denominators. Nothing about axis B can be
pooled with A and C until it is on the frozen source.

WHAT THIS DOES. Rebuilds the weekly within-group network metrics from the frozen narrow
export, matching the legacy `combined_1100_1600` definition:

  * two sampling hours per day, 11:00 and 16:00 UTC, so 14 timestamps in a full week
    (verified against the legacy file, whose n_timestamps is exactly 7 per scope and 14
    combined);
  * observed positions only -- carried positions are excluded;
  * cluster identity from `association_event_id`, which is shared by animals in the
    same spatial cluster in that hour;
  * split when a unit's observed animals occupy more than one cluster in a timestamp;
    multi-animal split when at least two of those clusters hold two or more animals;
  * association network over dyads co-observed at least once, edge weight = co-clustered
    timestamps / co-observed timestamps;
  * modularity of a greedy-modularity partition of that weighted graph.

It then re-runs the four axis-B findings that the framing rests on and prints them
beside the legacy values, so any that do not survive the cohort change are visible:

  1. the detectability ladder (modularity against collar count)
  2. the coverage-matched between-unit comparison
  3. the between-unit ICC on well-covered group-weeks
  4. week-to-week persistence

Outputs: outputs/general_structure_2026_09/phase4d_axis_b_frozen/
"""

from __future__ import annotations

import argparse
import json
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
LEGACY = Path(
    "outputs/canonical_within_group_density_modularity_ridges/"
    "canonical_within_group_weekly_network_metrics.csv"
)
OUT = Path("outputs/general_structure_2026_09/phase4d_axis_b_frozen")

SAMPLE_HOURS = (11, 16)
MIN_ANIMALS = 3       # below three nodes modularity is not meaningful
WELL_COVERED = 12
MATCHED_BAND = (12, 16)
SEED = 20260904


def build_weekly(narrow: Path) -> pd.DataFrame:
    import networkx as nx
    from networkx.algorithms.community import greedy_modularity_communities, modularity

    cols = ["window_start", "animal_id", "dynamic_social_unit",
            "association_event_id", "is_observed"]
    d = pd.read_csv(narrow, usecols=cols, parse_dates=["window_start"])
    d = d[d["is_observed"] & d["window_start"].dt.hour.isin(SAMPLE_HOURS)].copy()
    d["period_start"] = (d["window_start"]
                         - pd.to_timedelta(d["window_start"].dt.weekday, unit="D")
                         ).dt.normalize()
    print("  frozen source, sampled hours %s: %d animal-timestamps, %d units, "
          "%d weeks" % (SAMPLE_HOURS, len(d), d["dynamic_social_unit"].nunique(),
                        d["period_start"].nunique()))

    rows = []
    for (unit, week), g in d.groupby(["dynamic_social_unit", "period_start"], sort=True):
        animals = g["animal_id"].unique()
        if len(animals) < MIN_ANIMALS:
            continue
        stamps = sorted(g["window_start"].unique())

        # per timestamp: cluster structure of this unit's observed animals
        n_split = 0
        n_multi_split = 0
        largest_fracs = []
        for ts, h in g.groupby("window_start", sort=False):
            sizes = h.groupby("association_event_id").size()
            if len(sizes) > 1:
                n_split += 1
                if (sizes >= 2).sum() >= 2:
                    n_multi_split += 1
            largest_fracs.append(float(sizes.max() / sizes.sum()))

        # dyad co-observation and co-clustering across the week's timestamps
        by_ts = {ts: dict(zip(h["animal_id"], h["association_event_id"]))
                 for ts, h in g.groupby("window_start", sort=False)}
        co_obs: dict = {}
        co_clu: dict = {}
        for ts, mapping in by_ts.items():
            present = sorted(mapping)
            for a, b in combinations(present, 2):
                key = (a, b)
                co_obs[key] = co_obs.get(key, 0) + 1
                if mapping[a] == mapping[b]:
                    co_clu[key] = co_clu.get(key, 0) + 1
        if not co_obs:
            continue

        tot_obs = sum(co_obs.values())
        tot_clu = sum(co_clu.values())
        G = nx.Graph()
        G.add_nodes_from(animals)
        for key, n in co_obs.items():
            w = co_clu.get(key, 0) / n
            if w > 0:
                G.add_edge(key[0], key[1], weight=w)

        mod, ncomm, lcf = 0.0, 1, 1.0
        if G.number_of_edges() > 0:
            try:
                parts = list(greedy_modularity_communities(G, weight="weight"))
                if parts:
                    mod = float(modularity(G, parts, weight="weight"))
                    ncomm = len(parts)
                    lcf = float(max(len(p) for p in parts) / len(animals))
            except Exception:
                pass
        else:
            ncomm = len(animals)
            lcf = 1.0 / len(animals)

        rows.append({
            "dynamic_social_unit": unit,
            "period_start": week,
            "n_animals_observed": int(len(animals)),
            "n_timestamps": int(len(stamps)),
            "possible_dyads": int(len(animals) * (len(animals) - 1) / 2),
            "observed_dyads": int(len(co_obs)),
            "positive_edges": int(G.number_of_edges()),
            "association_density": round(tot_clu / tot_obs, 6) if tot_obs else 0.0,
            "modularity": round(mod, 6),
            "n_communities": int(ncomm),
            "largest_community_fraction": round(lcf, 6),
            "split_timestamp_fraction": round(n_split / len(stamps), 6),
            "multi_animal_split_timestamp_fraction": round(n_multi_split / len(stamps), 6),
            "mean_largest_cluster_fraction": round(float(np.mean(largest_fracs)), 6),
        })
    return pd.DataFrame(rows)


def findings(w: pd.DataFrame, tag: str) -> dict:
    w = w.copy()
    w["period_start"] = pd.to_datetime(w["period_start"])
    w["is_modular"] = w["modularity"] > 0.001
    w["had_split"] = w["split_timestamp_fraction"] > 0

    bands = [(1, 4), (5, 7), (8, 9), (10, 11), (12, 14), (15, 19), (20, 99)]
    ladder = []
    for lo, hi in bands:
        s = w[w["n_animals_observed"].between(lo, hi)]
        if s.empty:
            continue
        ladder.append({"collars_observed": "%d-%d" % (lo, hi),
                       "group_weeks": len(s), "units": s["dynamic_social_unit"].nunique(),
                       "pct_weeks_modular": round(100 * s["is_modular"].mean(), 1),
                       "pct_weeks_with_split": round(100 * s["had_split"].mean(), 1)})

    m = w[w["n_animals_observed"].between(*MATCHED_BAND)]
    matched = m.groupby("dynamic_social_unit").agg(
        group_weeks=("modularity", "size"),
        median_collars=("n_animals_observed", "median"),
        pct_weeks_modular=("is_modular", lambda x: round(100 * x.mean(), 1)),
        max_modularity=("modularity", "max"),
    ).reset_index()
    matched = matched[matched["group_weeks"] >= 5].sort_values(
        "pct_weeks_modular", ascending=False)

    wc = w[w["n_animals_observed"].ge(WELL_COVERED)].copy()
    icc = {"group_weeks": int(len(wc)), "units": int(wc["dynamic_social_unit"].nunique())}
    try:
        import statsmodels.formula.api as smf
        md = smf.mixedlm("modularity ~ n_animals_observed", wc,
                         groups=wc["dynamic_social_unit"]).fit(reml=True)
        vg, ve = float(md.cov_re.iloc[0, 0]), float(md.scale)
        icc.update({"between_unit_variance": round(vg, 6),
                    "residual_variance": round(ve, 6),
                    "icc_between_unit": round(vg / (vg + ve), 3),
                    "collar_pvalue": round(float(md.pvalues.get("n_animals_observed",
                                                                np.nan)), 4)})
    except Exception as exc:
        icc["error"] = str(exc)

    wc = wc.sort_values(["dynamic_social_unit", "period_start"])
    pers = []
    for unit, g in wc.groupby("dynamic_social_unit"):
        g = g.sort_values("period_start")
        gap = g["period_start"].diff().dt.days
        prev = g["modularity"].shift()
        ok = gap.eq(7) & prev.notna()
        if ok.sum() >= 20:
            pers.append({"dynamic_social_unit": unit,
                         "consecutive_week_pairs": int(ok.sum()),
                         "lag1_autocorrelation": round(float(np.corrcoef(
                             g.loc[ok, "modularity"], prev[ok])[0, 1]), 3)})

    return {"cohort": tag,
            "group_weeks": int(len(w)),
            "units": int(w["dynamic_social_unit"].nunique()),
            "window": [str(w["period_start"].min().date()),
                       str(w["period_start"].max().date())],
            "spearman_modularity_vs_collars": round(
                w["modularity"].corr(w["n_animals_observed"], method="spearman"), 3),
            "spearman_split_fraction_vs_collars": round(
                w["split_timestamp_fraction"].corr(w["n_animals_observed"],
                                                   method="spearman"), 3),
            "detectability_ladder": ladder,
            "coverage_matched": matched.to_dict("records"),
            "icc": icc,
            "persistence": pers}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--narrow", type=Path, default=NARROW)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Rebuilding weekly within-group metrics on the frozen cohort...")
    frozen = build_weekly(args.narrow)
    frozen.to_csv(args.output_dir / "weekly_network_metrics_frozen.csv", index=False)
    print("  -> %d group-weeks, %d units" % (len(frozen),
                                             frozen["dynamic_social_unit"].nunique()))

    legacy = pd.read_csv(LEGACY)
    legacy = legacy[legacy["scope"].eq("combined_1100_1600")]

    report = {
        "seed": SEED,
        "sample_hours_utc": list(SAMPLE_HOURS),
        "min_animals_for_a_group_week": MIN_ANIMALS,
        "frozen": findings(frozen, "frozen export (1,924,104 rows, 26 units, 2024-03-01+)"),
        "legacy": findings(legacy, "legacy source (1,703,133 rows, 21 units, 2025-01-01+)"),
    }

    f, l = report["frozen"], report["legacy"]
    print("\n" + "=" * 82)
    print("COHORT COMPARISON -- does each axis-B finding survive the frozen source?")
    print("=" * 82)
    print("\n%-34s %-28s %-28s" % ("", "FROZEN", "LEGACY"))
    print("%-34s %-28s %-28s" % ("group-weeks / units",
                                 "%d / %d" % (f["group_weeks"], f["units"]),
                                 "%d / %d" % (l["group_weeks"], l["units"])))
    print("%-34s %-28s %-28s" % ("window", " to ".join(f["window"]),
                                 " to ".join(l["window"])))
    print("%-34s %-28s %-28s" % ("Spearman(modularity, collars)",
                                 f["spearman_modularity_vs_collars"],
                                 l["spearman_modularity_vs_collars"]))
    print("%-34s %-28s %-28s" % ("Spearman(split frac, collars)",
                                 f["spearman_split_fraction_vs_collars"],
                                 l["spearman_split_fraction_vs_collars"]))
    print("%-34s %-28s %-28s" % ("ICC between-unit",
                                 f["icc"].get("icc_between_unit"),
                                 l["icc"].get("icc_between_unit")))
    print("%-34s %-28s %-28s" % ("  on group-weeks / units",
                                 "%s / %s" % (f["icc"]["group_weeks"], f["icc"]["units"]),
                                 "%s / %s" % (l["icc"]["group_weeks"], l["icc"]["units"])))
    print("%-34s %-28s %-28s" % ("  collar term p",
                                 f["icc"].get("collar_pvalue"),
                                 l["icc"].get("collar_pvalue")))

    print("\n--- detectability ladder ---")
    for k, blk in (("FROZEN", f), ("LEGACY", l)):
        print("\n%s:" % k)
        print(pd.DataFrame(blk["detectability_ladder"]).to_string(index=False))

    print("\n--- coverage-matched (%d-%d collars, >=5 group-weeks) ---" % MATCHED_BAND)
    for k, blk in (("FROZEN", f), ("LEGACY", l)):
        print("\n%s:" % k)
        print(pd.DataFrame(blk["coverage_matched"]).to_string(index=False))

    print("\n--- week-to-week persistence (>=20 consecutive pairs) ---")
    for k, blk in (("FROZEN", f), ("LEGACY", l)):
        print("\n%s:" % k)
        print(pd.DataFrame(blk["persistence"]).to_string(index=False)
              if blk["persistence"] else "  (none)")

    with open(args.output_dir / "cohort_comparison.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
