"""Temporal variation in social structure: stable group signatures, or wandering?

Figure 3 showed between-group composition varying at 4-9x a multinomial null. But it
pooled each group's events over the whole record, so that spread has two possible
sources and the figure cannot tell them apart:

  STABLE SIGNATURE   each group has a characteristic composition it holds over time
  TEMPORAL WANDER    each group moves around, and the pooled snapshots differ because
                     different groups were caught in different phases

This script separates them, using the same decomposition already applied to modularity
(where the answer was ICC 0.169 -- 85% of the variance was within group over time).

  1. VARIANCE DECOMPOSITION of each composition proportion across group-quarters:
     between-group against within-group-over-time, giving an ICC per proportion.
  2. PERSISTENCE: lag-1 autocorrelation of composition within group across consecutive
     quarters. A stable signature is autocorrelated; wandering is not.
  3. DO THE EXTREMES HOLD? For the groups Figure 3 singled out, composition quarter by
     quarter, so a reader can see whether the extreme is a signature or one phase.
  4. PARTNER TURNOVER: Jaccard similarity of a group's merge-partner set between
     consecutive quarters, against a null that redraws partners from the pool that was
     actually available in that quarter.
  5. RATE over time, with the group's own collar count in that quarter as an effort
     covariate -- within-group temporal comparison is exactly where the 121-fold effort
     ramp bites hardest.

Outputs: outputs/general_structure_2026_09/phase4f_temporal_variation/
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

EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/"
              "canonical_event_size_duration_all_events.csv")
WEEKLY = Path("outputs/canonical_within_group_density_modularity_ridges/"
              "canonical_within_group_weekly_network_metrics.csv")
OUT = Path("outputs/general_structure_2026_09/phase4f_temporal_variation")

FAMILIES = ["group_merge", "within_group", "individual"]
MERGE_TYPES = ["large_merge", "medium_partial_merge", "small_subset_merge"]
MIN_N = 20         # events needed for a group-quarter to carry a composition estimate
MIN_WINDOWS = 4    # quarters a group needs to enter the decomposition
N_NULL = 2000
SEED = 20260904


def group_quarter_table(ev: pd.DataFrame) -> pd.DataFrame:
    ev = ev.copy()
    ev["quarter"] = pd.to_datetime(ev["start_time"]).dt.to_period("Q").astype(str)
    rows = []
    ind = ev[ev["event_family"].isin(["individual", "within_group"])]
    rows.append(ind.assign(group=ind["origin_group"])[
        ["group", "quarter", "event_family", "event_type", "pair"]])
    mg = ev[ev["event_family"].eq("group_merge")]
    for col in ("group_a", "group_b"):
        other = "group_b" if col == "group_a" else "group_a"
        rows.append(mg.assign(group=mg[col], partner=mg[other])[
            ["group", "quarter", "event_family", "event_type", "pair", "partner"]])
    out = pd.concat(rows, ignore_index=True)
    return out[out["group"].notna()].copy()


def comp_by_window(tab: pd.DataFrame, keys: list[str], col: str) -> pd.DataFrame:
    ct = (tab.groupby(["group", "quarter", col]).size().unstack(fill_value=0)
          .reindex(columns=keys, fill_value=0))
    ct["n"] = ct.sum(axis=1)
    for k in keys:
        ct["p_" + k] = ct[k] / ct["n"]
    return ct.reset_index()


def logit(p, n):
    """Logit with a Haldane-Anscombe style correction so 0 and 1 are usable."""
    p = (p * n + 0.5) / (n + 1.0)
    return np.log(p / (1 - p))


def decompose(ct: pd.DataFrame, keys: list[str]) -> dict:
    """Between-group against within-group-over-time variance, per proportion."""
    use = ct[ct["n"] >= MIN_N].copy()
    counts = use.groupby("group")["quarter"].nunique()
    keep = counts[counts >= MIN_WINDOWS].index
    use = use[use["group"].isin(keep)]
    out = {"group_quarters": int(len(use)), "groups": int(use["group"].nunique()),
           "min_events_per_window": MIN_N, "min_windows_per_group": MIN_WINDOWS,
           "quarters": sorted(use["quarter"].unique().tolist()), "by_proportion": {}}
    try:
        import statsmodels.formula.api as smf
        for k in keys:
            f = use.copy()
            f["y"] = logit(f["p_" + k].to_numpy(), f["n"].to_numpy())
            md = smf.mixedlm("y ~ 1", f, groups=f["group"]).fit(reml=True)
            vg, ve = float(md.cov_re.iloc[0, 0]), float(md.scale)
            out["by_proportion"][k] = {
                "between_group_variance": round(vg, 4),
                "within_group_over_time_variance": round(ve, 4),
                "icc_between_group": round(vg / (vg + ve), 3),
                "interpretation": ("stable group signature"
                                   if vg / (vg + ve) >= 0.5 else
                                   "mostly temporal variation within group"),
            }
    except Exception as exc:  # pragma: no cover
        out["error"] = str(exc)
    return out


def anova_icc_and_permutation(ct: pd.DataFrame, keys: list[str],
                              rng: np.random.Generator) -> dict:
    """ICC by one-way ANOVA, plus a group-label permutation test.

    Two reasons this exists alongside the mixed model. It needs no optimiser, so it is
    immune to the boundary-convergence warnings a near-zero variance component
    produces. And the permutation is the RIGHT null for "do groups differ": it shuffles
    group labels across group-quarters, keeping the temporal variation intact. The
    multinomial null used for Figure 3 assumed each group drew from a fixed pooled
    distribution with no temporal structure, which is too weak -- it can call a spread
    significant when the spread is really quarter-to-quarter drift.
    """
    use = ct[ct["n"] >= MIN_N].copy()
    counts = use.groupby("group")["quarter"].nunique()
    use = use[use["group"].isin(counts[counts >= MIN_WINDOWS].index)]
    out = {"group_quarters": int(len(use)), "groups": int(use["group"].nunique()),
           "by_proportion": {}}
    for k in keys:
        y = use["p_" + k].to_numpy(dtype=float)
        g = use["group"].to_numpy()
        groups = np.unique(g)
        grand = y.mean()
        ns = np.array([(g == u).sum() for u in groups], dtype=float)
        means = np.array([y[g == u].mean() for u in groups])
        ssb = float((ns * (means - grand) ** 2).sum())
        dfb = len(groups) - 1
        ssw = float(sum(((y[g == u] - means[i]) ** 2).sum()
                        for i, u in enumerate(groups)))
        dfw = len(y) - len(groups)
        msb, msw = ssb / dfb, ssw / dfw
        n0 = (ns.sum() - (ns ** 2).sum() / ns.sum()) / dfb
        va = max(0.0, (msb - msw) / n0)
        icc = va / (va + msw) if (va + msw) > 0 else 0.0
        perm = np.empty(N_NULL)
        for b in range(N_NULL):
            gp = rng.permutation(g)
            m2 = np.array([y[gp == u].mean() for u in groups])
            n2 = np.array([(gp == u).sum() for u in groups], dtype=float)
            perm[b] = (n2 * (m2 - grand) ** 2).sum()
        out["by_proportion"][k] = {
            "icc_anova": round(icc, 3),
            "F": round(msb / msw, 2) if msw > 0 else None,
            "permutation_p": round(float((perm >= ssb).mean()), 4),
            "verdict": ("groups differ" if (perm >= ssb).mean() < 0.05
                        else "not distinguishable from shuffling group labels"),
        }
    return out


def within_group_swing(ct: pd.DataFrame, key: str) -> pd.DataFrame:
    """How far one group's own proportion travels across its observed quarters."""
    use = ct[ct["n"] >= MIN_N].copy()
    use["qi"] = pd.PeriodIndex(use["quarter"], freq="Q").astype("int64")
    rows = []
    for g, sub in use.groupby("group"):
        if len(sub) < 2:
            continue
        rows.append({"group": g, "windows": int(len(sub)),
                     "total_events": int(sub["n"].sum()),
                     "min": round(float(sub[key].min()), 3),
                     "max": round(float(sub[key].max()), 3),
                     "swing": round(float(sub[key].max() - sub[key].min()), 3)})
    return pd.DataFrame(rows).sort_values("swing", ascending=False)


def persistence(ct: pd.DataFrame, keys: list[str]) -> dict:
    """Lag-1 autocorrelation of composition within group, consecutive quarters only."""
    use = ct[ct["n"] >= MIN_N].copy()
    use["qi"] = pd.PeriodIndex(use["quarter"], freq="Q").astype("int64")
    out = {}
    for k in keys:
        pairs = []
        for g, sub in use.groupby("group"):
            sub = sub.sort_values("qi")
            prev_q = sub["qi"].shift()
            prev_p = sub["p_" + k].shift()
            ok = (sub["qi"] - prev_q).eq(1) & prev_p.notna()
            for a, b in zip(sub.loc[ok, "p_" + k], prev_p[ok]):
                pairs.append((a, b))
        if len(pairs) >= 8:
            a = np.array([p[0] for p in pairs])
            b = np.array([p[1] for p in pairs])
            out[k] = {"consecutive_quarter_pairs": len(pairs),
                      "lag1_autocorrelation": round(float(np.corrcoef(a, b)[0, 1]), 3)}
        else:
            out[k] = {"consecutive_quarter_pairs": len(pairs),
                      "lag1_autocorrelation": None}
    return out


def partner_turnover(tab: pd.DataFrame, rng: np.random.Generator) -> dict:
    """Jaccard of a group's merge-partner set between consecutive quarters."""
    mg = tab[tab["event_family"].eq("group_merge") & tab["partner"].notna()]
    sets = (mg.groupby(["group", "quarter"])["partner"]
              .apply(lambda s: frozenset(s.unique())).reset_index())
    sets["qi"] = pd.PeriodIndex(sets["quarter"], freq="Q").astype("int64")
    # the pool of partners any group merged with, per quarter
    pool = mg.groupby("quarter")["partner"].apply(lambda s: sorted(set(s))).to_dict()

    obs, nulls, rows = [], [], []
    for g, sub in sets.groupby("group"):
        sub = sub.sort_values("qi")
        for (q1, s1, i1), (q2, s2, i2) in zip(
                zip(sub["quarter"], sub["partner"], sub["qi"]),
                list(zip(sub["quarter"], sub["partner"], sub["qi"]))[1:]):
            if i2 - i1 != 1 or not s1 or not s2:
                continue
            j = len(s1 & s2) / len(s1 | s2)
            obs.append(j)
            rows.append({"group": g, "from": q1, "to": q2, "n_from": len(s1),
                         "n_to": len(s2), "shared": len(s1 & s2),
                         "jaccard": round(j, 3)})
            avail = [p for p in pool.get(q2, []) if p != g]
            if len(avail) >= len(s2):
                draws = []
                for _ in range(200):
                    pick = frozenset(rng.choice(avail, len(s2), replace=False))
                    draws.append(len(s1 & pick) / len(s1 | pick))
                nulls.append(float(np.mean(draws)))
    res = {"consecutive_quarter_pairs": len(obs)}
    if obs:
        res.update({
            "mean_jaccard": round(float(np.mean(obs)), 3),
            "median_jaccard": round(float(np.median(obs)), 3),
            "null_mean_jaccard": round(float(np.mean(nulls)), 3) if nulls else None,
            "ratio_to_null": (round(float(np.mean(obs) / np.mean(nulls)), 2)
                              if nulls and np.mean(nulls) > 0 else None),
        })
    return res, pd.DataFrame(rows)


def rate_over_time(tab: pd.DataFrame, weekly: pd.DataFrame) -> dict:
    """Events per observed group-week by quarter, with collar count as effort."""
    w = weekly[weekly["scope"].eq("combined_1100_1600")].copy()
    w["quarter"] = pd.to_datetime(w["period_start"]).dt.to_period("Q").astype(str)
    exp = w.groupby(["dynamic_social_unit", "quarter"]).agg(
        group_weeks=("modularity", "size"),
        median_collars=("n_animals_observed", "median")).reset_index().rename(
        columns={"dynamic_social_unit": "group"})
    cnt = tab.groupby(["group", "quarter"]).size().rename("events").reset_index()
    m = cnt.merge(exp, on=["group", "quarter"], how="inner")
    m = m[m["group_weeks"] >= 4]
    m["rate"] = m["events"] / m["group_weeks"]
    res = {"group_quarters": int(len(m)), "groups": int(m["group"].nunique())}
    try:
        import statsmodels.formula.api as smf
        m["log_rate"] = np.log(m["rate"].clip(lower=0.1))
        md = smf.mixedlm("log_rate ~ median_collars", m, groups=m["group"]).fit(reml=True)
        vg, ve = float(md.cov_re.iloc[0, 0]), float(md.scale)
        res.update({
            "model": "log(events per group-week) ~ collars + (1|group)",
            "between_group_variance": round(vg, 4),
            "within_group_over_time_variance": round(ve, 4),
            "icc_between_group": round(vg / (vg + ve), 3),
            "collar_coefficient": round(float(md.params.get("median_collars", np.nan)), 4),
            "collar_pvalue": round(float(md.pvalues.get("median_collars", np.nan)), 4),
            "spearman_rate_vs_collars": round(
                float(m["rate"].corr(m["median_collars"], method="spearman")), 3),
        })
    except Exception as exc:  # pragma: no cover
        res["error"] = str(exc)
    return res, m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    ev = pd.read_csv(EVENTS, parse_dates=["start_time"])
    weekly = pd.read_csv(WEEKLY)
    tab = group_quarter_table(ev)

    fam = comp_by_window(tab, FAMILIES, "event_family")
    mtab = tab[tab["event_family"].eq("group_merge")]
    mrg = comp_by_window(mtab, MERGE_TYPES, "event_type")
    fam.to_csv(args.output_dir / "composition_by_group_quarter_family.csv", index=False)
    mrg.to_csv(args.output_dir / "composition_by_group_quarter_merge.csv", index=False)

    turn, turn_rows = partner_turnover(tab, rng)
    turn_rows.to_csv(args.output_dir / "partner_turnover.csv", index=False)
    rates, rate_rows = rate_over_time(tab, weekly)
    rate_rows.to_csv(args.output_dir / "rate_by_group_quarter.csv", index=False)

    report = {
        "seed": SEED,
        "window": "calendar quarter",
        "source_caveat": "legacy hourly source filtered to 2025-01-01 onward; six "
                         "quarters 2025Q1 to 2026Q2; gap 7.3 applies",
        "family": {"decomposition": decompose(fam, FAMILIES),
                   "anova_permutation": anova_icc_and_permutation(fam, FAMILIES, rng),
                   "persistence": persistence(fam, FAMILIES)},
        "merge_scale": {"decomposition": decompose(mrg, MERGE_TYPES),
                        "anova_permutation": anova_icc_and_permutation(mrg, MERGE_TYPES, rng),
                        "persistence": persistence(mrg, MERGE_TYPES)},
        "partner_turnover": turn,
        "event_rate": rates,
    }
    with open(args.output_dir / "temporal_variation_report.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    # ---------------------------------------------------------------- print
    print("=" * 86)
    print("TEMPORAL VARIATION IN SOCIAL STRUCTURE")
    print("=" * 86)
    for name, keys in (("family", FAMILIES), ("merge_scale", MERGE_TYPES)):
        d = report[name]["decomposition"]
        p = report[name]["persistence"]
        print("\n--- %s: %d group-quarters, %d groups, quarters %s ---"
              % (name, d["group_quarters"], d["groups"],
                 " ".join(d["quarters"])))
        print("  %-24s %10s %10s %8s   %-34s %s"
              % ("proportion", "between", "within", "ICC", "reading", "lag-1 r"))
        for k in keys:
            b = d["by_proportion"].get(k, {})
            pr = p.get(k, {})
            print("  %-24s %10.4f %10.4f %8.3f   %-34s %s"
                  % (k, b.get("between_group_variance", float("nan")),
                     b.get("within_group_over_time_variance", float("nan")),
                     b.get("icc_between_group", float("nan")),
                     b.get("interpretation", "?"),
                     pr.get("lag1_autocorrelation")))
    print("\n--- partner turnover between consecutive quarters ---")
    print("  " + json.dumps(report["partner_turnover"]))
    print("\n--- event rate per group-week ---")
    print("  " + json.dumps(report["event_rate"]))
    swing = within_group_swing(mrg, "p_large_merge")
    swing.to_csv(args.output_dir / "within_group_swing_p_large_merge.csv", index=False)
    report["within_group_swing_p_large_merge"] = swing.to_dict("records")
    with open(args.output_dir / "temporal_variation_report.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("\n--- ANOVA ICC and group-label permutation (the correct null) ---")
    print("    the permutation keeps temporal variation and shuffles only group labels,")
    print("    so it asks whether groups differ once quarter-to-quarter drift is allowed")
    for name, keys in (("family", FAMILIES), ("merge_scale", MERGE_TYPES)):
        a = report[name]["anova_permutation"]
        print("  %s (%d group-quarters, %d groups)"
              % (name, a["group_quarters"], a["groups"]))
        for k in keys:
            v = a["by_proportion"][k]
            print("      %-24s ICC %.3f  F %5.2f  perm p %.4f   %s"
                  % (k, v["icc_anova"], v["F"], v["permutation_p"], v["verdict"]))

    print("\n--- within-group swing in p_large_merge across its own quarters ---")
    print(swing.head(12).to_string(index=False))

    print("\n--- do the extremes hold? composition by quarter ---")
    for g in ("PhantomWest", "FireOpal", "Jade", "LapisSplinter"):
        sub = mrg[mrg["group"].eq(g) & mrg["n"].ge(MIN_N)]
        if len(sub):
            print("  %-14s merge-scale p_large by quarter: %s"
                  % (g, "  ".join("%s %.2f (n=%d)" % (r.quarter, r.p_large_merge, r.n)
                                  for r in sub.itertuples())))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
