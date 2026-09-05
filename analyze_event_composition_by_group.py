"""Event-type composition per group, and whether groups genuinely differ.

Three questions:

  1. PROPORTIONS. What mix of event types does each group produce -- between-group
     merges, within-group splits, individual separations -- and how do the merges
     themselves break down by scale?
  2. VARIATION. Do groups differ in that mix by more than sampling would produce?
     Tested against a multinomial null in which every group draws from the population
     type distribution in proportion to its own observed exposure.
  3. CONFOUND. Does composition track collar coverage? A group with more collars sees
     more of everything, but not necessarily equally: a split needs two clusters to be
     visible while a merge needs only one animal from each side, so detection is
     type-dependent and could manufacture composition differences on its own.

Attribution: merges are dyadic and count for BOTH participating groups, so merge counts
sum to twice the event total; splits and individual separations are attributed to the
origin group. Exposure is observed group-weeks from the weekly network metrics on the
same source as the events.

Outputs: outputs/general_structure_2026_09/phase4e_event_composition/
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/"
              "canonical_event_size_duration_all_events.csv")
WEEKLY_LEGACY = Path("outputs/canonical_within_group_density_modularity_ridges/"
                     "canonical_within_group_weekly_network_metrics.csv")
OUT = Path("outputs/general_structure_2026_09/phase4e_event_composition")

FAMILIES = ["group_merge", "within_group", "individual"]
MERGE_TYPES = ["large_merge", "medium_partial_merge", "small_subset_merge"]
MIN_EVENTS = 30          # a group needs this many to carry a composition estimate
N_BOOT = 2000
N_NULL = 2000
SEED = 20260904


def group_event_table(ev: pd.DataFrame) -> pd.DataFrame:
    """One row per (group, event), with merges duplicated across both participants."""
    rows = []
    ind = ev[ev["event_family"].isin(["individual", "within_group"])]
    rows.append(ind.assign(group=ind["origin_group"])[
        ["group", "event_family", "event_type", "duration_hours"]])
    mg = ev[ev["event_family"].eq("group_merge")]
    for col in ("group_a", "group_b"):
        rows.append(mg.assign(group=mg[col])[
            ["group", "event_family", "event_type", "duration_hours"]])
    out = pd.concat(rows, ignore_index=True)
    return out[out["group"].notna()].copy()


def exposure(weekly: pd.DataFrame) -> pd.DataFrame:
    w = weekly[weekly["scope"].eq("combined_1100_1600")]
    return w.groupby("dynamic_social_unit").agg(
        group_weeks=("modularity", "size"),
        median_collars=("n_animals_observed", "median"),
        nominal_size=("origin_group_total_size", "first"),
    ).reset_index().rename(columns={"dynamic_social_unit": "group"})


def composition(tab: pd.DataFrame, keys: list[str], colname: str) -> pd.DataFrame:
    ct = (tab.groupby(["group", colname]).size().unstack(fill_value=0)
          .reindex(columns=keys, fill_value=0))
    ct["total"] = ct.sum(axis=1)
    for k in keys:
        ct["p_" + k] = ct[k] / ct["total"]
    return ct.reset_index()


def boot_props(tab: pd.DataFrame, keys: list[str], colname: str,
               rng: np.random.Generator) -> dict:
    """Percentile intervals on each group's proportions, resampling that group's events."""
    out = {}
    for g, sub in tab.groupby("group"):
        n = len(sub)
        if n < MIN_EVENTS:
            continue
        codes = pd.Categorical(sub[colname], categories=keys).codes
        codes = codes[codes >= 0]
        draws = np.empty((N_BOOT, len(keys)))
        for b in range(N_BOOT):
            pick = rng.integers(0, len(codes), len(codes))
            draws[b] = np.bincount(codes[pick], minlength=len(keys)) / len(codes)
        out[g] = {keys[i]: {"point": float((codes == i).mean()),
                            "lo": float(np.percentile(draws[:, i], 2.5)),
                            "hi": float(np.percentile(draws[:, i], 97.5))}
                  for i in range(len(keys))}
        out[g]["n"] = int(n)
    return out


def variation_vs_null(ct: pd.DataFrame, keys: list[str],
                      rng: np.random.Generator) -> dict:
    """Observed between-group spread in composition against a multinomial null.

    Under the null every group draws its own number of events from the pooled type
    distribution, so any spread is sampling noise. The statistic is the mean absolute
    deviation of a group's proportion from the pooled proportion, averaged over types
    and weighted by group event count.
    """
    use = ct[ct["total"] >= MIN_EVENTS].copy()
    counts = use[keys].to_numpy(dtype=float)
    totals = counts.sum(axis=1)
    pooled = counts.sum(axis=0) / counts.sum()

    def stat(mat):
        props = mat / mat.sum(axis=1, keepdims=True)
        return float(np.average(np.abs(props - pooled).mean(axis=1), weights=totals))

    obs = stat(counts)
    null = np.empty(N_NULL)
    for b in range(N_NULL):
        sim = np.array([rng.multinomial(int(t), pooled) for t in totals], dtype=float)
        null[b] = stat(sim)
    return {
        "groups_used": int(len(use)),
        "min_events_per_group": MIN_EVENTS,
        "pooled_proportions": {k: round(float(p), 4) for k, p in zip(keys, pooled)},
        "observed_statistic": round(obs, 5),
        "null_mean": round(float(null.mean()), 5),
        "null_95th": round(float(np.percentile(null, 95)), 5),
        "ratio_observed_to_null": round(obs / max(1e-9, null.mean()), 2),
        "p_value": round(float((null >= obs).mean()), 5),
    }


def coverage_confound(ct: pd.DataFrame, exp: pd.DataFrame, keys: list[str]) -> dict:
    # `ct` may already carry the exposure columns from an earlier merge; only join
    # them when they are missing, otherwise pandas suffixes them and the lookup fails.
    m = ct if "median_collars" in ct.columns else ct.merge(exp, on="group", how="left")
    m = m[(m["total"] >= MIN_EVENTS) & m["median_collars"].notna()]
    out = {"groups": int(len(m)), "spearman_with_median_collars": {}}
    for k in keys:
        out["spearman_with_median_collars"][k] = round(
            float(m["p_" + k].corr(m["median_collars"], method="spearman")), 3)
    out["spearman_total_events_vs_collars"] = round(
        float(m["total"].corr(m["median_collars"], method="spearman")), 3)
    out["spearman_total_events_vs_group_weeks"] = round(
        float(m["total"].corr(m["group_weeks"], method="spearman")), 3)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    ev = pd.read_csv(EVENTS)
    weekly = pd.read_csv(WEEKLY_LEGACY)
    exp = exposure(weekly)
    tab = group_event_table(ev)

    fam = composition(tab, FAMILIES, "event_family").merge(exp, on="group", how="left")
    fam["events_per_group_week"] = fam["total"] / fam["group_weeks"]
    fam = fam.sort_values("total", ascending=False)
    fam.to_csv(args.output_dir / "composition_by_family.csv", index=False)

    mtab = tab[tab["event_family"].eq("group_merge")]
    mrg = composition(mtab, MERGE_TYPES, "event_type").merge(exp, on="group", how="left")
    mrg = mrg.sort_values("total", ascending=False)
    mrg.to_csv(args.output_dir / "composition_by_merge_scale.csv", index=False)

    report = {
        "seed": SEED,
        "attribution": "merges counted for both participating groups; splits and "
                       "individual separations attributed to origin_group",
        "source_caveat": "events and exposure both from the legacy source filtered to "
                         "2025-01-01 onward; gap 7.3 applies",
        "totals": {
            "events": int(len(ev)),
            "group_event_rows_after_merge_duplication": int(len(tab)),
            "groups": int(tab["group"].nunique()),
        },
        "family": {
            "variation_vs_null": variation_vs_null(fam, FAMILIES, rng),
            "coverage_confound": coverage_confound(fam, exp, FAMILIES),
            "bootstrap": boot_props(tab, FAMILIES, "event_family", rng),
        },
        "merge_scale": {
            "variation_vs_null": variation_vs_null(mrg, MERGE_TYPES, rng),
            "coverage_confound": coverage_confound(mrg, exp, MERGE_TYPES),
            "bootstrap": boot_props(mtab, MERGE_TYPES, "event_type", rng),
        },
    }

    # THE CONTROL. The multinomial null assumes every group detects each event type
    # with the same probability, and coverage violates that: a within-group split needs
    # at least two clusters visible, a merge needs only one animal from each side. So
    # repeat the test inside collar strata, where detection is comparable.
    strata = [("all groups", 0, 99), ("collars >= 8", 8, 99), ("collars 7-14", 7, 14)]
    for name, ct, keys in (("family", fam, FAMILIES),
                           ("merge_scale", mrg, MERGE_TYPES)):
        ladder = []
        for lab, lo, hi in strata:
            sub = ct[ct["median_collars"].between(lo, hi)]
            if (sub["total"] >= MIN_EVENTS).sum() < 4:
                ladder.append({"stratum": lab, "skipped": "fewer than 4 usable groups"})
                continue
            v = variation_vs_null(sub, keys, rng)
            ladder.append({"stratum": lab, "groups": v["groups_used"],
                           "observed": v["observed_statistic"],
                           "null_mean": v["null_mean"],
                           "ratio": v["ratio_observed_to_null"],
                           "p": v["p_value"]})
        report[name]["coverage_stratified"] = ladder
        pd.DataFrame(ladder).to_csv(
            args.output_dir / ("variation_stratified_%s.csv" % name), index=False)
    with open(args.output_dir / "event_composition_report.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    # ---------------------------------------------------------------- print
    print("=" * 84)
    print("EVENT COMPOSITION BY GROUP")
    print("=" * 84)
    show = ["group", "total", "median_collars", "nominal_size", "group_weeks",
            "events_per_group_week", "p_group_merge", "p_within_group", "p_individual"]
    print("\n--- family composition (groups with >= %d events) ---" % MIN_EVENTS)
    print(fam.loc[fam["total"].ge(MIN_EVENTS), show].round(3).to_string(index=False))

    print("\n--- merge-scale composition within merges ---")
    show_m = ["group", "total", "median_collars", "p_large_merge",
              "p_medium_partial_merge", "p_small_subset_merge"]
    print(mrg.loc[mrg["total"].ge(MIN_EVENTS), show_m].round(3).to_string(index=False))

    for name in ("family", "merge_scale"):
        v = report[name]["variation_vs_null"]
        c = report[name]["coverage_confound"]
        print("\n--- %s: is between-group variation more than sampling? ---" % name)
        print("  pooled proportions: %s" % json.dumps(v["pooled_proportions"]))
        print("  observed spread %.5f vs null mean %.5f (95th %.5f) -> %.1fx, p = %.4f"
              % (v["observed_statistic"], v["null_mean"], v["null_95th"],
                 v["ratio_observed_to_null"], v["p_value"]))
        print("  coverage confound, Spearman of each proportion with median collars:")
        for k, r in c["spearman_with_median_collars"].items():
            print("      %-24s %+0.3f" % (k, r))
        print("      total events vs collars      %+0.3f" % c["spearman_total_events_vs_collars"])
        print("      total events vs group-weeks  %+0.3f" % c["spearman_total_events_vs_group_weeks"])

    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
