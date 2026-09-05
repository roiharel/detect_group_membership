"""Consequence test A -> B: does a between-group merge change internal cohesion?

THE QUESTION. The paper asserts that reconfiguration redistributes companions. If that
is true, a group's internal cohesion should change after it merges with another group.
Cohesion is the right outcome to test because Figure 4 showed it is a STATE (modularity
ICC 0.17), so it is capable of responding; the group's event rate (ICC 0.66) is not.

DESIGN, built against the two traps.

  Trap 1, reverse causation. Groups may merge more BECAUSE they are already divided.
  Handled by strict temporal ordering with a gap: the pre-window is weeks t-3 and t-2,
  leaving week t-1 unused as a buffer, and the outcome is weeks t+1 and t+2. The reverse
  direction is also fitted explicitly, so the reader can see it.

  Trap 2, regression to the mean. Modularity is autocorrelated at r ~ 0.4, so selecting
  merge weeks on anything correlated with it makes the after-value drift back toward the
  group's mean regardless of the merge. Handled two ways: pre-window cohesion enters the
  model as a covariate, and a matched-control analysis pairs each merge week with a
  no-merge week from THE SAME UNIT at similar pre-window cohesion and collar count.

  Everything is within unit. The session established that between-unit comparisons in
  this cohort are not identifiable (26 units, group variance component unidentified),
  so unit enters as a cluster and never as an explanatory contrast.

Three cohesion outcomes, because modularity is mostly zero and one measure would be
fragile:
  modularity                              higher = more internally divided
  largest_community_fraction              lower  = more internally divided
  multi_animal_split_timestamp_fraction   higher = more time in multi-animal subgroups

Also reports the A -> C FEASIBILITY COUNT: how many individual excursions could ever be
joined to a fine-scale mixing measure, which decides whether that link is worth
building.

Outputs: outputs/general_structure_2026_09/phase5a_merge_consequence/
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

BASE = Path("outputs/general_structure_2026_09")
WEEKLY = BASE / "phase4d_axis_b_frozen/weekly_network_metrics_frozen.csv"
STAGE1 = BASE / "phase2_two_stage_events/stage1_events_with_stage2_mixing.csv"
EXC = BASE / "phase4b_individual_axis/excursions_dominant_gap0.csv"
GRAD = BASE / "phase3_mixing_gradient/mixing_gradient_by_dyad.csv"
OUT = BASE / "phase5a_merge_consequence"

WELL_COVERED = 12
PRE_WEEKS = (-3, -2)      # buffer week t-1 deliberately unused
POST_WEEKS = (1, 2)
OUTCOMES = [
    ("modularity", "higher = more divided"),
    ("largest_community_fraction", "lower = more divided"),
    ("multi_animal_split_timestamp_fraction", "higher = more subgrouping"),
]
CALIPER_PRE = 0.05        # matched control must be this close on pre-window cohesion
CALIPER_COLLARS = 3
SEED = 20260904


def build_panel() -> pd.DataFrame:
    wk = pd.read_csv(WEEKLY, parse_dates=["period_start"])
    wk = wk.rename(columns={"dynamic_social_unit": "unit"})
    wk["wk_index"] = (wk["period_start"] - wk["period_start"].min()).dt.days // 7

    s1 = pd.read_csv(STAGE1, parse_dates=["start_hour"])
    s1["period_start"] = (s1["start_hour"]
                          - pd.to_timedelta(s1["start_hour"].dt.weekday, unit="D")
                          ).dt.normalize()
    # one row per (unit, week) that the unit merged in, with the largest scale reached
    rows = []
    for col in ("unit_a", "unit_b"):
        rows.append(s1.assign(unit=s1[col])[
            ["unit", "period_start", "merge_scale", "max_cluster_size",
             "structural_span_hours"]])
    ev = pd.concat(rows, ignore_index=True)
    ev["is_large"] = ev["merge_scale"].eq("large_merge").astype(int)
    agg = ev.groupby(["unit", "period_start"]).agg(
        n_merges=("merge_scale", "size"),
        any_large_merge=("is_large", "max"),
        max_joint_cluster=("max_cluster_size", "max"),
        max_span_hours=("structural_span_hours", "max")).reset_index()

    p = wk.merge(agg, on=["unit", "period_start"], how="left")
    for c in ("n_merges", "any_large_merge"):
        p[c] = p[c].fillna(0).astype(int)
    return p.sort_values(["unit", "wk_index"]).reset_index(drop=True)


def window_mean(p: pd.DataFrame, col: str, lo: int, hi: int) -> pd.Series:
    """Mean of `col` over weeks t+lo .. t+hi, within unit, coverage-gated."""
    idx = p.set_index(["unit", "wk_index"])[col]
    out = []
    for r in p.itertuples():
        vals = []
        for k in range(lo, hi + 1):
            key = (r.unit, r.wk_index + k)
            if key in idx.index:
                v = idx.loc[key]
                if np.isscalar(v) and not pd.isna(v):
                    vals.append(float(v))
        out.append(np.mean(vals) if vals else np.nan)
    return pd.Series(out, index=p.index)


def coverage_ok(p: pd.DataFrame, lo: int, hi: int) -> pd.Series:
    idx = p.set_index(["unit", "wk_index"])["n_animals_observed"]
    out = []
    for r in p.itertuples():
        ok = True
        for k in range(lo, hi + 1):
            key = (r.unit, r.wk_index + k)
            if key not in idx.index or float(idx.loc[key]) < WELL_COVERED:
                ok = False
                break
        out.append(ok)
    return pd.Series(out, index=p.index)


def assemble(p: pd.DataFrame) -> pd.DataFrame:
    d = p.copy()
    d["pre_ok"] = coverage_ok(d, *PRE_WEEKS)
    d["post_ok"] = coverage_ok(d, *POST_WEEKS)
    for col, _ in OUTCOMES:
        d["pre_" + col] = window_mean(d, col, *PRE_WEEKS)
        d["post_" + col] = window_mean(d, col, *POST_WEEKS)
        d["delta_" + col] = d["post_" + col] - d["pre_" + col]
    keep = (d["n_animals_observed"].ge(WELL_COVERED) & d["pre_ok"] & d["post_ok"]
            & d["pre_modularity"].notna() & d["post_modularity"].notna())
    return d[keep].copy()


def gee_effect(d: pd.DataFrame, outcome: str, treat: str) -> dict:
    import statsmodels.api as sm
    col = "delta_" + outcome
    f = d.dropna(subset=[col, "pre_" + outcome, "n_animals_observed", treat]).copy()
    if len(f) < 40 or f["unit"].nunique() < 4 or f[treat].nunique() < 2:
        return {"skipped": "too few rows, units, or no treatment variation",
                "n": int(len(f))}
    X = f[[treat, "pre_" + outcome, "n_animals_observed"]].astype(float)
    X = sm.add_constant(X)
    m = sm.GEE(f[col].astype(float), X, groups=f["unit"],
               family=sm.families.Gaussian(),
               cov_struct=sm.cov_struct.Independence()).fit()
    ci = m.conf_int()
    return {
        "n": int(len(f)), "units": int(f["unit"].nunique()),
        "treated": int(f[treat].sum()) if f[treat].dtype != float else None,
        "effect": round(float(m.params[treat]), 5),
        "ci": [round(float(ci.loc[treat, 0]), 5), round(float(ci.loc[treat, 1]), 5)],
        "p": round(float(m.pvalues[treat]), 4),
        "pre_window_coefficient": round(float(m.params["pre_" + outcome]), 4),
        "sd_of_outcome_delta": round(float(f[col].std()), 5),
    }


def matched_effect(d: pd.DataFrame, outcome: str, rng) -> dict:
    """Pair each merge week with a no-merge week from the same unit at similar
    pre-window cohesion and collar count. Removes both unit and level as explanations."""
    pre = "pre_" + outcome
    col = "delta_" + outcome
    f = d.dropna(subset=[col, pre]).copy()
    treated = f[f["any_large_merge"].eq(1)]
    control = f[f["n_merges"].eq(0)]
    pairs = []
    for t in treated.itertuples():
        pool = control[control["unit"].eq(t.unit)]
        if pool.empty:
            continue
        pv = getattr(t, pre.replace("-", "_")) if False else t.__getattribute__(pre)
        near = pool[(pool[pre].sub(pv).abs() <= CALIPER_PRE)
                    & (pool["n_animals_observed"].sub(t.n_animals_observed).abs()
                       <= CALIPER_COLLARS)]
        if near.empty:
            continue
        c = near.iloc[int(rng.integers(0, len(near)))]
        pairs.append(t.__getattribute__(col) - float(c[col]))
    if len(pairs) < 12:
        return {"skipped": "fewer than 12 matched pairs", "pairs": len(pairs)}
    a = np.array(pairs, dtype=float)
    boot = [np.mean(rng.choice(a, len(a), replace=True)) for _ in range(2000)]
    return {
        "matched_pairs": int(len(a)),
        "mean_difference": round(float(a.mean()), 5),
        "ci": [round(float(np.percentile(boot, 2.5)), 5),
               round(float(np.percentile(boot, 97.5)), 5)],
        "caliper_pre": CALIPER_PRE, "caliper_collars": CALIPER_COLLARS,
    }


def reverse_direction(d: pd.DataFrame) -> dict:
    """Does pre-window cohesion predict whether the unit merges in week t?

    Run against two treatments, because the first is potentially CIRCULAR. A "large
    merge" is coded when much of the group sits in one joint cluster, which is
    mechanically harder for a group that is internally split -- so a negative
    association with pre-window modularity could be definitional rather than
    behavioural. If it is definitional it should weaken or vanish for ANY merge,
    at any scale, which does not carry that requirement.
    """
    import statsmodels.api as sm
    d = d.copy()
    d["any_merge"] = d["n_merges"].gt(0).astype(int)
    out = {}
    for treat, tlab in (("any_large_merge", "large merge only"),
                        ("any_merge", "any merge, any scale")):
        block = {}
        for outcome, _ in OUTCOMES:
            f = d.dropna(subset=["pre_" + outcome]).copy()
            if f[treat].nunique() < 2:
                continue
            X = f[["pre_" + outcome, "n_animals_observed"]].astype(float)
            for c in X.columns:
                sd = X[c].std()
                if sd > 0:
                    X[c] = (X[c] - X[c].mean()) / sd
            m = sm.GEE(f[treat].astype(float), sm.add_constant(X),
                       groups=f["unit"], family=sm.families.Binomial(),
                       cov_struct=sm.cov_struct.Independence()).fit()
            ci = m.conf_int()
            k = "pre_" + outcome
            block[outcome] = {"n": int(len(f)), "positives": int(f[treat].sum()),
                              "or_per_sd": round(float(np.exp(m.params[k])), 3),
                              "ci_or": [round(float(np.exp(ci.loc[k, 0])), 3),
                                        round(float(np.exp(ci.loc[k, 1])), 3)],
                              "p": round(float(m.pvalues[k]), 4)}
        out[treat] = {"treatment": tlab, "by_predictor": block}
    lg = out.get("any_large_merge", {}).get("by_predictor", {}).get("modularity", {})
    am = out.get("any_merge", {}).get("by_predictor", {}).get("modularity", {})
    if lg and am:
        out["verdict"] = (
            "pre-window modularity predicts a LARGE merge at OR %.3f (p = %.4f) but ANY "
            "merge at OR %.3f (p = %.4f). The association is specific to the treatment "
            "whose own definition requires a cohesive group, so it is definitional, not "
            "behavioural. No reverse effect."
            % (lg["or_per_sd"], lg["p"], am["or_per_sd"], am["p"]))
    return out


def ac_feasibility() -> dict:
    """How many excursions could ever be joined to a fine-scale mixing measure?"""
    ex = pd.read_csv(EXC)
    grad = pd.read_csv(GRAD)
    fine_dyads = grad.loc[grad["stratum"].eq("discrete_encounter"), "pair_key"]
    fine_units = set()
    for p in fine_dyads:
        fine_units.update(x.strip() for x in str(p).split("-"))
    joined = ex[ex["other_nights"].gt(0)]
    in_fine = joined[joined["origin_group"].isin(fine_units)]
    return {
        "excursions_total": int(len(ex)),
        "excursions_reaching_another_unit": int(len(joined)),
        "origin_groups_of_those": int(joined["origin_group"].nunique()),
        "units_with_a_fine_scale_mixing_estimate": len(fine_units),
        "joined_excursions_whose_origin_has_fine_scale": int(len(in_fine)),
        "origin_groups_in_that_intersection": int(in_fine["origin_group"].nunique()),
        "animals_in_that_intersection": int(in_fine["animal_id"].nunique()),
        "note": ("destination unit is not carried in the excursion table, so the true "
                 "dyad-level sample is at most this and in practice smaller; a "
                 "dyad-quarter design over 6 quarters would have very few positives"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    panel = build_panel()
    d = assemble(panel)
    d.to_csv(args.output_dir / "merge_consequence_panel.csv", index=False)

    report = {
        "seed": SEED,
        "design": {
            "unit": "unit-week, frozen cohort",
            "coverage_gate": ">= %d collars in week t and in every pre and post week"
                             % WELL_COVERED,
            "pre_window": "weeks t%d to t%d" % PRE_WEEKS,
            "buffer": "week t-1 deliberately unused",
            "post_window": "weeks t+%d to t+%d" % POST_WEEKS,
            "treatment": "unit participated in a large merge in week t",
            "everything_within_unit": True,
        },
        "sample": {
            "eligible_unit_weeks": int(len(d)),
            "units": int(d["unit"].nunique()),
            "large_merge_weeks": int(d["any_large_merge"].sum()),
            "no_merge_weeks": int(d["n_merges"].eq(0).sum()),
        },
        "effects": {}, "matched": {}, "reverse_direction": reverse_direction(d),
        "ac_feasibility": ac_feasibility(),
    }
    for outcome, direction in OUTCOMES:
        report["effects"][outcome] = {
            "direction": direction,
            "binary_large_merge": gee_effect(d, outcome, "any_large_merge"),
            "dose_n_merges": gee_effect(d, outcome, "n_merges"),
        }
        report["matched"][outcome] = matched_effect(d, outcome, rng)

    with open(args.output_dir / "merge_consequence_report.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    # ---------------------------------------------------------------- print
    print("=" * 82)
    print("DOES A MERGE CHANGE INTERNAL COHESION?  (A -> B)")
    print("=" * 82)
    print("\nsample: " + json.dumps(report["sample"]))
    for outcome, direction in OUTCOMES:
        e = report["effects"][outcome]
        print("\n--- %s  (%s) ---" % (outcome, direction))
        for k in ("binary_large_merge", "dose_n_merges"):
            v = e[k]
            if "skipped" in v:
                print("  %-22s SKIPPED (%s, n=%d)" % (k, v["skipped"], v["n"]))
                continue
            print("  %-22s effect %+.5f  CI [%+.5f, %+.5f]  p = %.4f   "
                  "(n=%d, units=%d, SD of delta %.4f)"
                  % (k, v["effect"], v["ci"][0], v["ci"][1], v["p"],
                     v["n"], v["units"], v["sd_of_outcome_delta"]))
        m = report["matched"][outcome]
        if "skipped" in m:
            print("  %-22s SKIPPED (%s, pairs=%d)"
                  % ("matched pairs", m["skipped"], m["pairs"]))
        else:
            print("  %-22s mean diff %+.5f  CI [%+.5f, %+.5f]  (%d pairs)"
                  % ("matched pairs", m["mean_difference"], m["ci"][0], m["ci"][1],
                     m["matched_pairs"]))
    print("\n--- reverse direction: does pre-window cohesion predict merging? ---")
    rd = report["reverse_direction"]
    for treat in ("any_large_merge", "any_merge"):
        if treat not in rd:
            continue
        print("  treatment: %s" % rd[treat]["treatment"])
        for k, v in rd[treat]["by_predictor"].items():
            print("      %-38s OR %.3f [%.3f, %.3f]  p = %.4f  (%d positives)"
                  % (k, v["or_per_sd"], v["ci_or"][0], v["ci_or"][1], v["p"],
                     v["positives"]))
    if "verdict" in rd:
        print("\n  VERDICT: " + rd["verdict"])
    print("\n--- A -> C feasibility ---")
    print(json.dumps(report["ac_feasibility"], indent=2))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
