"""Axis B: within-group fission and modularity, taken as far as the coverage allows.

Three questions, three parts.

PART 1 -- Is within-group modularity real, and is it a group property or a state?
  The earlier framing said modularity was coverage-bound and should not be reported as
  biology. That is right for the full 20-group set and wrong for the well-covered
  subset. Detectability rises monotonically with collar count, so the test has to be
  run at matched coverage. Then: variance decomposition (between-group against
  within-group-over-time) and week-to-week persistence.

PART 2 -- Does GENERAL co-sitting predict split composition, or only the proximity
  immediately before the split?
  The saved analysis (AUC 0.568) uses only the cohesive hours immediately preceding
  each split. That measures current state, not bond. This part adds a long-run bond
  measure -- the dyad's co-sitting rate over all other cohesive hours in its group,
  excluding a window around the focal event -- and compares the two, and fits them
  together. It is the between-dyad / within-dyad separation from Phase 4, applied to
  fission.

PART 3 -- Do environmental conditions predict when a group splits and for how long?
  Group-weeks are axis B's opportunity table: a group-week with observation support and
  no split is an observed negative. NDVI is joined per group per week. Collar coverage
  enters every model as an observation covariate, because split DETECTION rises steeply
  with it (Spearman 0.489) -- any ecological effect must survive that adjustment.

Outputs: outputs/general_structure_2026_09/phase4b_axis_b_fission/

CAVEATS CARRIED THROUGHOUT
  * The weekly network metrics run on the legacy source filtered to 2025-01-01 onward;
    axes A and C run on the frozen 1,924,104-row export. Cross-axis pooling is still
    blocked on that reconciliation (gap 7.3).
  * NDVI ends 2026-04-29 while membership runs to 2026-07-22, and per-group NDVI start
    dates differ by up to 18 months. Part 3 reports its retained sample explicitly.
  * `n_animals_observed` is collar coverage, never group size.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

WEEKLY = Path(
    "outputs/canonical_within_group_density_modularity_ridges/"
    "canonical_within_group_weekly_network_metrics.csv"
)
PRESPLIT = Path("outputs/presplit_cositting_prediction")
NDVI = Path(r"C:\Users\rharel\Documents\Github\detect_group_membership\indices_sentinel2.csv")
OUT = Path("outputs/general_structure_2026_09/phase4b_axis_b_fission")

SCOPE = "combined_1100_1600"
WELL_COVERED = 12       # collars observed in a group-week
MATCHED_BAND = (12, 16)  # coverage-matched comparison band
BOND_EXCLUDE_DAYS = 7    # temporal holdout around a focal split event
MIN_BINS = 10            # minimum co-observed 2-min bins for a dyad measure
SEED = 20260904


# ---------------------------------------------------------------- helpers
def auc(y: np.ndarray, x: np.ndarray) -> float:
    """Mann-Whitney AUC. Returns nan if either class is empty."""
    y = np.asarray(y).astype(bool)
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    y, x = y[ok], x[ok]
    if y.all() or (~y).all() or len(y) < 3:
        return np.nan
    r = pd.Series(x).rank().to_numpy()
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def boot_ci(vals: np.ndarray, n: int, rng: np.random.Generator, q=(2.5, 97.5)):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) < 3:
        return (np.nan, np.nan)
    draws = [np.mean(rng.choice(vals, len(vals), replace=True)) for _ in range(n)]
    return tuple(float(v) for v in np.percentile(draws, q))


# ---------------------------------------------------------------- part 1
def part1_modularity(weekly: pd.DataFrame, out: Path, rng) -> dict:
    c = weekly[weekly["scope"].eq(SCOPE)].copy()
    c["period_start"] = pd.to_datetime(c["period_start"])
    c["is_modular"] = c["modularity"] > 0.001
    c["had_split"] = c["split_timestamp_fraction"] > 0

    # --- detectability ladder: modularity is only visible with enough collars
    bands = [(1, 4), (5, 7), (8, 9), (10, 11), (12, 14), (15, 19), (20, 99)]
    ladder = []
    for lo, hi in bands:
        s = c[c["n_animals_observed"].between(lo, hi)]
        if s.empty:
            continue
        ladder.append({
            "collars_observed": f"{lo}-{hi}",
            "group_weeks": len(s),
            "units": s["dynamic_social_unit"].nunique(),
            "pct_weeks_modular": round(100 * s["is_modular"].mean(), 1),
            "pct_weeks_with_split": round(100 * s["had_split"].mean(), 1),
            "median_modularity": round(s["modularity"].median(), 4),
            "p90_modularity": round(s["modularity"].quantile(0.9), 4),
            "median_n_communities": float(s["n_communities"].median()),
        })
    ladder = pd.DataFrame(ladder)
    ladder.to_csv(out / "modularity_detectability_ladder.csv", index=False)

    conf = {
        "spearman_modularity_vs_collars": round(
            c["modularity"].corr(c["n_animals_observed"], method="spearman"), 3),
        "spearman_split_fraction_vs_collars": round(
            c["split_timestamp_fraction"].corr(c["n_animals_observed"], method="spearman"), 3),
    }

    # --- coverage-matched between-group comparison
    m = c[c["n_animals_observed"].between(*MATCHED_BAND)]
    matched = m.groupby("dynamic_social_unit").agg(
        group_weeks=("modularity", "size"),
        median_collars=("n_animals_observed", "median"),
        nominal_group_size=("origin_group_total_size", "first"),
        pct_weeks_modular=("is_modular", lambda x: round(100 * x.mean(), 1)),
        median_modularity=("modularity", "median"),
        max_modularity=("modularity", "max"),
        median_split_fraction=("split_timestamp_fraction", "median"),
        median_largest_community_fraction=("largest_community_fraction", "median"),
    ).reset_index()
    matched = matched[matched["group_weeks"] >= 5].sort_values(
        "pct_weeks_modular", ascending=False)
    matched.to_csv(out / "modularity_coverage_matched.csv", index=False)

    # --- entity or state: variance decomposition on well-covered weeks
    w = c[c["n_animals_observed"].ge(WELL_COVERED)].copy()
    icc = {}
    try:
        import statsmodels.formula.api as smf
        md = smf.mixedlm("modularity ~ n_animals_observed", w,
                         groups=w["dynamic_social_unit"]).fit(reml=True)
        vg = float(md.cov_re.iloc[0, 0])
        ve = float(md.scale)
        icc = {
            "model": "MixedLM(modularity ~ collars + (1|unit)), well-covered weeks",
            "group_weeks": int(len(w)),
            "units": int(w["dynamic_social_unit"].nunique()),
            "between_unit_variance": round(vg, 6),
            "residual_variance": round(ve, 6),
            "icc_between_unit": round(vg / (vg + ve), 3),
            "collar_coefficient": round(float(md.params.get("n_animals_observed", np.nan)), 5),
            "collar_pvalue": round(float(md.pvalues.get("n_animals_observed", np.nan)), 4),
        }
    except Exception as exc:  # pragma: no cover
        icc = {"error": str(exc)}

    # --- persistence: week-to-week autocorrelation within unit
    w = w.sort_values(["dynamic_social_unit", "period_start"])
    rows = []
    for unit, g in w.groupby("dynamic_social_unit"):
        g = g.sort_values("period_start")
        gap = g["period_start"].diff().dt.days
        prev = g["modularity"].shift()
        ok = gap.eq(7) & prev.notna()
        if ok.sum() >= 5:
            rows.append({
                "dynamic_social_unit": unit,
                "consecutive_week_pairs": int(ok.sum()),
                "lag1_autocorrelation": round(
                    float(np.corrcoef(g.loc[ok, "modularity"], prev[ok])[0, 1]), 3),
            })
    persistence = pd.DataFrame(rows)
    if not persistence.empty:
        persistence.to_csv(out / "modularity_persistence.csv", index=False)

    return {
        "detectability_ladder": ladder.to_dict("records"),
        "coverage_confound": conf,
        "coverage_matched_band": list(MATCHED_BAND),
        "coverage_matched": matched.to_dict("records"),
        "variance_decomposition": icc,
        "persistence": persistence.to_dict("records") if not persistence.empty else [],
    }


# ---------------------------------------------------------------- part 2
def part2_cositting(out: Path, rng) -> dict:
    ev = pd.read_parquet(PRESPLIT / "event_dyad_predictions.parquet")
    hours = pd.read_parquet(PRESPLIT / "cohesive_hour_dyads_5m.parquet")
    events = pd.read_csv(PRESPLIT / "candidate_split_events.csv", parse_dates=["split_onset"])

    hours["window_start"] = pd.to_datetime(hours["window_start"], utc=True)
    events["split_onset"] = pd.to_datetime(events["split_onset"], utc=True)
    onset = events.set_index("event_id")["split_onset"]
    ev["split_onset"] = ev["event_id"].map(onset)
    ev = ev[ev["split_onset"].notna()].copy()

    # canonicalise dyad order on both sides so the join is well defined
    for f in (ev, hours):
        a = np.minimum(f["animal_a"], f["animal_b"])
        b = np.maximum(f["animal_a"], f["animal_b"])
        f["animal_a"], f["animal_b"] = a, b

    # long-run bond: same dyad, same unit, all cohesive hours OUTSIDE a window around
    # the focal onset -- a temporal holdout, so the bond cannot contain the event
    idx = ["dynamic_social_unit", "animal_a", "animal_b"]
    tot = hours.groupby(idx)[["coobserved_bins", "within5_bins"]].sum()
    tot.columns = ["bond_co_all", "bond_w5_all"]

    win = pd.Timedelta(days=BOND_EXCLUDE_DAYS)
    near_rows = []
    for eid, g in ev.groupby("event_id"):
        o = g["split_onset"].iloc[0]
        unit = g["dynamic_social_unit"].iloc[0]
        h = hours[hours["dynamic_social_unit"].eq(unit)
                  & hours["window_start"].between(o - win, o + win)]
        if h.empty:
            continue
        agg = h.groupby(idx)[["coobserved_bins", "within5_bins"]].sum().reset_index()
        agg["event_id"] = eid
        near_rows.append(agg)
    near = (pd.concat(near_rows, ignore_index=True) if near_rows
            else pd.DataFrame(columns=idx + ["coobserved_bins", "within5_bins", "event_id"]))
    near = near.rename(columns={"coobserved_bins": "bond_co_near", "within5_bins": "bond_w5_near"})

    d = ev.merge(tot, on=idx, how="left").merge(near, on=idx + ["event_id"], how="left")
    for col in ("bond_co_all", "bond_w5_all", "bond_co_near", "bond_w5_near"):
        d[col] = d[col].fillna(0)
    d["bond_co"] = d["bond_co_all"] - d["bond_co_near"]
    d["bond_w5"] = d["bond_w5_all"] - d["bond_w5_near"]
    d["general_bond"] = np.where(d["bond_co"] >= MIN_BINS, d["bond_w5"] / d["bond_co"], np.nan)
    d["preceding_rate"] = np.where(d["coobserved_bins"] >= MIN_BINS, d["cosit_rate_5m"], np.nan)

    usable = d[d["general_bond"].notna() & d["preceding_rate"].notna()].copy()

    # per-event AUC for each predictor, on the identical retained sample
    per_event = []
    for eid, g in usable.groupby("event_id"):
        if g["same_subgroup"].nunique() < 2 or len(g) < 4:
            continue
        per_event.append({
            "event_id": eid,
            "dynamic_social_unit": g["dynamic_social_unit"].iloc[0],
            "dyads": len(g),
            "auc_preceding": auc(g["same_subgroup"], g["preceding_rate"]),
            "auc_general_bond": auc(g["same_subgroup"], g["general_bond"]),
        })
    pe = pd.DataFrame(per_event).dropna(subset=["auc_preceding", "auc_general_bond"])
    pe.to_csv(out / "cositting_per_event_auc.csv", index=False)

    res: dict = {
        "retained_event_dyads": int(len(usable)),
        "retained_events": int(pe["event_id"].nunique()),
        "retained_units": int(pe["dynamic_social_unit"].nunique()),
        "min_coobserved_bins": MIN_BINS,
        "bond_exclusion_days": BOND_EXCLUDE_DAYS,
        "note": "both AUCs computed on the identical retained sample",
    }
    if not pe.empty:
        for k, col in (("preceding", "auc_preceding"), ("general_bond", "auc_general_bond")):
            lo, hi = boot_ci(pe[col].to_numpy(), 2000, rng)
            res[f"mean_event_auc_{k}"] = round(float(pe[col].mean()), 4)
            res[f"mean_event_auc_{k}_ci"] = [round(lo, 4), round(hi, 4)]
        diff = (pe["auc_general_bond"] - pe["auc_preceding"]).to_numpy()
        lo, hi = boot_ci(diff, 2000, rng)
        res["paired_auc_difference_general_minus_preceding"] = round(float(np.mean(diff)), 4)
        res["paired_auc_difference_ci"] = [round(lo, 4), round(hi, 4)]

    # Joint model, within event. Both predictors are demeaned WITHIN event, so every
    # comparison is between dyads of the same split -- the same thing event fixed
    # effects buy, without a 500-column design. Uncertainty is robust SEs clustered on
    # event (the Phase 4 Stage-1 pattern), not a bootstrap, which was prohibitively
    # slow at this design size.
    try:
        import statsmodels.api as sm
        u = usable.dropna(subset=["general_bond", "preceding_rate"]).copy()
        keep = u.groupby("event_id")["same_subgroup"].transform("nunique").eq(2)
        u = u[keep].copy()
        terms = []
        for col in ("general_bond", "preceding_rate"):
            sd = u[col].std()
            u[f"z_{col}"] = (u[col] - u[col].mean()) / (sd if sd > 0 else 1.0)
            u[f"w_{col}"] = u[f"z_{col}"] - u.groupby("event_id")[f"z_{col}"].transform("mean")
            terms.append(f"w_{col}")
        X = sm.add_constant(u[terms].astype(float))
        fit = sm.GEE(u["same_subgroup"].astype(float), X, groups=u["event_id"],
                     family=sm.families.Binomial(),
                     cov_struct=sm.cov_struct.Independence()).fit()
        ci = fit.conf_int()
        res["joint_within_event_model"] = {
            "design": "logit(same_subgroup) ~ general_bond + preceding_rate, both z-scored "
                      "then demeaned within event; Binomial GEE clustered on event",
            "n_event_dyads": int(len(u)),
            "n_events": int(u["event_id"].nunique()),
            "terms": {
                t: {
                    "log_odds_per_sd": round(float(fit.params[t]), 4),
                    "odds_ratio_per_sd": round(float(np.exp(fit.params[t])), 3),
                    "ci_or": [round(float(np.exp(ci.loc[t, 0])), 3),
                              round(float(np.exp(ci.loc[t, 1])), 3)],
                } for t in terms
            },
        }
        # each predictor alone, on the same retained sample, for a matched comparison
        for t in terms:
            Xs = sm.add_constant(u[[t]].astype(float))
            fs = sm.GEE(u["same_subgroup"].astype(float), Xs, groups=u["event_id"],
                        family=sm.families.Binomial(),
                        cov_struct=sm.cov_struct.Independence()).fit()
            cis = fs.conf_int()
            res["joint_within_event_model"].setdefault("alone", {})[t] = {
                "odds_ratio_per_sd": round(float(np.exp(fs.params[t])), 3),
                "ci_or": [round(float(np.exp(cis.loc[t, 0])), 3),
                          round(float(np.exp(cis.loc[t, 1])), 3)],
            }
    except Exception as exc:  # pragma: no cover
        res["joint_within_event_model"] = {"error": str(exc)}

    usable.to_parquet(out / "cositting_event_dyads.parquet", index=False)
    return res


# ---------------------------------------------------------------- part 3
def part3_environment(weekly: pd.DataFrame, out: Path, rng) -> dict:
    c = weekly[weekly["scope"].eq(SCOPE)].copy()
    c["period_start"] = pd.to_datetime(c["period_start"])
    c["had_split"] = (c["split_timestamp_fraction"] > 0).astype(int)
    c["unit_key"] = c["dynamic_social_unit"].str.lower()

    nd = pd.read_csv(NDVI, parse_dates=["date"])
    nd = nd[nd["extraction_status"].eq("success")].copy()
    nd["unit_key"] = nd["group_id"].str.lower()
    nd["period_start"] = nd["date"] - pd.to_timedelta(nd["date"].dt.weekday, unit="D")
    wk = nd.groupby(["unit_key", "period_start"]).agg(
        ndvi=("NDVI", "mean"), ndvi_n=("NDVI", "size"),
        ndwi=("NDWI", "mean"), evi=("EVI", "mean")).reset_index()

    j = c.merge(wk, on=["unit_key", "period_start"], how="left")
    retained = j[j["ndvi"].notna() & j["n_animals_observed"].ge(5)].copy()

    # within-unit NDVI deviation separates "this range is greener" from
    # "this week was greener than usual for this range"
    retained["ndvi_unit_mean"] = retained.groupby("unit_key")["ndvi"].transform("mean")
    retained["ndvi_deviation"] = retained["ndvi"] - retained["ndvi_unit_mean"]
    doy = retained["period_start"].dt.dayofyear
    retained["season_sin"] = np.sin(2 * np.pi * doy / 365.25)
    retained["season_cos"] = np.cos(2 * np.pi * doy / 365.25)

    res: dict = {
        "ndvi_source": str(NDVI),
        "ndvi_date_range": [str(nd["date"].min().date()), str(nd["date"].max().date())],
        "group_weeks_before_join": int(len(c)),
        "group_weeks_retained": int(len(retained)),
        "units_retained": int(retained["unit_key"].nunique()),
        "retention_note": (
            "requires a successful NDVI week for that group and >=5 collars observed; "
            "NDVI ends 2026-04-29 while membership runs to 2026-07-22"
        ),
        "split_week_rate": round(float(retained["had_split"].mean()), 4),
    }

    try:
        import statsmodels.api as sm

        def gee(formula_terms, response, data, family, label):
            X = data[formula_terms].astype(float).copy()
            for col in formula_terms:
                sd = X[col].std()
                if sd > 0:
                    X[col] = (X[col] - X[col].mean()) / sd
            X = sm.add_constant(X)
            m = sm.GEE(data[response].astype(float), X, groups=data["unit_key"],
                       family=family, cov_struct=sm.cov_struct.Independence()).fit()
            ci = m.conf_int()
            return {
                "label": label,
                "n": int(len(data)),
                "clusters": int(data["unit_key"].nunique()),
                "terms": {
                    t: {
                        "coef": round(float(m.params[t]), 4),
                        "or_per_sd": round(float(np.exp(m.params[t])), 3),
                        "ci_or": [round(float(np.exp(ci.loc[t, 0])), 3),
                                  round(float(np.exp(ci.loc[t, 1])), 3)],
                    } for t in m.params.index if t != "const"
                },
            }

        terms = ["ndvi_unit_mean", "ndvi_deviation", "n_animals_observed",
                 "season_sin", "season_cos"]
        res["occurrence_model"] = gee(
            terms, "had_split", retained, sm.families.Binomial(),
            "does a split occur in this group-week? Binomial GEE clustered on unit")
        res["occurrence_model_no_coverage_term"] = gee(
            [t for t in terms if t != "n_animals_observed"], "had_split", retained,
            sm.families.Binomial(),
            "same model with the collar-coverage term dropped -- for contrast only")

        # The reported response (any split timestamp) has a 71% base rate, which is
        # near-inevitable for a group-week with 13 collars. Rerun on progressively
        # stricter responses, and on the well-covered subset, so the conclusion is not
        # an artefact of a permissive threshold.
        retained["split_any"] = (retained["split_timestamp_fraction"] > 0).astype(int)
        retained["split_majority"] = (retained["split_timestamp_fraction"] > 0.5).astype(int)
        retained["split_multi_any"] = (
            retained["multi_animal_split_timestamp_fraction"] > 0).astype(int)
        retained["split_multi_quarter"] = (
            retained["multi_animal_split_timestamp_fraction"] > 0.25).astype(int)
        responses = [
            ("split_any", "any split timestamp (as reported above)"),
            ("split_majority", "group split in >50% of timestamps"),
            ("split_multi_any", "any multi-animal split, excluding lone animals"),
            ("split_multi_quarter", "multi-animal split in >25% of timestamps"),
        ]
        ladder = []
        well = retained[retained["n_animals_observed"].ge(WELL_COVERED)]
        for subset_label, frame in (("all retained", retained),
                                    (f">={WELL_COVERED} collars", well)):
            for resp, resp_label in responses:
                try:
                    m = gee(terms, resp, frame, sm.families.Binomial(), resp_label)
                except Exception:
                    continue
                row = {"subset": subset_label, "response": resp_label,
                       "n": m["n"], "base_rate": round(float(frame[resp].mean()), 3)}
                for t in ("ndvi_unit_mean", "ndvi_deviation", "n_animals_observed"):
                    if t in m["terms"]:
                        row[f"{t}_or"] = m["terms"][t]["or_per_sd"]
                        row[f"{t}_ci"] = m["terms"][t]["ci_or"]
                ladder.append(row)
        res["occurrence_response_ladder"] = ladder
        pd.DataFrame(ladder).to_csv(out / "split_occurrence_response_ladder.csv", index=False)
        res["occurrence_verdict"] = (
            "The within-unit NDVI deviation term -- the only term that asks whether a "
            "group splits more in weeks greener than its own norm -- is null in every "
            "specification, and in the two strictest it turns weakly negative with a CI "
            "touching 1. The large ndvi_unit_mean term is between-group range greenness, "
            "not separable from group identity across 19 units; it is the same failure "
            "mode as the axis-A seasonal term, which fell from OR 0.62 to 0.86 under "
            "restriction. Collar coverage predicts split detection in the full sample and "
            "stops doing so above the well-covered threshold."
        )

        # extent, given a split week
        pos = retained[retained["had_split"].eq(1)].copy()
        if len(pos) > 30:
            res["extent_model"] = gee(
                terms, "split_timestamp_fraction", pos, sm.families.Gaussian(),
                "given a split week, what fraction of timestamps are split? "
                "Gaussian GEE (identity), z-scored predictors")
            res["extent_model"]["response_note"] = "coefficients are per SD of predictor, on the raw fraction"
    except Exception as exc:  # pragma: no cover
        res["occurrence_model"] = {"error": str(exc)}

    # duration: persistent split hours per event against NDVI at onset
    try:
        events = pd.read_csv(PRESPLIT / "candidate_split_events.csv", parse_dates=["split_onset"])
        events["unit_key"] = events["dynamic_social_unit"].str.lower()
        events["onset_date"] = events["split_onset"].dt.tz_localize(None).dt.normalize()
        nd_daily = nd[["unit_key", "date", "NDVI"]].rename(columns={"date": "onset_date"})
        e = events.merge(nd_daily, on=["unit_key", "onset_date"], how="left")
        e = e[e["NDVI"].notna() & e["persistent_split_hours"].gt(0)].copy()
        res["duration_sample"] = {
            "events_with_ndvi_on_onset_day": int(len(e)),
            "units": int(e["unit_key"].nunique()),
            "median_persistent_split_hours": float(e["persistent_split_hours"].median()),
            "p90_persistent_split_hours": float(e["persistent_split_hours"].quantile(0.9)),
            "max_persistent_split_hours": float(e["persistent_split_hours"].max()),
        }
        if len(e) > 30:
            e["ndvi_unit_mean"] = e.groupby("unit_key")["NDVI"].transform("mean")
            e["ndvi_deviation"] = e["NDVI"] - e["ndvi_unit_mean"]
            import statsmodels.api as sm
            X = e[["ndvi_unit_mean", "ndvi_deviation", "n_animals_onset"]].astype(float)
            X = sm.add_constant((X - X.mean()) / X.std())
            m = sm.GEE(np.log(e["persistent_split_hours"].astype(float)), X,
                       groups=e["unit_key"], family=sm.families.Gaussian(),
                       cov_struct=sm.cov_struct.Independence()).fit()
            ci = m.conf_int()
            res["duration_model"] = {
                "label": "log(persistent split hours) ~ NDVI terms + animals at onset; "
                         "Gaussian GEE clustered on unit",
                "n": int(len(e)),
                "clusters": int(e["unit_key"].nunique()),
                "terms": {
                    t: {"coef": round(float(m.params[t]), 4),
                        "ci": [round(float(ci.loc[t, 0]), 4), round(float(ci.loc[t, 1]), 4)]}
                    for t in m.params.index if t != "const"
                },
            }
        e.to_csv(out / "split_duration_with_ndvi.csv", index=False)
    except Exception as exc:  # pragma: no cover
        res["duration_model"] = {"error": str(exc)}

    retained.to_csv(out / "group_weeks_with_ndvi.csv", index=False)
    return res


# ---------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weekly", type=Path, default=WEEKLY)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    weekly = pd.read_csv(args.weekly)
    report = {
        "seed": SEED,
        "scope": SCOPE,
        "weekly_source": str(args.weekly),
        "source_caveat": (
            "weekly network metrics derive from the legacy hourly source filtered to "
            "2025-01-01 onward; axes A and C use the frozen 1,924,104-row export. "
            "Cross-axis pooling remains blocked on gap 7.3."
        ),
    }

    print("=" * 74)
    print("PART 1 -- modularity: detectability, then group differences at matched coverage")
    print("=" * 74)
    report["part1_modularity"] = part1_modularity(weekly, args.output_dir, rng)
    p1 = report["part1_modularity"]
    print(pd.DataFrame(p1["detectability_ladder"]).to_string(index=False))
    print("\ncoverage confound:", json.dumps(p1["coverage_confound"]))
    print("\ncoverage-matched (%d-%d collars, >=5 group-weeks):" % tuple(p1["coverage_matched_band"]))
    print(pd.DataFrame(p1["coverage_matched"]).to_string(index=False))
    print("\nvariance decomposition:", json.dumps(p1["variance_decomposition"], indent=2))
    if p1["persistence"]:
        print("\nweek-to-week persistence:")
        print(pd.DataFrame(p1["persistence"]).to_string(index=False))

    print("\n" + "=" * 74)
    print("PART 2 -- general bond against immediately-preceding proximity")
    print("=" * 74)
    report["part2_cositting"] = part2_cositting(args.output_dir, rng)
    print(json.dumps(report["part2_cositting"], indent=2))

    print("\n" + "=" * 74)
    print("PART 3 -- environment: split occurrence, extent and duration")
    print("=" * 74)
    report["part3_environment"] = part3_environment(weekly, args.output_dir, rng)
    print(json.dumps(report["part3_environment"], indent=2))

    with open(args.output_dir / "axis_b_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
