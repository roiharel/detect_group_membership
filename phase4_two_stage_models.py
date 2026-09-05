"""Phase 4: fit the two stages with inference matched to the design.

Stage 1  encounter occurrence on the Phase 1 opportunity table.
Stage 2  mixing during an encounter on the Phase 2 event table.

Design decisions this script exists to enforce:

* Covariates are strictly preceding. Same-day distance is endogenous during an
  encounter, so the within-dyad distance term is a lag-1 deviation and the
  event-level distance term is measured in the days before the event starts.
* Between-dyad and within-dyad proximity are separate terms (a within-between
  decomposition), because "these two ranges overlap" and "they are unusually
  close today" are different claims.
* A group contributes ONE effect that applies whichever end of the dyad it sits
  on, rather than separate alphabetical group_a and group_b effects.
* Collar counts are named collar coverage and treated as observation
  covariates. They are not group size.
* Every model in a comparison is fitted on one identical retained sample.
* Stage 2 has 12 dyads, which cannot support between-dyad inference, so it is
  fitted as a within-dyad estimator with dyad fixed effects and uncertainty
  from a cluster bootstrap over encounters.

Nothing here is called a hierarchical Bayesian fit unless it is one. The mixed
model is mean-field variational Bayes and is reported as such.

Run from the project root:
    python phase4_two_stage_models.py
"""
from __future__ import annotations

import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
from statsmodels.genmod.cov_struct import Independence
from statsmodels.genmod.families import Binomial

PROJECT = Path(r"C:\Users\rharel\Documents\group_mebership")
UPSTREAM = Path(r"C:\Users\rharel\Documents\New project")
CANON = UPSTREAM / "outputs" / "canonical_robust_hourly_membership_shared_full_20260722"
GENERAL = PROJECT / "outputs" / "general_structure_2026_09"
OPPORTUNITY = GENERAL / "phase1_opportunity" / "opportunity_dyad_day.csv"
EVENTS = GENERAL / "phase2_two_stage_events" / "stage1_events_with_stage2_mixing.csv"
FINE_5M = (CANON / "canonical_5m_shared_history_shuffle_expectation"
           / "canonical_5m_shuffle_expectation_2min_rows.csv")
OUT_DIR = GENERAL / "phase4_two_stage_models"

N_BOOTSTRAP = 2000
RNG_SEED = 20260904
CV_FOLDS = 5
PRE_EVENT_DAYS = 7
PLAUSIBLE_OPPORTUNITY_M = 5000.0
VB_PRIOR_SCALES = [0.5, 1.0, 2.0]


def zscore(values: pd.Series) -> pd.Series:
    sd = values.std(ddof=0)
    return (values - values.mean()) / (sd if sd and np.isfinite(sd) else 1.0)


def logit(p: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    p = np.clip(p, eps, 1.0 - eps)
    return np.log(p / (1.0 - p))


# --------------------------------------------------------------------------- #
# Stage 1
# --------------------------------------------------------------------------- #

def build_stage1_sample() -> tuple[pd.DataFrame, dict[str, object]]:
    days = pd.read_csv(OPPORTUNITY, parse_dates=["period_start"])
    audit: dict[str, object] = {"opportunity_rows": int(len(days))}

    supported = days[days["state"].ne("insufficient_support")].copy()
    audit["dropped_insufficient_support"] = int(len(days) - len(supported))
    supported["encounter"] = supported["state"].eq("detected_encounter").astype(int)
    supported = supported.sort_values(["pair_key", "period_start"]).reset_index(drop=True)

    # strictly preceding history: cumulative encounters BEFORE this row
    grouped = supported.groupby("pair_key", observed=True)
    supported["prior_encounters"] = grouped["encounter"].cumsum() - supported["encounter"]

    # proximity, decomposed between dyads and within a dyad over time
    supported["log_med_dist"] = np.log1p(supported["median_centroid_distance_m"].clip(lower=0))
    supported["dyad_mean_log_dist"] = grouped["log_med_dist"].transform("mean")
    # lag-1 so the term cannot contain the encounter it is meant to predict
    supported["prev_log_med_dist"] = grouped["log_med_dist"].shift(1)
    supported["prev_days_gap"] = (
        supported["period_start"] - grouped["period_start"].shift(1)
    ).dt.days
    supported["within_dyad_prev_dist_dev"] = (
        supported["prev_log_med_dist"] - supported["dyad_mean_log_dist"]
    )

    # observation covariates -- collar coverage, NOT group size
    coverage_min = supported[["max_observed_a", "max_observed_b"]].min(axis=1)
    coverage_total = supported["max_observed_a"] + supported["max_observed_b"]
    supported["log_collar_coverage_min"] = np.log1p(coverage_min)
    supported["log_collar_coverage_total"] = np.log1p(coverage_total)
    supported["log_shared_hours_exact"] = np.log1p(supported["shared_hours_exact"])

    doy = supported["period_start"].dt.dayofyear
    supported["season_sin"] = np.sin(2.0 * np.pi * doy / 365.25)
    supported["season_cos"] = np.cos(2.0 * np.pi * doy / 365.25)

    fitted = supported.dropna(subset=["within_dyad_prev_dist_dev"]).copy()
    audit["dropped_no_previous_day"] = int(len(supported) - len(fitted))
    # a lagged term is only meaningful if the previous observation is recent
    fitted = fitted[fitted["prev_days_gap"].le(3)].copy()
    audit["dropped_stale_previous_day_gap_over_3"] = int(
        len(supported) - audit["dropped_no_previous_day"] - len(fitted)
    )

    fitted["z_between_dyad_dist"] = zscore(fitted["dyad_mean_log_dist"])
    fitted["z_within_dyad_prev_dist"] = zscore(fitted["within_dyad_prev_dist_dev"])
    fitted["z_prior_encounters"] = zscore(np.log1p(fitted["prior_encounters"]))
    fitted["z_collar_coverage_min"] = zscore(fitted["log_collar_coverage_min"])
    fitted["z_collar_coverage_total"] = zscore(fitted["log_collar_coverage_total"])
    fitted["z_shared_hours"] = zscore(fitted["log_shared_hours_exact"])

    audit["retained_rows"] = int(len(fitted))
    audit["retained_dyads"] = int(fitted["pair_key"].nunique())
    audit["retained_encounters"] = int(fitted["encounter"].sum())
    audit["retained_encounter_rate"] = round(float(fitted["encounter"].mean()), 5)
    return fitted, audit


def shared_group_matrix(rows: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """One column per group, set to 1 whichever end of the dyad the group is on."""
    groups = sorted(set(rows["unit_a"]) | set(rows["unit_b"]))
    data = {}
    for group in groups:
        data[f"grp_{group}"] = (
            rows["unit_a"].eq(group).astype(float) + rows["unit_b"].eq(group).astype(float)
        )
    return pd.DataFrame(data, index=rows.index), [f"grp_{g}" for g in groups]


STAGE1_BLOCKS = {
    "M0_proximity": ["z_between_dyad_dist", "z_within_dyad_prev_dist"],
    "M1_plus_history": ["z_prior_encounters"],
    "M2_plus_coverage": ["z_collar_coverage_min", "z_collar_coverage_total", "z_shared_hours"],
    "M3_plus_season": ["season_sin", "season_cos"],
}


def fit_stage1_gee(rows: pd.DataFrame, terms: list[str], group_cols: list[str],
                   label: str) -> tuple[pd.DataFrame, dict[str, object]]:
    exog = sm.add_constant(rows[terms + group_cols], has_constant="add")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = sm.GEE(rows["encounter"], exog, groups=rows["pair_key"],
                       family=Binomial(), cov_struct=Independence())
        result = model.fit(maxiter=100)
    coef = pd.DataFrame(
        {
            "model": label,
            "term": result.params.index,
            "estimate": result.params.to_numpy(),
            "std_error": result.bse.to_numpy(),
        }
    )
    coef["z"] = coef["estimate"] / coef["std_error"]
    coef["ci_low"] = coef["estimate"] - 1.96 * coef["std_error"]
    coef["ci_high"] = coef["estimate"] + 1.96 * coef["std_error"]
    coef["odds_ratio"] = np.exp(coef["estimate"])
    coef["is_group_effect"] = coef["term"].str.startswith("grp_")
    fitted_p = np.asarray(result.fittedvalues, dtype=float)
    info = {
        "model": label,
        "rows": int(len(rows)),
        "dyad_clusters": int(rows["pair_key"].nunique()),
        "encounters": int(rows["encounter"].sum()),
        "n_terms": int(exog.shape[1]),
        "in_sample_auc": round(float(roc_auc_score(rows["encounter"], fitted_p)), 4),
    }
    return coef, info


def stage1_cross_validate(rows: pd.DataFrame, terms: list[str],
                          group_cols: list[str]) -> tuple[dict[str, object], pd.DataFrame]:
    """Out-of-sample calibration with whole dyads held out."""
    splitter = GroupKFold(n_splits=CV_FOLDS)
    predictions = np.full(len(rows), np.nan)
    for train_idx, test_idx in splitter.split(rows, groups=rows["pair_key"]):
        train, test = rows.iloc[train_idx], rows.iloc[test_idx]
        exog_train = sm.add_constant(train[terms + group_cols], has_constant="add")
        exog_test = sm.add_constant(test[terms + group_cols], has_constant="add")
        exog_test = exog_test.reindex(columns=exog_train.columns, fill_value=0.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = sm.GEE(train["encounter"], exog_train, groups=train["pair_key"],
                         family=Binomial(), cov_struct=Independence()).fit(maxiter=100)
        predictions[test_idx] = fit.predict(exog_test)
    valid = np.isfinite(predictions)
    observed = rows["encounter"].to_numpy()[valid]
    predicted = predictions[valid]
    deciles = pd.qcut(predicted, 10, labels=False, duplicates="drop")
    calibration = (
        pd.DataFrame({"decile": deciles, "predicted": predicted, "observed": observed})
        .groupby("decile", observed=True)
        .agg(rows=("observed", "size"), mean_predicted=("predicted", "mean"),
             observed_rate=("observed", "mean"))
        .reset_index()
    )
    slope = sm.GLM(observed, sm.add_constant(logit(predicted)),
                   family=Binomial()).fit()
    summary = {
        "folds": CV_FOLDS,
        "held_out_by": "whole dyads",
        "out_of_sample_auc": round(float(roc_auc_score(observed, predicted)), 4),
        "calibration_intercept": round(float(slope.params[0]), 4),
        "calibration_slope": round(float(slope.params[1]), 4),
        "mean_predicted": round(float(predicted.mean()), 5),
        "observed_rate": round(float(observed.mean()), 5),
    }
    return summary, calibration


def fit_stage1_mixed(rows: pd.DataFrame, terms: list[str],
                     prior_scale: float) -> dict[str, object]:
    """Mean-field variational Bayes mixed GLM: dyad effect + shared group effect."""
    exog = sm.add_constant(rows[terms], has_constant="add").to_numpy(dtype=float)
    groups = sorted(set(rows["unit_a"]) | set(rows["unit_b"]))
    group_index = {g: i for i, g in enumerate(groups)}
    dyads = sorted(rows["pair_key"].unique())
    dyad_index = {d: i for i, d in enumerate(dyads)}

    n = len(rows)
    rows_idx = np.arange(n)
    a_col = rows["unit_a"].map(group_index).to_numpy()
    b_col = rows["unit_b"].map(group_index).to_numpy()
    group_design = sp.coo_matrix(
        (np.ones(2 * n), (np.concatenate([rows_idx, rows_idx]), np.concatenate([a_col, b_col]))),
        shape=(n, len(groups)),
    )
    dyad_design = sp.coo_matrix(
        (np.ones(n), (rows_idx, rows["pair_key"].map(dyad_index).to_numpy())),
        shape=(n, len(dyads)),
    )
    exog_vc = sp.hstack([group_design, dyad_design]).tocsr()
    ident = np.concatenate([np.zeros(len(groups), dtype=int), np.ones(len(dyads), dtype=int)])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = BinomialBayesMixedGLM(
            rows["encounter"].to_numpy(dtype=float), exog, exog_vc, ident,
            vcp_p=prior_scale, fe_p=2.0,
        )
        result = model.fit_vb(verbose=False)

    names = ["const"] + terms
    fixed = {
        name: {"posterior_mean": round(float(result.fe_mean[i]), 4),
               "posterior_sd": round(float(result.fe_sd[i]), 4)}
        for i, name in enumerate(names)
    }
    return {
        "prior_scale_vcp": prior_scale,
        "variance_components": {
            "shared_group_effect_log_sd": round(float(result.vcp_mean[0]), 4),
            "dyad_effect_log_sd": round(float(result.vcp_mean[1]), 4),
            "shared_group_effect_sd": round(float(np.exp(result.vcp_mean[0])), 4),
            "dyad_effect_sd": round(float(np.exp(result.vcp_mean[1])), 4),
        },
        "fixed_effects": fixed,
        "n_groups": len(groups),
        "n_dyads": len(dyads),
    }


# --------------------------------------------------------------------------- #
# Stage 2
# --------------------------------------------------------------------------- #

def build_stage2_sample(stage1_rows: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    events = pd.read_csv(EVENTS, parse_dates=["start_hour", "end_hour"])
    audit: dict[str, object] = {"stage1_events": int(len(events))}

    fine = pd.read_csv(
        FINE_5M,
        usecols=["bin_2min", "pair_key", "cross_edges", "total_edges",
                 "shuffle_expected_cross_edge_fraction"],
        parse_dates=["bin_2min"],
    )
    fine["pair_key"] = fine["pair_key"].astype(str)
    windows = events[["stage1_event_id", "pair_key", "start_hour", "end_hour"]].copy()
    windows["window_end_exclusive"] = windows["end_hour"] + pd.Timedelta(hours=1)

    chunks = []
    for pair, bins in fine.groupby("pair_key", observed=True):
        pair_windows = windows[windows["pair_key"].eq(pair)]
        if pair_windows.empty:
            continue
        matched = pd.merge_asof(
            bins.sort_values("bin_2min"),
            pair_windows.sort_values("start_hour")[
                ["start_hour", "window_end_exclusive", "stage1_event_id"]
            ],
            left_on="bin_2min", right_on="start_hour", direction="backward",
        )
        chunks.append(matched[matched["bin_2min"] < matched["window_end_exclusive"]])
    bins_in_events = pd.concat(chunks, ignore_index=True)

    # edge-weighted expectation, so the offset matches the estimand
    per_event = (
        bins_in_events.groupby("stage1_event_id", observed=True)
        .apply(
            lambda f: pd.Series(
                {
                    "cross_edges": float(f["cross_edges"].sum()),
                    "total_edges": float(f["total_edges"].sum()),
                    "expected_fraction": float(
                        np.average(f["shuffle_expected_cross_edge_fraction"],
                                   weights=f["total_edges"])
                    ) if f["total_edges"].sum() > 0 else np.nan,
                    "eligible_bins": int(f["bin_2min"].nunique()),
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )
    sample = events.merge(per_event, on="stage1_event_id", how="inner")
    sample = sample[sample["total_edges"].gt(0) & sample["expected_fraction"].notna()].copy()
    audit["events_with_finescale_edges"] = int(len(sample))

    # strictly preceding: encounters and typical distance before the event starts
    history = (
        stage1_rows[["pair_key", "period_start", "encounter", "median_centroid_distance_m"]]
        .sort_values(["pair_key", "period_start"])
    )
    pre_encounters, pre_distance = [], []
    for record in sample.itertuples(index=False):
        past = history[
            history["pair_key"].eq(record.pair_key)
            & (history["period_start"] < record.start_hour.floor("D"))
        ]
        pre_encounters.append(int(past["encounter"].sum()))
        window = past[past["period_start"]
                      >= record.start_hour.floor("D") - pd.Timedelta(days=PRE_EVENT_DAYS)]
        pre_distance.append(
            float(np.log1p(window["median_centroid_distance_m"].median()))
            if len(window) else np.nan
        )
    sample["prior_encounters"] = pre_encounters
    sample["pre_event_log_distance"] = pre_distance

    sample["observed_fraction"] = sample["cross_edges"] / sample["total_edges"]
    sample["logit_expected_offset"] = logit(sample["expected_fraction"].to_numpy())
    sample["log_structural_span"] = np.log(sample["structural_span_hours"].clip(lower=1.0))
    sample["log_supported_exposure"] = np.log(sample["eligible_bins"].clip(lower=1))
    sample["log_collar_coverage"] = np.log1p(
        sample[["mean_observed_a", "mean_observed_b"]].min(axis=1)
    )
    sample["is_sustained"] = sample["encounter_class"].eq("sustained_association").astype(float)

    sample = sample.dropna(subset=["pre_event_log_distance"]).copy()
    audit["dropped_no_preceding_distance"] = int(len(per_event) - len(sample))
    for col, name in [
        ("log_structural_span", "z_log_span"),
        ("log_supported_exposure", "z_log_exposure"),
        ("log_collar_coverage", "z_collar_coverage"),
        ("pre_event_log_distance", "z_pre_event_distance"),
    ]:
        sample[name] = zscore(sample[col])
    sample["z_prior_encounters"] = zscore(np.log1p(sample["prior_encounters"]))

    # within a dyad, prior_encounters is close to a time index, so carry an
    # explicit calendar-time term to test whether the history effect survives it
    sample["days_since_study_start"] = (
        sample["start_hour"] - sample["start_hour"].min()
    ).dt.total_seconds() / 86400.0
    sample["z_event_time"] = zscore(sample["days_since_study_start"])

    audit["collinearity_warnings"] = {
        "is_sustained_vs_z_log_span_pearson": round(float(
            sample["is_sustained"].corr(sample["z_log_span"])), 3),
        "z_prior_encounters_vs_z_event_time_pearson": round(float(
            sample["z_prior_encounters"].corr(sample["z_event_time"])), 3),
        "note": "a sustained association is DEFINED by its span, so the two must "
                "not enter one model together; within a dyad, prior encounters "
                "largely indexes time",
    }

    audit["retained_events"] = int(len(sample))
    audit["retained_dyads"] = int(sample["pair_key"].nunique())
    audit["retained_discrete"] = int((sample["encounter_class"] == "discrete_encounter").sum())
    audit["retained_sustained"] = int((sample["encounter_class"] == "sustained_association").sum())
    return sample, audit


STAGE2_TERMS = ["z_pre_event_distance", "z_prior_encounters", "z_log_span",
                "z_log_exposure", "z_collar_coverage"]


def stage2_design(rows: pd.DataFrame, terms: list[str], dyad_fixed: bool) -> pd.DataFrame:
    exog = rows[terms].copy()
    if dyad_fixed:
        dummies = pd.get_dummies(rows["pair_key"], prefix="dyad", drop_first=True, dtype=float)
        exog = pd.concat([exog, dummies], axis=1)
    return sm.add_constant(exog, has_constant="add")


def fit_stage2(rows: pd.DataFrame, terms: list[str], label: str,
               dyad_fixed: bool, rng: np.random.Generator) -> tuple[pd.DataFrame, dict[str, object]]:
    """Binomial GLM on cross-edge counts with the composition expectation as offset.

    The intercept is therefore the mean logit deficit relative to expectation.
    Uncertainty comes from a cluster bootstrap over encounters within dyads,
    because 12 dyads cannot support asymptotic between-dyad inference.
    """
    exog = stage2_design(rows, terms, dyad_fixed)
    endog = np.column_stack([rows["cross_edges"], rows["total_edges"] - rows["cross_edges"]])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = sm.GLM(endog, exog, family=Binomial(),
                     offset=rows["logit_offset_used"].to_numpy()).fit()

    keep = [c for c in exog.columns if not c.startswith("dyad_")]
    draws = {name: [] for name in keep}
    dyad_groups = {d: g for d, g in rows.groupby("pair_key", observed=True)}
    for _ in range(N_BOOTSTRAP):
        parts = []
        for frame in dyad_groups.values():
            picks = rng.integers(0, len(frame), size=len(frame))
            parts.append(frame.iloc[picks])
        sample = pd.concat(parts, ignore_index=True)
        try:
            b_exog = stage2_design(sample, terms, dyad_fixed).reindex(
                columns=exog.columns, fill_value=0.0
            )
            b_endog = np.column_stack(
                [sample["cross_edges"], sample["total_edges"] - sample["cross_edges"]]
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                b_fit = sm.GLM(b_endog, b_exog, family=Binomial(),
                               offset=sample["logit_offset_used"].to_numpy()).fit()
            for name in keep:
                draws[name].append(float(b_fit.params[name]))
        except Exception:
            continue

    records = []
    for name in keep:
        values = np.asarray(draws[name], dtype=float)
        values = values[np.isfinite(values)]
        records.append(
            {
                "model": label,
                "term": name,
                "estimate": float(fit.params[name]),
                "boot_low": float(np.percentile(values, 2.5)) if values.size else np.nan,
                "boot_high": float(np.percentile(values, 97.5)) if values.size else np.nan,
                "boot_draws": int(values.size),
                "model_se_unadjusted": float(fit.bse[name]),
            }
        )
    coef = pd.DataFrame(records)
    info = {
        "model": label,
        "events": int(len(rows)),
        "dyads": int(rows["pair_key"].nunique()),
        "dyad_fixed_effects": dyad_fixed,
        "total_edges": int(rows["total_edges"].sum()),
        "pearson_chi2_over_df": round(float(fit.pearson_chi2 / fit.df_resid), 2),
        "mean_observed_fraction": round(float(
            rows["cross_edges"].sum() / rows["total_edges"].sum()), 4),
        "mean_expected_fraction": round(float(np.average(
            rows["expected_fraction"], weights=rows["total_edges"])), 4),
    }
    return coef, info


# --------------------------------------------------------------------------- #

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    report: dict[str, object] = {
        "phase": "4 - two-stage models",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }

    # ---- Stage 1 --------------------------------------------------------- #
    print("building stage-1 sample ...")
    stage1, stage1_audit = build_stage1_sample()
    print(f"  {stage1_audit['retained_rows']:,} rows, "
          f"{stage1_audit['retained_dyads']} dyads, "
          f"{stage1_audit['retained_encounters']:,} encounters")
    report["stage1_sample"] = stage1_audit

    group_frame, group_cols = shared_group_matrix(stage1)
    stage1 = pd.concat([stage1, group_frame], axis=1)
    # one group is the reference so the shared-group columns are identifiable
    group_cols_fit = group_cols[1:]
    report["stage1_shared_group_effects"] = {
        "n_groups": len(group_cols),
        "reference_group": group_cols[0].replace("grp_", ""),
        "encoding": "one column per group, 1 whichever end of the dyad it occupies",
    }

    coefs, infos = [], []
    terms: list[str] = []
    for label, block in STAGE1_BLOCKS.items():
        terms = terms + block
        coef, info = fit_stage1_gee(stage1, terms, group_cols_fit, label)
        coefs.append(coef)
        infos.append(info)
        print(f"  {label}: AUC {info['in_sample_auc']}  ({info['n_terms']} terms)")
    stage1_coef = pd.concat(coefs, ignore_index=True)
    stage1_coef.to_csv(OUT_DIR / "stage1_coefficients.csv", index=False)
    report["stage1_models"] = infos
    report["stage1_matched_sample_note"] = (
        "all stage-1 models are fitted on the identical retained sample of "
        f"{stage1_audit['retained_rows']} rows, so nested comparisons are valid"
    )

    print("cross-validating stage 1 by held-out dyads ...")
    cv_summary, calibration = stage1_cross_validate(stage1, terms, group_cols_fit)
    calibration.to_csv(OUT_DIR / "stage1_calibration_deciles.csv", index=False)
    report["stage1_cross_validation"] = cv_summary
    print(f"  out-of-sample AUC {cv_summary['out_of_sample_auc']}, "
          f"calibration slope {cv_summary['calibration_slope']}")

    print("fitting stage-1 mixed model (variational Bayes) ...")
    mixed_runs = []
    for scale in VB_PRIOR_SCALES:
        try:
            mixed_runs.append(fit_stage1_mixed(stage1, terms, scale))
            print(f"  prior scale {scale}: done")
        except Exception as exc:  # noqa: BLE001
            mixed_runs.append({"prior_scale_vcp": scale, "error": str(exc)[:200]})
            print(f"  prior scale {scale}: failed ({str(exc)[:80]})")
    report["stage1_mixed_model"] = {
        "method": "mean-field variational Bayes mixed GLM "
                  "(statsmodels BinomialBayesMixedGLM.fit_vb); NOT MCMC, and no "
                  "convergence diagnostics of the MCMC kind apply",
        "random_effects": "shared group effect (26 levels) + dyad effect",
        "prior_sensitivity": mixed_runs,
    }

    print("stage-1 sensitivity: plausible-opportunity subset ...")
    plausible = stage1[stage1["prev_log_med_dist"].le(np.log1p(PLAUSIBLE_OPPORTUNITY_M))].copy()
    sens_coef, sens_info = fit_stage1_gee(plausible, terms, group_cols_fit, "S1_within_5km")
    sens_coef.to_csv(OUT_DIR / "stage1_sensitivity_within_5km.csv", index=False)
    sens_info["encounter_rate"] = round(float(plausible["encounter"].mean()), 4)
    report["stage1_sensitivity_within_5km"] = sens_info
    print(f"  {sens_info['rows']:,} rows, encounter rate {sens_info['encounter_rate']}")

    # ---- Stage 2 --------------------------------------------------------- #
    print("building stage-2 sample ...")
    stage2, stage2_audit = build_stage2_sample(stage1)
    print(f"  {stage2_audit['retained_events']} events, {stage2_audit['retained_dyads']} dyads")
    report["stage2_sample"] = stage2_audit
    stage2["logit_offset_used"] = stage2["logit_expected_offset"]
    stage2.to_csv(OUT_DIR / "stage2_event_sample.csv", index=False)

    stage2_coefs, stage2_infos = [], []
    discrete = stage2[stage2["encounter_class"].eq("discrete_encounter")].copy()

    # primary: discrete encounters only, within dyad
    coef, info = fit_stage2(discrete, STAGE2_TERMS, "D0_within_dyad", True, rng)
    stage2_coefs.append(coef); stage2_infos.append(info)
    # how much the estimate moves when dyad effects are dropped
    coef, info = fit_stage2(discrete, STAGE2_TERMS, "D1_no_dyad_effects", False, rng)
    stage2_coefs.append(coef); stage2_infos.append(info)
    # does the history term survive an explicit calendar-time term?
    coef, info = fit_stage2(discrete, STAGE2_TERMS + ["z_event_time"],
                            "D2_plus_calendar_time", True, rng)
    stage2_coefs.append(coef); stage2_infos.append(info)
    # stratum contrast, WITHOUT the span term that defines the stratum
    strata_terms = [t for t in STAGE2_TERMS if t != "z_log_span"] + ["is_sustained"]
    coef, info = fit_stage2(stage2, strata_terms, "A1_stratum_contrast_no_span", True, rng)
    stage2_coefs.append(coef); stage2_infos.append(info)
    # sustained events are by construction the highest-exposure ones, so exposure
    # absorbs the stratum too; this is the unconditioned within-dyad contrast
    coef, info = fit_stage2(stage2, ["is_sustained"], "A2_stratum_contrast_alone", True, rng)
    stage2_coefs.append(coef); stage2_infos.append(info)

    stage2_coef = pd.concat(stage2_coefs, ignore_index=True)
    stage2_coef.to_csv(OUT_DIR / "stage2_coefficients.csv", index=False)
    report["stage2_models"] = stage2_infos
    report["stage2_design_note"] = (
        "12 dyads cannot support between-dyad inference, so the primary stage-2 "
        "model carries dyad fixed effects and is a within-dyad estimator; "
        "D1 drops them only to show how much the estimate moves"
    )
    report["stage2_offset_note"] = (
        "the response is cross-group edges out of total edges with "
        "logit(edge-weighted composition expectation) as an offset, so the "
        "intercept is the mean logit deficit relative to expectation"
    )

    (OUT_DIR / "phase4_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print()
    show = ["model", "term", "estimate", "std_error", "ci_low", "ci_high", "odds_ratio"]
    print("STAGE 1 -- non-group terms")
    print(stage1_coef[~stage1_coef["is_group_effect"]][show].round(4).to_string(index=False))
    print()
    print("STAGE 1 -- within-5km sensitivity, non-group terms")
    print(sens_coef[~sens_coef["is_group_effect"]][show].round(4).to_string(index=False))
    print()
    print("STAGE 2")
    print(stage2_coef.round(4).to_string(index=False))
    print()
    print(json.dumps({k: v for k, v in report.items()
                      if k in {"stage1_cross_validation", "stage1_sensitivity_within_5km"}},
                     indent=2))
    print()
    print("stage-2 collinearity check:",
          json.dumps(stage2_audit["collinearity_warnings"], indent=2))
    print("stage-2 mixed-model variance components (stage 1):",
          json.dumps([r.get("variance_components") for r in mixed_runs], indent=2))


if __name__ == "__main__":
    main()
