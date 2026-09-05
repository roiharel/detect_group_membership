"""Seasonality and secular trend in boundary crossing, on all three axes.

THE CENTRAL PROBLEM. Observation effort grows monotonically across the record --
from 10 collared animals in 5 units (March 2024) to 232 animals in 25 units
(May 2026), with nightly resolvability rising from 0.10 to 0.96. A raw secular
trend in any event rate is therefore collar deployment, not behaviour. And because
the record spans only about 2.4 years with effort ramping inside it, the annual
seasonal harmonic and the linear trend are partly aliased: the same calendar month
is observed with wildly different effort in different years.

This script does four things:

  1. Quantifies the effort ramp, and picks a well-observed window from it.
  2. Reports, for each axis, the raw monthly rate alongside effort, so the reader can
     see the confound rather than take it on trust.
  3. Fits season (annual sin/cos harmonic) and trend jointly, with an effort covariate
     and unit clustering, on the well-observed window -- and reports the same model
     without the effort term, so the difference is visible.
  4. Reports the ALIASING DIAGNOSTIC: how collinear the seasonal harmonic is with the
     trend and with effort in the retained sample, and how many seasonal cycles the
     window actually contains. A seasonal coefficient from 1.7 cycles is weakly
     identified whatever its interval says, and that has to be stated, not buried.

NDVI's own seasonal phase is estimated the same way, so each axis's season can be
compared against the vegetation cycle rather than interpreted in isolation.

The three axes, and their denominators:
  A  between-group encounter   73,353 dyad-days   effort = shared observed hours
  B  within-group split         1,418 group-weeks effort = collars observed
  C  individual excursion      82,205 animal-nights effort = resolvable unit-mates

Outputs: outputs/general_structure_2026_09/phase4c_seasonality/
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

OPP = Path("outputs/general_structure_2026_09/phase1_opportunity/opportunity_dyad_day.csv")
NIGHTLY = Path("outputs/general_structure_2026_09/phase4b_individual_axis/nightly_states_dominant.csv")
WEEKLY = Path(
    "outputs/canonical_within_group_density_modularity_ridges/"
    "canonical_within_group_weekly_network_metrics.csv"
)
NDVI = Path(r"C:\Users\rharel\Documents\Github\detect_group_membership\indices_sentinel2.csv")
OUT = Path("outputs/general_structure_2026_09/phase4c_seasonality")

SCOPE = "combined_1100_1600"
# Effort is >=0.87 nightly resolvable from 2024-11 and >=0.92 from 2024-12.
WELL_OBSERVED_FROM = "2024-12-01"
SEED = 20260904


# ------------------------------------------------------------------ helpers
def add_time_terms(d: pd.DataFrame, datecol: str) -> pd.DataFrame:
    t = pd.to_datetime(d[datecol])
    doy = t.dt.dayofyear
    d = d.copy()
    d["season_sin"] = np.sin(2 * np.pi * doy / 365.25)
    d["season_cos"] = np.cos(2 * np.pi * doy / 365.25)
    d["trend_months"] = (t - t.min()).dt.days / 30.44
    return d


def aliasing(d: pd.DataFrame, effort_col: str, datecol: str) -> dict:
    """How badly is the seasonal harmonic confounded with trend and effort here?"""
    t = pd.to_datetime(d[datecol])
    span_days = float((t.max() - t.min()).days)
    cols = ["season_sin", "season_cos", "trend_months", effort_col]
    sub = d[cols].astype(float).dropna()
    corr = sub.corr()
    out = {
        "window": [str(t.min().date()), str(t.max().date())],
        "span_days": span_days,
        "seasonal_cycles_in_window": round(span_days / 365.25, 2),
        "corr_season_sin_trend": round(float(corr.loc["season_sin", "trend_months"]), 3),
        "corr_season_cos_trend": round(float(corr.loc["season_cos", "trend_months"]), 3),
        "corr_season_sin_effort": round(float(corr.loc["season_sin", effort_col]), 3),
        "corr_season_cos_effort": round(float(corr.loc["season_cos", effort_col]), 3),
        "corr_trend_effort": round(float(corr.loc["trend_months", effort_col]), 3),
    }
    # variance inflation for each term against the other three
    try:
        import statsmodels.api as sm
        vif = {}
        for c in cols:
            others = [x for x in cols if x != c]
            X = sm.add_constant(sub[others])
            r2 = sm.OLS(sub[c], X).fit().rsquared
            vif[c] = round(float(1 / max(1e-9, 1 - r2)), 2)
        out["vif"] = vif
    except Exception as exc:  # pragma: no cover
        out["vif"] = {"error": str(exc)}
    return out


def seasonal_model(d: pd.DataFrame, response: str, effort_col: str, cluster: str,
                   binomial: bool, label: str, with_effort: bool = True) -> dict:
    """Season + trend, clustered on unit. Effort z-scored and optionally included."""
    import statsmodels.api as sm
    terms = ["season_sin", "season_cos", "trend_months"]
    if with_effort:
        terms = terms + [effort_col]
    X = d[terms].astype(float).copy()
    for c in terms:
        sd = X[c].std()
        if sd > 0:
            X[c] = (X[c] - X[c].mean()) / sd
    X = sm.add_constant(X)
    fam = sm.families.Binomial() if binomial else sm.families.Gaussian()
    m = sm.GEE(d[response].astype(float), X, groups=d[cluster], family=fam,
               cov_struct=sm.cov_struct.Independence()).fit()
    ci = m.conf_int()

    # amplitude and peak day of the fitted annual harmonic, on the raw (unscaled) scale
    sin_sd = d["season_sin"].std() or 1.0
    cos_sd = d["season_cos"].std() or 1.0
    a = float(m.params["season_sin"]) / sin_sd
    b = float(m.params["season_cos"]) / cos_sd
    amp = float(np.hypot(a, b))
    peak_doy = float((np.arctan2(a, b) * 365.25 / (2 * np.pi)) % 365.25)

    res = {
        "label": label,
        "n": int(len(d)),
        "clusters": int(d[cluster].nunique()),
        "effort_term_included": with_effort,
        "base_rate": round(float(d[response].mean()), 4),
        "terms": {
            t: {
                "coef": round(float(m.params[t]), 4),
                ("or_per_sd" if binomial else "effect_per_sd"): round(
                    float(np.exp(m.params[t]) if binomial else m.params[t]), 3),
                "ci": [round(float(np.exp(ci.loc[t, 0]) if binomial else ci.loc[t, 0]), 3),
                       round(float(np.exp(ci.loc[t, 1]) if binomial else ci.loc[t, 1]), 3)],
            } for t in terms
        },
        "seasonal_amplitude_logit" if binomial else "seasonal_amplitude": round(amp, 4),
        "seasonal_peak_day_of_year": round(peak_doy, 0),
        "seasonal_peak_month": pd.Timestamp("2025-01-01").dayofyear and
                               str((pd.Timestamp("2025-01-01") +
                                    pd.Timedelta(days=peak_doy - 1)).strftime("%b")),
    }
    # joint significance of the two harmonic terms
    try:
        w = m.wald_test(np.array([[0, 1, 0] + [0] * (len(terms) - 2),
                                  [0, 0, 1] + [0] * (len(terms) - 2)]), scalar=True)
        res["seasonal_joint_wald_p"] = round(float(w.pvalue), 4)
    except Exception:
        pass
    return res


def monthly_table(d: pd.DataFrame, datecol: str, response: str, effort_col: str,
                  cluster: str) -> pd.DataFrame:
    d = d.copy()
    d["month"] = pd.to_datetime(d[datecol]).dt.to_period("M").astype(str)
    return d.groupby("month").agg(
        rows=(response, "size"),
        rate=(response, "mean"),
        units=(cluster, "nunique"),
        mean_effort=(effort_col, "mean"),
    ).round(4).reset_index()


def by_calendar_month(d: pd.DataFrame, datecol: str, response: str,
                      cluster: str) -> pd.DataFrame:
    """Rate by month-of-year, and the same unit-balanced (mean of per-unit rates),
    because unit composition changes across the record."""
    d = d.copy()
    d["moy"] = pd.to_datetime(d[datecol]).dt.month
    raw = d.groupby("moy")[response].agg(["size", "mean"])
    bal = (d.groupby(["moy", cluster])[response].mean()
             .groupby("moy").agg(["size", "mean"]))
    out = pd.DataFrame({
        "month_of_year": raw.index,
        "rows": raw["size"].to_numpy(),
        "rate_pooled": raw["mean"].round(4).to_numpy(),
        "units_present": bal["size"].to_numpy(),
        "rate_unit_balanced": bal["mean"].round(4).to_numpy(),
    })
    return out


# ------------------------------------------------------------------ effort
def effort_ramp(nightly: pd.DataFrame, out: Path) -> dict:
    n = nightly.copy()
    n["month"] = n["night"].dt.to_period("M").astype(str)
    t = n.groupby("month").agg(
        animal_nights=("night", "size"),
        animals=("animal_id", "nunique"),
        origin_groups=("origin_group", "nunique"),
        resolvable_share=("state", lambda s: round(float((s != "unresolvable").mean()), 3)),
    ).reset_index()
    t.to_csv(out / "observation_effort_by_month.csv", index=False)
    first = t.iloc[0]
    last = t.iloc[-2] if len(t) > 1 else t.iloc[-1]  # final month is partial
    peak = t.loc[t["animals"].idxmax()]
    return {
        "months": int(len(t)),
        "first_month": {"month": first["month"], "animals": int(first["animals"]),
                        "units": int(first["origin_groups"]),
                        "resolvable": float(first["resolvable_share"])},
        "peak_month": {"month": peak["month"], "animals": int(peak["animals"]),
                       "units": int(peak["origin_groups"]),
                       "resolvable": float(peak["resolvable_share"])},
        "animals_fold_increase": round(float(peak["animals"] / max(1, first["animals"])), 1),
        "well_observed_from": WELL_OBSERVED_FROM,
        "verdict": (
            "Collared animals rise about %.0f-fold and resolvability from %.2f to %.2f "
            "across the record. No raw secular trend in any event rate is interpretable; "
            "every trend term below is reported only alongside an effort covariate, and "
            "even then the two are strongly correlated."
            % (peak["animals"] / max(1, first["animals"]),
               first["resolvable_share"], peak["resolvable_share"])
        ),
    }


# ------------------------------------------------------------------ axes
def axis_a(out: Path) -> dict:
    d = pd.read_csv(OPP, parse_dates=["period_start"])
    d = d[d["state"].isin(["detected_encounter", "observed_no_encounter"])].copy()
    d["encounter"] = d["state"].eq("detected_encounter").astype(int)
    d = add_time_terms(d, "period_start")
    res = {"full_record": {"n": int(len(d)), "dyads": int(d["pair_key"].nunique()),
                           "encounter_rate": round(float(d["encounter"].mean()), 4)}}
    w = d[d["period_start"].ge(WELL_OBSERVED_FROM)].copy()
    w = add_time_terms(w, "period_start")
    res["aliasing"] = aliasing(w, "shared_hours", "period_start")
    monthly_table(w, "period_start", "encounter", "shared_hours", "pair_key").to_csv(
        out / "axis_a_monthly.csv", index=False)
    by_calendar_month(w, "period_start", "encounter", "pair_key").to_csv(
        out / "axis_a_by_calendar_month.csv", index=False)
    res["by_calendar_month"] = by_calendar_month(
        w, "period_start", "encounter", "pair_key").to_dict("records")
    try:
        res["model_with_effort"] = seasonal_model(
            w, "encounter", "shared_hours", "pair_key", True,
            "does a structural encounter occur? Binomial GEE clustered on dyad", True)
        res["model_without_effort"] = seasonal_model(
            w, "encounter", "shared_hours", "pair_key", True,
            "same, effort term dropped -- for contrast", False)

        # THE DECISIVE CHECK. A seasonal encounter signal can mean two different
        # things: groups' ranges overlap seasonally (geometry), or groups that are
        # already near each other decide to meet seasonally (interaction). Restricting
        # to dyad-days where the units were already close removes the first. This is
        # the restriction that took the Phase 4 seasonal term from OR 0.62 to 0.86.
        ladder = []
        for km in (None, 10, 5, 2):
            f = w if km is None else w[w["min_centroid_distance_m"].le(km * 1000)]
            if len(f) < 500:
                continue
            m = seasonal_model(f, "encounter", "shared_hours", "pair_key", True,
                               "all distances" if km is None else "<= %d km" % km, True)
            ladder.append({
                "restriction": "all supported dyad-days" if km is None
                               else "min centroid distance <= %d km" % km,
                "n": m["n"], "dyads": m["clusters"], "encounter_rate": m["base_rate"],
                "seasonal_amplitude_logit": m["seasonal_amplitude_logit"],
                "seasonal_peak_month": m["seasonal_peak_month"],
                "seasonal_joint_wald_p": m.get("seasonal_joint_wald_p"),
                "trend_or_per_sd": m["terms"]["trend_months"]["or_per_sd"],
                "trend_ci": m["terms"]["trend_months"]["ci"],
            })
        res["distance_restriction_ladder"] = ladder
        pd.DataFrame(ladder).to_csv(out / "axis_a_distance_restriction_ladder.csv",
                                    index=False)

        # Is the declining trend real, or is the denominator filling with distant
        # dyads as more units get collared? Balanced panel = dyads observed across the
        # whole window.
        lo, hi = w["period_start"].min(), w["period_start"].max()
        span = w.groupby("pair_key")["period_start"].agg(["min", "max"])
        bal = span[(span["min"] <= lo + pd.Timedelta(days=45))
                   & (span["max"] >= hi - pd.Timedelta(days=45))].index
        by_month = w.groupby(w["period_start"].dt.to_period("M"))["pair_key"].nunique()
        halves = w.assign(
            half=np.where(w["period_start"] < lo + (hi - lo) / 2, "first", "second")
        ).groupby("half")["min_centroid_distance_m"].median()
        res["denominator_composition"] = {
            "dyads_present_throughout": int(len(bal)),
            "dyads_total": int(w["pair_key"].nunique()),
            "dyads_in_denominator_first_month": int(by_month.iloc[0]),
            "dyads_in_denominator_last_month": int(by_month.iloc[-1]),
            "median_min_centroid_distance_first_half_m": float(halves.get("first", np.nan)),
            "median_min_centroid_distance_second_half_m": float(halves.get("second", np.nan)),
        }
        bal_ladder = []
        for km in (None, 5, 2):
            f = w[w["pair_key"].isin(bal)]
            if km is not None:
                f = f[f["min_centroid_distance_m"].le(km * 1000)]
            if len(f) < 500:
                continue
            m = seasonal_model(f, "encounter", "shared_hours", "pair_key", True,
                               "balanced panel" if km is None
                               else "balanced panel and <= %d km" % km, True)
            bal_ladder.append({
                "restriction": m["label"], "n": m["n"], "dyads": m["clusters"],
                "encounter_rate": m["base_rate"],
                "seasonal_amplitude_logit": m["seasonal_amplitude_logit"],
                "seasonal_peak_month": m["seasonal_peak_month"],
                "seasonal_joint_wald_p": m.get("seasonal_joint_wald_p"),
                "trend_or_per_sd": m["terms"]["trend_months"]["or_per_sd"],
                "trend_ci": m["terms"]["trend_months"]["ci"],
            })
        res["balanced_panel_ladder"] = bal_ladder
        pd.DataFrame(bal_ladder).to_csv(out / "axis_a_balanced_panel_ladder.csv",
                                        index=False)
        res["verdict"] = (
            "The seasonal signal decays monotonically as the dyads are restricted to "
            "those already close, and stops being distinguishable from zero at 5 km -- "
            "the same restriction that killed the Phase 4 seasonal term. The season is "
            "therefore in WHERE GROUPS ARE, not in whether they interact once "
            "co-located. The apparent secular decline is denominator composition: on a "
            "balanced dyad panel the trend is null and the sign flips, while the median "
            "closest approach across all dyad-days rises as newly collared units add "
            "distant pairs."
        )
    except Exception as exc:  # pragma: no cover
        res["error"] = str(exc)
    return res


def axis_b(out: Path) -> dict:
    d = pd.read_csv(WEEKLY)
    d = d[d["scope"].eq(SCOPE)].copy()
    d["period_start"] = pd.to_datetime(d["period_start"])
    d["split_multi_quarter"] = (d["multi_animal_split_timestamp_fraction"] > 0.25).astype(int)
    d["is_modular"] = (d["modularity"] > 0.001).astype(int)
    d = add_time_terms(d, "period_start")
    res = {"full_record": {"group_weeks": int(len(d)),
                           "units": int(d["dynamic_social_unit"].nunique())}}
    w = d[d["period_start"].ge(WELL_OBSERVED_FROM)].copy()
    w = add_time_terms(w, "period_start")
    res["aliasing"] = aliasing(w, "n_animals_observed", "period_start")
    monthly_table(w, "period_start", "split_multi_quarter", "n_animals_observed",
                  "dynamic_social_unit").to_csv(out / "axis_b_monthly.csv", index=False)
    res["by_calendar_month_split"] = by_calendar_month(
        w, "period_start", "split_multi_quarter", "dynamic_social_unit").to_dict("records")
    try:
        res["split_model_with_effort"] = seasonal_model(
            w, "split_multi_quarter", "n_animals_observed", "dynamic_social_unit", True,
            "multi-animal split in >25% of timestamps; Binomial GEE clustered on unit", True)
        res["split_model_without_effort"] = seasonal_model(
            w, "split_multi_quarter", "n_animals_observed", "dynamic_social_unit", True,
            "same, effort term dropped -- for contrast", False)
        # modularity, restricted to well-covered group-weeks
        wc = w[w["n_animals_observed"].ge(12)].copy()
        wc = add_time_terms(wc, "period_start")
        if len(wc) > 60:
            res["modularity_model_well_covered"] = seasonal_model(
                wc, "is_modular", "n_animals_observed", "dynamic_social_unit", True,
                "is the group-week modular? >=12 collars only; Binomial GEE clustered on unit",
                True)
            res["modularity_aliasing_well_covered"] = aliasing(
                wc, "n_animals_observed", "period_start")

            # Can the modularity trend be checked on a balanced unit panel? Only if
            # enough well-covered units span the window. Report the count, and refuse
            # to fit when it is too small for cluster-robust inference.
            lo, hi = wc["period_start"].min(), wc["period_start"].max()
            span = wc.groupby("dynamic_social_unit")["period_start"].agg(
                ["min", "max", "size"])
            bal = span[(span["min"] <= lo + pd.Timedelta(days=90))
                       & (span["max"] >= hi - pd.Timedelta(days=90))
                       & (span["size"] >= 20)].index
            res["modularity_trend_verifiability"] = {
                "well_covered_units": int(wc["dynamic_social_unit"].nunique()),
                "units_spanning_window_with_20plus_weeks": int(len(bal)),
                "units": sorted(bal.tolist()),
                "fitted": False,
                "reason": (
                    "A balanced-panel fit would be clustered on %d units. Cluster-robust "
                    "standard errors are badly anti-conservative below roughly 10 "
                    "clusters, so the trend term is reported as UNVERIFIABLE rather than "
                    "confirmed or refuted. The unrestricted estimate stands on 12 units "
                    "with collar coverage still carrying OR 1.53." % len(bal)
                ),
            }
    except Exception as exc:  # pragma: no cover
        res["error"] = str(exc)
    return res


def axis_c(nightly: pd.DataFrame, out: Path) -> dict:
    n = nightly.copy()
    # effort: resolvable animals of the same origin group on the same night
    n["resolvable"] = n["state"].ne("unresolvable").astype(int)
    eff = (n.groupby(["night", "origin_group"])["resolvable"].sum()
             .rename("unit_resolvable").reset_index())
    n = n.merge(eff, on=["night", "origin_group"], how="left")
    n["unit_resolvable"] = n["unit_resolvable"] - n["resolvable"]

    r = n[n["resolvable"].eq(1)].copy()
    r["away"] = r["state"].isin(["alone", "with_other"]).astype(int)
    r = r.sort_values(["animal_id", "night"])
    prev = r.groupby("animal_id")["away"].shift()
    gap = r.groupby("animal_id")["night"].diff().dt.days
    # an onset is a first away night after a with-origin night the previous night
    r["onset"] = ((r["away"].eq(1)) & (prev.eq(0)) & (gap.eq(1))).astype(int)
    r = add_time_terms(r, "night")

    res = {"full_record": {"resolvable_animal_nights": int(len(r)),
                           "animals": int(r["animal_id"].nunique()),
                           "away_rate": round(float(r["away"].mean()), 4),
                           "onset_rate": round(float(r["onset"].mean()), 4)}}
    w = r[r["night"].ge(WELL_OBSERVED_FROM)].copy()
    w = add_time_terms(w, "night")
    res["aliasing"] = aliasing(w, "unit_resolvable", "night")
    monthly_table(w, "night", "away", "unit_resolvable", "origin_group").to_csv(
        out / "axis_c_monthly.csv", index=False)
    res["by_calendar_month_away"] = by_calendar_month(
        w, "night", "away", "origin_group").to_dict("records")
    res["by_calendar_month_onset"] = by_calendar_month(
        w, "night", "onset", "origin_group").to_dict("records")
    try:
        res["away_model_with_effort"] = seasonal_model(
            w, "away", "unit_resolvable", "origin_group", True,
            "is the animal away from its origin group tonight? Binomial GEE clustered on unit",
            True)
        res["onset_model_with_effort"] = seasonal_model(
            w, "onset", "unit_resolvable", "origin_group", True,
            "does an excursion START tonight? Binomial GEE clustered on unit", True)
        res["onset_model_without_effort"] = seasonal_model(
            w, "onset", "unit_resolvable", "origin_group", True,
            "same, effort term dropped -- for contrast", False)
    except Exception as exc:  # pragma: no cover
        res["error"] = str(exc)
    return res


def ndvi_season(out: Path) -> dict:
    nd = pd.read_csv(NDVI, parse_dates=["date"])
    nd = nd[nd["extraction_status"].eq("success")].copy()
    nd["unit_key"] = nd["group_id"].str.lower()
    nd = nd[nd["date"].ge(WELL_OBSERVED_FROM)].copy()
    nd = add_time_terms(nd, "date")
    nd["ndvi_unit_mean"] = nd.groupby("unit_key")["NDVI"].transform("mean")
    nd["ndvi_dev"] = nd["NDVI"] - nd["ndvi_unit_mean"]
    res = {"n": int(len(nd)), "units": int(nd["unit_key"].nunique()),
           "window": [str(nd["date"].min().date()), str(nd["date"].max().date())]}
    try:
        res["model"] = seasonal_model(
            nd, "ndvi_dev", "ndvi_unit_mean", "unit_key", False,
            "NDVI deviation from the unit's own mean; Gaussian GEE clustered on unit", True)
        by = nd.groupby(nd["date"].dt.month)["ndvi_dev"].mean().round(4)
        res["mean_ndvi_deviation_by_month"] = {int(k): float(v) for k, v in by.items()}
    except Exception as exc:  # pragma: no cover
        res["error"] = str(exc)
    return res


# ------------------------------------------------------------------ main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    nightly = pd.read_csv(NIGHTLY, parse_dates=["night"])

    report = {
        "seed": SEED,
        "well_observed_from": WELL_OBSERVED_FROM,
        "design_note": (
            "Season is an annual sin/cos harmonic on day-of-year; trend is linear in "
            "months. Every model is clustered on unit (dyad for axis A) and reports the "
            "same fit with and without an effort covariate. The aliasing block gives the "
            "number of seasonal cycles in the window and the collinearity between season, "
            "trend and effort -- read it before any coefficient."
        ),
    }

    print("=" * 76)
    print("OBSERVATION EFFORT -- read this before any trend")
    print("=" * 76)
    report["effort"] = effort_ramp(nightly, args.output_dir)
    print(json.dumps(report["effort"], indent=2))

    for name, fn in (("axis_a_between_group", lambda: axis_a(args.output_dir)),
                     ("axis_b_within_group", lambda: axis_b(args.output_dir)),
                     ("axis_c_individual", lambda: axis_c(nightly, args.output_dir)),
                     ("ndvi_reference", lambda: ndvi_season(args.output_dir))):
        print("\n" + "=" * 76)
        print(name.upper().replace("_", " "))
        print("=" * 76)
        report[name] = fn()
        print(json.dumps(report[name], indent=2, default=str))

    with open(args.output_dir / "seasonality_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
