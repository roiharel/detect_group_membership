"""Is seasonality absent from DEPTH and DURATION too, or only from occurrence?

analyze_seasonality_and_trends.py tested one kind of response on each axis: whether a
crossing HAPPENS. That licensed the statement "opportunity is seasonal, the crossing
decision is not" -- but only for occurrence. Depth and duration were never tested, and
there is a specific reason to expect they might differ: the project's own central
finding is that frequency and depth have different predictors, and on two axes duration
is the thing that carries depth (axis A: log span is the only robust mixing predictor;
axis C: alone-only excursions last 1 night, joined 4, settlement 7+). If season acts on
duration, depth could be seasonal while occurrence is not.

This script tests, with the same annual harmonic and the same clustering:

  AXIS A   encounter duration (structural span)          1,705 Stage-1 events
           mixing depth given an encounter               751 events with 5 m bins
  AXIS B   split duration (persistent split hours)       860 candidate split events
  AXIS C   excursion duration (away-nights)              338 excursions
           excursion depth: reached another unit         338 excursions
           excursion depth: reached settlement           338 excursions

Every model is on the well-observed window (2024-12-01 onward), which holds about 1.6
seasonal cycles. Samples here are far smaller than the occurrence tables, so these are
weaker tests than the occurrence nulls -- a null below is a weak claim of absence and
is reported as such.

Outputs: outputs/general_structure_2026_09/phase4c_seasonality/depth_duration/
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
STAGE1 = BASE / "phase2_two_stage_events/stage1_events_with_stage2_mixing.csv"
EXCURSIONS = BASE / "phase4b_individual_axis/excursions_dominant_gap0.csv"
SPLITS = Path("outputs/presplit_cositting_prediction/candidate_split_events.csv")
OUT = BASE / "phase4c_seasonality/depth_duration"

WELL_OBSERVED_FROM = "2024-12-01"


def harmonic(d: pd.DataFrame, datecol: str) -> pd.DataFrame:
    d = d.copy()
    t = pd.to_datetime(d[datecol])
    if getattr(t.dtype, "tz", None) is not None:
        t = t.dt.tz_localize(None)
    d["_date"] = t
    doy = t.dt.dayofyear
    d["season_sin"] = np.sin(2 * np.pi * doy / 365.25)
    d["season_cos"] = np.cos(2 * np.pi * doy / 365.25)
    d["trend_months"] = (t - t.min()).dt.days / 30.44
    return d


def fit(d: pd.DataFrame, response: str, cluster: str, extra: list[str],
        family: str, label: str, offset: np.ndarray | None = None,
        weights: np.ndarray | None = None) -> dict:
    """Annual harmonic + trend + extra covariates, GEE clustered on `cluster`."""
    import statsmodels.api as sm
    terms = ["season_sin", "season_cos", "trend_months"] + list(extra)
    d = d.dropna(subset=[response] + terms).copy()
    if len(d) < 40 or d[cluster].nunique() < 3:
        return {"label": label, "n": int(len(d)),
                "clusters": int(d[cluster].nunique()) if len(d) else 0,
                "skipped": "fewer than 40 rows or fewer than 3 clusters"}
    X = d[terms].astype(float).copy()
    sds = {c: X[c].std() for c in terms}
    for c in terms:
        if sds[c] and sds[c] > 0:
            X[c] = (X[c] - X[c].mean()) / sds[c]
    X = sm.add_constant(X)
    fam = {"binomial": sm.families.Binomial(),
           "gaussian": sm.families.Gaussian()}[family]
    kw = {}
    if offset is not None:
        kw["offset"] = np.asarray(offset, dtype=float)[d.index.map(
            lambda i: True).values] if False else np.asarray(offset, dtype=float)
    m = sm.GEE(d[response].astype(float), X, groups=d[cluster], family=fam,
               cov_struct=sm.cov_struct.Independence(),
               weights=None if weights is None else np.asarray(weights, dtype=float),
               **kw).fit()
    ci = m.conf_int()
    a = m.params["season_sin"] / (sds["season_sin"] or 1)
    b = m.params["season_cos"] / (sds["season_cos"] or 1)
    amp = float(np.hypot(a, b))
    peak = float((np.arctan2(a, b) * 365.25 / (2 * np.pi)) % 365.25)
    R = np.zeros((2, len(terms) + 1))
    R[0, 1] = 1
    R[1, 2] = 1
    try:
        pval = round(float(m.wald_test(R, scalar=True).pvalue), 4)
    except Exception:
        pval = None
    binom = family == "binomial"
    return {
        "label": label,
        "n": int(len(d)),
        "clusters": int(d[cluster].nunique()),
        "response_mean": round(float(d[response].mean()), 4),
        "terms": {
            t: {
                ("or_per_sd" if binom else "effect_per_sd"): round(
                    float(np.exp(m.params[t]) if binom else m.params[t]), 3),
                "ci": [round(float(np.exp(ci.loc[t, 0]) if binom else ci.loc[t, 0]), 3),
                       round(float(np.exp(ci.loc[t, 1]) if binom else ci.loc[t, 1]), 3)],
            } for t in terms
        },
        "seasonal_amplitude": round(amp, 4),
        "seasonal_peak_month": (pd.Timestamp("2025-01-01")
                                + pd.Timedelta(days=peak - 1)).strftime("%b"),
        "seasonal_joint_wald_p": pval,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report: dict = {
        "well_observed_from": WELL_OBSERVED_FROM,
        "scope_note": (
            "Occurrence responses were tested in analyze_seasonality_and_trends.py. "
            "This script tests DEPTH and DURATION, which that script did not. Samples "
            "are one to two orders of magnitude smaller, so nulls here are weaker "
            "claims of absence than the occurrence nulls."
        ),
    }

    # ---------------------------------------------------------------- axis A
    a = pd.read_csv(STAGE1)
    a = harmonic(a, "start_hour")
    a = a[a["_date"].ge(WELL_OBSERVED_FROM)].copy()
    a["log_span"] = np.log(a["structural_span_hours"].clip(lower=1))
    axa: dict = {"stage1_events_in_window": int(len(a)),
                 "dyads": int(a["pair_key"].nunique())}
    axa["encounter_duration"] = fit(
        a, "log_span", "pair_key", ["observed_hour_fraction"], "gaussian",
        "log encounter structural span; Gaussian GEE clustered on dyad")

    # depth: cross-edge fraction against the composition expectation carried in the
    # product, on events with 5 m bins. Deficit on the logit scale, matching Phase 4.
    f = a[a["has_finescale_5m"].astype(str).str.lower().isin(["true", "1"])].copy()
    f = f[f["5m_total_edges"].gt(0) & f["5m_mean_shuffle_expected"].between(1e-4, 1 - 1e-4)]
    if len(f):
        obs = (f["5m_cross_edges"] / f["5m_total_edges"]).clip(1e-4, 1 - 1e-4)
        exp = f["5m_mean_shuffle_expected"].clip(1e-4, 1 - 1e-4)
        f["logit_deficit"] = np.log(obs / (1 - obs)) - np.log(exp / (1 - exp))
        axa["events_with_5m"] = int(len(f))
        axa["mixing_depth"] = fit(
            f, "logit_deficit", "pair_key", ["log_span"], "gaussian",
            "logit mixing deficit given an encounter; Gaussian GEE clustered on dyad")
        axa["mixing_depth_no_span_term"] = fit(
            f, "logit_deficit", "pair_key", [], "gaussian",
            "same, without the duration term -- season could act through duration")
        f.to_csv(args.output_dir / "axis_a_depth_events.csv", index=False)
    report["axis_a"] = axa

    # ---------------------------------------------------------------- axis B
    b = pd.read_csv(SPLITS)
    b = harmonic(b, "split_onset")
    b = b[b["_date"].ge(WELL_OBSERVED_FROM)].copy()
    b["log_hours"] = np.log(b["persistent_split_hours"].clip(lower=1))
    report["axis_b"] = {
        "split_events_in_window": int(len(b)),
        "units": int(b["dynamic_social_unit"].nunique()),
        "median_persistent_split_hours": float(b["persistent_split_hours"].median()),
        "split_duration": fit(
            b, "log_hours", "dynamic_social_unit", ["n_animals_onset"], "gaussian",
            "log persistent split hours; Gaussian GEE clustered on unit"),
        "source_caveat": (
            "candidate_split_events.csv derives from the presplit pipeline, not the "
            "frozen export; gap 7.3 applies"
        ),
    }

    # ---------------------------------------------------------------- axis C
    c = pd.read_csv(EXCURSIONS)
    c = harmonic(c, "start_night")
    c = c[c["_date"].ge(WELL_OBSERVED_FROM)].copy()
    c["log_away"] = np.log(c["away_nights"].clip(lower=1))
    c["reached_other"] = c["other_nights"].gt(0).astype(int)
    c["reached_settlement"] = c["settlement_candidate"].astype(str).str.lower().isin(
        ["true", "1"]).astype(int)
    report["axis_c"] = {
        "excursions_in_window": int(len(c)),
        "units": int(c["origin_group"].nunique()),
        "animals": int(c["animal_id"].nunique()),
        "excursion_duration": fit(
            c, "log_away", "origin_group", [], "gaussian",
            "log away-nights per excursion; Gaussian GEE clustered on unit"),
        "depth_reached_another_unit": fit(
            c, "reached_other", "origin_group", [], "binomial",
            "did the excursion reach another unit? Binomial GEE clustered on unit"),
        "depth_reached_settlement": fit(
            c, "reached_settlement", "origin_group", [], "binomial",
            "did the excursion reach settlement (>=7 joined nights)? "
            "Binomial GEE clustered on unit"),
    }

    with open(args.output_dir / "depth_duration_seasonality.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    # ---------------------------------------------------------------- print
    def show(block, keys):
        for k in keys:
            m = block.get(k)
            if not m:
                continue
            if "skipped" in m:
                print("  %-46s SKIPPED (%s)" % (k, m["skipped"]))
                continue
            print("  %-46s n=%5d clus=%3d  amplitude %.3f  peak %s  joint p = %s"
                  % (k, m["n"], m["clusters"], m["seasonal_amplitude"],
                     m["seasonal_peak_month"], m["seasonal_joint_wald_p"]))
            for t, v in m["terms"].items():
                kk = "or_per_sd" if "or_per_sd" in v else "effect_per_sd"
                print("        %-24s %-14s %.3f [%.3f, %.3f]"
                      % (t, kk, v[kk], v["ci"][0], v["ci"][1]))

    print("=" * 78)
    print("SEASONALITY IN DEPTH AND DURATION -- the responses occurrence tests missed")
    print("=" * 78)
    print("\nAXIS A  (%d Stage-1 events, %d dyads in window)"
          % (report["axis_a"]["stage1_events_in_window"], report["axis_a"]["dyads"]))
    show(report["axis_a"], ["encounter_duration", "mixing_depth",
                            "mixing_depth_no_span_term"])
    print("\nAXIS B  (%d split events, %d units)"
          % (report["axis_b"]["split_events_in_window"], report["axis_b"]["units"]))
    show(report["axis_b"], ["split_duration"])
    print("\nAXIS C  (%d excursions, %d animals, %d units)"
          % (report["axis_c"]["excursions_in_window"], report["axis_c"]["animals"],
             report["axis_c"]["units"]))
    show(report["axis_c"], ["excursion_duration", "depth_reached_another_unit",
                            "depth_reached_settlement"])
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
