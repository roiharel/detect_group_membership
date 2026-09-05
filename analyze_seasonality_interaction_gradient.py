"""Is the TYPE of interaction seasonal, across the whole radius gradient?

The earlier seasonality work collapsed "how they interacted" into one number: the
edge-weighted cross-group contact fraction at 5 m. That is a thin reading of the
question. When two groups occupy the same spatial cluster -- a scale of hundreds of
metres -- what happens inside it spans a gradient:

  ~10^2 m  cluster co-membership   how much of EACH group actually joins
     5 m   close association       cross-group contact rate, and whether the merged
                                   mass holds two origin communities or one
     2 m   near-contact            the same, at a radius where only real proximity
                                   survives
   ratio   2 m / 5 m               how much of the 5 m structure survives tightening

Each level is a different social question. A merge can be large but structurally
segregated, or small but thoroughly mixed. Testing only the 5 m cross-fraction cannot
distinguish those, so this script tests the annual harmonic against every level, plus
the within-merge network descriptors (modularity, entropy, balance) that say whether
the merged group behaves as one community or two.

POWER DECAYS DOWN THE GRADIENT, and that is itself a result:
    cluster scale   1,705 events, 68 dyads
    5 m             751 events, 12 dyads
    2 m             668 events, 12 dyads
    betweenness     1 dyad (Copper-Lilac), 77 weeks, 4 radii
So a null at the cluster scale is a much stronger claim than a null at 2 m.

Outputs: outputs/general_structure_2026_09/phase4c_seasonality/interaction_gradient/
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

STAGE1 = Path("outputs/general_structure_2026_09/phase2_two_stage_events/"
              "stage1_events_with_stage2_mixing.csv")
BROKER_WEEKLY = Path("outputs/brokerage_individuals_2026-09-03/weekly_brokerage_by_animal.csv")
BROKER_CONC = Path("outputs/brokerage_individuals_2026-09-03/broker_concentration_weekly.csv")
OUT = Path("outputs/general_structure_2026_09/phase4c_seasonality/interaction_gradient")

WELL_OBSERVED_FROM = "2024-12-01"
MERGE_SCALE_ORDER = {"small_subset_merge": 0, "medium_partial_merge": 1, "large_merge": 2}


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
        family: str, label: str) -> dict:
    import statsmodels.api as sm
    terms = ["season_sin", "season_cos", "trend_months"] + list(extra)
    f = d.dropna(subset=[response] + terms).copy()
    if len(f) < 40 or f[cluster].nunique() < 3:
        return {"label": label, "n": int(len(f)),
                "clusters": int(f[cluster].nunique()) if len(f) else 0,
                "skipped": "under 40 rows or under 3 clusters"}
    X = f[terms].astype(float).copy()
    sds = {c: X[c].std() for c in terms}
    for c in terms:
        if sds[c] and sds[c] > 0:
            X[c] = (X[c] - X[c].mean()) / sds[c]
    X = sm.add_constant(X)
    fam = {"binomial": sm.families.Binomial(),
           "gaussian": sm.families.Gaussian()}[family]
    m = sm.GEE(f[response].astype(float), X, groups=f[cluster], family=fam,
               cov_struct=sm.cov_struct.Independence()).fit()
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
        "n": int(len(f)),
        "clusters": int(f[cluster].nunique()),
        "response_mean": round(float(f[response].mean()), 4),
        "seasonal_amplitude": round(amp, 4),
        "seasonal_peak_month": (pd.Timestamp("2025-01-01")
                                + pd.Timedelta(days=peak - 1)).strftime("%b"),
        "seasonal_joint_wald_p": pval,
        "terms": {
            t: {("or_per_sd" if binom else "effect_per_sd"): round(
                    float(np.exp(m.params[t]) if binom else m.params[t]), 3),
                "ci": [round(float(np.exp(ci.loc[t, 0]) if binom else ci.loc[t, 0]), 3),
                       round(float(np.exp(ci.loc[t, 1]) if binom else ci.loc[t, 1]), 3)]}
            for t in terms
        },
    }


def logit(x: pd.Series, lo: float = 1e-4) -> pd.Series:
    p = x.clip(lo, 1 - lo)
    return np.log(p / (1 - p))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    a = pd.read_csv(STAGE1)
    a = harmonic(a, "start_hour")
    a = a[a["_date"].ge(WELL_OBSERVED_FROM)].copy()
    a["log_span"] = np.log(a["structural_span_hours"].clip(lower=1))
    a["mutual_participation"] = a[["max_frac_a", "max_frac_b"]].min(axis=1)
    a["logit_participation"] = logit(a["mutual_participation"])
    a["is_large_merge"] = a["merge_scale"].eq("large_merge").astype(int)
    a["merge_scale_ord"] = a["merge_scale"].map(MERGE_SCALE_ORDER)

    for r in ("5m", "2m"):
        obs = a[f"{r}_edge_weighted_cross_fraction"]
        exp = a["5m_mean_shuffle_expected"] if r == "5m" else a["2m_mean_shuffle_expected"]
        ok = a[f"{r}_total_edges"].gt(0) & exp.between(1e-4, 1 - 1e-4)
        a[f"{r}_logit_deficit"] = np.where(ok, logit(obs) - logit(exp), np.nan)
    # how much of the 5 m cross-contact survives tightening to 2 m
    both = a["5m_edge_weighted_cross_fraction"].gt(0) & a["2m_edge_weighted_cross_fraction"].notna()
    a["radius_ratio_2m_over_5m"] = np.where(
        both, a["2m_edge_weighted_cross_fraction"] / a["5m_edge_weighted_cross_fraction"], np.nan)
    a["log_radius_ratio"] = np.where(a["radius_ratio_2m_over_5m"].gt(0),
                                     np.log(a["radius_ratio_2m_over_5m"]), np.nan)

    a.to_csv(args.output_dir / "stage1_events_gradient.csv", index=False)

    report: dict = {
        "well_observed_from": WELL_OBSERVED_FROM,
        "design_note": (
            "Each level of the radius gradient is a separate social question. Power "
            "decays down the gradient (68 dyads at the cluster scale, 12 at 5 m and "
            "2 m, 1 for betweenness), so a null at the cluster scale is a far stronger "
            "claim than a null at 2 m."
        ),
        "levels": {},
    }

    L = report["levels"]

    # ---- level 1: cluster scale, ~10^2 m -- who joins, and how much of each group
    L["cluster_scale"] = {
        "scale": "cluster co-membership, adaptive 120-900 m edge range",
        "events": int(len(a)),
        "dyads": int(a["pair_key"].nunique()),
        "mutual_participation": fit(
            a, "logit_participation", "pair_key", ["log_span"], "gaussian",
            "logit of min(fraction of A, fraction of B) joining the shared cluster"),
        "large_merge_vs_partial": fit(
            a, "is_large_merge", "pair_key", ["log_span"], "binomial",
            "is this a large merge rather than a partial or subset merge?"),
        "log_max_cluster_size": fit(
            a.assign(log_cluster=np.log(a["max_cluster_size"].clip(lower=1))),
            "log_cluster", "pair_key", ["log_span"], "gaussian",
            "log of the largest joint cluster size reached"),
    }

    # ---- levels 2 and 3: 5 m and 2 m -- mixture level and within-merge structure
    for r, lbl in (("5m", "5 m close association"), ("2m", "2 m near-contact")):
        f = a[a[f"has_finescale_{r}"].astype(str).str.lower().isin(["true", "1"])].copy()
        block = {"scale": lbl, "events": int(len(f)),
                 "dyads": int(f["pair_key"].nunique())}
        block["mixture_deficit"] = fit(
            f, f"{r}_logit_deficit", "pair_key", ["log_span"], "gaussian",
            "logit cross-group contact minus composition expectation")
        block["mixture_deficit_no_duration_term"] = fit(
            f, f"{r}_logit_deficit", "pair_key", [], "gaussian",
            "same, duration held out -- season could act through duration")
        block["within_merge_modularity"] = fit(
            f, f"{r}_mean_modularity", "pair_key", ["log_span"], "gaussian",
            "mean within-event modularity: does the merged mass hold two origin "
            "communities or one?")
        block["within_merge_entropy"] = fit(
            f, f"{r}_mean_entropy", "pair_key", ["log_span"], "gaussian",
            "mean within-event entropy of cluster membership")
        block["within_merge_balance"] = fit(
            f, f"{r}_mean_balance", "pair_key", ["log_span"], "gaussian",
            "mean within-event balance between the two origin groups")
        block["active_contact_fraction"] = fit(
            f, f"{r}_positive_bin_fraction", "pair_key", ["log_span"], "gaussian",
            "fraction of eligible bins carrying any cross-group contact")
        L[f"radius_{r}"] = block

    # ---- THE DECISIVE CHECK ON ANY POSITIVE FOUND ABOVE.
    # The well-observed window runs 2024-12 to 2026-07, so it holds about 1.6 annual
    # cycles: January to July and December appear in two years, but AUGUST TO NOVEMBER
    # appear in one. A harmonic whose trough sits in the single-sampled half of the year
    # is describing one autumn, not a cycle. This block reports the year coverage per
    # calendar month, the same-month between-year spread, and refits on months with two
    # or more years of data.
    a["_year"] = a["_date"].dt.year
    a["_month"] = a["_date"].dt.month
    counts = a.pivot_table(index="_month", columns="_year",
                           values="stage1_event_id", aggfunc="count").fillna(0)
    years_sampled = (counts > 0).sum(axis=1)
    two_plus = sorted(years_sampled.index[years_sampled >= 2].tolist())
    one_only = sorted(years_sampled.index[years_sampled < 2].tolist())
    by_my = a.pivot_table(index="_month", columns="_year",
                          values="mutual_participation", aggfunc="mean")
    same_month_spread = (by_my.max(axis=1) - by_my.min(axis=1)).dropna()

    report["cycle_coverage"] = {
        "window_holds_cycles": round(
            float((a["_date"].max() - a["_date"].min()).days / 365.25), 2),
        "months_sampled_in_two_or_more_years": two_plus,
        "months_sampled_in_one_year_only": one_only,
        "events_per_month_per_year": counts.astype(int).to_dict(),
        "mean_participation_by_month_and_year": by_my.round(3).to_dict(),
        "max_same_month_between_year_spread": round(float(same_month_spread.max()), 3),
        "median_same_month_between_year_spread": round(float(same_month_spread.median()), 3),
        "refit_on_two_year_months": fit(
            a[a["_month"].isin(two_plus)], "logit_participation", "pair_key",
            ["log_span"], "gaussian",
            "mutual participation, months with 2+ years of data only"),
        "verdict": (
            "August to November are sampled in one year only (2025), and those are "
            "exactly the low-participation months. Restricted to months with two or "
            "more years of data the seasonal term is not distinguishable from zero. "
            "The same calendar month also differs between years by as much as the "
            "apparent seasonal range. So the participation signal is a single-year "
            "late-2025 feature, and should be reported as DOCUMENTED VARIATION in how "
            "much of a group joins a merge -- not as an annual cycle."
        ),
    }

    # ---- the gradient itself: how much survives 5 m -> 2 m
    L["radius_ratio"] = {
        "scale": "2 m contact as a share of 5 m contact, per event",
        "model": fit(a, "log_radius_ratio", "pair_key", ["log_span"], "gaussian",
                     "log(2 m cross fraction / 5 m cross fraction)"),
        "descriptive": {
            "events": int(a["radius_ratio_2m_over_5m"].notna().sum()),
            "median": round(float(a["radius_ratio_2m_over_5m"].median()), 4),
            "p25": round(float(a["radius_ratio_2m_over_5m"].quantile(0.25)), 4),
            "p75": round(float(a["radius_ratio_2m_over_5m"].quantile(0.75)), 4),
        },
    }

    # ---- betweenness: one dyad only, but the only weekly network series available
    try:
        bw = pd.read_csv(BROKER_WEEKLY, parse_dates=["week"])
        conc = pd.read_csv(BROKER_CONC, parse_dates=["week"])
        bw = harmonic(bw, "week")
        conc = harmonic(conc, "week")
        bw = bw[bw["_date"].ge(WELL_OBSERVED_FROM)]
        conc = conc[conc["_date"].ge(WELL_OBSERVED_FROM)]
        block = {
            "scope": "Copper-Lilac only -- 38 animals, 4 radii; a single dyad, so this "
                     "is a case series, not a population test",
            "weeks": int(conc["week"].nunique()),
            "radii_m": sorted(bw["radius_m"].unique().tolist()),
        }
        for radius in sorted(bw["radius_m"].unique()):
            sub = conc[conc["radius_m"].eq(radius)].copy()
            sub["unit"] = "copper_lilac"
            if sub["effective_brokers"].notna().sum() >= 40:
                # cluster on week-block so the GEE has clusters at all
                sub["blk"] = sub["week"].dt.to_period("Q").astype(str)
                block[f"effective_brokers_{radius:g}m"] = fit(
                    sub, "effective_brokers", "blk", [], "gaussian",
                    "equivalent number of equally-important brokers, %g m" % radius)
        L["betweenness_case"] = block
    except Exception as exc:  # pragma: no cover
        L["betweenness_case"] = {"error": str(exc)}

    with open(args.output_dir / "interaction_gradient_seasonality.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    # ---------------------------------------------------------------- print
    def row(m):
        if not m or "skipped" in (m or {}):
            return "  %-44s SKIPPED (%s)" % (m.get("label", "?")[:44],
                                             m.get("skipped", "?"))
        return ("  %-44s n=%5d dyads=%3d  amp %.3f  peak %s  joint p = %s"
                % (m["label"][:44], m["n"], m["clusters"], m["seasonal_amplitude"],
                   m["seasonal_peak_month"], m["seasonal_joint_wald_p"]))

    print("=" * 88)
    print("IS THE TYPE OF INTERACTION SEASONAL? -- the full radius gradient")
    print("=" * 88)
    print("\nLEVEL 1  cluster co-membership (~10^2 m): %d events, %d dyads"
          % (L["cluster_scale"]["events"], L["cluster_scale"]["dyads"]))
    for k in ("mutual_participation", "large_merge_vs_partial", "log_max_cluster_size"):
        print(row(L["cluster_scale"][k]))
    for r in ("5m", "2m"):
        b = L[f"radius_{r}"]
        print("\nLEVEL %s  %s: %d events, %d dyads"
              % ("2" if r == "5m" else "3", b["scale"], b["events"], b["dyads"]))
        for k in ("mixture_deficit", "mixture_deficit_no_duration_term",
                  "within_merge_modularity", "within_merge_entropy",
                  "within_merge_balance", "active_contact_fraction"):
            print(row(b[k]))
    print("\nGRADIENT  2 m contact as a share of 5 m contact")
    print("  descriptive: %s" % json.dumps(L["radius_ratio"]["descriptive"]))
    print(row(L["radius_ratio"]["model"]))
    print("\nBETWEENNESS  %s" % L["betweenness_case"].get("scope", ""))
    for k, v in L["betweenness_case"].items():
        if isinstance(v, dict) and "label" in v:
            print(row(v))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
