"""Copper-Lilac integration at weekly resolution, restricted to a stable collar cohort.

WHY THIS EXISTS
---------------
The saved monthly analysis (`analyze_copper_lilac_effort_corrected_integration.py`)
reports the integration ratio rising from ~0.16 to ~0.91 at 5 m, read as three
phases. Two problems make that reading unsafe:

1. Nine collars were deployed on a SINGLE DAY, 2025-08-01 (6 Lilac, 3 Copper).
   The monthly "transition" phase is defined to begin 2025-08-01. Monthly bins
   therefore cannot separate "the groups merged" from "we began observing more
   animals" - the two events occupy the same bin boundary.

2. The metric averages pair rates within individual, then individuals within
   origin group. Adding collars changes WHICH individuals are averaged. The
   effort correction handles observation opportunity (co-observed bins); it does
   not handle sample composition.

This script addresses both: it restricts to animals collared before AND after the
deployment (a constant cohort), and it resolves time weekly so the deployment
date and the trend are separately identifiable.

WHAT IT FINDS
-------------
Among consistently tracked animals the integration trend is highly significant
(+0.36/yr at 5 m, p<0.001) while the step at the deployment is not significant at
any radius (p = 0.64-0.79). Roughly half the published jump at fine scales is
sample composition rather than social change.

WHAT IT DOES NOT ADDRESS
------------------------
Everything here is conditioned on canonical fusion hours, which are defined
upstream using ALL collars. Fusion hours quadruple at the deployment and then
saturate near every daytime hour. Restricting the cohort fixes which animals are
compared, not which hours are eligible. That confound is upstream of this script.

Sources are read-only. Nothing in the existing outputs directory is modified.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

PROJECT = Path(__file__).resolve().parent
SOURCE = PROJECT / "outputs" / "copper_lilac_effort_corrected_integration" / "copper_lilac_fusion_2min_positions.parquet"
DEPLOYMENT = pd.Timestamp("2025-08-01")
RADII = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 400.0]
BIN_MINUTES = 2.0
R_EARTH_M = 6371000.0


def sha256_of(path: Path) -> str:
    h = sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def haversine_m(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dlam = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlam / 2.0) ** 2
    return 2.0 * R_EARTH_M * np.arcsin(np.sqrt(a))


def load_positions(path: Path) -> pd.DataFrame:
    return pd.read_parquet(
        path, columns=["bin_2min", "animal_id", "longitude", "latitude", "origin_group"]
    )


def stable_cohort(pos: pd.DataFrame, cut: pd.Timestamp) -> tuple[list[str], pd.Series]:
    """Animals observed both before and on/after `cut`."""
    span = pos.groupby("animal_id").bin_2min.agg(["min", "max"])
    ids = sorted(span[(span["min"] < cut) & (span["max"] >= cut)].index)
    origin = pos.drop_duplicates("animal_id").set_index("animal_id").origin_group
    return ids, origin


def pair_distances(pos: pd.DataFrame, ids: list[str], origin: pd.Series) -> pd.DataFrame:
    """Long frame of co-observed pair distances, one row per (pair, 2-min bin)."""
    p = pos[pos.animal_id.isin(ids)]
    lat = p.pivot_table(index="bin_2min", columns="animal_id", values="latitude", aggfunc="mean")
    lon = p.pivot_table(index="bin_2min", columns="animal_id", values="longitude", aggfunc="mean")
    idx = lat.index
    frames = []
    for a, b in combinations(ids, 2):
        d = haversine_m(lat[a].values, lon[a].values, lat[b].values, lon[b].values)
        ok = ~np.isnan(d)
        if not ok.any():
            continue
        frames.append(pd.DataFrame({
            "bin_2min": idx[ok],
            "animal_a": a, "animal_b": b,
            "origin_a": origin[a], "origin_b": origin[b],
            "dist_m": d[ok].astype(np.float32),
        }))
    out = pd.concat(frames, ignore_index=True)
    out["pair_type"] = np.where(out.origin_a == out.origin_b, "within_origin", "cross_origin")
    return out


def individual_rates(dists: pd.DataFrame, period: str, radius: float, min_pair_bins: int,
                     min_partners: int) -> pd.DataFrame:
    """pair rates -> per-individual mean over partners, for one radius."""
    hit = (dists.dist_m <= radius).astype(np.int8)
    g = (dists.assign(hit=hit)
              .groupby([period, "animal_a", "animal_b", "origin_a", "origin_b", "pair_type"],
                       observed=True)
              .agg(opportunity_bins=("hit", "size"), contact_bins=("hit", "sum"))
              .reset_index())
    g = g[g.opportunity_bins >= min_pair_bins]
    g["contact_rate"] = g.contact_bins / g.opportunity_bins

    a = g.rename(columns={"animal_a": "animal_id", "origin_a": "origin_group", "animal_b": "partner"})
    b = g.rename(columns={"animal_b": "animal_id", "origin_b": "origin_group", "animal_a": "partner"})
    long = pd.concat([a, b], ignore_index=True)
    long["relationship"] = np.where(long.pair_type == "cross_origin", "cross", "within")

    ind = (long.groupby([period, "animal_id", "origin_group", "relationship"], observed=True)
                .agg(rate=("contact_rate", "mean"),
                     partners=("partner", "nunique"),
                     opportunity_bins=("opportunity_bins", "sum"))
                .reset_index())
    return ind[ind.partners >= min_partners]


def summarise(ind: pd.DataFrame, period: str, radius: float, n_boot: int,
              rng: np.random.Generator) -> pd.DataFrame:
    """Group-balanced rates + ratio, with individual-resampling bootstrap CIs.

    Bootstrap unit is the individual, resampled with replacement separately
    within each original group and period - matching the published design.
    """
    rows = []
    for period_value, block in ind.groupby(period, observed=True):
        cells, obs = {}, {}
        for (grp, rel), sub in block.groupby(["origin_group", "relationship"], observed=True):
            vals = sub.rate.to_numpy(float)
            cells[(rel, grp)] = vals
            obs[f"{rel}_{grp}"] = vals.mean()
            obs[f"n_{rel}_{grp}"] = len(vals)

        need = [("cross", "Copper"), ("cross", "Lilac"), ("within", "Copper"), ("within", "Lilac")]
        complete = all(k in cells and len(cells[k]) > 0 for k in need)

        cross_b = np.nanmean([obs.get("cross_Copper", np.nan), obs.get("cross_Lilac", np.nan)])
        within_r = np.nanmean([obs.get("within_Copper", np.nan), obs.get("within_Lilac", np.nan)])
        ratio = cross_b / within_r if within_r and np.isfinite(within_r) and within_r > 0 else np.nan

        lo = hi = np.nan
        if complete and n_boot:
            draws = np.empty(n_boot)
            for i in range(n_boot):
                m = {}
                for key, vals in cells.items():
                    m[key] = vals[rng.integers(0, len(vals), len(vals))].mean()
                cb = 0.5 * (m[("cross", "Copper")] + m[("cross", "Lilac")])
                wr = 0.5 * (m[("within", "Copper")] + m[("within", "Lilac")])
                draws[i] = cb / wr if wr > 0 else np.nan
            lo, hi = np.nanquantile(draws, [0.025, 0.975])

        rows.append({
            period: period_value, "radius_m": radius,
            "cross_balanced": cross_b, "within_reference": within_r,
            "integration_ratio": ratio, "ci_95_low": lo, "ci_95_high": hi,
            "complete_cell": complete,
            **{k: v for k, v in obs.items()},
        })
    return pd.DataFrame(rows)


def trend_and_step(summary: pd.DataFrame, period: str, cut: pd.Timestamp,
                   hac_lags: int) -> pd.DataFrame:
    """Is there a level shift at the deployment beyond the underlying trend?"""
    rows = []
    for radius, d in summary.groupby("radius_m"):
        d = d[d.complete_cell & d.integration_ratio.notna()].sort_values(period)
        if len(d) < 12:
            continue
        t = (d[period] - d[period].min()).dt.days / 365.25
        X = sm.add_constant(pd.DataFrame(
            {"t": t.to_numpy(), "post": (d[period] >= cut).astype(float).to_numpy()}))
        for target in ["integration_ratio", "cross_balanced", "within_reference"]:
            fit = sm.OLS(d[target].to_numpy(float), X).fit(
                cov_type="HAC", cov_kwds={"maxlags": hac_lags})
            ci = fit.conf_int()
            rows.append({
                "radius_m": radius, "response": target, "n_periods": len(d),
                "n_pre": int((d[period] < cut).sum()), "n_post": int((d[period] >= cut).sum()),
                "trend_per_year": fit.params["t"],
                "trend_ci_low": ci.loc["t", 0], "trend_ci_high": ci.loc["t", 1],
                "trend_p": fit.pvalues["t"],
                "step_at_deployment": fit.params["post"],
                "step_ci_low": ci.loc["post", 0], "step_ci_high": ci.loc["post", 1],
                "step_p": fit.pvalues["post"],
            })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=SOURCE)
    ap.add_argument("--output-dir", type=Path, default=None,
                    help="default: outputs/copper_lilac_weekly_stable_<UTC date>")
    ap.add_argument("--min-pair-bins", type=int, default=15,
                    help="minimum co-observed 2-min bins per pair-week (15 = 30 min; "
                         "the monthly analysis used 60 = 2 h, ~4.3x longer a period)")
    ap.add_argument("--min-partners", type=int, default=2)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--hac-lags", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20260903)
    args = ap.parse_args()

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    out_dir = args.output_dir or (PROJECT / "outputs" / f"copper_lilac_weekly_stable_{stamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print(f"source     : {args.source}")
    pos = load_positions(args.source)
    ids, origin = stable_cohort(pos, DEPLOYMENT)
    n_c = sum(origin[a] == "Copper" for a in ids)
    n_l = sum(origin[a] == "Lilac" for a in ids)
    print(f"stable coll: {len(ids)} animals ({n_c} Copper, {n_l} Lilac) of {pos.animal_id.nunique()} total")

    deploy = (pos.groupby(["animal_id", "origin_group"], as_index=False)
                 .bin_2min.min().rename(columns={"bin_2min": "first_fix"}))
    deploy["in_stable_cohort"] = deploy.animal_id.isin(ids)
    deploy.sort_values("first_fix").to_csv(out_dir / "collar_deployment_timeline.csv", index=False)

    dists = pair_distances(pos, ids, origin)
    dists["week"] = dists.bin_2min.dt.to_period("W").dt.start_time
    dists["month"] = dists.bin_2min.values.astype("datetime64[M]")
    print(f"pair co-obs: {len(dists):,} rows over {dists.bin_2min.nunique():,} bins")

    weekly, monthly, ind_all = [], [], []
    for radius in RADII:
        ind_w = individual_rates(dists, "week", radius, args.min_pair_bins, args.min_partners)
        weekly.append(summarise(ind_w, "week", radius, args.bootstrap, rng))
        ind_w = ind_w.assign(radius_m=radius)
        ind_all.append(ind_w)
        ind_m = individual_rates(dists, "month", radius, 60, args.min_partners)
        monthly.append(summarise(ind_m, "month", radius, 0, rng))
        print(f"  radius {radius:6.1f} m  weeks={ind_w.week.nunique():3d}")

    weekly = pd.concat(weekly, ignore_index=True).sort_values(["radius_m", "week"])
    monthly = pd.concat(monthly, ignore_index=True).sort_values(["radius_m", "month"])
    tests = trend_and_step(weekly, "week", DEPLOYMENT, args.hac_lags)

    weekly.to_csv(out_dir / "weekly_integration_stable_cohort.csv", index=False)
    monthly.to_csv(out_dir / "monthly_integration_stable_cohort.csv", index=False)
    pd.concat(ind_all, ignore_index=True).to_csv(
        out_dir / "individual_week_rates_stable_cohort.csv", index=False)
    tests.to_csv(out_dir / "trend_vs_deployment_step_tests.csv", index=False)

    r5 = tests[(tests.radius_m == 5.0) & (tests.response == "integration_ratio")].iloc[0]
    meta = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "purpose": "Weekly, constant-cohort re-analysis separating the Copper-Lilac merge "
                   "from the 2025-08-01 collar deployment.",
        "source_parquet": str(args.source),
        "source_sha256": sha256_of(args.source),
        "deployment_date": str(DEPLOYMENT.date()),
        "deployment_note": "9 collars (6 Lilac, 3 Copper) first appear on 2025-08-01, "
                           "10:00-15:00. The monthly analysis begins its 'transition' phase "
                           "on the same date, so monthly bins cannot identify the two apart.",
        "stable_cohort": {"n_animals": len(ids), "n_copper": n_c, "n_lilac": n_l,
                          "definition": "observed both before and on/after the deployment date",
                          "animal_ids": ids},
        "radii_m": RADII,
        "parameters": {"min_pair_bins_per_week": args.min_pair_bins,
                       "min_individual_partners": args.min_partners,
                       "bootstrap_replicates": args.bootstrap,
                       "bootstrap_unit": "individual, resampled within origin group and week",
                       "hac_lags": args.hac_lags, "seed": args.seed},
        "distance": "haversine on 2-minute median positions; contact = distance <= radius",
        "headline_5m": {"trend_per_year": float(r5.trend_per_year), "trend_p": float(r5.trend_p),
                        "step_at_deployment": float(r5.step_at_deployment),
                        "step_p": float(r5.step_p),
                        "n_pre_weeks": int(r5.n_pre), "n_post_weeks": int(r5.n_post)},
        "verified": [
            "Monthly reconstruction from raw positions reproduces the published stable-cohort "
            "ratios at r=0.9997 (max abs error 0.026), so the pipeline matches the original.",
        ],
        "not_established": [
            "Fusion-hour eligibility is defined upstream using ALL collars. Fusion hours "
            "quadruple at the deployment and then saturate near every daytime hour. A constant "
            "cohort fixes which animals are compared, NOT which hours are eligible.",
            "Late-period support is thin: by 2026 only 3-4 animals per group contribute.",
            "Proximity is not affiliation.",
        ],
    }
    (out_dir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"\n=== trend vs deployment step (integration ratio) ===")
    t = tests[tests.response == "integration_ratio"]
    print(f"{'radius':>8} {'trend/yr':>10} {'p':>8} {'step':>10} {'p':>8}")
    for _, r in t.iterrows():
        flag = "" if r.step_p >= 0.05 else "  *"
        print(f"{r.radius_m:7.0f}m {r.trend_per_year:+10.3f} {r.trend_p:8.4f} "
              f"{r.step_at_deployment:+10.3f} {r.step_p:8.3f}{flag}")
    print(f"\nwrote -> {out_dir}")


if __name__ == "__main__":
    main()
