"""Is between-group variation in permeability real, or an artefact of collar coverage?

The claim this supports is not that the population behaves a certain way on
average, but that groups are NOT interchangeable units: permeability varies
widely across the population, and that variation is the result.

That claim only holds if the variation survives its obvious confound. Collar
coverage ranges 3%-42% across the 27 groups, and every permeability measure in
the project depends on detection, which depends on coverage. So for each measure
the question is the same: how much of the between-group spread is explained by
how well each group is watched, and how much remains?

Measures assembled here
  fusion rate per dyad      fusion hours / hours where the >=2-per-group quorum
                            was satisfiable  (all 248 dyads)
  split probability         per group, from the pipeline's own split-detection
                            evidence, conditional on classifiable hours
  scale sensitivity         per dyad, the share of hours on which DBSCAN(500 m)
                            and HDBSCAN disagree about whether the two groups
                            were together - how ambiguous that boundary is

Read-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT = Path(__file__).resolve().parent
OUT = PROJECT / "outputs"
DEMOG = OUT / "authoritative_group_names_2026-09-03" / "group_demographics.csv"
FUSION = OUT / "fusion_quorum_all_dyads_2026-09-03" / "fusion_rate_by_dyad.csv"
SPLIT = OUT / "fusion_quorum_all_dyads_2026-09-03" / "split_probability_rank.csv"
SPLIT_ALT = Path(r"C:\Users\rharel\Documents\New project\outputs"
                 r"\canonical_robust_hourly_membership_shared_full_20260722"
                 r"\split_detection_evidence\split_evidence_summary_by_dynamic_group.csv")
AGREE = OUT / "clustering_method_agreement_2026-09-04" / "group_pair_calls_by_hour.csv"


def norm(s):
    return s.astype(str).str.replace(" ", "", regex=False).str.lower()


def spread(v, label, unit=""):
    v = pd.Series(v).dropna()
    if len(v) < 3:
        return None
    lo, hi = v.min(), v.max()
    q1, q3 = v.quantile(.25), v.quantile(.75)
    cv = v.std() / v.mean() if v.mean() else np.nan
    fold = (hi / lo) if lo > 0 else np.inf
    print(f"  {label:<26} n={len(v):>3}  median {v.median():.3f}{unit}  "
          f"IQR {q1:.3f}-{q3:.3f}  range {lo:.3f}-{hi:.3f}  CV {cv:.2f}  "
          f"fold {'inf' if np.isinf(fold) else f'{fold:.0f}x'}")
    return {"n": int(len(v)), "median": float(v.median()), "q1": float(q1), "q3": float(q3),
            "min": float(lo), "max": float(hi), "cv": float(cv)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-opportunity", type=int, default=2000,
                    help="minimum opportunity hours for a dyad to be included")
    ap.add_argument("--output-dir", type=Path,
                    default=OUT / "permeability_variation_2026-09-04")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    report = {}

    demo = pd.read_csv(DEMOG)
    demo["key"] = norm(demo.group_id)
    cov = demo.set_index("key")
    print(f"demographics: {len(demo)} groups, coverage "
          f"{demo.collar_coverage_percent.min():.1f}%-{demo.collar_coverage_percent.max():.1f}%")

    # ---------- 1. fusion rate per dyad ----------
    print(f"\n{'=' * 92}\nHOW MUCH DO DYADS DIFFER?\n{'=' * 92}")
    f = pd.read_csv(FUSION)
    f = f[f.opportunity_hours >= a.min_opportunity].copy()
    f[["ga", "gb"]] = f.dyad.str.split(" - ", expand=True)
    for s, c in (("ga", "a"), ("gb", "b")):
        f[f"cov_{c}"] = norm(f[s]).map(cov.collar_coverage_percent)
        f[f"size_{c}"] = norm(f[s]).map(cov.group_size)
    f["cov_min"] = f[["cov_a", "cov_b"]].min(axis=1)
    f["cov_mean"] = f[["cov_a", "cov_b"]].mean(axis=1)
    f["size_total"] = f.size_a + f.size_b
    report["fusion_rate"] = spread(f.fusion_rate, "fusion rate per dyad")
    nz = f[f.fusion_rate > 0]
    print(f"  {'':<26} {len(nz)} of {len(f)} dyads show any fusion at all")

    # ---------- 2. split probability per group ----------
    sp = pd.read_csv(SPLIT) if SPLIT.exists() else None
    if sp is None or "split_prob" not in sp.columns:
        s2 = pd.read_csv(SPLIT_ALT)
        sp = s2.rename(columns={"dynamic_social_unit": "dynamic_social_unit",
                                "split_probability_given_classifiable": "split_prob"})
        sp = sp[sp.opportunity_hours >= 500]
    sp["key"] = norm(sp.dynamic_social_unit)
    sp["coverage"] = sp.key.map(cov.collar_coverage_percent)
    sp["group_size"] = sp.key.map(cov.group_size)
    report["split_prob"] = spread(sp.split_prob, "split probability per group")

    # ---------- 3. scale sensitivity per dyad ----------
    sens = None
    if AGREE.exists():
        g = pd.read_csv(AGREE)
        g["disagree"] = g.fixed_500 != g.hdbscan
        sens = (g.groupby(["group_a", "group_b"], as_index=False)
                  .agg(hours=("disagree", "size"), disagree=("disagree", "mean"),
                       together_500=("fixed_500", "mean")))
        sens = sens[sens.hours >= 500]
        report["scale_sensitivity"] = spread(sens.disagree, "scale sensitivity per dyad")

    # ---------- the confound ----------
    print(f"\n{'=' * 92}\nIS IT JUST COVERAGE?\n{'=' * 92}")
    tests = []
    def rho(x, y, label):
        d = pd.DataFrame({"x": x, "y": y}).dropna()
        if len(d) < 6:
            return
        r, p = spearmanr(d.x, d.y)
        tests.append({"comparison": label, "n": int(len(d)), "rho": float(r), "p": float(p)})
        flag = " *" if p < 0.05 else ""
        print(f"  {label:<52} n={len(d):>3}  rho {r:+.3f}  p={p:.3f}{flag}")

    rho(f.cov_min, f.fusion_rate, "dyad fusion rate  vs  lower coverage of the pair")
    rho(f.cov_mean, f.fusion_rate, "dyad fusion rate  vs  mean coverage of the pair")
    rho(f.size_total, f.fusion_rate, "dyad fusion rate  vs  combined demographic size")
    rho(sp.coverage, sp.split_prob, "group split probability  vs  its collar coverage")
    rho(sp.group_size, sp.split_prob, "group split probability  vs  its demographic size")
    if sens is not None:
        sens["cov_min"] = np.minimum(norm(sens.group_a).map(cov.collar_coverage_percent),
                                     norm(sens.group_b).map(cov.collar_coverage_percent))
        rho(sens.cov_min, sens.disagree, "scale sensitivity  vs  lower coverage of the pair")
        rho(sens.together_500, sens.disagree, "scale sensitivity  vs  how often they are together")
    report["confound_tests"] = tests

    # ---------- what the spread looks like once ranked ----------
    print(f"\n{'=' * 92}\nTHE VARIATION ITSELF\n{'=' * 92}")
    print("\n  fusion rate, dyads with the most opportunity:")
    t = f.nlargest(10, "fusion_rate")[["dyad", "fusion_rate", "opportunity_hours", "cov_min"]]
    for r in t.itertuples(index=False):
        bar = "#" * int(round(r.fusion_rate * 44))
        print(f"    {r.dyad:<26} {r.fusion_rate:>6.3f}  cov>={r.cov_min:>4.1f}%  {bar}")

    print("\n  split probability by group:")
    for r in sp.nlargest(12, "split_prob").itertuples(index=False):
        bar = "#" * int(round(r.split_prob * 110))
        cv = f"{r.coverage:>4.1f}%" if np.isfinite(r.coverage) else "   ?"
        print(f"    {r.dynamic_social_unit:<16} {r.split_prob:>6.3f}  cov {cv}  "
              f"size {r.group_size if np.isfinite(r.group_size) else '?':>5}  {bar}")

    if sens is not None:
        print("\n  scale sensitivity (DBSCAN-500 vs HDBSCAN disagreement):")
        for r in sens.nlargest(10, "disagree").itertuples(index=False):
            bar = "#" * int(round(r.disagree * 300))
            print(f"    {r.group_a + ' - ' + r.group_b:<30} {r.disagree:>6.3f}  {bar}")

    f.to_csv(a.output_dir / "dyad_fusion_with_coverage.csv", index=False)
    sp.to_csv(a.output_dir / "group_split_with_coverage.csv", index=False)
    if sens is not None:
        sens.to_csv(a.output_dir / "dyad_scale_sensitivity.csv", index=False)
    (a.output_dir / "permeability_variation.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
