"""Sleep-site membership from the EAS night-location table.

Supersedes the dusk/dawn proxy in analyze_lilac_range_and_sleep_sites.py. That
script had to bracket the night because the GPS collars record 03:00-16:00 UTC
only. EAS holds the proper product:

  processed/2025/gps/individual_night_locations.parquet
      date, animal_id, group_id, lat, lon, cluster, cluster_united
      -> one night location per animal-night, already assigned to a sleep-site
         cluster, so "same site" is the project's own definition rather than an
         arbitrary distance threshold.

  processed/2025/gps/dyad_dist_by_night.parquet
      night, id_i, id_j, distance  -> precomputed pairwise night distances.

Questions:
  1. Do the nine animals collared as Lilac on 2025-07-31 sleep at Lilac sites or
     Copper sites?
  2. 24AE04_6L7M has the largest Copper affinity of any Lilac-labelled animal.
     When did it change, and did the membership pipeline notice?

Read-only with respect to EAS.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
EAS_GPS = Path(r"Z:\baboon\working\data\processed\2025\gps")
NIGHT_LOC = EAS_GPS / "individual_night_locations.parquet"
DYAD_DIST = EAS_GPS / "dyad_dist_by_night.parquet"
MEMBERSHIP = Path(r"C:\Users\rharel\Documents\New project\outputs"
                  r"\canonical_robust_hourly_membership_shared_full_20260722"
                  r"\canonical_hourly_membership.parquet")
DEPLOYMENT = pd.Timestamp("2025-08-01")
NEW_LILAC = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
             "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"}
NEW_COPPER = {"25AB47_2ROS", "25AB49_1V3W", "25AB54_6F3G", "25AB55_OH7I"}


def cohort_of(animal: str, group: str) -> str:
    if animal in NEW_LILAC:
        return "Lilac_new"
    if animal in NEW_COPPER:
        return "Copper_new"
    return f"{group}_original"


def same_site_pairs(nl: pd.DataFrame, min_nights: int) -> pd.DataFrame:
    """Fraction of shared nights on which a pair occupied the same sleep site."""
    piv = nl.pivot_table(index="date", columns="animal_id", values="site",
                         aggfunc="first")
    animals = list(piv.columns)
    arr = piv.to_numpy(dtype=object)
    rows = []
    for i in range(len(animals)):
        ai = arr[:, i]
        for j in range(i + 1, len(animals)):
            aj = arr[:, j]
            ok = pd.notna(ai) & pd.notna(aj)
            n = int(ok.sum())
            if n < min_nights:
                continue
            rows.append({"animal_a": animals[i], "animal_b": animals[j],
                         "shared_nights": n,
                         "frac_same_site": float(np.mean(ai[ok] == aj[ok]))})
    return pd.DataFrame(rows)


def affinity(pairs: pd.DataFrame, cohorts: dict, value="frac_same_site") -> pd.DataFrame:
    a = pairs.rename(columns={"animal_a": "focal", "animal_b": "other"})
    b = pairs.rename(columns={"animal_b": "focal", "animal_a": "other"})
    long = pd.concat([a, b], ignore_index=True)
    long["other_cohort"] = long.other.map(cohorts)
    long["focal_cohort"] = long.focal.map(cohorts)
    g = (long[long.other_cohort.isin(["Lilac_original", "Copper_original"])]
         .groupby(["focal", "focal_cohort", "other_cohort"], as_index=False)[value].mean()
         .pivot_table(index=["focal", "focal_cohort"], columns="other_cohort", values=value)
         .reset_index())
    for c in ("Lilac_original", "Copper_original"):
        if c not in g:
            g[c] = np.nan
    g["lilac_minus_copper"] = g.Lilac_original - g.Copper_original
    g["closer_to"] = np.where(g.lilac_minus_copper > 0, "Lilac", "Copper")
    return g.sort_values(["focal_cohort", "lilac_minus_copper"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site-col", default="cluster_united",
                    choices=["cluster", "cluster_united"])
    ap.add_argument("--min-nights", type=int, default=20)
    ap.add_argument("--focal", default="24AE04_6L7M")
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "sleep_sites_eas_2026-09-03")
    a = ap.parse_args()
    if not NIGHT_LOC.exists():
        raise SystemExit(f"EAS not reachable: {NIGHT_LOC}")
    a.output_dir.mkdir(parents=True, exist_ok=True)

    nl = pd.read_parquet(NIGHT_LOC)
    for c in ("animal_id", "group_id", "cluster_united", "cluster"):
        nl[c] = nl[c].astype(str)
    nl["date"] = pd.to_datetime(nl["date"])
    print(f"night locations: {len(nl):,} animal-nights, {nl.animal_id.nunique()} animals, "
          f"{nl.date.min():%Y-%m-%d} to {nl.date.max():%Y-%m-%d}")
    print(f"sleep sites ({a.site_col}): {nl[a.site_col].nunique()} distinct")

    cl = nl[nl.group_id.isin(["Copper", "Lilac"])].copy()
    cl["site"] = cl[a.site_col]
    cohorts = {r.animal_id: cohort_of(r.animal_id, r.group_id)
               for r in cl.drop_duplicates("animal_id").itertuples(index=False)}
    print(f"Copper/Lilac: {len(cl):,} animal-nights, {cl.animal_id.nunique()} animals")
    print("  cohorts:", pd.Series(cohorts).value_counts().to_dict())

    post = cl[cl.date >= DEPLOYMENT]
    pairs = same_site_pairs(post, a.min_nights)
    pairs.to_csv(a.output_dir / "sleep_site_sharing_pairs.csv", index=False)

    p2 = pairs.copy()
    p2["cls"] = [" x ".join(sorted([cohorts.get(x, "?"), cohorts.get(y, "?")]))
                 for x, y in zip(p2.animal_a, p2.animal_b)]
    summ = (p2.groupby("cls", as_index=False)
              .agg(pairs=("frac_same_site", "size"), mean=("frac_same_site", "mean")))
    print(f"\n=== same sleep site, post-deployment (by {a.site_col}) ===")
    print(summ.sort_values("mean", ascending=False).round(4).to_string(index=False))

    aff = affinity(pairs, cohorts)
    aff.to_csv(a.output_dir / "sleep_site_affinity.csv", index=False)
    print(f"\n=== per-animal affinity ===")
    print(aff.round(4).to_string(index=False))
    for c in ("Lilac_new", "Copper_new", "Lilac_original", "Copper_original"):
        s = aff[aff.focal_cohort == c]
        if len(s):
            n = int((s.closer_to == "Lilac").sum())
            print(f"  {c:16s}: {n}/{len(s)} closer to LILAC sites, "
                  f"{len(s) - n}/{len(s)} closer to COPPER sites")

    # ---- focal case ----
    print(f"\n{'=' * 74}\nCASE: {a.focal}\n{'=' * 74}")
    f = cl[cl.animal_id == a.focal]
    print(f"night locations {f.date.min():%Y-%m-%d} to {f.date.max():%Y-%m-%d}, {len(f)} nights")
    lil = [x for x, c in cohorts.items() if c == "Lilac_original" and x != a.focal]
    cop = [x for x, c in cohorts.items() if c == "Copper_original" and x != a.focal]
    piv = cl.pivot_table(index="date", columns="animal_id", values="site", aggfunc="first")
    rows = []
    for night, r in piv.iterrows():
        if a.focal not in r or pd.isna(r.get(a.focal)):
            continue
        me = r[a.focal]
        rec = {"night": night}
        for lab, members in (("lilac", lil), ("copper", cop)):
            vals = [1.0 if r[o] == me else 0.0 for o in members
                    if o in r and pd.notna(r[o])]
            rec[f"with_{lab}"] = float(np.mean(vals)) if vals else np.nan
            rec[f"n_{lab}"] = len(vals)
        rows.append(rec)
    nightly = pd.DataFrame(rows).dropna(subset=["with_lilac", "with_copper"])
    nightly["month"] = nightly.night.values.astype("datetime64[M]")
    monthly = (nightly.groupby("month", as_index=False)
               .agg(nights=("night", "size"), with_lilac=("with_lilac", "mean"),
                    with_copper=("with_copper", "mean")))
    monthly["lilac_minus_copper"] = monthly.with_lilac - monthly.with_copper
    monthly.to_csv(a.output_dir / f"case_{a.focal}_monthly.csv", index=False)
    print(monthly.round(3).to_string(index=False))

    m = pd.read_parquet(MEMBERSHIP, columns=["window_start", "animal_id",
                                             "dynamic_social_unit", "dynamic_assignment"])
    fm = m[m.animal_id == a.focal]
    if len(fm):
        print(f"\npipeline said: dynamic_social_unit = "
              f"{sorted(fm.dynamic_social_unit.unique())}, "
              f"assignment = {sorted(fm.dynamic_assignment.unique())}")

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "source_night_locations": str(NIGHT_LOC),
        "source_mtime_utc": pd.Timestamp(NIGHT_LOC.stat().st_mtime, unit="s").isoformat(),
        "site_column": a.site_col, "min_shared_nights": a.min_nights,
        "window": "post-deployment (>= 2025-08-01) for the cohort comparison; "
                  "full record for the focal case",
        "note": "Same sleep site is the project's own cluster assignment, not a distance "
                "threshold, so it supersedes the dusk/dawn 100 m proxy.",
        "cohorts": cohorts,
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
