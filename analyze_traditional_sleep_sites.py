"""Traditional sleep sites: whose ground is each animal sleeping on, over time?

The earlier measure - "same site as Lilac animals" - loses meaning once the groups
fuse, because then everyone is at the same site by definition and the question
becomes circular. This anchors to the SITES instead.

METHOD
------
1. Find when Copper and Lilac started sharing sleep sites, from the data.
2. In the pre-fusion baseline, classify every sleep-site cluster by which group
   used it: Copper-traditional, Lilac-traditional, or shared.
3. For every month afterwards, ask what proportion of an animal's nights were
   spent on Copper-traditional ground versus Lilac-traditional ground.

Because site identity is fixed from the baseline, this stays interpretable after
the merge: an animal can be at a Lilac-traditional site whether or not any other
Lilac animal is there that night.

Source: EAS processed/2025/gps/individual_night_locations.parquet (cluster_united).
Read-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
NIGHT_LOC = Path(r"Z:\baboon\working\data\processed\2025\gps\individual_night_locations.parquet")
NEW_LILAC = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
             "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"}
NEW_COPPER = {"25AB47_2ROS", "25AB49_1V3W", "25AB54_6F3G", "25AB55_OH7I"}


def cohort_of(animal, group):
    if animal in NEW_LILAC:
        return "Lilac_new"
    if animal in NEW_COPPER:
        return "Copper_new"
    return f"{group}_original"


def classify_sites(base: pd.DataFrame, ratio: float, min_nights: int) -> pd.DataFrame:
    """Label each site by which group used it during the baseline."""
    use = (base.groupby(["site", "group_id"]).size().unstack(fill_value=0)
           .reindex(columns=["Copper", "Lilac"], fill_value=0))
    tot = use.sum(axis=0).replace(0, np.nan)
    share = use / tot                      # per-group share, so unequal collar counts cancel
    out = pd.DataFrame(index=use.index)
    out["nights_copper"] = use["Copper"]
    out["nights_lilac"] = use["Lilac"]
    out["share_copper"] = share["Copper"]
    out["share_lilac"] = share["Lilac"]
    eps = 1e-9
    r = (out.share_lilac + eps) / (out.share_copper + eps)
    out["site_type"] = np.where(out.nights_copper + out.nights_lilac < min_nights, "sparse",
                                np.where(r >= ratio, "Lilac_traditional",
                                         np.where(r <= 1 / ratio, "Copper_traditional",
                                                  "shared")))
    return out.reset_index()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site-col", default="cluster_united")
    ap.add_argument("--baseline-end", default="2025-05-01",
                    help="baseline = all nights before this date")
    ap.add_argument("--ratio", type=float, default=3.0,
                    help="a site is group-traditional when that group's usage share is "
                         "this many times the other's")
    ap.add_argument("--min-site-nights", type=int, default=10)
    ap.add_argument("--focal", default="24AE04_6L7M")
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "traditional_sleep_sites_2026-09-03")
    a = ap.parse_args()
    if not NIGHT_LOC.exists():
        raise SystemExit(f"EAS not reachable: {NIGHT_LOC}")
    a.output_dir.mkdir(parents=True, exist_ok=True)

    nl = pd.read_parquet(NIGHT_LOC)
    for c in ("animal_id", "group_id", a.site_col):
        nl[c] = nl[c].astype(str)
    nl["date"] = pd.to_datetime(nl["date"])
    cl = nl[nl.group_id.isin(["Copper", "Lilac"])].copy()
    cl["site"] = cl[a.site_col]
    cl["month"] = cl.date.values.astype("datetime64[M]")
    cl["cohort"] = [cohort_of(r.animal_id, r.group_id)
                    for r in cl.itertuples(index=False)]
    print(f"Copper/Lilac night locations: {len(cl):,} animal-nights, "
          f"{cl.animal_id.nunique()} animals, {cl.date.min():%Y-%m} to {cl.date.max():%Y-%m}")

    # ---- 1. when did the two groups start sharing sites? ----
    per_month = []
    for m, s in cl.groupby("month"):
        by_group = s.groupby(["date", "group_id"]).site.agg(
            lambda x: x.value_counts().idxmax())
        w = by_group.unstack()
        if {"Copper", "Lilac"} <= set(w.columns):
            ok = w[["Copper", "Lilac"]].dropna()
            if len(ok):
                per_month.append({"month": m, "nights": len(ok),
                                  "same_site": float((ok.Copper == ok.Lilac).mean())})
    fusion = pd.DataFrame(per_month)
    fusion.to_csv(a.output_dir / "group_site_sharing_monthly.csv", index=False)
    print("\n=== nights on which the two groups' MODAL site was the same ===")
    print(fusion.round(3).to_string(index=False))

    # ---- 2. classify sites from the pre-fusion baseline ----
    base = cl[cl.date < pd.Timestamp(a.baseline_end)]
    sites = classify_sites(base, a.ratio, a.min_site_nights)
    sites.to_csv(a.output_dir / "site_classification.csv", index=False)
    print(f"\n=== sites classified from baseline (< {a.baseline_end}) ===")
    print(sites.site_type.value_counts().to_string())
    top = sites[sites.site_type != "sparse"].nlargest(10, "nights_copper")
    print("\nlargest sites:")
    print(top[["site", "nights_copper", "nights_lilac", "site_type"]].to_string(index=False))

    # ---- 3. proportion of nights on each group's traditional ground ----
    stype = dict(zip(sites.site, sites.site_type))
    cl["site_type"] = cl.site.map(stype).fillna("novel")
    cl.to_csv(a.output_dir / "animal_nights_labelled.csv.gz", index=False,
              compression="gzip")

    def monthly_share(frame, by):
        g = (frame.groupby(by + ["month", "site_type"]).size()
             .unstack(fill_value=0))
        for c in ["Copper_traditional", "Lilac_traditional", "shared", "sparse", "novel"]:
            if c not in g:
                g[c] = 0
        tot = g.sum(axis=1).replace(0, np.nan)
        out = (g[["Copper_traditional", "Lilac_traditional", "shared"]]
               .div(tot, axis=0).reset_index())
        out["nights"] = tot.values
        return out

    coh = monthly_share(cl, ["cohort"])
    coh.to_csv(a.output_dir / "cohort_monthly_site_type.csv", index=False)
    foc = monthly_share(cl[cl.animal_id == a.focal], [])
    foc.to_csv(a.output_dir / f"case_{a.focal}_monthly_site_type.csv", index=False)
    print(f"\n=== {a.focal}: proportion of nights by site type ===")
    print(foc.round(3).to_string(index=False))

    print(f"\n=== cohort means, before vs from 2025-08 ===")
    coh["era"] = np.where(coh.month < pd.Timestamp("2025-08-01"), "pre", "post")
    print(coh.groupby(["cohort", "era"])[["Copper_traditional", "Lilac_traditional", "shared"]]
          .mean().round(3).to_string())

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "source": str(NIGHT_LOC), "site_column": a.site_col,
        "baseline_end": a.baseline_end, "traditional_ratio": a.ratio,
        "min_site_nights": a.min_site_nights,
        "site_type_counts": sites.site_type.value_counts().to_dict(),
        "method": "Site usage is converted to a within-group share before comparison, so a "
                  "group with more collars does not appear to own more sites.",
        "caveat": "Sites are classified once, from the baseline. A site that only came into "
                  "use later is 'novel' and is excluded from the three plotted proportions.",
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
