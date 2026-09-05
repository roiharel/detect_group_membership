"""Do the nine animals collared as Lilac on 2025-07-31 range and sleep with Lilac,
with Copper, or with neither?

WHY THIS IS INDEPENDENT EVIDENCE
--------------------------------
Everything established so far about that cohort comes from fine-scale contact
inside canonical Copper-Lilac FUSION hours - an hour set that is itself
collar-dependent and saturates after the deployment. Home range and sleeping site
are measured from the full GPS record: all hours, no fusion conditioning, no
quorum. If they agree with the contact result, the conclusion no longer rests on
the fusion pipeline at all.

Sleeping site is also the strongest natural group signature available. Baboon
groups return to a limited set of sleeping sites, and a group's members sleep
together nightly whether or not they mixed at close range that day. Two animals
that sleep in the same place are in the same group in a way that daytime
proximity cannot establish on its own.

MEASURES
--------
  sleep-site sharing   For each pair, the fraction of nights BOTH were located in
                       which their night positions fall within `--sleep-radius`.
                       Effort-corrected by construction: the denominator is
                       shared nights, not calendar nights.

  range overlap        Bhattacharyya coefficient between the two animals' gridded
                       occupancy distributions. 1 = identical use of space,
                       0 = disjoint. Insensitive to differences in fix counts.

Each new-cohort animal is then scored against the original Lilac animals and
against the original Copper animals, and the two are compared.

Read-only.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

PROJECT = Path(__file__).resolve().parent
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
DEPLOYMENT = pd.Timestamp("2025-08-01")
NEW_LILAC = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
             "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"}
NEW_COPPER = {"25AB47_2ROS", "25AB49_1V3W", "25AB54_6F3G", "25AB55_OH7I"}
R_EARTH_M = 6371000.0
# These collars record 03:00-16:00 UTC only (06:00-19:00 local) - there are NO
# overnight fixes, which is why the membership pipeline carries a night position
# rather than measuring one. A true sleeping site is therefore not observable.
# The best available proxy brackets the night: the last fixes before dark and the
# first after dawn. Baboons are at or beside the roost at both moments, so two
# animals that agree at BOTH ends spent the night in the same place.
DUSK_FROM_UTC = 15      # 18:00 local, last hour of recording
DAWN_TO_UTC = 4         # 07:00 local, first full hour of recording


def haversine_m(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2.0) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2.0) ** 2)
    return 2.0 * R_EARTH_M * np.arcsin(np.sqrt(a))


def load_gps(path: Path, groups: list[str]) -> pd.DataFrame:
    f = pq.ParquetFile(path)
    cols = ["animal_id", "timestamp", "location.long", "location.lat", "group_id"]
    keep = []
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=cols)
        t = t.filter(pc.is_in(t["group_id"], value_set=pd.array(groups).__arrow_array__()))
        if t.num_rows:
            keep.append(t.to_pandas())
    d = pd.concat(keep, ignore_index=True)
    d = d.rename(columns={"location.long": "lon", "location.lat": "lat",
                          "group_id": "origin_group"})
    d["timestamp"] = pd.to_datetime(d["timestamp"], utc=True).dt.tz_localize(None)
    # parquet hands these back as categoricals; grouping on a categorical builds the
    # full cartesian product of levels, which blows up the aggregation
    for c in ("animal_id", "origin_group"):
        d[c] = d[c].astype(str)
    return d.dropna(subset=["lat", "lon"])


def sleep_sites(d: pd.DataFrame) -> pd.DataFrame:
    """Roost proxy per animal-night: dusk position of day D and dawn of day D+1.

    Returned separately so they can be required to agree, rather than averaged
    into a single point that would hide disagreement.
    """
    h = d.timestamp.dt.hour
    dusk = d[h >= DUSK_FROM_UTC].copy()
    dusk["night"] = dusk.timestamp.dt.floor("D")
    dawn = d[h <= DAWN_TO_UTC].copy()
    # a dawn fix belongs to the night that just ended, i.e. the previous day
    dawn["night"] = dawn.timestamp.dt.floor("D") - pd.Timedelta(days=1)

    out = []
    for lab, s in (("dusk", dusk), ("dawn", dawn)):
        g = (s.groupby(["animal_id", "origin_group", "night"], as_index=False, observed=True)
              .agg(lat=("lat", "median"), lon=("lon", "median"), fixes=("lat", "size")))
        out.append(g[g.fixes >= 2].assign(edge=lab))
    return pd.concat(out, ignore_index=True)


def sleep_sharing(sites: pd.DataFrame, radius_m: float, min_nights: int) -> pd.DataFrame:
    """Pairwise roost sharing, requiring agreement at BOTH ends of the night."""
    piv = {}
    for edge in ("dusk", "dawn"):
        s = sites[sites.edge == edge]
        piv[edge] = (s.pivot_table(index="night", columns="animal_id", values="lat"),
                     s.pivot_table(index="night", columns="animal_id", values="lon"))
    animals = sorted(set(piv["dusk"][0].columns) & set(piv["dawn"][0].columns))
    rows = []
    for a, b in combinations(animals, 2):
        per_edge, counts = {}, {}
        for edge in ("dusk", "dawn"):
            la, lo = piv[edge]
            if a not in la.columns or b not in la.columns:
                continue
            xa, ya = la[a].values, lo[a].values
            xb, yb = la[b].values, lo[b].values
            ok = ~(np.isnan(xa) | np.isnan(xb))
            if ok.sum() == 0:
                continue
            dist = haversine_m(xa[ok], ya[ok], xb[ok], yb[ok])
            per_edge[edge] = (la.index[ok], dist)
            counts[edge] = int(ok.sum())
        if len(per_edge) < 2 or min(counts.values()) < min_nights:
            continue
        # a shared roost requires being together at dusk AND again at dawn
        du = pd.Series(per_edge["dusk"][1] <= radius_m, index=per_edge["dusk"][0])
        da = pd.Series(per_edge["dawn"][1] <= radius_m, index=per_edge["dawn"][0])
        both = du.align(da, join="inner")
        n = len(both[0])
        if n < min_nights:
            continue
        rows.append({
            "animal_a": a, "animal_b": b, "shared_nights": n,
            "frac_same_site": float((both[0] & both[1]).mean()),
            "frac_dusk_only": float(both[0].mean()),
            "median_night_distance_m": float(np.median(per_edge["dusk"][1])),
        })
    return pd.DataFrame(rows)


def range_overlap(d: pd.DataFrame, cell_m: float, min_fixes: int) -> pd.DataFrame:
    lat0 = d.lat.mean()
    mx = 111320.0 * np.cos(np.radians(lat0))
    d = d.assign(gx=np.floor(d.lon * mx / cell_m).astype(int),
                 gy=np.floor(d.lat * 110540.0 / cell_m).astype(int))
    occ = {}
    for aid, s in d.groupby("animal_id"):
        if len(s) < min_fixes:
            continue
        c = s.groupby(["gx", "gy"]).size()
        occ[aid] = c / c.sum()
    rows = []
    for a, b in combinations(sorted(occ), 2):
        pa, pb = occ[a].align(occ[b], join="inner", fill_value=0.0)
        bc = float(np.sum(np.sqrt(pa.values * pb.values)))
        rows.append({"animal_a": a, "animal_b": b, "range_overlap": bc})
    return pd.DataFrame(rows)


def affinity(pairs: pd.DataFrame, value: str, cohorts: dict) -> pd.DataFrame:
    """For each animal, mean `value` against original Lilac vs original Copper."""
    a = pairs.rename(columns={"animal_a": "focal", "animal_b": "other"})
    b = pairs.rename(columns={"animal_b": "focal", "animal_a": "other"})
    long = pd.concat([a, b], ignore_index=True)
    long["other_cohort"] = long.other.map(cohorts)
    long["focal_cohort"] = long.focal.map(cohorts)
    g = (long[long.other_cohort.isin(["Lilac_original", "Copper_original"])]
         .groupby(["focal", "focal_cohort", "other_cohort"], as_index=False)[value].mean()
         .pivot_table(index=["focal", "focal_cohort"], columns="other_cohort",
                      values=value).reset_index())
    for c in ("Lilac_original", "Copper_original"):
        if c not in g:
            g[c] = np.nan
    g["lilac_minus_copper"] = g.Lilac_original - g.Copper_original
    g["closer_to"] = np.where(g.lilac_minus_copper > 0, "Lilac", "Copper")
    return g.sort_values(["focal_cohort", "lilac_minus_copper"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gps", type=Path, default=GPS)
    ap.add_argument("--sleep-radius", type=float, default=100.0)
    ap.add_argument("--min-shared-nights", type=int, default=15)
    ap.add_argument("--cell-m", type=float, default=250.0)
    ap.add_argument("--min-fixes", type=int, default=2000)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "lilac_range_sleep_sites_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    print("loading GPS for Copper and Lilac (full record, all hours)...")
    d = load_gps(a.gps, ["Copper", "Lilac"])
    print(f"  {len(d):,} fixes, {d.animal_id.nunique()} animals, "
          f"{d.timestamp.min():%Y-%m-%d} to {d.timestamp.max():%Y-%m-%d}")

    cohorts = {}
    for aid, grp in d.drop_duplicates("animal_id")[["animal_id", "origin_group"]].values:
        if aid in NEW_LILAC:
            cohorts[aid] = "Lilac_new"
        elif aid in NEW_COPPER:
            cohorts[aid] = "Copper_new"
        else:
            cohorts[aid] = f"{grp}_original"
    print("  cohorts:", pd.Series(cohorts).value_counts().to_dict())

    post = d[d.timestamp >= DEPLOYMENT]
    sites = sleep_sites(post)
    print(f"\nsleep sites: {len(sites):,} animal-nights, {sites.animal_id.nunique()} animals "
          f"(post-deployment only)")

    share = sleep_sharing(sites, a.sleep_radius, a.min_shared_nights)
    share.to_csv(a.output_dir / "sleep_site_sharing_pairs.csv", index=False)
    rng = range_overlap(post, a.cell_m, a.min_fixes)
    rng.to_csv(a.output_dir / "range_overlap_pairs.csv", index=False)

    for lab, pairs, val in [("SLEEPING SITE (fraction of shared nights within "
                             f"{a.sleep_radius:g} m)", share, "frac_same_site"),
                            (f"HOME RANGE OVERLAP (Bhattacharyya, {a.cell_m:g} m cells)",
                             rng, "range_overlap")]:
        print(f"\n{'=' * 76}\n{lab}\n{'=' * 76}")
        a2 = pairs.rename(columns={"animal_a": "x", "animal_b": "y"})
        a2["ca"] = a2.x.map(cohorts); a2["cb"] = a2.y.map(cohorts)
        a2["cls"] = [" x ".join(sorted([p, q])) for p, q in zip(a2.ca, a2.cb)]
        summ = a2.groupby("cls", as_index=False).agg(pairs=(val, "size"), mean=(val, "mean"))
        print(summ.sort_values("mean", ascending=False).round(4).to_string(index=False))

        aff = affinity(pairs, val, cohorts)
        aff.to_csv(a.output_dir / f"affinity_{val}.csv", index=False)
        print(f"\n-- per-animal affinity: {val} --")
        print(aff.round(4).to_string(index=False))
        for c in ("Lilac_new", "Copper_new", "Lilac_original", "Copper_original"):
            s = aff[aff.focal_cohort == c]
            if len(s):
                n_lil = int((s.closer_to == "Lilac").sum())
                print(f"  {c}: {n_lil}/{len(s)} closer to original LILAC, "
                      f"{len(s) - n_lil}/{len(s)} closer to original COPPER")

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "gps_source": str(a.gps), "window": "post-deployment only (>= 2025-08-01)",
        "roost_proxy": f"dusk >= {DUSK_FROM_UTC}:00 UTC and dawn <= {DAWN_TO_UTC}:00 UTC "
                       f"(18:00 and 07:00 local); collars record 03:00-16:00 UTC only, so no "
                       f"true overnight position exists and a shared roost is scored only when "
                       f"the pair is together at BOTH ends of the night",
        "sleep_radius_m": a.sleep_radius, "min_shared_nights": a.min_shared_nights,
        "range_cell_m": a.cell_m, "min_fixes_for_range": a.min_fixes,
        "cohorts": cohorts,
        "why_independent": "Uses the full GPS record, not the fusion-hour extract, so it does "
                           "not inherit the collar-dependent fusion quorum or its saturation.",
        "caveats": [
            "Sleeping-site sharing at 100 m cannot distinguish one sleeping tree from an "
            "adjacent one; it indexes site fidelity, not literal co-sleeping.",
            "Bhattacharyya overlap is symmetric and ignores time - two animals using the same "
            "ground in different months score as overlapping.",
        ],
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
