"""Case study: 24AE04_6L7M, a Lilac-origin animal that roosts with Copper.

In the roost analysis this animal has the largest Copper affinity of any in the
study (-0.47), despite being labelled Lilac. This script asks when that happened,
whether it was ever with Lilac, and what the pipeline made of it.

Roost proxy: median of the last 15 minutes of recording each day (dusk) and the
first 15 minutes of the next (dawn); the collars stop at 16:00 UTC so no true
overnight fix exists.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

PROJECT = Path(__file__).resolve().parent
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
MEMBERSHIP = Path(r"C:\Users\rharel\Documents\New project\outputs"
                  r"\canonical_robust_hourly_membership_shared_full_20260722"
                  r"\canonical_hourly_membership.parquet")
LOCAL_META = Path(r"C:\Users\rharel\Documents\Github\MBRP_basic\data\metadata_animal_name.csv")
R_EARTH_M = 6371000.0
NEW_LILAC = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
             "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"}
NEW_COPPER = {"25AB47_2ROS", "25AB49_1V3W", "25AB54_6F3G", "25AB55_OH7I"}


def haversine_m(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = (np.sin((p2 - p1) / 2.0) ** 2
         + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2.0) ** 2)
    return 2.0 * R_EARTH_M * np.arcsin(np.sqrt(a))


def load_gps(groups):
    f = pq.ParquetFile(GPS)
    cols = ["animal_id", "timestamp", "location.long", "location.lat", "group_id"]
    keep = []
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=cols)
        t = t.filter(pc.is_in(t["group_id"], value_set=pd.array(groups).__arrow_array__()))
        if t.num_rows:
            keep.append(t.to_pandas())
    d = pd.concat(keep, ignore_index=True).rename(
        columns={"location.long": "lon", "location.lat": "lat", "group_id": "origin_group"})
    d["timestamp"] = pd.to_datetime(d["timestamp"], utc=True).dt.tz_localize(None)
    for c in ("animal_id", "origin_group"):
        d[c] = d[c].astype(str)
    return d.dropna(subset=["lat", "lon"])


def roost_positions(d: pd.DataFrame, window_min: int = 15) -> pd.DataFrame:
    """Median of the last `window_min` of each day, and the first of the next."""
    d = d.sort_values(["animal_id", "timestamp"])
    d["day"] = d.timestamp.dt.floor("D")
    out = []
    for lab, keyfun in (("dusk", "last"), ("dawn", "first")):
        ext = d.groupby(["animal_id", "day"]).timestamp.transform(keyfun)
        if lab == "dusk":
            sel = d[d.timestamp >= ext - pd.Timedelta(minutes=window_min)].copy()
            sel["night"] = sel.day
        else:
            sel = d[d.timestamp <= ext + pd.Timedelta(minutes=window_min)].copy()
            sel["night"] = sel.day - pd.Timedelta(days=1)
        g = (sel.groupby(["animal_id", "night"], as_index=False, observed=True)
                .agg(lat=("lat", "median"), lon=("lon", "median"), n=("lat", "size")))
        out.append(g.assign(edge=lab))
    return pd.concat(out, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--focal", default="24AE04_6L7M")
    ap.add_argument("--radius", type=float, default=100.0)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "case_24AE04_transfer_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(LOCAL_META, dtype=str).drop_duplicates("animal-id")
    row = meta[meta["animal-id"] == a.focal]
    print(f"=== {a.focal} ===")
    if len(row):
        r = row.iloc[0]
        for c in ["animal-group-id", "animal-sex", "animal-comments", "deploy-on-date",
                  "deploy-off-date", "animal-mass", "animal-offspring", "tag-id"]:
            if c in row.columns:
                print(f"  {c:22s} {r.get(c)}")

    d = load_gps(["Copper", "Lilac"])
    print(f"\nGPS: {d[d.animal_id == a.focal].timestamp.min():%Y-%m-%d} to "
          f"{d[d.animal_id == a.focal].timestamp.max():%Y-%m-%d}, "
          f"{(d.animal_id == a.focal).sum():,} fixes")

    cohort = {}
    for aid, grp in d.drop_duplicates("animal_id")[["animal_id", "origin_group"]].values:
        cohort[aid] = ("Lilac_new" if aid in NEW_LILAC else
                       "Copper_new" if aid in NEW_COPPER else f"{grp}_original")

    roost = roost_positions(d)
    roost = roost[roost.n >= 2]
    piv = {e: (roost[roost.edge == e].pivot_table(index="night", columns="animal_id", values="lat"),
               roost[roost.edge == e].pivot_table(index="night", columns="animal_id", values="lon"))
           for e in ("dusk", "dawn")}

    lil = [x for x, c in cohort.items() if c == "Lilac_original" and x != a.focal]
    cop = [x for x, c in cohort.items() if c == "Copper_original" and x != a.focal]
    def pos(edge, night, animal):
        la, lo = piv[edge]
        if night not in la.index or animal not in la.columns:
            return None
        x, y = la.at[night, animal], lo.at[night, animal]
        return None if (np.isnan(x) or np.isnan(y)) else (x, y)

    rows = []
    for night in piv["dusk"][0].index:
        # a night with no focal position tells us nothing - it must be dropped, not
        # scored as "not together", which would fabricate zeros outside the collar's life
        fd, fw = pos("dusk", night, a.focal), pos("dawn", night, a.focal)
        if fd is None or fw is None:
            continue
        rec = {"night": night}
        for lab, members in (("lilac", lil), ("copper", cop)):
            vals = []
            for other in members:
                od, ow = pos("dusk", night, other), pos("dawn", night, other)
                if od is None or ow is None:
                    continue          # partner unobserved: exclude, do not score 0
                together = (haversine_m(fd[0], fd[1], od[0], od[1]) <= a.radius
                            and haversine_m(fw[0], fw[1], ow[0], ow[1]) <= a.radius)
                vals.append(1.0 if together else 0.0)
            rec[f"frac_with_{lab}"] = float(np.mean(vals)) if vals else np.nan
            rec[f"n_{lab}"] = len(vals)
        rows.append(rec)
    nightly = pd.DataFrame(rows).dropna(subset=["frac_with_lilac", "frac_with_copper"])
    nightly["month"] = nightly.night.values.astype("datetime64[M]")
    nightly.to_csv(a.output_dir / "nightly_roost_affinity.csv", index=False)

    monthly = (nightly.groupby("month", as_index=False)
               .agg(nights=("night", "size"), with_lilac=("frac_with_lilac", "mean"),
                    with_copper=("frac_with_copper", "mean")))
    monthly["lilac_minus_copper"] = monthly.with_lilac - monthly.with_copper
    monthly.to_csv(a.output_dir / "monthly_roost_affinity.csv", index=False)
    print(f"\n=== monthly roost affinity (within {a.radius:g} m at dusk AND dawn) ===")
    print(monthly.round(3).to_string(index=False))

    m = pd.read_parquet(MEMBERSHIP, columns=["window_start", "animal_id", "origin_group",
                                             "dynamic_social_unit", "dynamic_assignment"])
    f = m[m.animal_id == a.focal].copy()
    if len(f):
        f["month"] = pd.to_datetime(f.window_start).values.astype("datetime64[M]")
        print(f"\n=== what the pipeline assigned ===")
        print(f.groupby(["month", "dynamic_social_unit"], observed=True)
               .size().unstack(fill_value=0).to_string())
        print("\ndynamic_assignment values:")
        print(f.dynamic_assignment.value_counts().to_string())

    switch = monthly[monthly.lilac_minus_copper < 0]
    first_neg = switch.month.min() if len(switch) else None
    (a.output_dir / "case_summary.json").write_text(json.dumps({
        "focal": a.focal,
        "roost_radius_m": a.radius,
        "roost_proxy": "median of last 15 min of recording (dusk) and first 15 min of the "
                       "following day (dawn); a shared roost requires both",
        "first_month_favouring_copper": str(first_neg) if first_neg is not None else None,
        "monthly": monthly.to_dict("records"),
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
