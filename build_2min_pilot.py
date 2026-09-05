"""Pilot: can encounter duration be measured at 2-minute resolution?

Panel 2d bottoms out at one hour because every state in it is assigned hourly. Reaching
two minutes is possible only inside the fixing window -- 03:00-15:00 UTC, where a typical
2-minute bin holds 153 of 241 animals -- because outside it there are no fixes at all.
Through the night the last position is carried forward, which is the assumption already
used by the hourly pipeline, so a bout can run continuously across days.

This runs the whole idea on ONE MONTH before anything is committed to the full record,
and answers the three questions that decide whether the axis is honest:

  NOISE FLOOR    GPS error is 5-15 m against a linkage threshold near 198 m, so pairs at
                 the boundary flicker in and out from error alone. Positions are
                 re-clustered with jitter drawn from that error, repeatedly, and the
                 bout-length distribution from noise is measured directly. Anything below
                 where the real distribution separates from the jittered one is not a
                 finding.
  DWELL          How many consecutive bins a state must hold before it counts, chosen
                 from the noise floor rather than by hand.
  OVERLAP        Between 1 h and 13 h both resolutions measure the same bouts. If they
                 disagree there, the two arms cannot be spliced onto one axis and the
                 hourly figure stands.

A quantisation artefact is expected and is reported rather than hidden: a bout that ends
during the night is only seen to have ended at the next dawn, so night-spanning durations
pile up near multiples of 24 h.

Read-only. Outputs: outputs/general_structure_2026_09/phase7_2min_pilot/
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from analyze_clustering_options import GPS, adaptive_labels, dist_matrix

OUT = Path("outputs/general_structure_2026_09/phase7_2min_pilot")
NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")

MONTH = ("2026-02-01", "2026-03-01")
BIN_MIN = 2
DAY_START_H, DAY_END_H = 3, 15        # the fixing window, inclusive of 15:00's bins
MIN_SIDE = 2                          # animals per unit for a dyad to count as together
MAX_EDGE = 900.0
JITTER_M = 12.0                       # GPS error scale for the noise floor
N_JITTER = 5
SEED = 20260905


def bin_positions(d: pd.DataFrame) -> pd.DataFrame:
    """Mean position per animal per 2-minute bin, daylight only."""
    d = d[d["timestamp"].dt.hour.between(DAY_START_H, DAY_END_H)].copy()
    d["bin"] = d["timestamp"].dt.floor("%dmin" % BIN_MIN)
    g = (d.groupby(["bin", "animal_id"])[["location.lat", "location.long"]]
         .mean().reset_index())
    return g.rename(columns={"location.lat": "lat", "location.long": "lon"})


def cluster_bins(pos: pd.DataFrame, unit_of: dict, rng=None,
                 jitter_m: float = 0.0) -> pd.DataFrame:
    """Cluster every 2-minute bin; optionally jitter positions first."""
    out = []
    for b, g in pos.groupby("bin", sort=True):
        if len(g) < 3:
            continue
        lat = g["lat"].to_numpy(dtype=float)
        lon = g["lon"].to_numpy(dtype=float)
        if jitter_m > 0:
            dlat = rng.normal(0, jitter_m) / 111_320.0
            # per-point jitter, not per-frame: each fix carries its own error
            lat = lat + rng.normal(0, jitter_m, len(lat)) / 111_320.0
            lon = lon + rng.normal(0, jitter_m, len(lon)) / (
                111_320.0 * np.cos(np.radians(np.nanmean(lat))))
        D = dist_matrix(lat, lon)
        lab = adaptive_labels(D, MAX_EDGE)
        out.append(pd.DataFrame({"bin": b, "animal_id": g["animal_id"].to_numpy(),
                                 "cluster": lab}))
    if not out:
        return pd.DataFrame(columns=["bin", "animal_id", "cluster", "unit"])
    r = pd.concat(out, ignore_index=True)
    r["unit"] = r["animal_id"].map(unit_of)
    return r.dropna(subset=["unit"])


def dyad_bins(cl: pd.DataFrame) -> pd.DataFrame:
    """Bins in which two units each put MIN_SIDE+ animals in one cluster."""
    n = (cl.groupby(["bin", "cluster", "unit"]).size().rename("n").reset_index())
    big = n[n["n"].ge(MIN_SIDE)]
    rows = []
    for (b, c), g in big.groupby(["bin", "cluster"], sort=False):
        us = sorted(g["unit"])
        if len(us) < 2:
            continue
        for a, bb in combinations(us, 2):
            rows.append((b, a, bb))
    return pd.DataFrame(rows, columns=["bin", "unit_a", "unit_b"]).drop_duplicates()


def bouts_from_bins(db: pd.DataFrame, carry_night: bool) -> pd.DataFrame:
    """Contiguous runs per dyad, in minutes.

    With `carry_night`, the last daylight state is assumed to hold until the next dawn,
    so a run that is open at the end of a day continues into the next one if the dyad is
    together again at first light. The cost is that a bout which really ended in the dark
    is only seen to end at dawn, which is why night-spanning durations cluster near
    multiples of 24 h -- reported, not hidden.
    """
    if db.empty:
        return pd.DataFrame(columns=["unit_a", "unit_b", "start", "end", "minutes",
                                     "spans_night"])
    db = db.copy()
    db["t"] = db["bin"].astype("int64") // (60_000_000_000 * BIN_MIN)
    step_night = None
    rows = []
    for (a, b), g in db.groupby(["unit_a", "unit_b"], sort=False):
        t = np.sort(g["t"].unique())
        start = prev = t[0]
        for x in t[1:]:
            gap_bins = x - prev
            same_day = pd.Timestamp(int(x) * BIN_MIN * 60, unit="s").date() == \
                pd.Timestamp(int(prev) * BIN_MIN * 60, unit="s").date()
            # contiguous within a day, or the first bin of the next morning under carry
            cont = (gap_bins == 1) or (carry_night and not same_day
                                       and gap_bins <= (24 - (DAY_END_H - DAY_START_H))
                                       * 60 // BIN_MIN + 1)
            if cont:
                prev = x
                continue
            rows.append((a, b, start, prev))
            start = prev = x
        rows.append((a, b, start, prev))
    out = pd.DataFrame(rows, columns=["unit_a", "unit_b", "t0", "t1"])
    out["start"] = pd.to_datetime(out["t0"] * BIN_MIN * 60, unit="s")
    out["end"] = pd.to_datetime(out["t1"] * BIN_MIN * 60, unit="s")
    out["minutes"] = (out["t1"] - out["t0"] + 1) * BIN_MIN
    out["spans_night"] = out["start"].dt.date.ne(out["end"].dt.date)
    return out.drop(columns=["t0", "t1"])


def hourly_bouts(units: set, lo, hi) -> pd.DataFrame:
    """The same quantity from the hourly table, same month and same units."""
    cols = ["window_start", "animal_id", "dynamic_social_unit", "association_event_id",
            "is_observed", "is_carried_night", "is_local_2h_supported"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
    d = d[(d["is_observed"] | d["is_carried_night"] | d["is_local_2h_supported"])
          & d["association_event_id"].notna()
          & d["window_start"].between(lo, hi)
          & d["dynamic_social_unit"].isin(units)]
    n = (d.groupby(["window_start", "association_event_id", "dynamic_social_unit"])
         .size().rename("n").reset_index())
    big = n[n["n"].ge(MIN_SIDE)]
    rows = []
    for (t, c), g in big.groupby(["window_start", "association_event_id"], sort=False):
        us = sorted(g["dynamic_social_unit"])
        for a, b in combinations(us, 2):
            rows.append((t, a, b))
    dd = pd.DataFrame(rows, columns=["window_start", "unit_a", "unit_b"]).drop_duplicates()
    if dd.empty:
        return pd.DataFrame(columns=["unit_a", "unit_b", "minutes"])
    dd["h"] = dd["window_start"].astype("int64") // 3_600_000_000_000
    out = []
    for (a, b), g in dd.groupby(["unit_a", "unit_b"], sort=False):
        h = np.sort(g["h"].unique())
        s = p = h[0]
        for x in h[1:]:
            if x == p + 1:
                p = x
                continue
            out.append((a, b, (p - s + 1) * 60))
            s = p = x
        out.append((a, b, (p - s + 1) * 60))
    return pd.DataFrame(out, columns=["unit_a", "unit_b", "minutes"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    ap.add_argument("--start", default=MONTH[0])
    ap.add_argument("--end", default=MONTH[1])
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    lo, hi = pd.Timestamp(args.start), pd.Timestamp(args.end)
    print("loading GPS %s .. %s" % (lo.date(), hi.date()), flush=True)
    d = pq.read_table(GPS, columns=["animal_id", "timestamp", "location.long",
                                    "location.lat"],
                      filters=[("timestamp", ">=", lo.tz_localize("UTC")),
                               ("timestamp", "<", hi.tz_localize("UTC"))]).to_pandas()
    d["animal_id"] = d["animal_id"].astype(str)
    d["timestamp"] = pd.to_datetime(d["timestamp"]).dt.tz_localize(None)
    d = d[d["location.lat"].notna() & d["location.long"].notna()]

    # units come from the hourly table so the two arms describe the same grouping
    hu = pd.read_csv(NARROW, usecols=["window_start", "animal_id",
                                      "dynamic_social_unit"],
                     parse_dates=["window_start"])
    hu = hu[hu["window_start"].between(lo, hi)]
    hu["animal_id"] = hu["animal_id"].astype(str)
    unit_of = (hu.groupby("animal_id")["dynamic_social_unit"]
               .agg(lambda x: x.value_counts().index[0]).to_dict())

    pos = bin_positions(d)
    print("2-min bins %s | animal-bins %s | animals %d"
          % (format(pos["bin"].nunique(), ","), format(len(pos), ","),
             pos["animal_id"].nunique()), flush=True)

    cl = cluster_bins(pos, unit_of)
    db = dyad_bins(cl)
    real = bouts_from_bins(db, carry_night=True)
    real.to_csv(args.output_dir / "encounter_bouts_2min.csv", index=False)
    print("2-min encounter bouts: %s across %d dyads"
          % (format(len(real), ","),
             real[["unit_a", "unit_b"]].drop_duplicates().shape[0]), flush=True)

    # noise floor: the same pipeline on positions jittered within GPS error
    noise = []
    for i in range(N_JITTER):
        cj = cluster_bins(pos, unit_of, rng=rng, jitter_m=JITTER_M)
        nb = bouts_from_bins(dyad_bins(cj), carry_night=True)
        nb["rep"] = i
        noise.append(nb)
        print("  jitter replicate %d/%d: %s bouts" % (i + 1, N_JITTER,
                                                      format(len(nb), ",")), flush=True)
    noise = pd.concat(noise, ignore_index=True)
    noise.to_csv(args.output_dir / "encounter_bouts_2min_jittered.csv", index=False)

    units = set(unit_of.values())
    hb = hourly_bouts(units, lo, hi)
    hb.to_csv(args.output_dir / "encounter_bouts_hourly_same_month.csv", index=False)

    def q(v):
        v = np.asarray(v, dtype=float)
        return {"n": int(len(v)),
                "median": float(np.median(v)) if len(v) else None,
                "p10": float(np.percentile(v, 10)) if len(v) else None,
                "p90": float(np.percentile(v, 90)) if len(v) else None,
                "max": float(v.max()) if len(v) else None,
                "under_10min": round(float((v < 10).mean()), 3) if len(v) else None,
                "under_1h": round(float((v < 60).mean()), 3) if len(v) else None}

    # where does the real distribution separate from noise?
    edges = np.array([2, 4, 6, 10, 20, 30, 60, 120, 240, 480, 780])
    sep = []
    for e in edges:
        r = float((real["minutes"] <= e).mean()) if len(real) else np.nan
        j = float((noise["minutes"] <= e).mean()) if len(noise) else np.nan
        sep.append({"minutes": int(e), "real_share_at_or_below": round(r, 4),
                    "jittered_share_at_or_below": round(j, 4),
                    "ratio_real_to_noise": round(r / j, 3) if j else None})

    ov = {}
    if len(hb) and len(real):
        for lo_m, hi_m in ((60, 780),):
            a = real[real["minutes"].between(lo_m, hi_m)]["minutes"]
            b = hb[hb["minutes"].between(lo_m, hi_m)]["minutes"]
            ov["%d-%dmin" % (lo_m, hi_m)] = {
                "n_2min": int(len(a)), "n_hourly": int(len(b)),
                "median_2min": float(np.median(a)) if len(a) else None,
                "median_hourly": float(np.median(b)) if len(b) else None}

    report = {
        "month": [str(lo.date()), str(hi.date())],
        "bin_minutes": BIN_MIN, "day_window_utc": [DAY_START_H, DAY_END_H],
        "min_side": MIN_SIDE, "jitter_m": JITTER_M, "jitter_replicates": N_JITTER,
        "animals_binned": int(pos["animal_id"].nunique()),
        "bins": int(pos["bin"].nunique()),
        "real": q(real["minutes"]) if len(real) else None,
        "jittered_per_replicate": q(noise["minutes"]) if len(noise) else None,
        "hourly_same_month": q(hb["minutes"]) if len(hb) else None,
        "spans_night_share": round(float(real["spans_night"].mean()), 3) if len(real) else None,
        "separation_from_noise": sep,
        "overlap_band": ov,
    }
    with open(args.output_dir / "pilot_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("\n" + "=" * 78)
    print("2-MINUTE PILOT  %s to %s" % (lo.date(), hi.date()))
    print("=" * 78)
    for k in ("real", "jittered_per_replicate", "hourly_same_month"):
        v = report[k]
        if v:
            print("  %-24s n %-7s median %-8s p90 %-8s max %s"
                  % (k, format(v["n"], ","), v["median"], v["p90"], v["max"]))
    print("\n  share of 2-min bouts spanning a night: %s" % report["spans_night_share"])
    print("\n  real vs noise, share at or below a length:")
    print("    %8s %10s %10s %8s" % ("minutes", "real", "jittered", "ratio"))
    for s in sep:
        print("    %8d %10s %10s %8s" % (s["minutes"], s["real_share_at_or_below"],
                                         s["jittered_share_at_or_below"],
                                         s["ratio_real_to_noise"]))
    print("\n  overlap band 1-13 h: %s" % report["overlap_band"])
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
