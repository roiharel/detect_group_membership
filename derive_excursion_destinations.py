"""Resolve where each individual excursion actually went.

The excursion table records that an animal was away from its origin unit and whether it
was alone or with another unit, but not WHICH unit. Every figure that needed the
destination has resolved it ad hoc, and the last one left 5 of 15 away-nights
unidentified. This derives it once, properly, with an explicit resolution rate.

METHOD. For each away-night of each excursion, take the association clusters the focal
animal occupied that night, find the animals from other units sharing those clusters,
and assign the modal other unit. A night with no co-clustered animal from another unit
is left `unresolved` rather than guessed -- that happens when the animal is with an
uncollared group, or when the other unit's collars were not reporting.

Two products:
  excursion_night_destinations.csv   one row per away-night, with its destination
  excursion_destinations.csv         one row per excursion: modal destination, the set
                                     of units visited, and how much of it resolved
  transfer_edges.csv                 origin -> destination edge list for the network

Outputs: outputs/general_structure_2026_09/phase4g_excursion_destinations/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
BASE = Path("outputs/general_structure_2026_09/phase4b_individual_axis")
OUT = Path("outputs/general_structure_2026_09/phase4g_excursion_destinations")

AWAY_WITH_OTHERS = {"other", "mixed_without_origin_unit"}
RULES = {"dominant": "excursions_dominant_gap0.csv",
         "any_away": "excursions_any_away_gap0.csv"}


def load_hourly() -> pd.DataFrame:
    cols = ["window_start", "animal_id", "origin_group", "dynamic_social_unit",
            "social_context", "association_event_id", "is_observed",
            "is_carried_night", "is_local_2h_supported"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
# A row's state is KNOWN if it has a position from a real fix -- exactly at the
# hour, borrowed from a fix within 60 min (`local_2h`), or carried across the
# night. Omitting local_2h while accepting carried_night was indefensible and
# became visible only when coverage improved: at 17:00 the 2026-09-05 build
# supplies 106,147 local_2h rows where the frozen build had carried_night, so the
# old predicate silently discarded almost the whole hour. Interpolated rows are
# deliberately NOT known: their position is inferred, not observed.
    d = d[d["is_observed"] | d["is_carried_night"]
          | d["is_local_2h_supported"]].copy()
    d["night"] = (d["window_start"] - pd.Timedelta(hours=10)).dt.normalize()
    return d


def night_destinations(d: pd.DataFrame) -> pd.DataFrame:
    """For every away animal-night, the modal other unit sharing the focal's clusters."""
    away = d[d["social_context"].isin(AWAY_WITH_OTHERS)]
    if away.empty:
        return pd.DataFrame(columns=["animal_id", "origin_group", "night",
                                     "destination", "support_hours"])
    # the (hour, cluster) slots each away animal occupied
    slots = away[["animal_id", "origin_group", "night", "window_start",
                  "association_event_id"]].drop_duplicates()
    # candidate co-members: everyone else in those slots, from a different unit
    pool = d[["window_start", "association_event_id", "animal_id",
              "dynamic_social_unit"]].rename(
        columns={"animal_id": "other_animal",
                 "dynamic_social_unit": "other_unit"})
    j = slots.merge(pool, on=["window_start", "association_event_id"], how="left")
    j = j[j["other_animal"].ne(j["animal_id"])
          & j["other_unit"].notna()
          & j["other_unit"].ne(j["origin_group"])]
    if j.empty:
        res = pd.DataFrame(columns=["animal_id", "origin_group", "night",
                                    "destination", "support_hours"])
    else:
        cnt = (j.groupby(["animal_id", "origin_group", "night", "other_unit"])
               .size().rename("hours").reset_index())
        idx = cnt.groupby(["animal_id", "origin_group", "night"])["hours"].idxmax()
        res = (cnt.loc[idx]
               .rename(columns={"other_unit": "destination", "hours": "support_hours"})
               .reset_index(drop=True))
    # every away-night, resolved or not
    allnights = (away.groupby(["animal_id", "origin_group", "night"]).size()
                 .rename("away_hours").reset_index())
    out = allnights.merge(res, on=["animal_id", "origin_group", "night"], how="left")
    out["destination"] = out["destination"].fillna("unresolved")
    out["support_hours"] = out["support_hours"].fillna(0).astype(int)
    return out


def per_excursion(exc: pd.DataFrame, nd: pd.DataFrame) -> pd.DataFrame:
    """Attach destinations to each excursion by matching its away-night span."""
    nd = nd.copy()
    nd["night"] = pd.to_datetime(nd["night"])
    rows = []
    for r in exc.itertuples():
        s, e = pd.Timestamp(r.start_night), pd.Timestamp(r.end_night)
        sub = nd[nd["animal_id"].eq(r.animal_id) & nd["night"].between(s, e)]
        resolved = sub[sub["destination"].ne("unresolved")]
        counts = resolved["destination"].value_counts()
        rows.append({
            "animal_id": r.animal_id,
            "origin_group": r.origin_group,
            "start_night": s.date(),
            "end_night": e.date(),
            "away_nights": int(r.away_nights),
            "alone_nights": int(r.alone_nights),
            "other_nights": int(r.other_nights),
            "joined_nights_seen": int(len(sub)),
            "joined_nights_resolved": int(len(resolved)),
            "destination": counts.index[0] if len(counts) else None,
            "n_destinations": int(len(counts)),
            "destinations": ";".join(counts.index) if len(counts) else "",
            "settlement_candidate": bool(r.settlement_candidate),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    d = load_hourly()
    nd = night_destinations(d)
    nd.to_csv(args.output_dir / "excursion_night_destinations.csv", index=False)

    report = {
        "source": str(NARROW),
        "method": ("modal other unit among animals sharing the focal animal's "
                   "association clusters on that night; nights with no co-clustered "
                   "animal from another unit are left unresolved"),
        "away_animal_nights": int(len(nd)),
        "resolved": int((nd["destination"].ne("unresolved")).sum()),
        "resolution_rate": round(float((nd["destination"].ne("unresolved")).mean()), 4),
        "distinct_destinations": int(nd.loc[nd["destination"].ne("unresolved"),
                                            "destination"].nunique()),
        "by_rule": {},
    }

    for rule, fname in RULES.items():
        exc = pd.read_csv(BASE / fname, parse_dates=["start_night", "end_night"])
        pe = per_excursion(exc, nd)
        pe.to_csv(args.output_dir / ("excursion_destinations_%s.csv" % rule),
                  index=False)
        joined = pe[pe["other_nights"].gt(0)]
        with_dest = joined[joined["destination"].notna()]
        edges = (with_dest.groupby(["origin_group", "destination"])
                 .agg(excursions=("animal_id", "size"),
                      animals=("animal_id", "nunique"),
                      away_nights=("away_nights", "sum"),
                      settlements=("settlement_candidate", "sum")).reset_index())
        edges.to_csv(args.output_dir / ("transfer_edges_%s.csv" % rule), index=False)
        report["by_rule"][rule] = {
            "excursions": int(len(pe)),
            "excursions_reaching_another_unit": int(len(joined)),
            "with_a_named_destination": int(len(with_dest)),
            "destination_rate": round(float(len(with_dest) / max(1, len(joined))), 4),
            "multi_destination_excursions": int((with_dest["n_destinations"] > 1).sum()),
            "edges": int(len(edges)),
            "origin_units": int(edges["origin_group"].nunique()),
            "destination_units": int(edges["destination"].nunique()),
            "units_in_network": int(len(set(edges["origin_group"])
                                        | set(edges["destination"]))),
            "animals": int(with_dest["animal_id"].nunique()),
        }

    with open(args.output_dir / "destination_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("=" * 78)
    print("EXCURSION DESTINATIONS")
    print("=" * 78)
    print("away animal-nights %d, resolved %d (%.1f%%), %d distinct destinations"
          % (report["away_animal_nights"], report["resolved"],
             100 * report["resolution_rate"], report["distinct_destinations"]))
    for rule, b in report["by_rule"].items():
        print("\n--- %s rule ---" % rule)
        for k, v in b.items():
            print("    %-34s %s" % (k, v))
    top = (pd.read_csv(args.output_dir / "transfer_edges_dominant.csv")
           .sort_values("excursions", ascending=False).head(12))
    print("\ntop transfer edges (dominant rule):")
    print(top.to_string(index=False))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
