"""Derive the three things the level-2 figure cannot get from the summary tables.

1. ENCOUNTER BOUTS AT HOURLY RESOLUTION. The dyad-day opportunity table can only say a
   bout lasted N days, which floors every encounter at 24 h and hides the short ones.
   Two units sharing an association cluster is an hourly fact in the narrow export, so
   the elapsed duration of an encounter is derivable directly: contiguous hours in which
   both units have two or more animals in the same cluster. Gaps of up to GAP_H hours
   are bridged, because 02:00 is a known coverage hole and would otherwise cut every
   bout in two every night.

2. THE PER-DYAD LEDGER OF CROSSINGS. A dyad can carry group encounters and single-animal
   crossings at the same time, and the same dyad can carry several single-animal
   crossings years apart. Those are separate dispersal events, not one merge and not one
   edge -- so each dyad gets a count of group-encounter bouts, of single-animal
   episodes, and of the distinct animals behind them.

3. TWO LILAC WEEKS. The association network in the week Lilac is most modular and in a
   matched week where it is not, with communities, so the figure can show what a
   modularity score of 0.45 versus 0.00 actually looks like.

Outputs: outputs/general_structure_2026_09/phase4h_level2_inputs/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
BASE = Path("outputs/general_structure_2026_09")
FROZEN = BASE / "phase4d_axis_b_frozen/weekly_network_metrics_frozen.csv"
OPP = BASE / "phase1_opportunity/opportunity_dyad_day.csv"
NIGHTS = BASE / "phase4g_excursion_destinations/excursion_night_destinations.csv"
OUT = BASE / "phase4h_level2_inputs"

MIN_SIDE = 2          # animals per unit for a dyad-hour to count as a group encounter
GAP_H = 2             # hours of missing coverage bridged inside one bout
EXAMPLE_UNIT = "Lilac"
WELL_COVERED = 12
MIN_EDGE = 0.05     # association rate below this is not drawn in the inset
AWAY = ("other", "isolated", "mixed_without_origin_unit")
# four weeks of ONE unit, spanning its range, so every example ties to one time
# series: the top, a three-community week, a weakly modular one, and a flat one
EXAMPLE_WEEKS = (("Lilac", 0.45), ("Lilac", 0.275), ("Lilac", 0.14), ("Lilac", 0.0))
SPRING_K = 1.6        # spring-layout optimal distance; see lilac_weeks()
SEED = 20260904


def load_hourly() -> pd.DataFrame:
    cols = ["window_start", "animal_id", "dynamic_social_unit", "association_event_id",
            "social_context", "origin_group", "is_observed", "is_carried_night",
            "is_local_2h_supported"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
# A row's state is KNOWN if it has a position from a real fix -- exactly at the
# hour, borrowed from a fix within 60 min (`local_2h`), or carried across the
# night. Omitting local_2h while accepting carried_night was indefensible and
# became visible only when coverage improved: at 17:00 the 2026-09-05 build
# supplies 106,147 local_2h rows where the frozen build had carried_night, so the
# old predicate silently discarded almost the whole hour. Interpolated rows are
# deliberately NOT known: their position is inferred, not observed.
    known = (d["is_observed"] | d["is_carried_night"]
             | d["is_local_2h_supported"])
    d = d[known & d["association_event_id"].notna()]
    return d.drop(columns=["is_observed", "is_carried_night",
                           "is_local_2h_supported"])


def runs(sorted_ints, gap):
    """(start, end) of each run in a sorted int array, bridging gaps of <= `gap`."""
    out = []
    if len(sorted_ints) == 0:
        return out
    start = prev = sorted_ints[0]
    for x in sorted_ints[1:]:
        if x - prev <= gap + 1:
            prev = x
        else:
            out.append((start, prev))
            start = prev = x
    out.append((start, prev))
    return out


def dyad_hours(d: pd.DataFrame) -> pd.DataFrame:
    """Every (hour, unit-pair) where both units had MIN_SIDE+ animals in one cluster."""
    n = (d.groupby(["window_start", "association_event_id", "dynamic_social_unit"])
         .size().rename("n").reset_index())
    big = n[n["n"].ge(MIN_SIDE)]
    # a cluster-hour only matters if two or more units clear the threshold in it
    k = big.groupby(["window_start", "association_event_id"])["dynamic_social_unit"]
    multi = k.transform("size").ge(2)
    big = big[multi]
    rows = []
    for (ws, ev), g in big.groupby(["window_start", "association_event_id"], sort=False):
        us = sorted(g["dynamic_social_unit"])
        for i in range(len(us)):
            for j in range(i + 1, len(us)):
                rows.append((ws, us[i], us[j]))
    out = pd.DataFrame(rows, columns=["window_start", "unit_a", "unit_b"])
    return out.drop_duplicates()


def encounter_bouts_hourly(dh: pd.DataFrame) -> pd.DataFrame:
    """Contiguous hours of group encounter, per dyad, in elapsed hours."""
    if dh.empty:
        return pd.DataFrame(columns=["unit_a", "unit_b", "start", "end", "hours"])
    epoch = dh["window_start"].astype("int64") // 3_600_000_000_000
    dh = dh.assign(hr=epoch)
    rows = []
    for (a, b), g in dh.groupby(["unit_a", "unit_b"], sort=False):
        hrs = np.sort(g["hr"].unique())
        for s, e in runs(hrs, GAP_H):
            rows.append({"unit_a": a, "unit_b": b,
                         "start": pd.Timestamp(s * 3600, unit="s"),
                         "end": pd.Timestamp(e * 3600, unit="s"),
                         "hours": float(e - s + 1)})
    return pd.DataFrame(rows)


def split_bouts(d: pd.DataFrame) -> pd.DataFrame:
    """Within-group split at hourly resolution: the hourly face of modularity.

    Weekly modularity is a week-scale summary and cannot say how long a group stayed
    divided. The hourly fact underneath it is simpler: in this hour, do the unit's
    animals occupy two or more association clusters, each holding two or more of them?
    Restricted to unit-hours with WELL_COVERED or more animals present, for the same
    reason the weekly metric is -- below that a cluster count is a coverage statistic.
    """
    per = (d.groupby(["window_start", "dynamic_social_unit", "association_event_id"])
           .size().rename("n").reset_index())
    tot = (per.groupby(["window_start", "dynamic_social_unit"])["n"].sum()
           .rename("unit_n").reset_index())
    big = per[per["n"].ge(2)]
    k = (big.groupby(["window_start", "dynamic_social_unit"])["association_event_id"]
         .nunique().rename("clusters").reset_index())
    st = tot.merge(k, on=["window_start", "dynamic_social_unit"], how="left")
    st["clusters"] = st["clusters"].fillna(0)
    st = st[st["unit_n"].ge(WELL_COVERED)]
    st = st[st["clusters"].ge(2)]
    if st.empty:
        return pd.DataFrame(columns=["unit", "start", "end", "hours"])
    st = st.assign(hr=st["window_start"].astype("int64") // 3_600_000_000_000)
    rows = []
    for u, g in st.groupby("dynamic_social_unit", sort=False):
        for a, b in runs(np.sort(g["hr"].unique()), GAP_H):
            rows.append({"unit": u, "start": pd.Timestamp(a * 3600, unit="s"),
                         "end": pd.Timestamp(b * 3600, unit="s"),
                         "hours": float(b - a + 1)})
    return pd.DataFrame(rows)


def excursion_bouts(d: pd.DataFrame) -> pd.DataFrame:
    """Individual excursions at hourly resolution: contiguous hours away from origin.

    The excursion table counts away-NIGHTS, which floors every trip at 24 h and cannot
    see an animal that leaves and returns inside a day. `social_context` is assigned
    hourly, so the elapsed hours of a trip are already in the export.
    """
    a = d[d["social_context"].isin(AWAY)]
    if a.empty:
        return pd.DataFrame(columns=["animal_id", "origin_group", "start", "end",
                                     "hours", "alone_hours"])
    a = a.assign(hr=a["window_start"].astype("int64") // 3_600_000_000_000)
    rows = []
    for (an, og), g in a.groupby(["animal_id", "origin_group"], sort=False):
        hrs = np.sort(g["hr"].unique())
        alone = set(g.loc[g["social_context"].eq("isolated"), "hr"])
        for s, e in runs(hrs, GAP_H):
            rows.append({"animal_id": an, "origin_group": og,
                         "start": pd.Timestamp(s * 3600, unit="s"),
                         "end": pd.Timestamp(e * 3600, unit="s"),
                         "hours": float(e - s + 1),
                         "alone_hours": float(sum(1 for h in range(s, e + 1)
                                                  if h in alone))})
    return pd.DataFrame(rows)


def crossing_episodes() -> pd.DataFrame:
    """Single-animal crossings, one row per animal per contiguous stay in a unit.

    The same dyad can carry several of these months or years apart. They are separate
    dispersal events and must never collapse into one edge, or into a group merge.
    """
    nd = pd.read_csv(NIGHTS, parse_dates=["night"])
    nd = nd[nd["destination"].ne("unresolved")]
    epoch = nd["night"].astype("int64") // 86_400_000_000_000
    nd = nd.assign(day=epoch)
    rows = []
    for (an, og, de), g in nd.groupby(["animal_id", "origin_group", "destination"],
                                      sort=False):
        days = np.sort(g["day"].unique())
        for s, e in runs(days, 1):
            rows.append({"animal_id": an, "origin_group": og, "destination": de,
                         "start_night": pd.Timestamp(s * 86400, unit="s").date(),
                         "end_night": pd.Timestamp(e * 86400, unit="s").date(),
                         "nights": int(e - s + 1)})
    return pd.DataFrame(rows)


def dyad_ledger(bouts: pd.DataFrame, cross: pd.DataFrame) -> pd.DataFrame:
    """One row per unit pair: group encounters and single-animal crossings side by side."""
    def key(a, b):
        return "|".join(sorted((a, b)))
    b = bouts.assign(pair=[key(a, c) for a, c in zip(bouts["unit_a"], bouts["unit_b"])])
    c = cross.assign(pair=[key(a, d) for a, d in
                           zip(cross["origin_group"], cross["destination"])])
    gb = b.groupby("pair").agg(encounter_bouts=("hours", "size"),
                               encounter_hours=("hours", "sum"),
                               longest_bout_h=("hours", "max"))
    gc = c.groupby("pair").agg(crossing_episodes=("nights", "size"),
                               crossing_nights=("nights", "sum"),
                               crossing_animals=("animal_id", "nunique"),
                               longest_stay_nights=("nights", "max"))
    led = gb.join(gc, how="outer").fillna(0).reset_index()
    led[["unit_a", "unit_b"]] = led["pair"].str.split("|", expand=True)
    led["kind"] = np.where(led["encounter_bouts"].gt(0) & led["crossing_episodes"].gt(0),
                           "both",
                           np.where(led["encounter_bouts"].gt(0), "encounter_only",
                                    "crossing_only"))
    return led.sort_values(["kind", "encounter_hours"], ascending=[True, False])


def overlap_check(bouts: pd.DataFrame, cross: pd.DataFrame) -> dict:
    """How much of the group-encounter record sits inside a crossing's stay.

    If a dyad's encounters all happen while one of its animals is living in the other
    unit, the two chords in the figure are not independent evidence for that dyad.
    """
    def key(a, b):
        return "|".join(sorted((a, b)))
    b = bouts.assign(pair=[key(a, c) for a, c in zip(bouts["unit_a"], bouts["unit_b"])])
    c = cross.assign(pair=[key(a, d) for a, d in
                           zip(cross["origin_group"], cross["destination"])])
    per, tot_h, ins_h = [], 0.0, 0.0
    for p in sorted(set(b["pair"]) & set(c["pair"])):
        gb, gc = b[b["pair"].eq(p)], c[c["pair"].eq(p)]
        inside = np.zeros(len(gb), bool)
        for r in gc.itertuples():
            s = pd.Timestamp(r.start_night)
            e = pd.Timestamp(r.end_night) + pd.Timedelta(days=1)
            inside |= (gb["start"] < e).to_numpy() & (gb["end"] >= s).to_numpy()
        per.append({"pair": p, "bouts": int(len(gb)),
                    "bouts_inside_a_stay": int(inside.sum()),
                    "hours": float(gb["hours"].sum()),
                    "hours_inside_a_stay": float(gb.loc[inside, "hours"].sum()),
                    "share_hours_inside": round(
                        float(gb.loc[inside, "hours"].sum() / gb["hours"].sum()), 3)})
        tot_h += float(gb["hours"].sum())
        ins_h += float(gb.loc[inside, "hours"].sum())
    df = pd.DataFrame(per)
    return {"shared_dyads": len(df),
            "encounter_hours_on_shared_dyads": round(tot_h, 1),
            "hours_inside_a_stay": round(ins_h, 1),
            "share_inside": round(ins_h / tot_h, 3) if tot_h else None,
            "dyads_over_half_inside": int((df["share_hours_inside"] > 0.5).sum())
            if len(df) else 0,
            "per_dyad": df.to_dict("records")}


def pick_week(g: pd.DataFrame, target: float, ref_n=None):
    """The well-covered week whose modularity is closest to `target`.

    At target 0 the tie is broken on collar count, so the flat example is matched to the
    modular one on coverage and cannot be dismissed as a thinner week.
    """
    g = g.copy()
    g["d"] = (g["modularity"] - target).abs()
    if target <= 0.001 and ref_n is not None:
        g["dn"] = (g["n_animals_observed"] - ref_n).abs()
        return g.sort_values(["d", "dn", "period_start"]).iloc[0]
    return g.sort_values(["d", "period_start"]).iloc[0]


def example_networks(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Association networks for a graded set of weeks, so a Q value gets a picture."""
    wk = pd.read_csv(FROZEN, parse_dates=["period_start"])
    rows, meta = [], {}
    for unit, target in EXAMPLE_WEEKS:
        g = wk[wk["dynamic_social_unit"].eq(unit)
               & wk["n_animals_observed"].ge(WELL_COVERED)].copy()
        ref = g.loc[g["modularity"].idxmax(), "n_animals_observed"]
        wkrow = pick_week(g, target, ref_n=ref)
        key = "%s_%s" % (unit, str(target).replace(".", "p"))
        tag = target
        EXAMPLE_UNIT_LOCAL = unit
        s = wkrow["period_start"]
        e = s + pd.Timedelta(days=7)
        sub = d[d["window_start"].ge(s) & d["window_start"].lt(e)
                & d["dynamic_social_unit"].eq(EXAMPLE_UNIT_LOCAL)]
        animals = sorted(sub["animal_id"].unique())
        idx = {a: i for i, a in enumerate(animals)}
        present = sub.groupby("animal_id")["window_start"].nunique()
        # co-cluster count per dyad, over hours in which both animals were present
        pairs: dict = {}
        for _, cl in sub.groupby(["window_start", "association_event_id"], sort=False):
            aa = sorted(cl["animal_id"].unique())
            for i in range(len(aa)):
                for j in range(i + 1, len(aa)):
                    pairs[(aa[i], aa[j])] = pairs.get((aa[i], aa[j]), 0) + 1
        both = sub[["window_start", "animal_id"]].drop_duplicates()
        hours_of = both.groupby("animal_id")["window_start"].apply(set).to_dict()
        G = nx.Graph()
        G.add_nodes_from(animals)
        for (a, b), cnt in pairs.items():
            denom = len(hours_of[a] & hours_of[b])
            if denom == 0:
                continue
            w = cnt / denom
            if w >= MIN_EDGE:
                G.add_edge(a, b, weight=w)
        comms = list(nx.algorithms.community.greedy_modularity_communities(
            G, weight="weight")) if G.number_of_edges() else [set(animals)]
        comm_of = {a: i for i, c in enumerate(comms) for a in c}
        # a plain spring layout collapses each community to a point, because within
        # a community the association weights are all near 1; a larger optimal
        # distance keeps the members apart enough to be countable at thumbnail size
        pos = nx.spring_layout(G, weight="weight", k=SPRING_K, seed=SEED,
                               iterations=400)
        for a in animals:
            rows.append({"week": key, "unit": unit, "level": tag,
                         "period_start": s.date(), "animal_id": a,
                         "node": idx[a], "community": comm_of.get(a, 0),
                         "x": float(pos[a][0]), "y": float(pos[a][1]),
                         "hours_present": int(present.get(a, 0)),
                         "kind": "node", "weight": np.nan})
        for a, b, dta in G.edges(data=True):
            rows.append({"week": key, "unit": unit, "level": tag,
                         "period_start": s.date(), "animal_id": a,
                         "node": idx[a], "community": idx[b], "x": np.nan, "y": np.nan,
                         "hours_present": np.nan, "kind": "edge",
                         "weight": float(dta["weight"])})
        meta[key] = {"unit": unit, "level": tag, "period_start": str(s.date()),
                     "modularity": round(float(wkrow["modularity"]), 3),
                     "n_animals": int(len(animals)),
                     "n_communities": len(comms),
                     "largest_community_fraction": round(
                         max(len(c) for c in comms) / len(animals), 3),
                     "edges_drawn": G.number_of_edges()}
    return pd.DataFrame(rows), meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    ap.add_argument("--networks-only", action="store_true",
                    help="rebuild only the example networks, reusing the saved bout "
                         "tables; the hourly state passes are the slow part")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    d = load_hourly()
    print("hourly rows with a cluster: %d" % len(d))

    if args.networks_only:
        net, meta = example_networks(d)
        net.to_csv(args.output_dir / "example_week_networks.csv", index=False)
        print("example networks only:")
        for k, v in meta.items():
            print("    %-22s %s" % (k, v))
        return

    dh = dyad_hours(d)
    bouts = encounter_bouts_hourly(dh)
    bouts.to_csv(args.output_dir / "encounter_bouts_hourly.csv", index=False)

    cross = crossing_episodes()
    cross.to_csv(args.output_dir / "crossing_episodes.csv", index=False)

    led = dyad_ledger(bouts, cross)
    led.to_csv(args.output_dir / "dyad_ledger.csv", index=False)

    ov = overlap_check(bouts, cross)
    splits = split_bouts(d)
    splits.to_csv(args.output_dir / "split_bouts_hourly.csv", index=False)
    exc = excursion_bouts(d)
    exc.to_csv(args.output_dir / "excursion_bouts_hourly.csv", index=False)
    net, meta = example_networks(d)
    net.to_csv(args.output_dir / "example_week_networks.csv", index=False)

    report = {
        "min_side": MIN_SIDE, "gap_bridged_hours": GAP_H,
        "dyad_hours": int(len(dh)),
        "encounter_bouts": int(len(bouts)),
        "encounter_dyads": int(led["encounter_bouts"].gt(0).sum()),
        "bout_hours": {k: round(float(v), 1) for k, v in
                       bouts["hours"].describe().to_dict().items()},
        "bouts_under_a_day": int((bouts["hours"] < 24).sum()),
        "bouts_under_6h": int((bouts["hours"] < 6).sum()),
        "crossing_episodes": int(len(cross)),
        "crossing_animals": int(cross["animal_id"].nunique()),
        "crossing_dyads": int(led["crossing_episodes"].gt(0).sum()),
        "dyads_with_repeat_crossings": int((led["crossing_episodes"] > 1).sum()),
        "dyads_by_kind": led["kind"].value_counts().to_dict(),
        "overlap": {k: v for k, v in ov.items() if k != "per_dyad"},
        "split_bouts": int(len(splits)),
        "split_units": int(splits["unit"].nunique()) if len(splits) else 0,
        "split_hours_median": float(splits["hours"].median()) if len(splits) else None,
        "split_under_a_day": int((splits["hours"] < 24).sum()) if len(splits) else 0,
        "excursion_bouts": int(len(exc)),
        "excursion_animals": int(exc["animal_id"].nunique()) if len(exc) else 0,
        "excursion_hours_median": float(exc["hours"].median()) if len(exc) else None,
        "excursion_under_a_day": int((exc["hours"] < 24).sum()) if len(exc) else 0,
        "example_weeks": meta,
    }
    with open(args.output_dir / "level2_inputs_report.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)
    pd.DataFrame(ov["per_dyad"]).to_csv(
        args.output_dir / "encounter_crossing_overlap.csv", index=False)

    print("=" * 78)
    print("LEVEL-2 INPUTS")
    print("=" * 78)
    for k, v in report.items():
        if k in ("bout_hours", "overlap", "example_weeks", "dyads_by_kind"):
            print("  %s:" % k)
            for kk, vv in v.items():
                print("      %-30s %s" % (kk, vv))
        else:
            print("  %-32s %s" % (k, v))
    print("\ntop dyads by encounter hours:")
    print(led.head(10).to_string(index=False))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
