"""Individual-axis diagnostic: merge isolation and dispersal into one state space,
and quantify how much of the raw "isolated" record is a collar-coverage artefact.

Isolation is not a separate phenomenon from dispersal. An animal away from its origin
group is either alone in the landscape or with another group, and those are graded
depths of one individual-level process. This script builds that four-state ledger and
tests the identifiability rule that makes "alone" a real observation rather than a
statement about which collars happened to be reporting.

Identifiability rule (the "located reference cluster"): an animal-hour is scored
`alone` only if
  1. the focal animal's own position is observed, not carried; and
  2. at least MIN_REF_CLUSTER collared animals of its origin group are observed in the
     same hour, co-clustered with each other, in a cluster the focal animal is not in.

Condition 2 is what removes the sparse-collar cases. If a group has one collar, there
is never a reference cluster. If it has two and they separate, neither side reaches
size 2, so neither animal is scored alone -- both go to `unresolvable`. The same rule
applies to `with_other`, so "it left" is only asserted when the origin group was
locatable somewhere else at the time.

Outputs: outputs/general_structure_2026_09/phase4b_individual_axis/

Source: the frozen narrow export (1,924,104 animal-hours, 350 animals, 26 origin
groups, 2024-03-01 to 2026-07-22).

NOTE ON THE NIGHT DEFINITION. The nightly aggregation here uses a simple
window_start - 10h calendar cut, not the canonical 16:00-06:00 biological night used by
export_canonical_nightly_membership.py. Nightly counts below are therefore a prototype
and are not interchangeable with the canonical 81,695 animal-nights. Reconciling the two
is part of Phase 4b proper.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
OUT = Path("outputs/general_structure_2026_09/phase4b_individual_axis")

# social_context values that mean "in a cluster containing the origin group".
# mixed_with_origin_present is a between-group merge with the origin group present, so
# the animal has not left; that belongs to axis A, not here.
WITH_ORIGIN = {"with_origin", "mixed_with_origin_present"}
# away from the origin group but with some other unit
AWAY_WITH_OTHERS = {"other", "mixed_without_origin_unit"}

MIN_REF_CLUSTER = 2
SUSTAINED_NIGHTS = 7  # matches the membership builder's dynamic-reassignment rule

USECOLS = [
    "window_start",
    "animal_id",
    "origin_group",
    "dynamic_social_unit",
    "social_context",
    "association_event_id",
    "is_observed",
    "is_carried_night", "is_local_2h_supported",
]


def load_hourly(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, usecols=USECOLS, parse_dates=["window_start"])
    d["obs"] = d["is_observed"].astype(int)
    # A night hour is not re-fixed, but the pipeline carries the evening state forward
    # (is_carried_night), and carried rows retain social_context, association_event_id
    # and temp_group_size. So there are TWO different questions and they need two
    # different denominators -- see two_tier_ledger().
# A row's state is KNOWN if it has a position from a real fix -- exactly at the
# hour, borrowed from a fix within 60 min (`local_2h`), or carried across the
# night. Omitting local_2h while accepting carried_night was indefensible and
# became visible only when coverage improved: at 17:00 the 2026-09-05 build
# supplies 106,147 local_2h rows where the frozen build had carried_night, so the
# old predicate silently discarded almost the whole hour. Interpolated rows are
# deliberately NOT known: their position is inferred, not observed.
    d["known"] = (d["is_observed"] | d["is_carried_night"]
                  | d["is_local_2h_supported"]).astype(int)
    n_obs = (
        d.groupby(["window_start", "origin_group"])["obs"].sum().rename("n_obs_origin")
    )
    n_known = (
        d.groupby(["window_start", "origin_group"])["known"].sum().rename("n_known_origin")
    )
    d = d.join(n_obs, on=["window_start", "origin_group"])
    d = d.join(n_known, on=["window_start", "origin_group"])
    d["partners_origin"] = d["n_obs_origin"] - d["obs"]
    d["partners_known"] = d["n_known_origin"] - d["known"]
    return d


def add_reference_cluster(d: pd.DataFrame) -> pd.DataFrame:
    """Size of the largest same-origin-group cluster the focal animal is NOT in.

    Clusters are association events, counted over observed animals only. For each
    hour x origin group we take the two largest same-unit cluster sizes; an animal
    sitting in the largest one is compared against the second largest, everyone else
    against the largest.
    """
    obs = d.loc[d["is_observed"], ["window_start", "origin_group", "association_event_id", "animal_id"]]
    sizes = (
        obs.groupby(["window_start", "origin_group", "association_event_id"])
        .size()
        .rename("own_k")
        .reset_index()
        .sort_values(["window_start", "origin_group", "own_k"], ascending=[True, True, False])
    )
    sizes["rank"] = sizes.groupby(["window_start", "origin_group"]).cumcount()
    wide = sizes.pivot_table(
        index=["window_start", "origin_group"], columns="rank", values="own_k", aggfunc="first"
    )
    second = wide[1] if 1 in wide.columns else pd.Series(index=wide.index, dtype=float)
    hour_unit = pd.concat([wide[0].rename("c1"), second.rename("c2")], axis=1)
    hour_unit["c2"] = hour_unit["c2"].fillna(0)
    hour_unit = hour_unit.reset_index()

    keys = ["window_start", "origin_group", "association_event_id"]
    d = d.merge(sizes[keys + ["own_k"]], on=keys, how="left")
    d = d.merge(hour_unit, on=["window_start", "origin_group"], how="left")
    # animals in the largest cluster are compared against the runner-up
    focal_is_largest = d["own_k"].eq(d["c1"]) & d["own_k"].gt(d["c2"])
    d["ref_elsewhere"] = np.where(focal_is_largest, d["c2"], d["c1"])
    d["ref_elsewhere"] = d["ref_elsewhere"].fillna(-1)
    return d


def assign_state(d: pd.DataFrame, min_ref: int = MIN_REF_CLUSTER) -> pd.DataFrame:
    resolvable = d["is_observed"] & d["ref_elsewhere"].ge(min_ref)
    d["resolvable"] = resolvable
    d["state"] = np.select(
        [
            d["social_context"].isin(WITH_ORIGIN) & d["is_observed"],
            resolvable & d["social_context"].eq("isolated"),
            resolvable & d["social_context"].isin(AWAY_WITH_OTHERS),
        ],
        ["with_origin", "alone", "with_other"],
        default="unresolvable",
    )
    return d


def two_tier_ledger(d: pd.DataFrame, min_ref: int = MIN_REF_CLUSTER) -> dict:
    """Two different questions, two different denominators.

    The single "44.2% unresolvable" figure conflated them, and it read as far more
    alarming than it should.

    TIER 1 -- IS THE STATE KNOWN? A night hour is not re-fixed, but the evening state is
    carried forward and carried rows keep social_context, association_event_id and
    temp_group_size. So the animal's social context IS known overnight, under the
    assumption that it stays at its sleep site and the interaction persists. This is the
    right denominator for anything about DURATION: how long an animal was alone, how long
    an encounter lasted, whether a state spanned a night.

    TIER 2 -- IS THERE INDEPENDENT EVIDENCE? A carried hour adds no new fix. Its located
    reference cluster is inherited from the last daytime observation, not confirmed
    again. This is the right denominator for anything that COUNTS observations: effective
    sample size, model degrees of freedom, precision.

    Both tiers use the same located-reference-cluster rule; they differ only in whether a
    carried hour counts. The sparse-collar conclusion rests on 1- and 2-collar contexts,
    not on night gaps, so it is unaffected by which tier is used -- reported here to make
    that explicit rather than leaving it to be assumed.
    """
    n = len(d)
    away = d["social_context"].eq("isolated") | d["social_context"].isin(AWAY_WITH_OTHERS)

    known_res = d["known"].eq(1) & d["ref_elsewhere"].ge(min_ref)
    obs_res = d["is_observed"] & d["ref_elsewhere"].ge(min_ref)

    def block(label, base, resolvable, question, use_for):
        # `with_origin` needs only that the animal's own state is available -- being
        # WITH the group is read directly. The located-reference-cluster rule is
        # required only to assert the animal is AWAY from its group, which is the claim
        # sparse collars can manufacture.
        st = np.select(
            [
                d["social_context"].isin(WITH_ORIGIN) & base,
                resolvable & d["social_context"].eq("isolated"),
                resolvable & d["social_context"].isin(AWAY_WITH_OTHERS),
            ],
            ["with_origin", "alone", "with_other"],
            default="unresolvable",
        )
        vc = pd.Series(st).value_counts()
        return {
            "tier": label,
            "question": question,
            "use_for": use_for,
            "animal_hours": int(n),
            "with_origin": int(vc.get("with_origin", 0)),
            "alone": int(vc.get("alone", 0)),
            "with_other": int(vc.get("with_other", 0)),
            "unresolvable": int(vc.get("unresolvable", 0)),
            "unresolvable_share": round(float(vc.get("unresolvable", 0)) / n, 4),
            "away_hours_retained": int(((st == "alone") | (st == "with_other")).sum()),
            "away_hours_raw": int(away.sum()),
            "share_of_raw_away_retained": round(
                float(((st == "alone") | (st == "with_other")).sum()) / max(1, int(away.sum())), 4
            ),
        }

    t1 = block(
        "state known",
        d["known"].eq(1),
        known_res,
        "is the animal's social context known in this hour?",
        "durations, spans, whether a state crossed a night",
    )
    t2 = block(
        "independent evidence",
        d["is_observed"],
        obs_res,
        "is there a fresh fix confirming it in this hour?",
        "effective sample size, precision, model degrees of freedom",
    )
    return {
        "min_ref_cluster": min_ref,
        "tier_1_state_known": t1,
        "tier_2_independent_evidence": t2,
        "carried_hours_that_only_tier_1_retains": int(
            (known_res & ~obs_res).sum()
        ),
        "verdict": (
            "State is known for %.1f%% of animal-hours and independently confirmed for "
            "%.1f%%. The earlier single figure of %.1f%% unresolvable was the tier-2 "
            "number reported as though it were tier 1, which understated what the record "
            "knows. Use tier 1 for durations and tier 2 for anything that counts "
            "observations."
            % (100 * (1 - t1["unresolvable_share"]),
               100 * (1 - t2["unresolvable_share"]),
               100 * t2["unresolvable_share"])
        ),
    }


def isolation_rule_ladder(d: pd.DataFrame) -> pd.DataFrame:
    """How many raw 'isolated' animal-hours survive each candidate support rule."""
    iso = d[d["social_context"].eq("isolated")]
    obs = iso[iso["is_observed"]]
    rules = [
        ("raw isolated animal-hours", iso, None),
        ("focal position observed", obs, None),
        ("+ >=1 observed group-mate anywhere", obs, obs["partners_origin"].ge(1)),
        ("+ >=3 observed group-mates anywhere", obs, obs["partners_origin"].ge(3)),
        ("+ located reference cluster >=2", obs, obs["ref_elsewhere"].ge(2)),
        ("+ located reference cluster >=3", obs, obs["ref_elsewhere"].ge(3)),
        ("+ located reference cluster >=4", obs, obs["ref_elsewhere"].ge(4)),
    ]
    rows = []
    for label, frame, mask in rules:
        sel = frame if mask is None else frame[mask]
        rows.append(
            {
                "rule": label,
                "animal_hours": len(sel),
                "share_of_raw": round(len(sel) / len(iso), 4),
                "animals": sel["animal_id"].nunique(),
                "origin_groups": sel["origin_group"].nunique(),
            }
        )
    return pd.DataFrame(rows)


def sparse_collar_cases(d: pd.DataFrame) -> pd.DataFrame:
    """The cases the rule exists to remove: 1- and 2-collar observation contexts."""
    iso = d[d["social_context"].eq("isolated") & d["is_observed"]]
    rows = []
    for n_obs, label in [(1, "only the focal animal observed"), (2, "exactly 2 animals observed")]:
        sel = iso[iso["n_obs_origin"].eq(n_obs)]
        passing = sel[sel["ref_elsewhere"].ge(MIN_REF_CLUSTER)]
        rows.append(
            {
                "observation_context": label,
                "n_obs_origin": n_obs,
                "animal_hours": len(sel),
                "animals": sel["animal_id"].nunique(),
                "origin_groups": sel["origin_group"].nunique(),
                "passing_rule": len(passing),
                "share_passing": round(len(passing) / len(sel), 4) if len(sel) else 0.0,
            }
        )
    return pd.DataFrame(rows)


def per_group_survival(d: pd.DataFrame) -> pd.DataFrame:
    iso = d[d["social_context"].eq("isolated") & d["is_observed"]].copy()
    iso["supported"] = iso["ref_elsewhere"].ge(MIN_REF_CLUSTER)
    t = iso.groupby("origin_group").agg(
        observed_isolated_hours=("supported", "size"),
        supported_hours=("supported", "sum"),
        median_observed_group_animals=("n_obs_origin", "median"),
    )
    t["pct_surviving"] = (100 * t["supported_hours"] / t["observed_isolated_hours"]).round(1)
    return t.sort_values("observed_isolated_hours", ascending=False).reset_index()


def nightly_states(d: pd.DataFrame, rule: str) -> pd.DataFrame:
    """Collapse hours to nights under one of two explicit rules.

    'any_away'  -- the night is away if ANY hour was resolvably away (permissive)
    'dominant'  -- the night takes its most frequent hourly state (conservative)

    The two are reported side by side because the choice moves the excursion count by
    an order of magnitude. It must be pinned, not defaulted.
    """
    d = d.copy()
    d["night"] = (d["window_start"] - pd.Timedelta(hours=10)).dt.normalize()

    if rule == "any_away":
        def pick(s: pd.Series) -> str:
            counts = s.value_counts()
            alone = counts.get("alone", 0)
            other = counts.get("with_other", 0)
            if alone or other:
                return "alone" if alone >= other else "with_other"
            if "with_origin" in counts.index:
                return "with_origin"
            return "unresolvable"
    elif rule == "dominant":
        order = ["alone", "with_other", "with_origin", "unresolvable"]

        def pick(s: pd.Series) -> str:
            counts = s.value_counts()
            top = counts.max()
            for k in order:
                if counts.get(k, 0) == top:
                    return k
            return "unresolvable"
    else:
        raise ValueError(rule)

    n = (
        d.groupby(["animal_id", "origin_group", "night"])["state"]
        .apply(pick)
        .reset_index(name="state")
        .sort_values(["animal_id", "night"])
    )
    return n


def build_excursions(nightly: pd.DataFrame, gap_nights: int) -> pd.DataFrame:
    """Runs of away-from-origin nights, bridging up to `gap_nights` unresolvable nights.

    A `with_origin` night always ends a run -- the animal came back. Only
    `unresolvable` nights are bridged, and only up to the tolerance.
    """
    out = []
    for animal, g in nightly.groupby("animal_id", sort=False):
        g = g.sort_values("night")
        cur = None
        pending_gap = 0
        for row in g.itertuples():
            if row.state in ("alone", "with_other"):
                if cur is None:
                    cur = {
                        "animal_id": animal,
                        "origin_group": row.origin_group,
                        "start_night": row.night,
                        "end_night": row.night,
                        "alone_nights": 0,
                        "other_nights": 0,
                        "bridged_unresolvable": 0,
                    }
                cur["end_night"] = row.night
                cur["bridged_unresolvable"] += pending_gap
                pending_gap = 0
                if row.state == "alone":
                    cur["alone_nights"] += 1
                else:
                    cur["other_nights"] += 1
            elif row.state == "unresolvable" and cur is not None:
                pending_gap += 1
                if pending_gap > gap_nights:
                    out.append(cur)
                    cur = None
                    pending_gap = 0
            else:
                if cur is not None:
                    out.append(cur)
                    cur = None
                pending_gap = 0
        if cur is not None:
            out.append(cur)

    if not out:
        return pd.DataFrame(
            columns=[
                "animal_id", "origin_group", "start_night", "end_night", "alone_nights",
                "other_nights", "away_nights", "span_nights", "depth", "settlement_candidate",
            ]
        )

    e = pd.DataFrame(out)
    e["away_nights"] = e["alone_nights"] + e["other_nights"]
    e["span_nights"] = (e["end_night"] - e["start_night"]).dt.days + 1
    e["depth"] = np.select(
        [e["other_nights"].eq(0), e["alone_nights"].eq(0)],
        ["alone_only", "joined_only"],
        default="alone_and_joined",
    )
    e["settlement_candidate"] = e["other_nights"].ge(SUSTAINED_NIGHTS)
    return e


def summarise_excursions(e: pd.DataFrame) -> dict:
    if e.empty:
        return {"excursions": 0}
    return {
        "excursions": int(len(e)),
        "animals": int(e["animal_id"].nunique()),
        "origin_groups": int(e["origin_group"].nunique()),
        "median_away_nights": float(e["away_nights"].median()),
        "p90_away_nights": float(e["away_nights"].quantile(0.9)),
        "max_away_nights": int(e["away_nights"].max()),
        "single_night_share": round(float((e["away_nights"] == 1).mean()), 4),
        "depth_alone_only": int((e["depth"] == "alone_only").sum()),
        "depth_joined_only": int((e["depth"] == "joined_only").sum()),
        "depth_alone_and_joined": int((e["depth"] == "alone_and_joined").sum()),
        "median_nights_alone_only": float(e.loc[e["depth"] == "alone_only", "away_nights"].median()),
        "median_nights_joined_only": float(e.loc[e["depth"] == "joined_only", "away_nights"].median()),
        "settlement_candidates": int(e["settlement_candidate"].sum()),
        "settlement_candidate_animals": int(e.loc[e["settlement_candidate"], "animal_id"].nunique()),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--narrow", type=Path, default=NARROW)
    ap.add_argument("--output-dir", type=Path, default=OUT)
    ap.add_argument("--min-ref-cluster", type=int, default=MIN_REF_CLUSTER)
    ap.add_argument(
        "--gap-nights", type=int, nargs="+", default=[0, 1, 2, 3, 7, 14],
        help="unresolvable-night tolerances to report; the first is the primary",
    )
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    d = load_hourly(args.narrow)
    d = add_reference_cluster(d)
    d = assign_state(d, args.min_ref_cluster)

    # --- the identifiability diagnostic -------------------------------------
    ladder = isolation_rule_ladder(d)
    sparse = sparse_collar_cases(d)
    per_group = per_group_survival(d)
    ladder.to_csv(args.output_dir / "isolation_rule_ladder.csv", index=False)
    sparse.to_csv(args.output_dir / "sparse_collar_cases.csv", index=False)
    per_group.to_csv(args.output_dir / "isolation_survival_by_group.csv", index=False)

    # --- the four-state ledger ---------------------------------------------
    hourly_states = (
        d["state"].value_counts().rename("animal_hours").rename_axis("state").reset_index()
    )
    hourly_states["share"] = (hourly_states["animal_hours"] / len(d)).round(4)
    hourly_states["animals"] = hourly_states["state"].map(
        lambda s: d.loc[d["state"].eq(s), "animal_id"].nunique()
    )
    hourly_states.to_csv(args.output_dir / "hourly_state_ledger.csv", index=False)

    unresolved = d[d["state"].eq("unresolvable")]
    unresolvable_reasons = {
        "focal_position_carried_not_observed": int((~unresolved["is_observed"]).sum()),
        "observed_but_no_located_reference_cluster": int(
            (unresolved["is_observed"] & unresolved["ref_elsewhere"].lt(args.min_ref_cluster)).sum()
        ),
        "of_which_mixed_tie_or_unclear": int(
            (
                unresolved["is_observed"]
                & unresolved["ref_elsewhere"].lt(args.min_ref_cluster)
                & unresolved["social_context"].eq("mixed_tie_or_unclear")
            ).sum()
        ),
    }

    # --- excursions, under both nightly rules and every gap tolerance -------
    report: dict = {
        "source": str(args.narrow),
        "animal_hours": int(len(d)),
        "animals": int(d["animal_id"].nunique()),
        "origin_groups": int(d["origin_group"].nunique()),
        "min_ref_cluster": args.min_ref_cluster,
        "sustained_nights": SUSTAINED_NIGHTS,
        "night_definition": "window_start - 10h calendar cut; NOT the canonical 16:00-06:00 night",
        "unresolvable_reasons": unresolvable_reasons,
        "two_tier_ledger": two_tier_ledger(d, args.min_ref_cluster),
        "nightly_rules": {},
    }

    for rule in ("any_away", "dominant"):
        nightly = nightly_states(d, rule)
        counts = nightly["state"].value_counts()
        block = {
            "animal_nights": int(len(nightly)),
            "state_counts": {k: int(counts.get(k, 0)) for k in
                             ["with_origin", "with_other", "alone", "unresolvable"]},
            "resolvable_share": round(1 - counts.get("unresolvable", 0) / len(nightly), 4),
            "gap_sensitivity": {},
        }
        nightly.to_csv(args.output_dir / f"nightly_states_{rule}.csv", index=False)
        for gap in args.gap_nights:
            e = build_excursions(nightly, gap)
            block["gap_sensitivity"][str(gap)] = summarise_excursions(e)
            if gap == args.gap_nights[0]:
                e.to_csv(args.output_dir / f"excursions_{rule}_gap{gap}.csv", index=False)
        report["nightly_rules"][rule] = block

    with open(args.output_dir / "individual_axis_report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, default=str)

    # --- console summary ----------------------------------------------------
    print("=== isolation support ladder (raw isolated animal-hours = %d) ==="
          % int(ladder.loc[0, "animal_hours"]))
    print(ladder.to_string(index=False))
    print("\n=== the cases the rule exists to remove ===")
    print(sparse.to_string(index=False))
    print("\n=== four-state hourly ledger ===")
    print(hourly_states.to_string(index=False))
    print("\nunresolvable breakdown:", json.dumps(unresolvable_reasons, indent=2))

    tt = report["two_tier_ledger"]
    print("\n=== two tiers: what is KNOWN vs what is INDEPENDENTLY CONFIRMED ===")
    print("  %-22s %12s %8s %11s %13s %9s"
          % ("tier", "with_origin", "alone", "with_other", "unresolvable", "unres. %"))
    for k in ("tier_1_state_known", "tier_2_independent_evidence"):
        b = tt[k]
        print("  %-22s %12d %8d %11d %13d %8.1f%%"
              % (b["tier"], b["with_origin"], b["alone"], b["with_other"],
                 b["unresolvable"], 100 * b["unresolvable_share"]))
    print()
    for k in ("tier_1_state_known", "tier_2_independent_evidence"):
        b = tt[k]
        print("  %-22s away-hours retained %6d of %6d raw (%.1f%%)"
              % (b["tier"], b["away_hours_retained"], b["away_hours_raw"],
                 100 * b["share_of_raw_away_retained"]))
        print("  %-22s use for: %s" % ("", b["use_for"]))
    print("\n  hours only tier 1 retains: %d"
          % tt["carried_hours_that_only_tier_1_retains"])
    print("\n  " + tt["verdict"])
    print("\n=== excursions ===")
    for rule, block in report["nightly_rules"].items():
        primary = block["gap_sensitivity"][str(args.gap_nights[0])]
        print("\nnightly rule %-9s  animal-nights %d, resolvable %.1f%%"
              % (rule, block["animal_nights"], 100 * block["resolvable_share"]))
        print("  states:", block["state_counts"])
        print("  primary (gap=%d): %d excursions, %d animals, %d groups; median %.0f away-nights"
              % (args.gap_nights[0], primary["excursions"], primary["animals"],
                 primary["origin_groups"], primary["median_away_nights"]))
        print("  depth: alone-only %d (median %.0f nights) | joined-only %d (median %.0f nights) | both %d"
              % (primary["depth_alone_only"], primary["median_nights_alone_only"],
                 primary["depth_joined_only"], primary["median_nights_joined_only"],
                 primary["depth_alone_and_joined"]))
        print("  settlement candidates (>=%d nights with another unit): %d in %d animals"
              % (SUSTAINED_NIGHTS, primary["settlement_candidates"],
                 primary["settlement_candidate_animals"]))
        print("  gap sensitivity (excursions): "
              + ", ".join("gap %s -> %d" % (g, v["excursions"])
                          for g, v in block["gap_sensitivity"].items()))
    print("\nwrote %s" % args.output_dir)


if __name__ == "__main__":
    main()
