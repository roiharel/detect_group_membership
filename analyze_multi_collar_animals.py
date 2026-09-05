"""Animals with more than one collar deployment: does a swap hide a transfer?

Reading the deployment metadata with drop_duplicates("animal-id") silently keeps
one row per animal, so a re-collared animal looks like a single continuous
deployment. Four Copper/Lilac animals have two rows. 24AE04_6L7M turned out to
have transferred groups 60 days BEFORE its swap; this checks the other three the
same way.

For each animal: full deployment history, then weekly sleep-site agreement with
original-cohort Copper animals versus original-cohort Lilac animals, from the EAS
night-location table. A transfer shows as a sustained sign change; a swap
artefact would coincide with the swap date and usually a gap in the record.

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
EAS_META = Path(r"Z:\baboon\working\data\processed\2025\acc\movebank_metadata.csv")
REF_META = Path(r"Z:\baboon\working\data\processed\2025\metadata"
                r"\Baboons MBRP Mpala Kenya-reference-data.csv")
LOCAL_META = Path(r"C:\Users\rharel\Documents\Github\MBRP_basic\data\metadata_animal_name.csv")
MEMBERSHIP = Path(r"C:\Users\rharel\Documents\New project\outputs"
                  r"\canonical_robust_hourly_membership_shared_full_20260722"
                  r"\canonical_hourly_membership.parquet")
NEW_COLLARS = {"25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
               "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ",
               "25AB47_2ROS", "25AB49_1V3W", "25AB54_6F3G", "25AB55_OH7I"}


def deployments() -> tuple[pd.DataFrame, set[str]]:
    """All deployment rows from every source, NOT deduplicated.

    The sources disagree about which animals were re-collared, so the set of
    multi-deployment animals is the UNION across them; dates are taken from
    whichever source actually has parseable ones (the EAS movebank export has
    its date column mangled to '00:00.0').
    """
    frames, multi = [], set()
    for p in (REF_META, LOCAL_META, EAS_META):
        if not p.exists():
            print(f"  (unavailable: {p.name})")
            continue
        m = pd.read_csv(p, dtype=str, low_memory=False)
        m.columns = [c.replace("-", "_") for c in m.columns]
        m["source"] = p.name
        n = m.groupby("animal_id").size()
        multi |= set(n[n > 1].index)
        frames.append(m)
        print(f"  {p.name}: {len(m)} rows, {int((n > 1).sum())} animals with >1 deployment")
    if not frames:
        raise SystemExit("no deployment metadata reachable")
    allm = pd.concat(frames, ignore_index=True)
    allm["on"] = pd.to_datetime(allm.get("deploy_on_date"), errors="coerce")
    allm["off"] = pd.to_datetime(allm.get("deploy_off_date"), errors="coerce")
    return allm, multi


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-weeks", type=int, default=3)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "multi_collar_animals_2026-09-03")
    a = ap.parse_args()
    if not NIGHT_LOC.exists():
        raise SystemExit(f"EAS not reachable: {NIGHT_LOC}")
    a.output_dir.mkdir(parents=True, exist_ok=True)

    nl = pd.read_parquet(NIGHT_LOC)
    for c in ("animal_id", "group_id", "cluster_united"):
        nl[c] = nl[c].astype(str)
    nl["date"] = pd.to_datetime(nl["date"])
    cl = nl[nl.group_id.isin(["Copper", "Lilac"])].copy()

    print("deployment metadata sources:")
    meta, multi_all = deployments()
    multi = sorted(multi_all & set(cl.animal_id.unique()))
    print(f"\nCopper/Lilac animals with >1 deployment (union of sources): "
          f"{len(multi)} -> {multi}\n")

    piv = cl.pivot_table(index="date", columns="animal_id", values="cluster_united",
                         aggfunc="first")
    grp = cl.drop_duplicates("animal_id").set_index("animal_id").group_id
    orig_lil = [x for x in piv.columns if grp.get(x) == "Lilac" and x not in NEW_COLLARS]
    orig_cop = [x for x in piv.columns if grp.get(x) == "Copper" and x not in NEW_COLLARS]

    summary, series = [], []
    for animal in multi:
        rows = (meta[meta.animal_id == animal]
                .dropna(subset=["on"])
                .drop_duplicates(subset=["tag_id", "on"])
                .sort_values("on"))
        swaps = list(rows.on.iloc[1:]) if len(rows) > 1 else []

        print("=" * 78)
        sex = rows.animal_sex.dropna().iloc[0] if rows.animal_sex.notna().any() else "?"
        age = rows.animal_comments.dropna().iloc[0] if rows.animal_comments.notna().any() else "?"
        print(f"{animal}   origin {grp.get(animal)}   sex {sex}   age {age}")
        for r in rows.itertuples(index=False):
            off = f"{r.off:%Y-%m-%d}" if pd.notna(r.off) else "(open)"
            print(f"   tag {str(r.tag_id):>6}   {r.on:%Y-%m-%d} -> {off}")

        if animal not in piv.columns:
            print("   no night locations\n")
            continue
        recs = []
        for night, r in piv.iterrows():
            if pd.isna(r.get(animal)):
                continue
            me = r[animal]
            v = {"night": night}
            for lab, mem in (("lilac", orig_lil), ("copper", orig_cop)):
                vals = [1.0 if r[o] == me else 0.0 for o in mem
                        if o != animal and o in r and pd.notna(r[o])]
                v[lab] = np.mean(vals) if vals else np.nan
            recs.append(v)
        d = pd.DataFrame(recs).dropna()
        if d.empty:
            print("   no comparable nights\n")
            continue
        d["week"] = d.night.dt.to_period("W").dt.start_time
        w = (d.groupby("week", as_index=False)[["lilac", "copper"]].mean()
             .assign(animal_id=animal))
        w["diff"] = w.lilac - w.copper
        series.append(w)

        gaps = d.night.diff().dt.days
        own = "lilac" if grp.get(animal) == "Lilac" else "copper"
        other = "copper" if own == "lilac" else "lilac"
        # sustained sign change: 4+ consecutive weeks favouring the other group
        favour_other = (w["diff"] < 0) if own == "lilac" else (w["diff"] > 0)
        runs, best, start = 0, 0, None
        best_start = None
        for wk, fl in zip(w.week, favour_other):
            if fl:
                runs += 1
                start = wk if runs == 1 else start
                if runs > best:
                    best, best_start = runs, start
            else:
                runs = 0
        print(f"   nights {len(d)}   weeks {len(w)}   largest record gap "
              f"{gaps.max():.0f} d")
        print(f"   mean same-site with own group {w[own].mean():.3f}, "
              f"with other {w[other].mean():.3f}")
        print(f"   longest run favouring the OTHER group: {best} weeks"
              + (f", from {best_start:%Y-%m-%d}" if best_start is not None else ""))
        verdict = "TRANSFER LIKELY" if best >= 4 else "no sustained switch"
        if best >= 4 and swaps:
            days = min(abs((best_start - s).days) for s in swaps)
            verdict += f" - starts {days} d from nearest collar swap"
        print(f"   -> {verdict}\n")
        summary.append({"animal_id": animal, "origin": grp.get(animal),
                        "sex": sex, "age": age,
                        "n_deployments": len(rows),
                        "swap_dates": ";".join(f"{s:%Y-%m-%d}" for s in swaps),
                        "nights": len(d), "max_gap_days": float(gaps.max()),
                        "mean_own": float(w[own].mean()),
                        "mean_other": float(w[other].mean()),
                        "longest_run_favouring_other_weeks": int(best),
                        "run_start": str(best_start.date()) if best_start is not None else None,
                        "verdict": verdict})

    out = pd.DataFrame(summary)
    out.to_csv(a.output_dir / "multi_collar_summary.csv", index=False)
    if series:
        pd.concat(series, ignore_index=True).to_csv(
            a.output_dir / "weekly_site_agreement.csv", index=False)

    m = pd.read_parquet(MEMBERSHIP, columns=["animal_id", "dynamic_social_unit",
                                             "dynamic_assignment"])
    print("=" * 78)
    print("what the membership pipeline assigned each of them:")
    for animal in multi:
        s = m[m.animal_id == animal]
        if len(s):
            print(f"  {animal:16s} units={sorted(s.dynamic_social_unit.unique())}  "
                  f"assignment={sorted(s.dynamic_assignment.unique())}")

    (a.output_dir / "metadata.json").write_text(json.dumps({
        "night_locations": str(NIGHT_LOC), "metadata_source": meta["source"].iloc[0],
        "animals_checked": multi,
        "transfer_rule": "4+ consecutive weeks in which sleep-site agreement favours the "
                         "non-origin group",
        "why": "drop_duplicates on animal-id hides second deployments; a re-collared animal "
               "can look continuous, and a real transfer near a swap can be mistaken for a "
               "re-identification artefact (or vice versa).",
        "summary": summary,
    }, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
