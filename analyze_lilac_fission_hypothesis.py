"""Did Lilac fission, and did the July-2025 collaring sample the two daughter units?

BACKGROUND
----------
Nine collars deployed 2025-07-31 are labelled Lilac in every available source
(EAS movebank metadata 2026-08-17, EAS reference data 2025-12-05, local copy
2026-04-24). The server actively corrects group names - 16 animals changed since
April - so the absence of a correction here is informative: the field team
considers these animals Lilac.

Yet post-deployment they form a contact block distinct from the seven older
Lilac collars: cross-cohort contact is 0.56x within-cohort, every one of 16
animals prefers its own cohort, and the same-day Copper deployment shows no such
structure (ratio 1.00). EAS demographics put Lilac at 60 animals with 13 collars.

The fission hypothesis: Lilac is not one social unit, and the collaring campaign
put ~9 collars on one daughter unit while ~7 remained on the other.

The canonical pipeline cannot express this. `dynamic_social_unit` reassigns an
animal to another NAMED group; it has no mechanism to split a group into two
daughter units, so all nine remain "Lilac" by construction.

WHAT THIS TESTS
---------------
1. Co-membership: how often do two Lilac animals occupy the same hourly spatial
   cluster (temp_group_id)? Within the old cohort, within the new cohort, and
   between them. Uses ALL hours, not only Copper-Lilac fusion hours.
2. Persistence: is the separation stable across months, or transient?
3. Context: Lilac's split probability and within-group modularity against every
   other group, from the pipeline's own split-detection outputs.
4. Control: Copper, which received collars the same day.

Read-only with respect to existing outputs.
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent
RUN = Path(r"C:\Users\rharel\Documents\New project\outputs"
           r"\canonical_robust_hourly_membership_shared_full_20260722")
MEMBERSHIP = RUN / "canonical_hourly_membership.parquet"
SPLIT_MONTH = RUN / "split_detection_evidence" / "split_evidence_summary_by_dynamic_group_month.csv"
MODULARITY = RUN / "within_split_modularity_gt7_collars" / "shared_full_within_split_modularity_by_group_month.csv"
DEPLOYMENT = pd.Timestamp("2025-08-01")
# Cohort membership is defined by the deployment batch, NOT by a date cutoff on
# first observation: the Lilac collars were deployed 2025-07-31 13:00 and the
# Copper ones 2025-08-01 13:00, so any cutoff near the deployment splits the
# batch. These lists come from the EAS movebank metadata deploy-on-date.
DEPLOY_BATCH = {
    "Lilac": ["25AB37_4X2Y", "25AB38_9Z5A", "25AB39_3B7C", "25AB41_5F9G", "25AB42_2H6I",
              "25AB43_0J4K", "25AB44_7L1M", "25AB45_3N8O", "25AB46_9PSQ"],
    "Copper": ["25AB47_2ROS", "25AB49_1V3W", "25AB54_6F3G", "25AB55_OH7I"],
}
# Animals first seen after this are later deployments (Nov 2025, Mar 2026) and are
# excluded so the control is the same-week batch, not every later collar.
PRE_BATCH_CUTOFF = pd.Timestamp("2025-07-01")


def co_membership(m: pd.DataFrame, animals: list[str]) -> pd.DataFrame:
    """Per pair-month: hours co-observed, hours in the same hourly cluster."""
    d = m[m.animal_id.isin(animals)][["window_start", "animal_id", "temp_group_id"]]
    piv = d.pivot_table(index="window_start", columns="animal_id", values="temp_group_id",
                        aggfunc="first")
    idx = piv.index
    month = idx.values.astype("datetime64[M]")
    rows = []
    for a, b in combinations([c for c in animals if c in piv.columns], 2):
        va, vb = piv[a].values, piv[b].values
        both = pd.notna(va) & pd.notna(vb)
        if not both.any():
            continue
        same = both & (va == vb)
        df = pd.DataFrame({"month": month[both], "same": same[both]})
        g = df.groupby("month", as_index=False).agg(hours=("same", "size"), together=("same", "sum"))
        g["animal_a"], g["animal_b"] = a, b
        rows.append(g)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    out["co_membership_rate"] = out.together / out.hours
    return out


def classify(pairs: pd.DataFrame, cohort: pd.Series) -> pd.DataFrame:
    p = pairs.copy()
    p["ca"] = p.animal_a.map(cohort)
    p["cb"] = p.animal_b.map(cohort)
    p["pair_class"] = np.where(
        p.ca == p.cb, np.where(p.ca == "original", "within original", "within new"), "between cohorts")
    return p


def pooled(p: pd.DataFrame, min_hours: int) -> pd.DataFrame:
    d = p[p.hours >= min_hours]
    g = (d.groupby("pair_class", as_index=False)
           .agg(pairs=("co_membership_rate", "size"), hours=("hours", "sum"),
                together=("together", "sum")))
    g["co_membership_rate"] = g.together / g.hours
    return g


PAIR_RATES = (PROJECT / "outputs" / "copper_lilac_effort_corrected_integration"
              / "copper_lilac_pair_month_contact_rates.csv")


def fine_scale_separation(radius: float = 5.0, min_bins: int = 60) -> pd.DataFrame:
    """Within-Lilac 5 m contact split by collar cohort, per month."""
    p = pd.read_csv(PAIR_RATES, parse_dates=["month"])
    d = p[(p.radius_m == radius) & (p.pair_type == "within_lilac")].copy()
    batch = set(DEPLOY_BATCH["Lilac"])
    d["ca"] = np.where(d.animal_a.isin(batch), "new", "old")
    d["cb"] = np.where(d.animal_b.isin(batch), "new", "old")
    d["cls"] = np.where(d.ca == d.cb,
                        np.where(d.ca == "new", "within new", "within old"), "between")
    g = (d[d.opportunity_bins >= min_bins].groupby(["month", "cls"], as_index=False)
         .agg(opp=("opportunity_bins", "sum"), con=("contact_bins", "sum")))
    g["rate"] = g.con / g.opp
    w = g.pivot(index="month", columns="cls", values="rate").reindex(
        columns=["within old", "within new", "between"])
    w["separation_ratio"] = w["between"] / w[["within old", "within new"]].mean(axis=1)
    return w.reset_index()


def plot_fission(fine: pd.DataFrame, cluster: pd.DataFrame, mod: pd.DataFrame,
                 out_path: Path) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    A, B, ACC = "#2a78d6", "#eb6834", "#eb6834"
    INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
    GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
        "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.labelcolor": INK2, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False})

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.8))
    fig.subplots_adjust(top=0.71, bottom=0.15, left=0.055, right=0.985, wspace=0.27)

    def dress(ax, t, s, yl):
        ax.set_title(t, loc="left", fontweight="bold", color=INK, pad=26, fontsize=10.5)
        ax.text(0, 1.035, s, transform=ax.transAxes, fontsize=8.2, color=MUTED)
        ax.set_ylabel(yl, fontsize=8.5)
        ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8, length=3)
        ax.axvline(DEPLOYMENT, color=ACC, lw=1.3, ls=(0, (4, 3)), zorder=2)
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=4))
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y-%m"))

    ax = axes[0]
    if len(mod):
        ax.plot(mod.month, mod.split_probability_given_gt7, color=B, lw=2, marker="o",
                ms=4, mec=SURF, mew=0.8, zorder=3)
    ax.set_ylim(0, 0.68)
    ax.text(DEPLOYMENT, 0.655, " collars deployed", color=ACC, fontsize=7.8, va="top")
    dress(ax, "A  Lilac fragments, then re-coheres",
          "pipeline split probability, >7 collars present", "split probability")

    ax = axes[1]
    f = fine.dropna(subset=["separation_ratio"])
    ax.plot(f.month, f["within old"], color=A, lw=1.8, marker="o", ms=4, mec=SURF,
            mew=0.8, zorder=3, label="within original cohort")
    ax.plot(f.month, f["within new"], color=B, lw=1.8, marker="s", ms=4, mec=SURF,
            mew=0.8, zorder=3, label="within new cohort")
    ax.plot(f.month, f["between"], color=INK2, lw=1.6, ls=(0, (3, 2)), zorder=3,
            label="between cohorts")
    ax.legend(loc="upper right", fontsize=8, labelcolor=INK2)
    ax.set_ylim(0, None)
    dress(ax, "B  The two cohorts converge", "within-Lilac 5 m contact rate", "contact rate")

    ax = axes[2]
    ax.plot(f.month, f.separation_ratio, color=B, lw=2.2, marker="o", ms=4.5,
            mec=SURF, mew=0.9, zorder=4, label="5 m contact")
    c = cluster.dropna(subset=["separation_ratio"])
    ax.plot(c.month, c.separation_ratio, color=A, lw=1.8, marker="s", ms=4,
            mec=SURF, mew=0.8, zorder=3, label="hourly cluster (120-900 m)")
    ax.axhline(1.0, color=AXIS, lw=1.0, ls=(0, (2, 2)), zorder=1)
    ax.text(f.month.max(), 1.03, "no separation", fontsize=7.5, color=MUTED, ha="right")
    ax.set_ylim(0, 1.18)
    ax.legend(loc="lower right", fontsize=8, labelcolor=INK2)
    dress(ax, "C  Separation closes at both scales",
          "between-cohort / within-cohort contact", "separation ratio")

    fig.suptitle("Lilac was fragmented when the July-2025 collars went on, and re-cohered",
                 x=0.055, y=0.965, ha="left", va="top", fontsize=13,
                 fontweight="bold", color=INK)
    fig.text(0.055, 0.885,
             "Nine collars deployed 2025-07-31 landed on a Lilac subgroup that was socially distinct from the seven already collared. Over the following\n"
             "year the two cohorts converged at every spatial scale. This is a transient fragmentation, not a permanent fission and not a naming error.",
             ha="left", va="top", fontsize=8.6, color=INK2, linespacing=1.5)
    fig.savefig(out_path, dpi=200, facecolor=SURF)
    print(f"wrote {out_path.name}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-pair-hours", type=int, default=20)
    ap.add_argument("--output-dir", type=Path,
                    default=PROJECT / "outputs" / "lilac_fission_hypothesis_2026-09-03")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)

    m = pd.read_parquet(MEMBERSHIP, columns=["window_start", "animal_id", "origin_group",
                                             "temp_group_id"])
    m["window_start"] = pd.to_datetime(m["window_start"])
    first = m.groupby("animal_id").window_start.min()

    report, monthly_all = {}, []
    for group in ["Lilac", "Copper"]:
        present = set(m[m.origin_group == group].animal_id.unique())
        batch = [i for i in DEPLOY_BATCH[group] if i in present]
        original = [i for i in present
                    if i not in DEPLOY_BATCH[group] and first[i] < PRE_BATCH_CUTOFF]
        ids = sorted(batch + original)
        cohort = pd.Series({i: ("new_Aug2025" if i in batch else "original") for i in ids})
        skipped = len(present) - len(ids)
        if skipped:
            print(f"[{group}] excluding {skipped} animals first seen after "
                  f"{PRE_BATCH_CUTOFF.date()} that are not in the Jul/Aug-2025 batch")
        pairs = co_membership(m, ids)
        if pairs.empty:
            continue
        pairs = classify(pairs, cohort)
        pairs["group"] = group
        monthly_all.append(pairs)

        post = pairs[pairs.month >= DEPLOYMENT]
        pre = pairs[pairs.month < DEPLOYMENT]
        gp, gq = pooled(post, a.min_pair_hours), pooled(pre, a.min_pair_hours)

        n_old = int((cohort == "original").sum())
        n_new = int((cohort == "new_Aug2025").sum())
        print(f"\n{'='*72}\n{group}: {len(ids)} animals ({n_old} original, {n_new} new)\n{'='*72}")
        print("POST-deployment hourly co-membership (same spatial cluster):")
        print(gp.to_string(index=False))
        if len(pre):
            print("\nPRE-deployment (original animals only):")
            print(gq.to_string(index=False))

        rates = gp.set_index("pair_class").co_membership_rate.to_dict()
        wo, wn, bt = rates.get("within original"), rates.get("within new"), rates.get("between cohorts")
        if wo and wn and bt:
            sep = bt / np.mean([wo, wn])
            report[f"{group}_within_original"] = float(wo)
            report[f"{group}_within_new"] = float(wn)
            report[f"{group}_between"] = float(bt)
            report[f"{group}_separation_ratio"] = float(sep)
            print(f"\n  between / mean(within) = {sep:.2f}"
                  f"   ({'SEPARATED' if sep < 0.7 else 'mixed'})")
        report[f"{group}_n_original"] = n_old
        report[f"{group}_n_new"] = n_new

        # persistence: monthly separation ratio
        mm = (post[post.hours >= a.min_pair_hours]
              .groupby(["month", "pair_class"], as_index=False)
              .agg(hours=("hours", "sum"), together=("together", "sum")))
        mm["rate"] = mm.together / mm.hours
        w = mm.pivot(index="month", columns="pair_class", values="rate")
        if {"within original", "within new", "between cohorts"} <= set(w.columns):
            w["separation_ratio"] = w["between cohorts"] / w[["within original", "within new"]].mean(axis=1)
            w.to_csv(a.output_dir / f"{group.lower()}_monthly_separation.csv")
            print(f"\n  monthly separation ratio, {group}:")
            print(w[["within original", "within new", "between cohorts",
                     "separation_ratio"]].round(3).to_string())
            report[f"{group}_months_separated"] = int((w.separation_ratio < 0.7).sum())
            report[f"{group}_months_scored"] = int(w.separation_ratio.notna().sum())

    if monthly_all:
        pd.concat(monthly_all, ignore_index=True).to_csv(
            a.output_dir / "pair_month_co_membership.csv", index=False)

    # context: pipeline's own split evidence, Lilac vs everyone
    if SPLIT_MONTH.exists():
        s = pd.read_csv(SPLIT_MONTH)
        s = s[s.opportunity_hours >= 50]
        rank = (s.groupby("dynamic_social_unit", as_index=False)
                 .agg(months=("month", "nunique"),
                      split_prob=("split_probability_given_classifiable", "mean"))
                 .sort_values("split_prob", ascending=False))
        rank.to_csv(a.output_dir / "split_probability_rank.csv", index=False)
        print(f"\n{'='*72}\nPipeline split probability (mean over months, >=50 opportunity hours)")
        print(rank.head(12).round(3).to_string(index=False))
        if "Lilac" in set(rank.dynamic_social_unit):
            pos = int(rank.reset_index(drop=True).query("dynamic_social_unit=='Lilac'").index[0]) + 1
            report["lilac_split_prob_rank"] = f"{pos} of {len(rank)}"
            report["lilac_split_prob"] = float(
                rank.loc[rank.dynamic_social_unit == "Lilac", "split_prob"].iloc[0])

    if MODULARITY.exists():
        mod = pd.read_csv(MODULARITY, parse_dates=["month"])
        lil = mod[mod.dynamic_social_unit == "Lilac"]
        if len(lil):
            print(f"\nLilac within-split modularity by month:")
            print(lil[["month", "eligible_hours_gt7_collars", "split_probability_given_gt7",
                       "mean_collars_present", "modularity", "n_communities"]]
                  .round(3).to_string(index=False))
            lil.to_csv(a.output_dir / "lilac_within_split_modularity.csv", index=False)
            report["lilac_mean_modularity"] = float(lil.modularity.mean())
            report["lilac_modal_n_communities"] = int(lil.n_communities.mode().iloc[0])
            others = mod[mod.dynamic_social_unit != "Lilac"]
            report["other_groups_mean_modularity"] = float(others.modularity.mean())

    fine = fine_scale_separation()
    fine.to_csv(a.output_dir / "lilac_fine_scale_cohort_separation.csv", index=False)
    fpost = fine[fine.month >= DEPLOYMENT].dropna(subset=["separation_ratio"])
    early = fpost[fpost.month < "2026-01-01"].separation_ratio.mean()
    late = fpost[fpost.month >= "2026-01-01"].separation_ratio.mean()
    report["lilac_5m_separation_aug_dec_2025"] = float(early)
    report["lilac_5m_separation_2026"] = float(late)
    print(f"\n=== within-Lilac 5 m separation between collar cohorts ===")
    print(fine.dropna(subset=["separation_ratio"]).round(4).to_string(index=False))
    print(f"\n  Aug-Dec 2025 {early:.3f}  ->  2026 {late:.3f}   ({late - early:+.3f})")

    report["interpretation"] = (
        "Lilac was internally fragmented when the 2025-07-31 collars were deployed. The nine "
        "new collars landed on a subgroup that was socially distinct at fine scale (5 m "
        "separation ratio ~0.43) and partly separate even at hourly-cluster scale. Over the "
        "following year both scales converged (5 m 0.44 -> 0.65; cluster 0.74 -> 1.00) and the "
        "pipeline's own split probability for Lilac fell from 0.57 back to ~0.02 while collar "
        "count stayed high, so the decay is not a detection effect. This is a TRANSIENT "
        "fragmentation that resolved - not a permanent fission, and not a naming error. "
        "Consequence for the Copper-Lilac result: 'Lilac' was not one social unit for part of "
        "the study, so the within-Lilac reference is not a stable yardstick across it. The "
        "stable-cohort analysis is unaffected because it uses only the original cohort.")

    lil_mod = pd.DataFrame()
    if MODULARITY.exists():
        mod_all = pd.read_csv(MODULARITY, parse_dates=["month"])
        lil_mod = mod_all[mod_all.dynamic_social_unit == "Lilac"].sort_values("month")
    cluster_sep = pd.DataFrame()
    p = a.output_dir / "lilac_monthly_separation.csv"
    if p.exists():
        cluster_sep = pd.read_csv(p, parse_dates=["month"])
    plot_fission(fine, cluster_sep, lil_mod, a.output_dir / "lilac_fission_hypothesis.png")

    (a.output_dir / "lilac_fission_hypothesis.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(f"\nwrote -> {a.output_dir}")


if __name__ == "__main__":
    main()
