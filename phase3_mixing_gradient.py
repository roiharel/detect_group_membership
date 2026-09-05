"""Phase 3: the general mixing gradient across dyads.

One observation-aware metric applied to every dyad with fine-scale support,
weighted by dyad rather than by bin, with uncertainty clustered on stage-1
encounters. Discrete encounters and sustained associations are estimated as
separate strata so a fused dyad cannot stand in for an encounter.

Two estimators are reported side by side, because their weighting differs:

  edge-weighted deficit   sum(cross_edges)/sum(total_edges) minus the
                          edge-weighted composition expectation
  bin-weighted deficit    mean over bins of (observed - expected), the
                          weighting used by the saved dyad summary

A common-denominator radius comparison is also reported: the fraction of the
5 m eligible bin set that carries at least one cross-group contact at 5 m and
at 2 m. This is not conditioned on a contact existing at the tighter radius.

Run from the project root:
    python phase3_mixing_gradient.py
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(r"C:\Users\rharel\Documents\group_mebership")
UPSTREAM = Path(r"C:\Users\rharel\Documents\New project")
CANON = UPSTREAM / "outputs" / "canonical_robust_hourly_membership_shared_full_20260722"
EVENTS = (PROJECT / "outputs" / "general_structure_2026_09" / "phase2_two_stage_events"
          / "stage1_events_with_stage2_mixing.csv")
FINE = {
    "5m": CANON / "canonical_5m_shared_history_shuffle_expectation" / "canonical_5m_shuffle_expectation_2min_rows.csv",
    "2m": CANON / "canonical_2m_shared_history_shuffle_expectation" / "canonical_5m_shuffle_expectation_2min_rows.csv",
}
OUT_DIR = PROJECT / "outputs" / "general_structure_2026_09" / "phase3_mixing_gradient"

N_BOOTSTRAP = 2000
RNG_SEED = 20260903
MIN_EVENTS_FOR_INTERVAL = 3

GROUP_COLOURS = {
    "Bronze": "#8C6239", "Chartreuse": "#A8C020", "Copper": "#B06B33", "Emerald": "#2E8B57",
    "Jade": "#00A86B", "Lapis": "#26619C", "LapisSplinter": "#6E9BC5", "Lilac": "#9B87C4",
    "Magenta": "#B03A76", "Maroon": "#7B2D3B", "Periwinkle": "#8FA8E8", "Purple": "#6A4C9C",
    "Sapphire": "#1560BD",
}


def load_bins_with_events(radius: str, events: pd.DataFrame) -> pd.DataFrame:
    """Assign every eligible fine-scale bin to its stage-1 encounter."""
    usecols = ["bin_2min", "pair_key", "cross_edges", "total_edges",
               "observed_cross_edge_fraction", "shuffle_expected_cross_edge_fraction"]
    fine = pd.read_csv(FINE[radius], usecols=usecols, parse_dates=["bin_2min"])
    fine["pair_key"] = fine["pair_key"].astype(str)

    windows = events[["stage1_event_id", "pair_key", "encounter_class", "start_hour", "end_hour"]].copy()
    windows["start_hour"] = pd.to_datetime(windows["start_hour"])
    windows["end_hour"] = pd.to_datetime(windows["end_hour"])
    windows["window_end_exclusive"] = windows["end_hour"] + pd.Timedelta(hours=1)

    chunks: list[pd.DataFrame] = []
    for pair, bins in fine.groupby("pair_key", observed=True):
        pair_windows = windows[windows["pair_key"].eq(pair)]
        if pair_windows.empty:
            continue
        matched = pd.merge_asof(
            bins.sort_values("bin_2min"),
            pair_windows.sort_values("start_hour")[
                ["start_hour", "window_end_exclusive", "stage1_event_id", "encounter_class"]
            ],
            left_on="bin_2min",
            right_on="start_hour",
            direction="backward",
        )
        chunks.append(matched[matched["bin_2min"] < matched["window_end_exclusive"]])
    matched = pd.concat(chunks, ignore_index=True)
    matched["bin_deficit"] = (
        matched["observed_cross_edge_fraction"] - matched["shuffle_expected_cross_edge_fraction"]
    )
    matched["has_cross_contact"] = matched["cross_edges"].gt(0)
    return matched


def edge_weighted_deficit(frame: pd.DataFrame) -> float:
    total = frame["total_edges"].sum()
    if total <= 0:
        return float("nan")
    observed = frame["cross_edges"].sum() / total
    expected = float(
        np.average(frame["shuffle_expected_cross_edge_fraction"], weights=frame["total_edges"])
    )
    return observed - expected


def bootstrap_deficit(frame: pd.DataFrame, rng: np.random.Generator) -> tuple[float, float, int]:
    """Percentile interval from resampling stage-1 encounters within the dyad."""
    groups = [g for _, g in frame.groupby("stage1_event_id", observed=True)]
    n = len(groups)
    if n < MIN_EVENTS_FOR_INTERVAL:
        return float("nan"), float("nan"), n
    draws = np.empty(N_BOOTSTRAP)
    for i in range(N_BOOTSTRAP):
        picks = rng.integers(0, n, size=n)
        sample = pd.concat([groups[p] for p in picks], ignore_index=True)
        draws[i] = edge_weighted_deficit(sample)
    finite = draws[np.isfinite(draws)]
    if finite.size == 0:
        return float("nan"), float("nan"), n
    return float(np.percentile(finite, 2.5)), float(np.percentile(finite, 97.5)), n


def summarise(bins: pd.DataFrame, rng: np.random.Generator, label: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for pair, frame in bins.groupby("pair_key", observed=True):
        low, high, n_events = bootstrap_deficit(frame, rng)
        total_edges = float(frame["total_edges"].sum())
        rows.append(
            {
                "pair_key": pair,
                "stratum": label,
                "stage1_events": n_events,
                "eligible_bins": int(frame["bin_2min"].nunique()),
                "contact_positive_bins": int(frame["has_cross_contact"].sum()),
                "contact_positive_bin_fraction": float(frame["has_cross_contact"].mean()),
                "cross_edges": int(frame["cross_edges"].sum()),
                "total_edges": int(total_edges),
                "observed_cross_edge_fraction": (
                    float(frame["cross_edges"].sum() / total_edges) if total_edges > 0 else np.nan
                ),
                "expected_cross_edge_fraction": (
                    float(np.average(frame["shuffle_expected_cross_edge_fraction"],
                                     weights=frame["total_edges"])) if total_edges > 0 else np.nan
                ),
                "edge_weighted_deficit": edge_weighted_deficit(frame),
                "edge_weighted_deficit_low": low,
                "edge_weighted_deficit_high": high,
                "bin_weighted_deficit": float(frame["bin_deficit"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values("edge_weighted_deficit")


def plot_gradient(summary: pd.DataFrame, sustained: pd.DataFrame, out_path: Path) -> None:
    """Ordered per-dyad deficits with support, plus the fused estimate where it exists."""
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D

    frame = summary.sort_values("edge_weighted_deficit").reset_index(drop=True)
    height = 1.1 + 0.42 * len(frame)
    fig, (ax, ax_support) = plt.subplots(
        1, 2, figsize=(12.0, height), gridspec_kw={"width_ratios": [3.0, 1.0]}, sharey=True
    )
    positions = np.arange(len(frame))
    sustained_lookup = dict(zip(sustained["pair_key"], sustained["edge_weighted_deficit"]))

    for pos, row in zip(positions, frame.itertuples(index=False)):
        left, right = row.pair_key.split(" - ")
        cmap = LinearSegmentedColormap.from_list(
            row.pair_key,
            [GROUP_COLOURS.get(left, "#7A7A7A"), GROUP_COLOURS.get(right, "#7A7A7A")],
        )
        value = row.edge_weighted_deficit
        # one bar per dyad, shaded from the first group's colour to the second's
        ax.imshow(
            np.linspace(0.0, 1.0, 256).reshape(1, -1),
            cmap=cmap, aspect="auto", interpolation="bilinear",
            extent=(value, 0.0, pos - 0.31, pos + 0.31), zorder=2,
        )
        if np.isfinite(row.edge_weighted_deficit_low):
            ax.plot(
                [row.edge_weighted_deficit_low, row.edge_weighted_deficit_high],
                [pos, pos], color="#1B1F1D", linewidth=1.5, solid_capstyle="butt", zorder=4,
            )
        fused = sustained_lookup.get(row.pair_key)
        if fused is not None and np.isfinite(fused):
            ax.plot([fused], [pos], marker="D", markersize=7, color="#FFFFFF",
                    markeredgecolor="#1B1F1D", markeredgewidth=1.4, zorder=5)

    ax.set_xlim(min(frame["edge_weighted_deficit"].min(), -0.52) - 0.02, 0.02)
    ax.axvline(0.0, color="#1B1F1D", linewidth=1.0, zorder=3)
    ax.set_yticks(positions)
    ax.set_yticklabels(frame["pair_key"], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Cross-group contact minus composition expectation (5 m, edge-weighted)")
    ax.set_title("Shared space rarely means close mixing", loc="left", fontsize=12, fontweight="bold")
    ax.grid(axis="x", color="#D6DAD6", linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.legend(
        handles=[
            Line2D([0], [0], color="#1B1F1D", linewidth=1.5,
                   label="95% bootstrap interval over encounters"),
            Line2D([0], [0], marker="D", linestyle="none", markersize=7, color="#FFFFFF",
                   markeredgecolor="#1B1F1D", markeredgewidth=1.4,
                   label="same dyad while in sustained association"),
        ],
        loc="lower left", frameon=False, fontsize=8.5,
    )

    ax_support.barh(positions, frame["eligible_bins"], height=0.62, color="#B7C2BC")
    ax_support.set_xscale("log")
    ax_support.set_xlabel("Eligible 2-min bins (log)")
    ax_support.set_title("Support", loc="left", fontsize=10)
    ax_support.grid(axis="x", color="#D6DAD6", linewidth=0.6)
    ax_support.set_axisbelow(True)
    for spine in ["top", "right", "left"]:
        ax_support.spines[spine].set_visible(False)
    for pos, (bins_n, events_n) in enumerate(zip(frame["eligible_bins"], frame["stage1_events"])):
        ax_support.annotate(f"{int(bins_n):,} / {int(events_n)} ev", (bins_n, pos), va="center",
                            fontsize=7.5, color="#3B423D", xytext=(4, 0), textcoords="offset points")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    events = pd.read_csv(EVENTS, parse_dates=["start_hour", "end_hour"])

    print("assigning fine-scale bins to encounters ...")
    bins_5m = load_bins_with_events("5m", events)
    bins_2m = load_bins_with_events("2m", events)
    print(f"  5m bins={len(bins_5m):,}  2m bins={len(bins_2m):,}")

    discrete = bins_5m[bins_5m["encounter_class"].eq("discrete_encounter")]
    sustained = bins_5m[bins_5m["encounter_class"].eq("sustained_association")]

    print("bootstrapping per-dyad deficits ...")
    summary_discrete = summarise(discrete, rng, "discrete_encounter")
    summary_sustained = summarise(sustained, rng, "sustained_association")
    summary_all = summarise(bins_5m, rng, "all_encounters")
    summary = pd.concat([summary_discrete, summary_sustained, summary_all], ignore_index=True)
    summary.to_csv(OUT_DIR / "mixing_gradient_by_dyad.csv", index=False)

    # -- common-denominator radius comparison ---------------------------------
    contact_2m = (
        bins_2m.loc[bins_2m["has_cross_contact"], ["pair_key", "bin_2min"]]
        .drop_duplicates()
        .assign(cross_contact_2m=True)
    )
    nested = bins_5m.merge(contact_2m, on=["pair_key", "bin_2min"], how="left")
    nested["cross_contact_2m"] = nested["cross_contact_2m"].fillna(False).astype(bool)
    radius_table = (
        nested.groupby(["pair_key", "encounter_class"], observed=True)
        .agg(
            eligible_bins_5m=("bin_2min", "nunique"),
            contact_rate_5m=("has_cross_contact", "mean"),
            contact_rate_2m=("cross_contact_2m", "mean"),
        )
        .reset_index()
    )
    radius_table["ratio_2m_to_5m"] = np.where(
        radius_table["contact_rate_5m"] > 0,
        radius_table["contact_rate_2m"] / radius_table["contact_rate_5m"],
        np.nan,
    )
    radius_table.sort_values(["encounter_class", "contact_rate_5m"], ascending=[True, False]).to_csv(
        OUT_DIR / "radius_nested_contact_rates.csv", index=False
    )

    # -- pooled estimates: dyad-weighted versus bin-weighted ------------------
    def pooled(frame_bins: pd.DataFrame, per_dyad: pd.DataFrame, label: str) -> dict[str, object]:
        values = per_dyad["edge_weighted_deficit"].dropna().to_numpy()
        draws = np.array([
            np.mean(rng.choice(values, size=values.size, replace=True)) for _ in range(N_BOOTSTRAP)
        ]) if values.size > 1 else np.array([np.nan])
        return {
            "stratum": label,
            "dyads": int(values.size),
            "dyad_weighted_mean_deficit": float(np.mean(values)) if values.size else np.nan,
            "dyad_weighted_low": float(np.percentile(draws, 2.5)) if values.size > 1 else np.nan,
            "dyad_weighted_high": float(np.percentile(draws, 97.5)) if values.size > 1 else np.nan,
            "dyads_below_expectation": int((values < 0).sum()),
            "bin_weighted_pooled_deficit": float(frame_bins["bin_deficit"].mean()),
            "edge_weighted_pooled_deficit": edge_weighted_deficit(frame_bins),
            "largest_dyad_share_of_bins": float(
                frame_bins.groupby("pair_key", observed=True)["bin_2min"].nunique().max()
                / frame_bins["bin_2min"].nunique()
            ),
            "largest_dyad_by_bins": str(
                frame_bins.groupby("pair_key", observed=True)["bin_2min"].nunique().idxmax()
            ),
        }

    pooled_rows = [
        pooled(discrete, summary_discrete, "discrete_encounter"),
        pooled(sustained, summary_sustained, "sustained_association"),
        pooled(bins_5m, summary_all, "all_encounters"),
    ]
    pooled_frame = pd.DataFrame(pooled_rows)
    pooled_frame.to_csv(OUT_DIR / "mixing_gradient_pooled.csv", index=False)

    # -- figure: discrete-encounter stratum, with the fused estimate overlaid --
    figure_path = OUT_DIR / "mixing_gradient_5m.png"
    plot_gradient(summary_discrete, summary_sustained, figure_path)

    metadata = {
        "phase": "3 - general mixing gradient",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "event_source": str(EVENTS),
        "finescale_sources": {k: str(v) for k, v in FINE.items()},
        "reference": "composition-preserving shuffle expectation carried in the fine-scale product",
        "primary_estimator": "edge-weighted deficit: sum(cross_edges)/sum(total_edges) minus the "
                             "edge-weighted expectation",
        "uncertainty": f"percentile bootstrap over stage-1 encounters within a dyad, "
                       f"{N_BOOTSTRAP} draws, seed {RNG_SEED}; no interval below "
                       f"{MIN_EVENTS_FOR_INTERVAL} encounters",
        "pooling": "dyad-weighted; the bin-weighted pooled value is reported alongside so the "
                   "difference in weighting is visible",
        "strata": {
            "discrete_encounter": int(summary_discrete["pair_key"].nunique()),
            "sustained_association": int(summary_sustained["pair_key"].nunique()),
        },
        "pooled": pooled_rows,
        "dyads_with_finescale_support": int(summary_all["pair_key"].nunique()),
        "structural_dyads_without_finescale_support": int(
            events["pair_key"].nunique() - summary_all["pair_key"].nunique()
        ),
        "outputs": {
            "by_dyad": str(OUT_DIR / "mixing_gradient_by_dyad.csv"),
            "pooled": str(OUT_DIR / "mixing_gradient_pooled.csv"),
            "radius_nested": str(OUT_DIR / "radius_nested_contact_rates.csv"),
            "figure": str(figure_path),
        },
    }
    (OUT_DIR / "phase3_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    show = ["pair_key", "stage1_events", "eligible_bins", "observed_cross_edge_fraction",
            "expected_cross_edge_fraction", "edge_weighted_deficit",
            "edge_weighted_deficit_low", "edge_weighted_deficit_high", "bin_weighted_deficit"]
    print()
    print("DISCRETE ENCOUNTERS")
    print(summary_discrete[show].round(4).to_string(index=False))
    print()
    print("SUSTAINED ASSOCIATIONS")
    print(summary_sustained[show].round(4).to_string(index=False))
    print()
    print(pooled_frame.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
