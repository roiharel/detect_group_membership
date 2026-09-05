from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\rharel\Documents\group_mebership")
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
# The origin map and the fusion hours MUST come from the same builder run. Until
# 2026-09-03 the origin map was read from `..._local_2h_support` while the fusion
# hours came from `..._shared_full_20260722` - two different runs, differing by
# 3 Copper/Lilac animals. The local_2h_support run also carries corrupted
# animal_id values ("26AB52_????", and "8DOE" where the other run has "8D0E"),
# so the shared_full run is the correct single source. No origin label conflicts
# between the two runs, so this repointing changes cohort coverage, not labels.
CANONICAL_RUN = Path(
    r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_shared_full_20260722"
)
CANONICAL = CANONICAL_RUN / "canonical_hourly_membership.parquet"
SHARED_HISTORY = CANONICAL_RUN / "canonical_5m_shared_history" / "canonical_5m_merge_membership_rows.csv"
PHASES = ROOT / "outputs/copper_lilac_phase_alluvial_from_202403/copper_lilac_data_driven_phases.csv"
OUT = ROOT / "outputs/copper_lilac_effort_corrected_integration"

GROUPS = ("Copper", "Lilac")
PAIR_KEY = "Copper - Lilac"
EARTH_RADIUS_M = 6_371_000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Effort-corrected Copper-Lilac mixing during canonical shared/fusion hours."
    )
    parser.add_argument("--gps", type=Path, default=GPS)
    parser.add_argument("--canonical", type=Path, default=CANONICAL)
    parser.add_argument("--shared-history", type=Path, default=SHARED_HISTORY)
    parser.add_argument("--phases", type=Path, default=PHASES)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    parser.add_argument("--radii", default="1,2,5,10,20,50,100,200,400")
    parser.add_argument("--primary-radius", type=float, default=5.0)
    parser.add_argument(
        "--min-pair-bins",
        type=int,
        default=60,
        help="Minimum simultaneous 2-min bins per pair-month (60 = 2 hours).",
    )
    parser.add_argument("--min-individual-partners", type=int, default=2)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--refresh-gps-cache", action="store_true")
    return parser.parse_args()


def daytime_mask(series: pd.Series) -> pd.Series:
    time = pd.to_datetime(series)
    hour = time.dt.hour + time.dt.minute / 60.0 + time.dt.second / 3600.0
    return hour.between(3.0, 16.0, inclusive="both")


def haversine_vector_m(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    lat1r = np.radians(lat1.astype(float))
    lon1r = np.radians(lon1.astype(float))
    lat2r = np.radians(lat2.astype(float))
    lon2r = np.radians(lon2.astype(float))
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2.0) ** 2
    return 2.0 * EARTH_RADIUS_M * np.arctan2(np.sqrt(a), np.sqrt(np.maximum(0.0, 1.0 - a)))


def load_origin_map(path: Path) -> pd.DataFrame:
    origin = pd.read_parquet(path, columns=["animal_id", "origin_group"])
    origin = origin[origin["origin_group"].isin(GROUPS)].drop_duplicates()
    conflicts = origin.groupby("animal_id", observed=True)["origin_group"].nunique()
    if conflicts.gt(1).any():
        raise ValueError("At least one animal has more than one original group label.")
    return origin


def load_shared_membership(path: Path, origin: pd.DataFrame) -> pd.DataFrame:
    columns = ["timestamp", "association_cluster_id", "pair_key", "animal_id"]
    membership = pd.read_csv(path, usecols=columns, parse_dates=["timestamp"])
    membership = membership[membership["pair_key"].eq(PAIR_KEY)].copy()
    membership = membership.merge(origin, on="animal_id", how="inner", validate="many_to_one")
    membership = membership[daytime_mask(membership["timestamp"])].copy()
    membership = membership.rename(columns={"timestamp": "hour"})
    membership = membership.drop_duplicates(["hour", "association_cluster_id", "animal_id"])
    if membership.empty:
        raise ValueError("No daytime Copper-Lilac shared/fusion membership hours were found.")
    return membership.sort_values(["hour", "association_cluster_id", "animal_id"]).reset_index(drop=True)


def load_gps_positions(
    path: Path,
    membership: pd.DataFrame,
    cache: Path,
    refresh: bool,
) -> pd.DataFrame:
    if cache.exists() and not refresh:
        return pd.read_parquet(cache)

    animals = sorted(membership["animal_id"].unique().tolist())
    start = membership["hour"].min().tz_localize("UTC")
    end = (membership["hour"].max() + pd.Timedelta(hours=1)).tz_localize("UTC")
    gps = pd.read_parquet(
        path,
        columns=["animal_id", "timestamp", "location.long", "location.lat"],
        filters=[("animal_id", "in", animals), ("timestamp", ">=", start), ("timestamp", "<", end)],
    )
    gps["timestamp"] = pd.to_datetime(gps["timestamp"], utc=True).dt.tz_localize(None)
    gps = gps[
        gps["animal_id"].isin(animals)
        & daytime_mask(gps["timestamp"])
        & gps["location.long"].notna()
        & gps["location.lat"].notna()
    ].copy()
    gps["hour"] = gps["timestamp"].dt.floor("h")
    gps["bin_2min"] = gps["timestamp"].dt.floor("2min")
    gps = (
        gps.groupby(["hour", "bin_2min", "animal_id"], observed=True, as_index=False)
        .agg(longitude=("location.long", "median"), latitude=("location.lat", "median"))
    )
    positions = gps.merge(membership, on=["hour", "animal_id"], how="inner", validate="many_to_many")
    positions = positions.drop_duplicates(["hour", "bin_2min", "association_cluster_id", "animal_id"])
    positions["month"] = positions["bin_2min"].dt.to_period("M").dt.to_timestamp()
    cache.parent.mkdir(parents=True, exist_ok=True)
    positions.to_parquet(cache, index=False)
    return positions


def build_pair_month_counts(positions: pd.DataFrame, radii: tuple[float, ...]) -> pd.DataFrame:
    outputs: list[pd.DataFrame] = []
    keys = ["hour", "bin_2min", "association_cluster_id"]
    for month, frame in positions.groupby("month", observed=True, sort=True):
        left = frame[keys + ["animal_id", "origin_group", "longitude", "latitude"]].rename(
            columns={
                "animal_id": "animal_a",
                "origin_group": "origin_a",
                "longitude": "longitude_a",
                "latitude": "latitude_a",
            }
        )
        right = frame[keys + ["animal_id", "origin_group", "longitude", "latitude"]].rename(
            columns={
                "animal_id": "animal_b",
                "origin_group": "origin_b",
                "longitude": "longitude_b",
                "latitude": "latitude_b",
            }
        )
        pairs = left.merge(right, on=keys, how="inner")
        pairs = pairs[pairs["animal_a"].astype(str).lt(pairs["animal_b"].astype(str))].copy()
        if pairs.empty:
            continue
        pairs["distance_m"] = haversine_vector_m(
            pairs["latitude_a"].to_numpy(),
            pairs["longitude_a"].to_numpy(),
            pairs["latitude_b"].to_numpy(),
            pairs["longitude_b"].to_numpy(),
        )
        pairs["pair_type"] = np.where(
            pairs["origin_a"].ne(pairs["origin_b"]),
            "cross_origin",
            np.where(pairs["origin_a"].eq("Copper"), "within_copper", "within_lilac"),
        )
        group_cols = ["animal_a", "origin_a", "animal_b", "origin_b", "pair_type"]
        for radius in radii:
            work = pairs[group_cols].copy()
            work["contact"] = pairs["distance_m"].le(radius).astype(np.int64)
            counts = (
                work.groupby(group_cols, observed=True)
                .agg(opportunity_bins=("contact", "size"), contact_bins=("contact", "sum"))
                .reset_index()
            )
            counts["month"] = pd.Timestamp(month)
            counts["radius_m"] = float(radius)
            outputs.append(counts)
    if not outputs:
        return pd.DataFrame()
    result = pd.concat(outputs, ignore_index=True)
    result["coobserved_hours"] = result["opportunity_bins"] * (2.0 / 60.0)
    result["contact_rate"] = result["contact_bins"] / result["opportunity_bins"]
    return result.sort_values(["radius_m", "month", "pair_type", "animal_a", "animal_b"]).reset_index(drop=True)


def build_individual_month_rates(pair_rates: pd.DataFrame, min_pair_bins: int, min_partners: int) -> pd.DataFrame:
    valid = pair_rates[pair_rates["opportunity_bins"].ge(min_pair_bins)].copy()
    a = valid.rename(columns={"animal_a": "animal_id", "origin_a": "origin_group", "animal_b": "partner_id"})
    b = valid.rename(columns={"animal_b": "animal_id", "origin_b": "origin_group", "animal_a": "partner_id"})
    columns = [
        "month",
        "radius_m",
        "animal_id",
        "origin_group",
        "partner_id",
        "pair_type",
        "opportunity_bins",
        "contact_bins",
        "contact_rate",
    ]
    long = pd.concat([a[columns], b[columns]], ignore_index=True)
    long["relationship"] = np.where(long["pair_type"].eq("cross_origin"), "cross", "within")
    individual = (
        long.groupby(["month", "radius_m", "animal_id", "origin_group", "relationship"], observed=True)
        .agg(
            partners=("partner_id", "nunique"),
            pair_equal_contact_rate=("contact_rate", "mean"),
            opportunity_bins=("opportunity_bins", "sum"),
            contact_bins=("contact_bins", "sum"),
        )
        .reset_index()
    )
    individual["exposure_pooled_contact_rate"] = individual["contact_bins"] / individual["opportunity_bins"]
    individual["coobserved_hours"] = individual["opportunity_bins"] * (2.0 / 60.0)
    individual = individual[individual["partners"].ge(min_partners)].copy()
    return individual.sort_values(["radius_m", "month", "relationship", "origin_group", "animal_id"])


def resample_mean(values: np.ndarray, rng: np.random.Generator) -> float:
    return float(np.mean(rng.choice(values, size=len(values), replace=True)))


def summarize_monthly(
    individual: pd.DataFrame,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for (month, radius), frame in individual.groupby(["month", "radius_m"], observed=True, sort=True):
        subsets = {
            "cross_copper": frame[(frame["relationship"].eq("cross")) & (frame["origin_group"].eq("Copper"))],
            "cross_lilac": frame[(frame["relationship"].eq("cross")) & (frame["origin_group"].eq("Lilac"))],
            "within_copper": frame[(frame["relationship"].eq("within")) & (frame["origin_group"].eq("Copper"))],
            "within_lilac": frame[(frame["relationship"].eq("within")) & (frame["origin_group"].eq("Lilac"))],
        }
        if any(x.empty for x in subsets.values()):
            continue
        values = {k: v["pair_equal_contact_rate"].to_numpy(float) for k, v in subsets.items()}
        observed = {
            "cross_copper": float(np.mean(values["cross_copper"])),
            "cross_lilac": float(np.mean(values["cross_lilac"])),
            "within_copper": float(np.mean(values["within_copper"])),
            "within_lilac": float(np.mean(values["within_lilac"])),
        }
        observed["cross_balanced"] = 0.5 * (observed["cross_copper"] + observed["cross_lilac"])
        observed["within_reference"] = 0.5 * (observed["within_copper"] + observed["within_lilac"])
        observed["integration_ratio"] = observed["cross_balanced"] / observed["within_reference"] if observed["within_reference"] > 0 else np.nan

        draws = {name: [] for name in observed}
        for _ in range(replicates):
            draw = {name: resample_mean(vals, rng) for name, vals in values.items()}
            draw["cross_balanced"] = 0.5 * (draw["cross_copper"] + draw["cross_lilac"])
            draw["within_reference"] = 0.5 * (draw["within_copper"] + draw["within_lilac"])
            draw["integration_ratio"] = draw["cross_balanced"] / draw["within_reference"] if draw["within_reference"] > 0 else np.nan
            for name, value in draw.items():
                draws[name].append(value)

        for metric, estimate in observed.items():
            array = np.asarray(draws[metric], dtype=float)
            rows.append(
                {
                    "month": pd.Timestamp(month),
                    "radius_m": float(radius),
                    "metric": metric,
                    "estimate": estimate,
                    "ci_95_low": float(np.nanquantile(array, 0.025)),
                    "ci_95_high": float(np.nanquantile(array, 0.975)),
                    "n_cross_copper_animals": len(values["cross_copper"]),
                    "n_cross_lilac_animals": len(values["cross_lilac"]),
                    "n_within_copper_animals": len(values["within_copper"]),
                    "n_within_lilac_animals": len(values["within_lilac"]),
                }
            )
    return pd.DataFrame(rows).sort_values(["radius_m", "month", "metric"]).reset_index(drop=True)


def coverage_summary(positions: pd.DataFrame, pair_rates: pd.DataFrame) -> pd.DataFrame:
    animal = (
        positions.groupby(["month", "origin_group"], observed=True)
        .agg(
            tracked_animals=("animal_id", "nunique"),
            animal_2min_bins=("bin_2min", "size"),
            fusion_hours=("hour", "nunique"),
        )
        .reset_index()
    )
    pair = (
        pair_rates.groupby(["month", "radius_m", "pair_type"], observed=True)
        .agg(
            observed_pairs=("contact_rate", "size"),
            pair_opportunity_bins=("opportunity_bins", "sum"),
            pair_contact_bins=("contact_bins", "sum"),
        )
        .reset_index()
    )
    return pair.merge(animal, on="month", how="left")


def add_phase(summary: pd.DataFrame, phase_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    phases = pd.read_csv(phase_path, parse_dates=["start_month", "end_month"])
    out = summary.copy()
    out["phase"] = pd.NA
    for phase in phases.itertuples(index=False):
        mask = out["month"].between(phase.start_month, phase.end_month)
        out.loc[mask, "phase"] = phase.phase
    return out, phases


def plot_primary(summary: pd.DataFrame, phases: pd.DataFrame, radius: float, output: Path) -> None:
    data = summary[summary["radius_m"].eq(radius)].copy()
    full_months = pd.date_range(data["month"].min(), data["month"].max(), freq="MS")
    contact_metrics = [
        ("cross_balanced", "Cross-origin", "#7b3294"),
        ("within_copper", "Within Copper", "#b87333"),
        ("within_lilac", "Within Lilac", "#8e63a9"),
    ]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9), sharex=True, gridspec_kw={"height_ratios": [2.0, 1.25]})
    phase_colors = ["#d9d9d9", "#cfe8e6", "#ecd9c6", "#d9c8e8"]
    for ax in (ax1, ax2):
        for color, phase in zip(phase_colors, phases.itertuples(index=False)):
            ax.axvspan(phase.start_month, phase.end_month + pd.offsets.MonthEnd(1), color=color, alpha=0.18, lw=0)
    for phase in phases.itertuples(index=False):
        midpoint = phase.start_month + (phase.end_month - phase.start_month) / 2
        ax1.text(midpoint, 0.985, str(phase.phase), transform=ax1.get_xaxis_transform(), ha="center", va="top", fontsize=9)

    for metric, label, color in contact_metrics:
        frame = data[data["metric"].eq(metric)].sort_values("month")
        if frame.empty:
            continue
        frame = frame.set_index("month").reindex(full_months).rename_axis("month").reset_index()
        ax1.plot(frame["month"], frame["estimate"] * 100, marker="o", lw=2, ms=4, label=label, color=color)
        ax1.fill_between(
            frame["month"],
            frame["ci_95_low"] * 100,
            frame["ci_95_high"] * 100,
            color=color,
            alpha=0.12,
        )
    ratio = data[data["metric"].eq("integration_ratio")].sort_values("month")
    ratio = ratio.set_index("month").reindex(full_months).rename_axis("month").reset_index()
    ax2.plot(ratio["month"], ratio["estimate"], marker="o", lw=2, ms=4, color="#2f6f9f")
    ax2.fill_between(ratio["month"], ratio["ci_95_low"], ratio["ci_95_high"], color="#2f6f9f", alpha=0.15)
    ax2.axhline(1.0, color="#444444", ls="--", lw=1.2)

    ax1.set_ylabel(f"Pair-equal contact rate at {radius:g} m (%)")
    ax1.set_title("Copper–Lilac mixing during canonical fusion hours\nEffort corrected by simultaneous 2-minute pair observation")
    ax1.legend(frameon=False, ncol=3, loc="upper left")
    ax2.set_ylabel("Cross / within contact-rate ratio")
    ax2.set_xlabel("Calendar month")
    ax2.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax2.tick_params(axis="x", rotation=45)
    for ax in (ax1, ax2):
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_phase_summary(summary: pd.DataFrame) -> pd.DataFrame:
    metrics = ["cross_balanced", "within_reference", "within_copper", "within_lilac", "integration_ratio"]
    data = summary[summary["metric"].isin(metrics) & summary["phase"].notna()].copy()
    return (
        data.groupby(["radius_m", "phase", "metric"], observed=True)["estimate"]
        .mean()
        .unstack("metric")
        .reset_index()
        .sort_values(["radius_m", "phase"])
        .reset_index(drop=True)
    )


def plot_scale_summary(phase_summary: pd.DataFrame, output: Path) -> None:
    phase_order = ["early mixing", "transition", "high merge"]
    styles = {
        "early mixing": ("o", "#5e3c99"),
        "transition": ("s", "#e66101"),
        "high merge": ("^", "#1b9e77"),
    }
    fig, ax = plt.subplots(figsize=(10, 6))
    for phase in phase_order:
        frame = phase_summary[phase_summary["phase"].eq(phase)].sort_values("radius_m")
        if frame.empty:
            continue
        marker, color = styles[phase]
        ax.plot(
            frame["radius_m"], frame["integration_ratio"],
            marker=marker, ms=6, lw=2.2, color=color, label=phase,
        )
    ax.axhline(1.0, color="#444444", ls="--", lw=1.2)
    ax.set_xscale("log")
    radii = sorted(phase_summary["radius_m"].unique())
    ax.set_xticks(radii)
    ax.set_xticklabels([f"{x:g}" for x in radii])
    ax.set_xlabel("Proximity radius (m; log scale)")
    ax.set_ylabel("Mean monthly cross / within contact-rate ratio")
    ax.set_ylim(bottom=0)
    ax.set_title("Copper–Lilac integration depends on spatial scale")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    radii = tuple(float(x.strip()) for x in args.radii.split(",") if x.strip())
    if args.primary_radius not in radii:
        raise ValueError("--primary-radius must be included in --radii")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    origin = load_origin_map(args.canonical)
    membership = load_shared_membership(args.shared_history, origin)
    positions = load_gps_positions(
        args.gps,
        membership,
        args.output_dir / "copper_lilac_fusion_2min_positions.parquet",
        args.refresh_gps_cache,
    )
    pair_rates = build_pair_month_counts(positions, radii)
    if pair_rates.empty:
        raise ValueError("No pairwise co-observation opportunities were generated.")
    individual = build_individual_month_rates(pair_rates, args.min_pair_bins, args.min_individual_partners)
    summary = summarize_monthly(individual, args.bootstrap_replicates, args.seed)
    summary, phases = add_phase(summary, args.phases)
    coverage = coverage_summary(positions, pair_rates)

    pair_rates.to_csv(args.output_dir / "copper_lilac_pair_month_contact_rates.csv", index=False)
    individual.to_csv(args.output_dir / "copper_lilac_individual_month_contact_rates.csv", index=False)
    summary.to_csv(args.output_dir / "copper_lilac_monthly_integration_summary.csv", index=False)
    coverage.to_csv(args.output_dir / "copper_lilac_monthly_tracking_coverage.csv", index=False)
    phase_summary = build_phase_summary(summary)
    phase_summary.to_csv(args.output_dir / "copper_lilac_phase_integration_by_radius.csv", index=False)
    for radius in radii:
        plot_primary(
            summary,
            phases,
            radius,
            args.output_dir / f"copper_lilac_effort_corrected_integration_{radius:g}m.png",
        )
    plot_primary(summary, phases, args.primary_radius, args.output_dir / "copper_lilac_effort_corrected_integration.png")
    plot_scale_summary(phase_summary, args.output_dir / "copper_lilac_integration_scale_summary.png")

    metadata = {
        "pair": PAIR_KEY,
        "groups_are_stable_original_groups": True,
        "fusion_context": "canonical between-group shared/fusion hours for the Copper-Lilac dynamic-unit pair",
        "daytime_window_utc": "03:00-16:00 inclusive",
        "radii_m": radii,
        "primary_radius_m": args.primary_radius,
        "pair_opportunity": "both animals have valid GPS positions in the same 2-minute bin and canonical fusion cluster",
        "pair_contact_rate": "contact bins within the selected radius / pairwise co-observed bins",
        "aggregation": "pair rates averaged within individual, individuals averaged equally within original group, original groups averaged equally",
        "within_reference": "equal mean of within-Copper and within-Lilac individual-average pair rates",
        "integration_ratio": "group-balanced cross-origin rate / within-reference rate; 1 means equal cross- and within-origin contact density",
        "min_pair_bins_per_month": args.min_pair_bins,
        "min_individual_partners_per_month": args.min_individual_partners,
        "bootstrap_unit": "individual, resampled separately within original group and month",
        "bootstrap_replicates": args.bootstrap_replicates,
        "fusion_membership_hours": int(membership["hour"].nunique()),
        "fusion_membership_animals": int(membership["animal_id"].nunique()),
        "gps_position_rows": int(len(positions)),
        "pair_month_rows": int(len(pair_rates)),
        "caution": "Inference is limited to tracked animals and canonical fusion hours; uncollared contacts are not recoverable.",
    }
    (args.output_dir / "copper_lilac_effort_corrected_integration_metadata.json").write_text(
        json.dumps(metadata, indent=2, default=str), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, default=str))
    print("\n5 m integration ratio")
    print(
        summary[
            summary["radius_m"].eq(args.primary_radius) & summary["metric"].eq("integration_ratio")
        ][["month", "estimate", "ci_95_low", "ci_95_high", "phase"]].to_string(index=False)
    )


if __name__ == "__main__":
    main()
