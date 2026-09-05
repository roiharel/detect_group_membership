from __future__ import annotations

import argparse
import html
import json
import warnings
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.genmod.families import Binomial, Gaussian
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.tools.sm_exceptions import DomainWarning


BASE = Path(r"C:\Users\rharel\Documents\group_mebership")
DEFAULT_DAILY_DIR = BASE / "outputs" / "dynamic_social_unit_merge_gamm" / "daily_interaction_gamm"
DEFAULT_RAW_INTERACTIONS = (
    BASE
    / "outputs"
    / "dynamic_social_unit_merge_gamm"
    / "group_merge_mixing_dynamics"
    / "bigmerge_dynamic_social_unit_2min_vs_hourly_5m_no_copper_lilac_2min_metric_rows.csv"
)
DEFAULT_OUT_DIR = BASE / "outputs" / "dynamic_social_unit_merge_gamm" / "daily_interaction_hurdle"
DEFAULT_PROXIMITY_STATUS = (
    BASE / "outputs" / "dynamic_social_unit_merge_gamm" / "proximity_status_dynamic_social_unit.parquet"
)
DEFAULT_VEDBA_DIR = Path(
    r"\\10.126.19.90\EAS_shared\baboon\working\data\processed\2025\acc\vedba"
)

EARTH_RADIUS_M = 6_371_000.0
SPLINE_DF_DISTANCE = 5
SPLINE_DF_DISTANCE_SIMPLE = 4
SPLINE_DF_NDVI = 4
MIN_MODEL_ROWS = 20
MIN_GROUP_CONTROL_ROWS = 5
DAY_START_HOUR = 3
DAY_END_HOUR = 16
EVENT_MAX_GAP_HOURS = 14.0
EVENT_BIN_WIDTH_MINUTES = 2.0
BAYESIAN_DRAWS = 4000
BAYESIAN_PRIOR_SHAPE = 2.0
BAYESIAN_PRIOR_SCALE = 1.0


def normalize_group_name(value: object) -> str:
    return "".join(str(value).split()).lower()


def sorted_pair(group_a: object, group_b: object) -> tuple[str, str, str]:
    a = str(group_a)
    b = str(group_b)
    left, right = sorted([a, b])
    return left, right, f"{left} - {right}"


def in_daytime_window(series: pd.Series) -> pd.Series:
    times = pd.to_datetime(series)
    hours = times.dt.hour + times.dt.minute / 60.0 + times.dt.second / 3600.0
    return hours.between(DAY_START_HOUR, DAY_END_HOUR, inclusive="both")


def week_start(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series).dt.to_period("W-SUN").dt.start_time


def zscore(series: pd.Series) -> pd.Series:
    sd = series.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(0.0, index=series.index)
    return (series - series.mean()) / sd


def haversine_m(lat1: pd.Series, lon1: pd.Series, lat2: pd.Series, lon2: pd.Series) -> pd.Series:
    lat1_rad = np.radians(lat1.astype(float))
    lon1_rad = np.radians(lon1.astype(float))
    lat2_rad = np.radians(lat2.astype(float))
    lon2_rad = np.radians(lon2.astype(float))
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(a))


def load_candidate_dyad_days(daily_dir: Path) -> pd.DataFrame:
    centroids = pd.read_csv(daily_dir / "weekly_group_centroids.csv", parse_dates=["period_start"])
    centroids = centroids.dropna(
        subset=["period_start", "dynamic_social_unit", "centroid_latitude", "centroid_longitude"]
    ).copy()
    rows: list[dict[str, object]] = []
    for period_start, group in centroids.groupby("period_start", observed=True):
        group = group.sort_values("dynamic_social_unit")
        for (_, left), (_, right) in combinations(group.iterrows(), 2):
            group_a, group_b, pair_key = sorted_pair(left["dynamic_social_unit"], right["dynamic_social_unit"])
            if group_a == str(left["dynamic_social_unit"]):
                a = left
                b = right
            else:
                a = right
                b = left
            rows.append(
                {
                    "period_start": period_start,
                    "group_a": group_a,
                    "group_b": group_b,
                    "pair_key": pair_key,
                    "group_a_key": normalize_group_name(group_a),
                    "group_b_key": normalize_group_name(group_b),
                    "group_a_centroid_latitude": a["centroid_latitude"],
                    "group_a_centroid_longitude": a["centroid_longitude"],
                    "group_b_centroid_latitude": b["centroid_latitude"],
                    "group_b_centroid_longitude": b["centroid_longitude"],
                    "group_a_position_rows": a.get("n_position_rows", np.nan),
                    "group_b_position_rows": b.get("n_position_rows", np.nan),
                    "group_a_animals": a.get("n_animals", np.nan),
                    "group_b_animals": b.get("n_animals", np.nan),
                }
            )
    candidates = pd.DataFrame(rows)
    candidates["daily_centroid_distance_m"] = haversine_m(
        candidates["group_a_centroid_latitude"],
        candidates["group_a_centroid_longitude"],
        candidates["group_b_centroid_latitude"],
        candidates["group_b_centroid_longitude"],
    )
    pair_range = (
        candidates.groupby("pair_key", observed=True)["daily_centroid_distance_m"]
        .agg(
            range_mean_centroid_distance_m="mean",
            range_median_centroid_distance_m="median",
            range_overlap_days="size",
        )
        .reset_index()
    )
    candidates = candidates.merge(pair_range, on="pair_key", how="left")
    return candidates


# --- demographic group size -------------------------------------------------
# `group_a_animals` / `group_b_animals` are COLLARED ANIMAL COUNTS, not group
# sizes. Until 2026-09-03 their sum was named `group_size_total` and plotted as
# "Total group size", which made collar coverage look like a demographic effect.
# Independent demographic sizes now exist, sourced from EAS
# (metadata/GS_collars_demographics.xlsx) via build_authoritative_group_names.py.
# Coverage ranges 3%-42% across groups, so the two quantities are far from
# interchangeable. Both are retained: demographic size as the size covariate,
# collar count as an observation covariate.
DEMOGRAPHICS_PATH = (
    Path(__file__).resolve().parent
    / "outputs" / "authoritative_group_names_2026-09-03" / "group_demographics.csv"
)
# Which size variable enters the models. "collars" reproduces every run before
# 2026-09-03; "demographic" uses the independent group sizes. Set from --size-source.
SIZE_SOURCE = "collars"


def load_group_demographics(path: Path) -> pd.DataFrame | None:
    """group_key -> demographic group_size and collar count, or None if absent."""
    if not path.exists():
        print(f"NOTE: demographic sizes not found at {path}; "
              f"falling back to collar counts. Run build_authoritative_group_names.py "
              f"with EAS reachable to create it.")
        return None
    d = pd.read_csv(path)
    d["group_key"] = d["group_id"].map(normalize_group_name)
    d = d.drop_duplicates("group_key").set_index("group_key")
    return d[["group_size", "no_collars", "collar_coverage_percent"]]


def add_demographic_size_columns(rows: pd.DataFrame, demo: pd.DataFrame | None) -> pd.DataFrame:
    """Attach demographic sizes alongside the existing collar-count columns."""
    rows["collar_count_total"] = rows["group_size_total"]
    rows["log1p_collar_count_total_z"] = rows["log1p_group_size_total_z"]
    if demo is None:
        rows["demographic_group_size_total"] = np.nan
        return rows
    for side in ("a", "b"):
        key = rows[f"group_{side}"].map(normalize_group_name)
        rows[f"group_{side}_demographic_size"] = key.map(demo["group_size"])
    a = pd.to_numeric(rows["group_a_demographic_size"], errors="coerce")
    b = pd.to_numeric(rows["group_b_demographic_size"], errors="coerce")
    rows["demographic_group_size_total"] = a + b
    rows["demographic_group_size_abs_diff"] = (a - b).abs()
    rows["log1p_demographic_group_size_total"] = np.log1p(rows["demographic_group_size_total"])
    rows["log1p_demographic_group_size_abs_diff"] = np.log1p(rows["demographic_group_size_abs_diff"])
    rows["log1p_demographic_group_size_total_z"] = zscore(rows["log1p_demographic_group_size_total"])
    rows["log1p_demographic_group_size_abs_diff_z"] = zscore(rows["log1p_demographic_group_size_abs_diff"])
    missing = rows["demographic_group_size_total"].isna().mean()
    if missing:
        print(f"NOTE: {missing:.1%} of rows lack a demographic size for one or both groups.")
    return rows


def size_terms(source: str) -> tuple[str, str]:
    """Column names for (total size, size difference) under the chosen source."""
    if source == "demographic":
        return ("log1p_demographic_group_size_total_z",
                "log1p_demographic_group_size_abs_diff_z")
    return "log1p_group_size_total_z", "log1p_group_size_abs_diff_z"


def load_daytime_candidate_dyad_days(proximity_path: Path) -> pd.DataFrame:
    columns = [
        "animal_id",
        "timestamp",
        "median_latitude",
        "median_longitude",
        "dynamic_social_unit",
    ]
    data = pd.read_parquet(proximity_path, columns=columns)
    data["timestamp"] = pd.to_datetime(data["timestamp"])
    data = data[
        in_daytime_window(data["timestamp"])
        & data["dynamic_social_unit"].notna()
        & data["median_latitude"].notna()
        & data["median_longitude"].notna()
    ].copy()
    data["period_start"] = data["timestamp"].dt.floor("D")
    data["dynamic_social_unit"] = data["dynamic_social_unit"].astype(str)
    centroids = (
        data.groupby(["period_start", "dynamic_social_unit"], observed=True)
        .agg(
            centroid_latitude=("median_latitude", "mean"),
            centroid_longitude=("median_longitude", "mean"),
            n_position_rows=("animal_id", "size"),
            n_animals=("animal_id", "nunique"),
            first_timestamp=("timestamp", "min"),
            last_timestamp=("timestamp", "max"),
        )
        .reset_index()
    )
    centroids["group_key"] = centroids["dynamic_social_unit"].map(normalize_group_name)

    rows: list[dict[str, object]] = []
    for period_start, group in centroids.groupby("period_start", observed=True):
        group = group.sort_values("dynamic_social_unit")
        for (_, left), (_, right) in combinations(group.iterrows(), 2):
            group_a, group_b, pair_key = sorted_pair(left["dynamic_social_unit"], right["dynamic_social_unit"])
            if group_a == str(left["dynamic_social_unit"]):
                a = left
                b = right
            else:
                a = right
                b = left
            rows.append(
                {
                    "period_start": period_start,
                    "group_a": group_a,
                    "group_b": group_b,
                    "pair_key": pair_key,
                    "group_a_key": normalize_group_name(group_a),
                    "group_b_key": normalize_group_name(group_b),
                    "group_a_centroid_latitude": a["centroid_latitude"],
                    "group_a_centroid_longitude": a["centroid_longitude"],
                    "group_b_centroid_latitude": b["centroid_latitude"],
                    "group_b_centroid_longitude": b["centroid_longitude"],
                    "group_a_position_rows": a["n_position_rows"],
                    "group_b_position_rows": b["n_position_rows"],
                    "group_a_animals": a["n_animals"],
                    "group_b_animals": b["n_animals"],
                    "group_a_first_daytime_timestamp": a["first_timestamp"],
                    "group_a_last_daytime_timestamp": a["last_timestamp"],
                    "group_b_first_daytime_timestamp": b["first_timestamp"],
                    "group_b_last_daytime_timestamp": b["last_timestamp"],
                }
            )
    candidates = pd.DataFrame(rows)
    candidates["daily_centroid_distance_m"] = haversine_m(
        candidates["group_a_centroid_latitude"],
        candidates["group_a_centroid_longitude"],
        candidates["group_b_centroid_latitude"],
        candidates["group_b_centroid_longitude"],
    )
    pair_range = (
        candidates.groupby("pair_key", observed=True)["daily_centroid_distance_m"]
        .agg(
            range_mean_centroid_distance_m="mean",
            range_median_centroid_distance_m="median",
            range_overlap_days="size",
        )
        .reset_index()
    )
    candidates = candidates.merge(pair_range, on="pair_key", how="left")
    return candidates


def load_daily_ndvi(daily_dir: Path) -> pd.DataFrame:
    path = daily_dir / "weekly_group_ndvi.csv"
    if not path.exists():
        return pd.DataFrame()
    ndvi = pd.read_csv(path, parse_dates=["period_start", "weekly_ndvi_first_date", "weekly_ndvi_last_date"])
    return ndvi


def add_daily_ndvi(rows: pd.DataFrame, ndvi: pd.DataFrame) -> pd.DataFrame:
    if ndvi.empty:
        return rows.copy()
    left = ndvi.rename(
        columns={
            "group_key": "group_a_key",
            "weekly_ndvi_mean": "group_a_daily_ndvi_mean",
            "weekly_ndvi_median": "group_a_daily_ndvi_median",
            "weekly_ndvi_n": "group_a_daily_ndvi_n",
            "weekly_ndvi_first_date": "group_a_daily_ndvi_first_date",
            "weekly_ndvi_last_date": "group_a_daily_ndvi_last_date",
        }
    )
    right = ndvi.rename(
        columns={
            "group_key": "group_b_key",
            "weekly_ndvi_mean": "group_b_daily_ndvi_mean",
            "weekly_ndvi_median": "group_b_daily_ndvi_median",
            "weekly_ndvi_n": "group_b_daily_ndvi_n",
            "weekly_ndvi_first_date": "group_b_daily_ndvi_first_date",
            "weekly_ndvi_last_date": "group_b_daily_ndvi_last_date",
        }
    )
    out = rows.merge(left, on=["period_start", "group_a_key"], how="left")
    out = out.merge(right, on=["period_start", "group_b_key"], how="left")
    out["dyad_daily_ndvi_mean"] = out[["group_a_daily_ndvi_mean", "group_b_daily_ndvi_mean"]].mean(axis=1)
    out["dyad_daily_ndvi_abs_diff"] = (
        out["group_a_daily_ndvi_mean"] - out["group_b_daily_ndvi_mean"]
    ).abs()
    return out


def load_positive_interactions(raw_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(raw_path, parse_dates=["timestamp", "bin_2min"])
    raw = raw.dropna(subset=["bin_2min", "group_a", "group_b", "pair_key"]).copy()
    raw = raw[in_daytime_window(raw["bin_2min"])].copy()
    raw["period_start"] = raw["bin_2min"].dt.floor("D")
    pair_parts = raw.apply(lambda row: sorted_pair(row["group_a"], row["group_b"]), axis=1)
    raw["group_a"] = [parts[0] for parts in pair_parts]
    raw["group_b"] = [parts[1] for parts in pair_parts]
    raw["pair_key"] = [parts[2] for parts in pair_parts]
    raw["cross_edges"] = pd.to_numeric(raw["cross_edges"], errors="coerce").fillna(0)
    raw["total_edges"] = pd.to_numeric(raw["total_edges"], errors="coerce").fillna(0)
    raw["has_cross_edge"] = raw["cross_edges"].gt(0)

    positive_bins = (
        raw.loc[raw["has_cross_edge"], ["period_start", "pair_key", "bin_2min"]]
        .drop_duplicates()
        .groupby(["period_start", "pair_key"], observed=True)
        .agg(positive_2min_bins=("bin_2min", "size"), first_positive_bin=("bin_2min", "min"), last_positive_bin=("bin_2min", "max"))
        .reset_index()
    )
    exposure = (
        raw.groupby(["period_start", "pair_key"], observed=True)
        .agg(
            eligible_2min_bins=("bin_2min", "nunique"),
            raw_metric_rows=("bin_2min", "size"),
            total_cross_edges=("cross_edges", "sum"),
            total_candidate_edges=("total_edges", "sum"),
            mean_pair_n=("pair_n", "mean"),
            mean_cluster_size_total=("cluster_size_total", "mean"),
            mean_group_a_fraction_present=("group_a_fraction_present", "mean"),
            mean_group_b_fraction_present=("group_b_fraction_present", "mean"),
            mean_min_group_fraction_present=("min_group_fraction_present", "mean"),
            mean_composition_entropy_norm=("composition_entropy_norm", "mean"),
            mean_pair_balance=("pair_balance", "mean"),
        )
        .reset_index()
    )
    daily = exposure.merge(positive_bins, on=["period_start", "pair_key"], how="left")
    daily["positive_2min_bins"] = daily["positive_2min_bins"].fillna(0).astype(int)
    daily["positive_duration_hours"] = daily["positive_2min_bins"] * (2.0 / 60.0)
    daily["any_interaction"] = daily["positive_2min_bins"].gt(0).astype(int)
    daily["edge_interaction_probability"] = np.where(
        daily["total_candidate_edges"].gt(0),
        daily["total_cross_edges"] / daily["total_candidate_edges"],
        np.nan,
    )

    episode_rows = build_positive_episodes(raw)
    return daily, episode_rows


def build_positive_episodes(raw: pd.DataFrame, max_gap_hours: float = EVENT_MAX_GAP_HOURS) -> pd.DataFrame:
    bin_cols = [
        "pair_key",
        "group_a",
        "group_b",
        "bin_2min",
        "merge_episode_id",
        "cluster_n_groups",
        "cluster_size_total",
        "pair_n",
        "cluster_groups",
        "cross_edges",
        "total_edges",
        "cross_edge_fraction",
        "edge_modularity_q",
        "composition_entropy_norm",
        "pair_balance",
    ]
    bins = raw.loc[raw["has_cross_edge"], [col for col in bin_cols if col in raw.columns]].drop_duplicates()
    if bins.empty:
        return pd.DataFrame()
    chunks = []
    episode_index = 0
    for pair_key, group in bins.sort_values(["pair_key", "bin_2min"]).groupby("pair_key", observed=True):
        group = group.sort_values("bin_2min").copy()
        gap_hours = group["bin_2min"].diff().dt.total_seconds() / 3600.0
        local_episode = (gap_hours.isna() | gap_hours.gt(max_gap_hours)).cumsum()
        for _, episode in group.groupby(local_episode, sort=True):
            episode_index += 1
            start = episode["bin_2min"].min()
            end = episode["bin_2min"].max()
            end_exclusive = end + pd.Timedelta(minutes=EVENT_BIN_WIDTH_MINUTES)
            active_contact_hours = float(episode["bin_2min"].nunique() * (EVENT_BIN_WIDTH_MINUTES / 60.0))
            max_cluster_n_groups = (
                float(pd.to_numeric(episode.get("cluster_n_groups"), errors="coerce").max())
                if "cluster_n_groups" in episode
                else np.nan
            )
            total_cross_edges = (
                float(pd.to_numeric(episode.get("cross_edges"), errors="coerce").fillna(0).sum())
                if "cross_edges" in episode
                else np.nan
            )
            total_edges = (
                float(pd.to_numeric(episode.get("total_edges"), errors="coerce").fillna(0).sum())
                if "total_edges" in episode
                else np.nan
            )
            integration_5m_fraction = (
                total_cross_edges / total_edges
                if np.isfinite(total_cross_edges) and np.isfinite(total_edges) and total_edges > 0
                else np.nan
            )
            merge_size_class = "large_merge" if max_cluster_n_groups >= 3 else "small_merge"
            chunks.append(
                {
                    "interaction_episode_id": f"INTERACTION_{episode_index:05d}",
                    "pair_key": pair_key,
                    "group_a": episode["group_a"].iloc[0],
                    "group_b": episode["group_b"].iloc[0],
                    "episode_start": start,
                    "episode_end": end,
                    "episode_end_exclusive": end_exclusive,
                    "episode_start_day": start.floor("D"),
                    "episode_end_day": end.floor("D"),
                    "positive_2min_bins": int(episode["bin_2min"].nunique()),
                    "active_contact_hours": active_contact_hours,
                    "total_cross_group_5m_edges": total_cross_edges,
                    "total_candidate_5m_edges": total_edges,
                    "integration_5m_fraction": integration_5m_fraction,
                    "duration_hours": float((end_exclusive - start).total_seconds() / 3600.0),
                    "duration_days": float((end_exclusive - start).total_seconds() / 86_400.0),
                    "max_gap_allowed_hours": max_gap_hours,
                    "merge_size_class": merge_size_class,
                    "max_cluster_n_groups": max_cluster_n_groups,
                    "median_cluster_n_groups": float(pd.to_numeric(episode.get("cluster_n_groups"), errors="coerce").median())
                    if "cluster_n_groups" in episode
                    else np.nan,
                    "median_cluster_size_total": float(
                        pd.to_numeric(episode.get("cluster_size_total"), errors="coerce").median()
                    )
                    if "cluster_size_total" in episode
                    else np.nan,
                    "median_pair_n": float(pd.to_numeric(episode.get("pair_n"), errors="coerce").median())
                    if "pair_n" in episode
                    else np.nan,
                    "mean_cross_edge_fraction": float(
                        pd.to_numeric(episode.get("cross_edge_fraction"), errors="coerce").mean()
                    )
                    if "cross_edge_fraction" in episode
                    else np.nan,
                    "mean_edge_modularity_q": float(
                        pd.to_numeric(episode.get("edge_modularity_q"), errors="coerce").mean()
                    )
                    if "edge_modularity_q" in episode
                    else np.nan,
                    "mean_composition_entropy_norm": float(
                        pd.to_numeric(episode.get("composition_entropy_norm"), errors="coerce").mean()
                    )
                    if "composition_entropy_norm" in episode
                    else np.nan,
                    "mean_pair_balance": float(pd.to_numeric(episode.get("pair_balance"), errors="coerce").mean())
                    if "pair_balance" in episode
                    else np.nan,
                    "merge_episode_ids": " | ".join(sorted(set(map(str, episode.get("merge_episode_id", []))))),
                    "episode_cluster_groups": " | ".join(sorted(set(map(str, episode.get("cluster_groups", []))))),
                }
            )
    return pd.DataFrame(chunks)


def build_model_rows(daily_dir: Path, raw_path: Path, proximity_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = load_daytime_candidate_dyad_days(proximity_path)
    ndvi = load_daily_ndvi(daily_dir)
    interaction_daily, episodes = load_positive_interactions(raw_path)
    rows = candidates.merge(interaction_daily, on=["period_start", "pair_key"], how="left")
    fill_zero = [
        "eligible_2min_bins",
        "raw_metric_rows",
        "total_cross_edges",
        "total_candidate_edges",
        "positive_2min_bins",
        "positive_duration_hours",
        "any_interaction",
    ]
    for col in fill_zero:
        rows[col] = rows[col].fillna(0)
    rows["any_interaction"] = rows["any_interaction"].astype(int)
    rows = rows.sort_values(["pair_key", "period_start"]).copy()
    rows["prior_pair_interaction_days"] = (
        rows.groupby("pair_key", observed=True)["any_interaction"]
        .cumsum()
        - rows["any_interaction"]
    )
    rows["prior_pair_interaction_duration_hours"] = (
        rows.groupby("pair_key", observed=True)["positive_duration_hours"]
        .cumsum()
        - rows["positive_duration_hours"]
    )
    rows["integration_5m_fraction"] = np.where(
        rows["total_candidate_edges"].gt(0),
        rows["total_cross_edges"] / rows["total_candidate_edges"],
        np.nan,
    )
    rows = add_daily_ndvi(rows, ndvi)
    rows["log1p_range_mean_centroid_distance_m"] = np.log1p(
        rows["range_mean_centroid_distance_m"].clip(lower=0)
    )
    rows["log1p_range_mean_centroid_distance_m_z"] = zscore(
        rows["log1p_range_mean_centroid_distance_m"]
    )
    rows["log1p_positive_duration_hours"] = np.log1p(rows["positive_duration_hours"])
    rows["group_size_total"] = (
        pd.to_numeric(rows["group_a_animals"], errors="coerce")
        + pd.to_numeric(rows["group_b_animals"], errors="coerce")
    )
    rows["log1p_group_size_total"] = np.log1p(rows["group_size_total"])
    rows["group_size_abs_diff"] = (
        pd.to_numeric(rows["group_a_animals"], errors="coerce")
        - pd.to_numeric(rows["group_b_animals"], errors="coerce")
    ).abs()
    rows["log1p_group_size_abs_diff"] = np.log1p(rows["group_size_abs_diff"])
    rows["log1p_group_size_total_z"] = zscore(rows["log1p_group_size_total"])
    rows["log1p_group_size_abs_diff_z"] = zscore(rows["log1p_group_size_abs_diff"])
    rows = add_demographic_size_columns(rows, load_group_demographics(DEMOGRAPHICS_PATH))
    rows["log1p_prior_pair_interaction_days"] = np.log1p(rows["prior_pair_interaction_days"])
    rows["log1p_prior_pair_interaction_days_z"] = zscore(rows["log1p_prior_pair_interaction_days"])
    rows["log1p_prior_pair_interaction_duration_hours"] = np.log1p(
        rows["prior_pair_interaction_duration_hours"]
    )
    rows["log1p_prior_pair_interaction_duration_hours_z"] = zscore(
        rows["log1p_prior_pair_interaction_duration_hours"]
    )
    if "dyad_daily_ndvi_mean" in rows:
        rows["dyad_daily_ndvi_mean_z"] = zscore(rows["dyad_daily_ndvi_mean"])
        rows["dyad_daily_ndvi_abs_diff_z"] = zscore(rows["dyad_daily_ndvi_abs_diff"])
    rows["group_a_control"] = collapse_sparse_levels(rows, "group_a")
    rows["group_b_control"] = collapse_sparse_levels(rows, "group_b")

    if not episodes.empty:
        episode_predictors = rows[
            [
                "period_start",
                "pair_key",
                "daily_centroid_distance_m",
                "range_mean_centroid_distance_m",
                "range_median_centroid_distance_m",
                "range_overlap_days",
                "dyad_daily_ndvi_mean",
                "dyad_daily_ndvi_abs_diff",
                "group_a_animals",
                "group_b_animals",
            ]
        ].copy()
        episodes = episodes.merge(
            episode_predictors,
            left_on=["episode_start_day", "pair_key"],
            right_on=["period_start", "pair_key"],
            how="left",
        )
    return rows, episodes


def aggregate_to_weekly_rows(daily_rows: pd.DataFrame) -> pd.DataFrame:
    rows = daily_rows.copy()
    rows["week_start"] = week_start(rows["period_start"])
    agg = (
        rows.groupby(["week_start", "pair_key"], observed=True)
        .agg(
            group_a=("group_a", "first"),
            group_b=("group_b", "first"),
            group_a_key=("group_a_key", "first"),
            group_b_key=("group_b_key", "first"),
            daily_centroid_distance_m=("daily_centroid_distance_m", "mean"),
            range_mean_centroid_distance_m=("range_mean_centroid_distance_m", "first"),
            range_median_centroid_distance_m=("range_median_centroid_distance_m", "first"),
            range_overlap_days=("range_overlap_days", "first"),
            group_a_animals=("group_a_animals", "mean"),
            group_b_animals=("group_b_animals", "mean"),
            any_interaction=("any_interaction", "max"),
            positive_duration_hours=("positive_duration_hours", "sum"),
            positive_2min_bins=("positive_2min_bins", "sum"),
            eligible_2min_bins=("eligible_2min_bins", "sum"),
            total_cross_edges=("total_cross_edges", "sum"),
            total_candidate_edges=("total_candidate_edges", "sum"),
            raw_metric_rows=("raw_metric_rows", "sum"),
            dyad_daily_ndvi_mean=("dyad_daily_ndvi_mean", "mean"),
            dyad_daily_ndvi_abs_diff=("dyad_daily_ndvi_abs_diff", "mean"),
            group_a_control=("group_a_control", "first"),
            group_b_control=("group_b_control", "first"),
            n_observed_days=("period_start", "nunique"),
        )
        .reset_index()
        .rename(columns={"week_start": "period_start"})
    )
    agg = agg.sort_values(["pair_key", "period_start"]).copy()
    agg["prior_pair_interaction_days"] = (
        agg.groupby("pair_key", observed=True)["any_interaction"].cumsum()
        - agg["any_interaction"]
    )
    agg["prior_pair_interaction_duration_hours"] = (
        agg.groupby("pair_key", observed=True)["positive_duration_hours"].cumsum()
        - agg["positive_duration_hours"]
    )
    agg["integration_5m_fraction"] = np.where(
        agg["total_candidate_edges"].gt(0),
        agg["total_cross_edges"] / agg["total_candidate_edges"],
        np.nan,
    )
    agg["edge_interaction_probability"] = agg["integration_5m_fraction"]
    agg["log1p_range_mean_centroid_distance_m"] = np.log1p(
        agg["range_mean_centroid_distance_m"].clip(lower=0)
    )
    agg["log1p_range_mean_centroid_distance_m_z"] = zscore(
        agg["log1p_range_mean_centroid_distance_m"]
    )
    agg["log1p_positive_duration_hours"] = np.log1p(agg["positive_duration_hours"])
    agg["group_size_total"] = (
        pd.to_numeric(agg["group_a_animals"], errors="coerce")
        + pd.to_numeric(agg["group_b_animals"], errors="coerce")
    )
    agg["log1p_group_size_total"] = np.log1p(agg["group_size_total"])
    agg["group_size_abs_diff"] = (
        pd.to_numeric(agg["group_a_animals"], errors="coerce")
        - pd.to_numeric(agg["group_b_animals"], errors="coerce")
    ).abs()
    agg["log1p_group_size_abs_diff"] = np.log1p(agg["group_size_abs_diff"])
    agg["log1p_group_size_total_z"] = zscore(agg["log1p_group_size_total"])
    agg["log1p_group_size_abs_diff_z"] = zscore(agg["log1p_group_size_abs_diff"])
    agg["log1p_prior_pair_interaction_days"] = np.log1p(agg["prior_pair_interaction_days"])
    agg["log1p_prior_pair_interaction_days_z"] = zscore(agg["log1p_prior_pair_interaction_days"])
    agg["log1p_prior_pair_interaction_duration_hours"] = np.log1p(
        agg["prior_pair_interaction_duration_hours"]
    )
    agg["log1p_prior_pair_interaction_duration_hours_z"] = zscore(
        agg["log1p_prior_pair_interaction_duration_hours"]
    )
    agg["dyad_daily_ndvi_mean_z"] = zscore(agg["dyad_daily_ndvi_mean"])
    agg["dyad_daily_ndvi_abs_diff_z"] = zscore(agg["dyad_daily_ndvi_abs_diff"])
    return agg


def collapse_sparse_levels(rows: pd.DataFrame, col: str) -> pd.Series:
    counts = rows[col].value_counts()
    keep = set(counts[counts >= MIN_GROUP_CONTROL_ROWS].index)
    suffix = "a" if col.endswith("_a") else "b"
    return rows[col].where(rows[col].isin(keep), f"other_group_{suffix}")


def build_formula(response: str, include_ndvi: bool) -> str:
    distance_term = (
        "log1p_range_mean_centroid_distance_m_z"
        if response in {"any_interaction", "integration_5m_fraction"}
        else f"bs(log1p_range_mean_centroid_distance_m_z, df={SPLINE_DF_DISTANCE_SIMPLE}, degree=3, include_intercept=False)"
    )
    size_total, size_diff = size_terms(SIZE_SOURCE)
    terms = [
        distance_term,
        size_total,
        size_diff,
        "log1p_prior_pair_interaction_days_z",
    ]
    if SIZE_SOURCE == "demographic":
        # keep collar coverage in the model as an observation covariate, so a
        # demographic size effect is not just re-expressing sampling effort
        terms.append("log1p_collar_count_total_z")
    if include_ndvi:
        terms.extend(["dyad_daily_ndvi_mean_z", "dyad_daily_ndvi_abs_diff_z"])
    return f"{response} ~ " + " + ".join(terms)


def fit_gee(rows: pd.DataFrame, response: str, include_ndvi: bool, family: object) -> tuple[object, str, str]:
    formula = build_formula(response, include_ndvi)
    size_total, size_diff = size_terms(SIZE_SOURCE)
    cols = [
        "pair_key",
        response,
        "log1p_range_mean_centroid_distance_m_z",
        size_total,
        size_diff,
        "log1p_prior_pair_interaction_days_z",
    ]
    if SIZE_SOURCE == "demographic":
        cols.append("log1p_collar_count_total_z")
    if response == "integration_5m_fraction":
        cols += ["total_candidate_edges"]
    if include_ndvi:
        cols += ["dyad_daily_ndvi_mean_z", "dyad_daily_ndvi_abs_diff_z"]
    model_df = rows[cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if response == "integration_5m_fraction":
        model_df = model_df[model_df["total_candidate_edges"].gt(0)].copy()
    if len(model_df) < MIN_MODEL_ROWS or model_df["pair_key"].nunique() < 2:
        raise ValueError(f"Not enough rows to fit {response}.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DomainWarning)
        if response in {"any_interaction", "integration_5m_fraction"}:
            fit_kwargs = {
                "maxiter": 200,
                "cov_type": "cluster",
                "cov_kwds": {"groups": model_df["pair_key"]},
            }
            if response == "integration_5m_fraction":
                model = sm.GLM.from_formula(
                    formula,
                    data=model_df,
                    family=family,
                    freq_weights=model_df["total_candidate_edges"],
                )
            else:
                model = sm.GLM.from_formula(formula, data=model_df, family=family)
            result = model.fit(
                maxiter=200,
                cov_type="cluster",
                cov_kwds={"groups": model_df["pair_key"]},
            )
            if not np.isfinite(result.params).all():
                raise ValueError(f"{response} GLM returned non-finite coefficients.")
            note = (
                "Binomial GLM with dyad-clustered robust standard errors"
                if response == "any_interaction"
                else "Binomial edge-fraction GLM with total 5 m candidate edges as frequency weights and dyad-clustered robust standard errors"
            )
            return result, formula, note
        try:
            model = GEE.from_formula(
                formula,
                groups="pair_key",
                cov_struct=Exchangeable(),
                data=model_df,
                family=family,
            )
            result = model.fit(maxiter=200)
            return result, formula, "GEE with dyad-clustered exchangeable working correlation"
        except Exception as exc:
            glm = sm.GLM.from_formula(formula, data=model_df, family=family)
            result = glm.fit(maxiter=200)
            return result, formula, f"GLM fallback after GEE error: {exc}"


def coefficient_table(result: object, model_name: str, response: str) -> pd.DataFrame:
    params = result.params
    conf = result.conf_int()
    out = pd.DataFrame(
        {
            "model": model_name,
            "response": response,
            "term": params.index,
            "estimate": params.values,
            "std_error": result.bse.reindex(params.index).values,
            "ci_low": conf.loc[params.index, 0].values,
            "ci_high": conf.loc[params.index, 1].values,
            "p_value": result.pvalues.reindex(params.index).values,
        }
    )
    if hasattr(result, "prior_mean"):
        out["prior_mean"] = result.prior_mean.reindex(params.index).values
    if hasattr(result, "prior_sd"):
        out["prior_sd"] = result.prior_sd.reindex(params.index).values
    return out


# Labels must say what the variable IS. Terms built from `group_*_animals` are
# counts of COLLARED animals, so they are labelled "Collared animals", not
# "Total group size" - the previous labels asserted a demographic effect the
# data could not support. Demographic terms carry the group-size label.
EVENT_EFFECT_LABELS = {
    "log1p_centroid_distance_median_1wk_before_after_m_z": "Centroid distance",
    "log1p_prior_meetings_z": "Prior meetings",
    "log1p_group_size_total_1wk_before_after_z": "Collared animals (both groups)",
    "log1p_group_size_abs_diff_1wk_before_after_z": "Collar count difference",
    "dyad_ndvi_mean_1wk_before_after_z": "NDVI",
    "log1p_event_gps_movement_mean_m_per_h_z": "GPS movement",
}

MEETING_EFFECT_LABELS = {
    "log1p_range_mean_centroid_distance_m_z": "Centroid distance",
    "log1p_prior_pair_interaction_days_z": "Prior meetings",
    "log1p_group_size_total_z": "Collared animals (both groups)",
    "log1p_group_size_abs_diff_z": "Collar count difference",
    "log1p_collar_count_total_z": "Collared animals (both groups)",
    "log1p_demographic_group_size_total_z": "Total group size (demographic)",
    "log1p_demographic_group_size_abs_diff_z": "Group size difference (demographic)",
    "dyad_daily_ndvi_mean_z": "NDVI",
}

EVENT_DISTANCE_WINDOWS = {
    "before_after": {
        "term": "log1p_centroid_distance_median_1wk_before_after_m_z",
        "raw": "centroid_distance_median_1wk_before_after_m",
        "log": "log1p_centroid_distance_median_1wk_before_after_m",
        "label": "Median centroid distance in +/- 1 week event window (m)",
    },
    "before": {
        "term": "log1p_centroid_distance_median_1wk_before_m_z",
        "raw": "centroid_distance_median_1wk_before_m",
        "log": "log1p_centroid_distance_median_1wk_before_m",
        "label": "Median centroid distance in 1 week before event (m)",
    },
    "after": {
        "term": "log1p_centroid_distance_median_1wk_after_m_z",
        "raw": "centroid_distance_median_1wk_after_m",
        "log": "log1p_centroid_distance_median_1wk_after_m",
        "label": "Median centroid distance in 1 week after event (m)",
    },
}


def fixed_effects_for_plot(coef: pd.DataFrame, label_map: dict[str, str]) -> pd.DataFrame:
    fixed = coef[coef["term"].isin(label_map)].copy()
    fixed["effect_label"] = fixed["term"].map(label_map)
    fixed["abs_estimate"] = fixed["estimate"].abs()
    return fixed


def plot_event_effect_sizes(
    coef: pd.DataFrame,
    out_dir: Path,
    filename: str,
    title: str,
    xlabel: str,
    color: str,
    label_map: dict[str, str] | None = None,
) -> Path:
    label_map = EVENT_EFFECT_LABELS if label_map is None else label_map
    fixed = fixed_effects_for_plot(coef, label_map)
    order = list(label_map.values())
    fixed["effect_label"] = pd.Categorical(fixed["effect_label"], categories=order, ordered=True)
    fixed = fixed.sort_values("effect_label")
    fig_height = max(3.6, 0.52 * len(fixed) + 1.4)
    fig, ax = plt.subplots(figsize=(8.6, fig_height))
    y = np.arange(len(fixed))
    ax.axvline(0, color="#333333", linewidth=1.0, alpha=0.65)
    ax.errorbar(
        fixed["estimate"],
        y,
        xerr=[
            fixed["estimate"] - fixed["ci_low"],
            fixed["ci_high"] - fixed["estimate"],
        ],
        fmt="o",
        color=color,
        ecolor=color,
        elinewidth=2.0,
        capsize=4,
        markersize=7,
        alpha=0.92,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(fixed["effect_label"])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(True, axis="x", alpha=0.22)
    fig.tight_layout()
    path = out_dir / filename
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_event_effect_size_comparison(
    duration_coef: pd.DataFrame,
    integration_coef: pd.DataFrame,
    out_dir: Path,
) -> Path:
    duration = fixed_effects_for_plot(duration_coef, EVENT_EFFECT_LABELS)
    integration = fixed_effects_for_plot(integration_coef, EVENT_EFFECT_LABELS)
    combined = pd.concat(
        [
            duration.assign(response_label="Duration"),
            integration.assign(response_label="5 m integration"),
        ],
        ignore_index=True,
    )
    order = list(EVENT_EFFECT_LABELS.values())
    combined["effect_label"] = pd.Categorical(combined["effect_label"], categories=order, ordered=True)
    combined = combined.sort_values(["effect_label", "response_label"])
    y_base = np.arange(len(order))
    offsets = {"Duration": -0.14, "5 m integration": 0.14}
    colors = {"Duration": "#4c78a8", "5 m integration": "#59a14f"}
    fig, ax = plt.subplots(figsize=(9.2, max(3.8, 0.58 * len(order) + 1.6)))
    ax.axvline(0, color="#333333", linewidth=1.0, alpha=0.65)
    for response_label, sub in combined.groupby("response_label", observed=True):
        positions = np.array([order.index(label) for label in sub["effect_label"].astype(str)]) + offsets[response_label]
        ax.errorbar(
            sub["estimate"],
            positions,
            xerr=[
                sub["estimate"] - sub["ci_low"],
                sub["ci_high"] - sub["estimate"],
            ],
            fmt="o",
            color=colors[response_label],
            ecolor=colors[response_label],
            elinewidth=2.0,
            capsize=4,
            markersize=7,
            alpha=0.92,
            label=response_label,
        )
    ax.set_yticks(y_base)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Standardized coefficient on model link scale")
    ax.set_title("Bayesian event model effect-size comparison")
    ax.grid(True, axis="x", alpha=0.22)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    path = out_dir / "event_bayesian_effect_size_comparison.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_three_part_effect_size_comparison(
    meeting_coef: pd.DataFrame,
    duration_coef: pd.DataFrame,
    integration_coef: pd.DataFrame,
    out_dir: Path,
) -> Path:
    meeting = fixed_effects_for_plot(meeting_coef, MEETING_EFFECT_LABELS)
    duration = fixed_effects_for_plot(duration_coef, EVENT_EFFECT_LABELS)
    integration = fixed_effects_for_plot(integration_coef, EVENT_EFFECT_LABELS)
    combined = pd.concat(
        [
            meeting.assign(response_label="Meeting probability"),
            duration.assign(response_label="Duration"),
            integration.assign(response_label="5 m integration"),
        ],
        ignore_index=True,
    )
    order = list(EVENT_EFFECT_LABELS.values())
    combined["effect_label"] = pd.Categorical(combined["effect_label"], categories=order, ordered=True)
    combined = combined.sort_values(["effect_label", "response_label"])
    y_base = np.arange(len(order))
    offsets = {"Meeting probability": -0.22, "Duration": 0.0, "5 m integration": 0.22}
    colors = {"Meeting probability": "#f28e2b", "Duration": "#4c78a8", "5 m integration": "#59a14f"}
    fig, ax = plt.subplots(figsize=(9.8, max(4.0, 0.62 * len(order) + 1.7)))
    ax.axvline(0, color="#333333", linewidth=1.0, alpha=0.65)
    for response_label, sub in combined.groupby("response_label", observed=True):
        positions = np.array([order.index(label) for label in sub["effect_label"].astype(str)]) + offsets[response_label]
        ax.errorbar(
            sub["estimate"],
            positions,
            xerr=[
                sub["estimate"] - sub["ci_low"],
                sub["ci_high"] - sub["estimate"],
            ],
            fmt="o",
            color=colors[response_label],
            ecolor=colors[response_label],
            elinewidth=2.0,
            capsize=4,
            markersize=7,
            alpha=0.92,
            label=response_label,
        )
    ax.set_yticks(y_base)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Standardized coefficient on model link scale")
    ax.set_title("Hierarchical three-part model effect-size comparison")
    ax.grid(True, axis="x", alpha=0.22)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    path = out_dir / "hierarchical_three_part_effect_size_comparison.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def summary_row(rows: pd.DataFrame, model_name: str, response: str, formula: str, note: str, result: object) -> dict[str, object]:
    model_rows = rows.dropna(subset=[response, "log1p_range_mean_centroid_distance_m_z"]).copy()
    if response == "log1p_positive_duration_hours":
        model_rows = model_rows[model_rows["any_interaction"].eq(1)].copy()
    if response == "integration_5m_fraction":
        model_rows = model_rows[model_rows["total_candidate_edges"].gt(0)].copy()
    return {
        "model": model_name,
        "response": response,
        "fit_note": note,
        "formula": formula,
        "n_dyad_days": int(len(model_rows)),
        "n_dyads": int(model_rows["pair_key"].nunique()),
        "n_positive_dyad_days": int(rows["any_interaction"].sum()),
        "overall_daily_interaction_probability": float(rows["any_interaction"].mean()),
        "total_positive_duration_hours": float(rows["positive_duration_hours"].sum()),
        "median_positive_duration_hours": float(rows.loc[rows["any_interaction"].eq(1), "positive_duration_hours"].median()),
        "overall_5m_integration_fraction": float(
            rows["total_cross_edges"].sum() / rows["total_candidate_edges"].sum()
        ),
        "n_eligible_5m_dyad_days": int(rows["total_candidate_edges"].gt(0).sum()),
        "range_mean_distance_min_m": float(model_rows["range_mean_centroid_distance_m"].min()),
        "range_mean_distance_median_m": float(model_rows["range_mean_centroid_distance_m"].median()),
        "range_mean_distance_max_m": float(model_rows["range_mean_centroid_distance_m"].max()),
        "qic": float(getattr(result, "qic", lambda: [np.nan])()[0]) if hasattr(result, "qic") else np.nan,
    }


def prediction_grid(
    rows: pd.DataFrame,
    result: object,
    model_name: str,
    response: str,
    include_ndvi: bool,
    predictor: str,
) -> pd.DataFrame:
    positive = rows[rows["any_interaction"].eq(1)]
    reference = positive if response == "log1p_positive_duration_hours" and not positive.empty else rows
    reference_cols = [
        "range_mean_centroid_distance_m",
        "log1p_range_mean_centroid_distance_m_z",
        "prior_pair_interaction_days",
        "log1p_group_size_total_z",
        "log1p_group_size_abs_diff_z",
        "log1p_prior_pair_interaction_days_z",
        "group_a",
        "group_b",
    ]
    if include_ndvi:
        reference_cols += ["dyad_daily_ndvi_mean_z", "dyad_daily_ndvi_abs_diff_z"]
    reference = reference.replace([np.inf, -np.inf], np.nan).dropna(subset=reference_cols).copy()
    n_grid = 140
    dist_mean = rows["log1p_range_mean_centroid_distance_m"].mean()
    dist_sd = rows["log1p_range_mean_centroid_distance_m"].std(ddof=0)
    if not np.isfinite(dist_sd) or dist_sd == 0:
        dist_sd = 1.0
    prior_mean = rows["log1p_prior_pair_interaction_days"].mean()
    prior_sd = rows["log1p_prior_pair_interaction_days"].std(ddof=0)
    if not np.isfinite(prior_sd) or prior_sd == 0:
        prior_sd = 1.0
    group_a = reference["group_a"].mode().iloc[0]
    group_b = reference["group_b"].mode().iloc[0]

    reference_distance = float(reference["range_mean_centroid_distance_m"].median())
    reference_prior = float(reference["prior_pair_interaction_days"].median())
    grid = pd.DataFrame(
        {
            "range_mean_centroid_distance_m": np.repeat(reference_distance, n_grid),
            "log1p_range_mean_centroid_distance_m_z": np.repeat(
                (np.log1p(reference_distance) - dist_mean) / dist_sd,
                n_grid,
            ),
            "prior_pair_interaction_days": np.repeat(reference_prior, n_grid),
            "log1p_prior_pair_interaction_days_z": np.repeat(
                (np.log1p(reference_prior) - prior_mean) / prior_sd,
                n_grid,
            ),
            "log1p_group_size_total_z": 0.0,
            "log1p_group_size_abs_diff_z": 0.0,
            "group_a": group_a,
            "group_b": group_b,
            "group_a_control": collapse_value(rows, "group_a", group_a),
            "group_b_control": collapse_value(rows, "group_b", group_b),
            "pair_key": f"{group_a} - {group_b}",
        }
    )
    predictor_labels = {
        "distance": "Average daytime centroid distance (m)",
        "ndvi_mean": "Dyad mean NDVI",
        "ndvi_abs_diff": "Dyad absolute NDVI difference",
        "history": "Prior positive weeks",
    }
    if include_ndvi:
        for col in ["dyad_daily_ndvi_mean", "dyad_daily_ndvi_abs_diff"]:
            mean = rows[col].mean()
            sd = rows[col].std(ddof=0)
            if not np.isfinite(sd) or sd == 0:
                sd = 1.0
            grid[f"{col}_z"] = 0.0
            grid[col] = mean

    if predictor == "distance":
        values = np.linspace(
            reference["range_mean_centroid_distance_m"].quantile(0.01),
            reference["range_mean_centroid_distance_m"].quantile(0.99),
            n_grid,
        )
        grid["range_mean_centroid_distance_m"] = values
        grid["log1p_range_mean_centroid_distance_m_z"] = (np.log1p(values) - dist_mean) / dist_sd
    elif predictor == "history":
        upper = max(float(reference["prior_pair_interaction_days"].quantile(0.99)), 1.0)
        values = np.linspace(0.0, upper, n_grid)
        grid["prior_pair_interaction_days"] = values
        grid["log1p_prior_pair_interaction_days_z"] = (np.log1p(values) - prior_mean) / prior_sd
    elif predictor in {"ndvi_mean", "ndvi_abs_diff"}:
        if not include_ndvi:
            return pd.DataFrame()
        col = "dyad_daily_ndvi_mean" if predictor == "ndvi_mean" else "dyad_daily_ndvi_abs_diff"
        values = np.linspace(reference[col].quantile(0.01), reference[col].quantile(0.99), n_grid)
        mean = rows[col].mean()
        sd = rows[col].std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        grid[col] = values
        grid[f"{col}_z"] = (values - mean) / sd
    else:
        raise ValueError(f"Unknown predictor grid: {predictor}")
    grid["predictor"] = predictor
    grid["predictor_value"] = values
    grid["predictor_label"] = predictor_labels[predictor]

    design = patsy.build_design_matrices([result.model.data.design_info], grid, return_type="dataframe")[0]
    pred_linear = np.asarray(design @ result.params)
    cov = result.cov_params()
    pred_var = np.einsum("ij,jk,ik->i", design.to_numpy(), cov, design.to_numpy())
    pred_se = np.sqrt(np.clip(pred_var, 0, np.inf))
    if response in {"any_interaction", "integration_5m_fraction"}:
        pred = 1 / (1 + np.exp(-pred_linear))
        low = 1 / (1 + np.exp(-(pred_linear - norm.ppf(0.975) * pred_se)))
        high = 1 / (1 + np.exp(-(pred_linear + norm.ppf(0.975) * pred_se)))
        scale = "probability" if response == "any_interaction" else "5m_integration_fraction"
    else:
        pred = np.expm1(pred_linear)
        low = np.expm1(pred_linear - norm.ppf(0.975) * pred_se)
        high = np.expm1(pred_linear + norm.ppf(0.975) * pred_se)
        scale = "hours"
    out = grid[
        [
            "predictor",
            "predictor_value",
            "predictor_label",
            "range_mean_centroid_distance_m",
            "prior_pair_interaction_days",
            "log1p_group_size_abs_diff_z",
            "group_a",
            "group_b",
            "pair_key",
        ]
    ].copy()
    if include_ndvi:
        out["dyad_daily_ndvi_mean"] = grid["dyad_daily_ndvi_mean"]
        out["dyad_daily_ndvi_abs_diff"] = grid["dyad_daily_ndvi_abs_diff"]
    out["predicted"] = pred
    out["ci_low"] = low
    out["ci_high"] = high
    out["model"] = model_name
    out["response"] = response
    out["scale"] = scale
    return out


def collapse_value(rows: pd.DataFrame, col: str, value: str) -> str:
    counts = rows[col].value_counts()
    suffix = "a" if col.endswith("_a") else "b"
    return value if counts.get(value, 0) >= MIN_GROUP_CONTROL_ROWS else f"other_group_{suffix}"


def plot_predictions(rows: pd.DataFrame, predictions: pd.DataFrame, out_dir: Path) -> list[Path]:
    paths = []
    response_specs = [
        ("any_interaction", "Weekly probability of any 5 m cross-group interaction", "probability"),
        ("log1p_positive_duration_hours", "Positive interaction duration (hours/week)", "duration"),
        ("integration_5m_fraction", "5 m integration: cross-group edge fraction", "5m_integration"),
    ]
    predictor_specs = [
        ("distance", "range_mean_centroid_distance_m", "Average daytime centroid distance for group pair (m, log scale)", True),
        ("ndvi_mean", "dyad_daily_ndvi_mean", "Dyad mean NDVI", False),
        ("ndvi_abs_diff", "dyad_daily_ndvi_abs_diff", "Absolute difference in group NDVI", False),
        ("history", "prior_pair_interaction_days", "Prior positive weeks", False),
    ]
    colors = {
        "any_interaction": "#4c78a8",
        "log1p_positive_duration_hours": "#59a14f",
        "integration_5m_fraction": "#af7aa1",
    }
    rng = np.random.default_rng(42)
    for response, ylabel, response_slug in response_specs:
        for predictor, x_col, xlabel, use_log_x in predictor_specs:
            if x_col not in rows.columns:
                continue
            pred = predictions[
                predictions["response"].eq(response)
                & predictions["predictor"].eq(predictor)
            ]
            if pred.empty:
                continue
            filename = f"daily_hurdle_{response_slug}_{predictor}.png"
            if predictor == "distance":
                legacy_names = {
                    "any_interaction": "daily_hurdle_probability_distance.png",
                    "log1p_positive_duration_hours": "daily_hurdle_duration_distance.png",
                    "integration_5m_fraction": "daily_hurdle_5m_integration_distance.png",
                }
                filename = legacy_names[response]

            fig, ax = plt.subplots(figsize=(9.5, 5.8))
            if response == "any_interaction":
                plot_df = rows.dropna(subset=[x_col]).copy()
                y = plot_df["any_interaction"]
                ax.scatter(
                    plot_df[x_col],
                    y + rng.normal(0, 0.015, size=len(plot_df)),
                    s=12,
                    alpha=0.18,
                    color=colors[response],
                    edgecolor="none",
                )
                ax.set_ylim(-0.08, 1.08)
            elif response == "log1p_positive_duration_hours":
                plot_df = rows[rows["any_interaction"].eq(1)].dropna(subset=[x_col]).copy()
                ax.scatter(
                    plot_df[x_col],
                    plot_df["positive_duration_hours"],
                    s=22,
                    alpha=0.32,
                    color=colors[response],
                    edgecolor="none",
                )
            else:
                plot_df = rows[rows["total_candidate_edges"].gt(0)].dropna(subset=[x_col]).copy()
                sizes = np.clip(np.sqrt(plot_df["total_candidate_edges"]) * 2.2, 16, 95)
                ax.scatter(
                    plot_df[x_col],
                    plot_df["integration_5m_fraction"],
                    s=sizes,
                    alpha=0.35,
                    color=colors[response],
                    edgecolor="none",
                )
                ax.set_ylim(-0.04, 1.04)
            for model_name, sub in pred.groupby("model", observed=True):
                sub = sub.sort_values("predictor_value")
                ax.plot(sub["predictor_value"], sub["predicted"], linewidth=2.2, label=model_name)
                ax.fill_between(sub["predictor_value"], sub["ci_low"], sub["ci_high"], alpha=0.16)
            if use_log_x:
                ax.set_xscale("log")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend(frameon=False)
            ax.grid(True, alpha=0.22)
            fig.tight_layout()
            path = out_dir / filename
            fig.savefig(path, dpi=220)
            plt.close(fig)
            paths.append(path)
    return paths


def build_event_duration_rows(episodes: pd.DataFrame, daily_rows: pd.DataFrame) -> pd.DataFrame:
    if episodes.empty:
        return pd.DataFrame()
    events = episodes.copy()
    for col in ["episode_start", "episode_end", "episode_end_exclusive", "episode_start_day", "episode_end_day"]:
        if col in events.columns:
            events[col] = pd.to_datetime(events[col])
    days = daily_rows.copy()
    days["period_start"] = pd.to_datetime(days["period_start"]).dt.floor("D")

    events = events.sort_values(["pair_key", "episode_start"]).copy()
    events["prior_meetings"] = events.groupby("pair_key", observed=True).cumcount()
    events["prior_meeting_duration_hours"] = (
        events.groupby("pair_key", observed=True)["duration_hours"].cumsum()
        - events["duration_hours"]
    )

    context_rows = []
    for idx, event in events.iterrows():
        start_day = pd.to_datetime(event["episode_start_day"]).floor("D")
        end_day = pd.to_datetime(event["episode_end_day"]).floor("D")
        before_start = start_day - pd.Timedelta(days=7)
        after_end = end_day + pd.Timedelta(days=7)
        pair_days = days[
            days["pair_key"].eq(event["pair_key"])
            & days["period_start"].between(before_start, after_end, inclusive="both")
        ].copy()
        before_days = pair_days[pair_days["period_start"].lt(start_day)]
        after_days = pair_days[pair_days["period_start"].gt(end_day)]
        context_rows.append(
            {
                "_event_index": idx,
                "centroid_distance_median_1wk_before_after_m": pair_days["daily_centroid_distance_m"].median(),
                "centroid_distance_median_1wk_before_m": before_days["daily_centroid_distance_m"].median(),
                "centroid_distance_median_1wk_after_m": after_days["daily_centroid_distance_m"].median(),
                "centroid_distance_n_days_1wk_before_after": int(pair_days["daily_centroid_distance_m"].notna().sum()),
                "dyad_ndvi_mean_1wk_before_after": pair_days["dyad_daily_ndvi_mean"].mean(),
                "dyad_ndvi_abs_diff_1wk_before_after": pair_days["dyad_daily_ndvi_abs_diff"].mean(),
                "group_a_animals_1wk_before_after": pair_days["group_a_animals"].mean(),
                "group_b_animals_1wk_before_after": pair_days["group_b_animals"].mean(),
                "group_size_total_1wk_before_after": (
                    pair_days["group_a_animals"].mean() + pair_days["group_b_animals"].mean()
                ),
                "group_size_abs_diff_1wk_before_after": (
                    pair_days["group_a_animals"] - pair_days["group_b_animals"]
                ).abs().mean(),
            }
        )
    context = pd.DataFrame(context_rows).set_index("_event_index")
    events = events.join(context, how="left")
    events["log1p_duration_hours"] = np.log1p(events["duration_hours"])
    for raw_col in [
        "centroid_distance_median_1wk_before_after_m",
        "centroid_distance_median_1wk_before_m",
        "centroid_distance_median_1wk_after_m",
    ]:
        events[f"log1p_{raw_col}"] = np.log1p(events[raw_col].clip(lower=0))
    events["log1p_prior_meetings"] = np.log1p(events["prior_meetings"])
    events["log1p_prior_meeting_duration_hours"] = np.log1p(events["prior_meeting_duration_hours"])
    events["log1p_group_size_total_1wk_before_after"] = np.log1p(
        pd.to_numeric(events["group_size_total_1wk_before_after"], errors="coerce")
    )
    events["log1p_group_size_abs_diff_1wk_before_after"] = np.log1p(
        pd.to_numeric(events["group_size_abs_diff_1wk_before_after"], errors="coerce")
    )
    if "integration_5m_fraction" in events.columns:
        integration = pd.to_numeric(events["integration_5m_fraction"], errors="coerce").clip(1e-4, 1 - 1e-4)
        events["logit_integration_5m_fraction"] = np.log(integration / (1 - integration))
    for metric_col in [
        "mean_cross_edge_fraction",
        "mean_composition_entropy_norm",
        "mean_pair_balance",
    ]:
        if metric_col in events.columns:
            metric = pd.to_numeric(events[metric_col], errors="coerce").clip(1e-4, 1 - 1e-4)
            events[f"logit_{metric_col}"] = np.log(metric / (1 - metric))
    for col in [
        "log1p_centroid_distance_median_1wk_before_after_m",
        "log1p_centroid_distance_median_1wk_before_m",
        "log1p_centroid_distance_median_1wk_after_m",
        "log1p_prior_meetings",
        "log1p_prior_meeting_duration_hours",
        "dyad_ndvi_mean_1wk_before_after",
        "dyad_ndvi_abs_diff_1wk_before_after",
        "log1p_group_size_total_1wk_before_after",
        "log1p_group_size_abs_diff_1wk_before_after",
    ]:
        events[f"{col}_z"] = zscore(events[col])
    return events


def build_event_bayesian_formula(
    response_col: str,
    include_ndvi: bool,
    distance_window: str = "before_after",
    extra_terms: list[str] | None = None,
) -> str:
    distance_term = EVENT_DISTANCE_WINDOWS[distance_window]["term"]
    terms = [
        distance_term,
        "log1p_prior_meetings_z",
        "log1p_group_size_total_1wk_before_after_z",
        "log1p_group_size_abs_diff_1wk_before_after_z",
    ]
    if include_ndvi:
        terms.append("dyad_ndvi_mean_1wk_before_after_z")
    if extra_terms:
        terms.extend(extra_terms)
    fixed = response_col + " ~ " + " + ".join(terms)
    random = " + (1|group_a) + (1|group_b) + (1|pair_key)"
    return fixed + random


def build_event_duration_formula(include_ndvi: bool, distance_window: str = "before_after") -> str:
    return build_event_bayesian_formula("log1p_duration_hours", include_ndvi, distance_window=distance_window)


class BayesianLinearMixedResult:
    def __init__(
        self,
        term_names: list[str],
        fixed_terms: list[str],
        posterior_draws: np.ndarray,
        sigma2_draws: np.ndarray,
        formula: str,
        prior_mean: np.ndarray,
        prior_sd: np.ndarray,
    ) -> None:
        self.term_names = term_names
        self.fixed_terms = fixed_terms
        self.posterior_draws = posterior_draws
        self.sigma2_draws = sigma2_draws
        self.formula = formula
        self.prior_mean = pd.Series(prior_mean, index=term_names)
        self.prior_sd = pd.Series(prior_sd, index=term_names)
        self.params = pd.Series(posterior_draws.mean(axis=0), index=term_names)
        self.bse = pd.Series(posterior_draws.std(axis=0, ddof=1), index=term_names)
        prob_pos = (posterior_draws > 0).mean(axis=0)
        self.pvalues = pd.Series(2 * np.minimum(prob_pos, 1 - prob_pos), index=term_names)

    def conf_int(self) -> pd.DataFrame:
        ci = np.quantile(self.posterior_draws, [0.025, 0.975], axis=0).T
        return pd.DataFrame(ci, index=self.term_names, columns=[0, 1])

    def fixed_draws(self) -> np.ndarray:
        fixed_idx = [self.term_names.index(term) for term in self.fixed_terms]
        return self.posterior_draws[:, fixed_idx]


def sigmoid(values: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-np.clip(values, -35, 35)))


class BayesianLogisticMixedResult:
    def __init__(
        self,
        term_names: list[str],
        fixed_terms: list[str],
        posterior_draws: np.ndarray,
        formula: str,
        prior_mean: np.ndarray,
        prior_sd: np.ndarray,
    ) -> None:
        self.term_names = term_names
        self.fixed_terms = fixed_terms
        self.posterior_draws = posterior_draws
        self.formula = formula
        self.prior_mean = pd.Series(prior_mean, index=term_names)
        self.prior_sd = pd.Series(prior_sd, index=term_names)
        self.params = pd.Series(posterior_draws.mean(axis=0), index=term_names)
        self.bse = pd.Series(posterior_draws.std(axis=0, ddof=1), index=term_names)
        prob_pos = (posterior_draws > 0).mean(axis=0)
        self.pvalues = pd.Series(2 * np.minimum(prob_pos, 1 - prob_pos), index=term_names)

    def conf_int(self) -> pd.DataFrame:
        ci = np.quantile(self.posterior_draws, [0.025, 0.975], axis=0).T
        return pd.DataFrame(ci, index=self.term_names, columns=[0, 1])

    def fixed_draws(self) -> np.ndarray:
        fixed_idx = [self.term_names.index(term) for term in self.fixed_terms]
        return self.posterior_draws[:, fixed_idx]


def build_meeting_probability_formula(include_ndvi: bool, hierarchical: bool) -> str:
    terms = [
        "log1p_range_mean_centroid_distance_m_z",
        "log1p_prior_pair_interaction_days_z",
        "log1p_group_size_total_z",
        "log1p_group_size_abs_diff_z",
    ]
    if include_ndvi:
        terms.append("dyad_daily_ndvi_mean_z")
    formula = "any_interaction ~ " + " + ".join(terms)
    if hierarchical:
        formula += " + (1|group_a) + (1|group_b) + (1|pair_key)"
    return formula


def make_bayesian_meeting_design(
    model_df: pd.DataFrame,
    include_ndvi: bool,
    hierarchical: bool,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str], np.ndarray, np.ndarray]:
    fixed_terms = [
        "Intercept",
        "log1p_range_mean_centroid_distance_m_z",
        "log1p_prior_pair_interaction_days_z",
        "log1p_group_size_total_z",
        "log1p_group_size_abs_diff_z",
    ]
    fixed_arrays = [
        np.ones(len(model_df)),
        model_df["log1p_range_mean_centroid_distance_m_z"].to_numpy(float),
        model_df["log1p_prior_pair_interaction_days_z"].to_numpy(float),
        model_df["log1p_group_size_total_z"].to_numpy(float),
        model_df["log1p_group_size_abs_diff_z"].to_numpy(float),
    ]
    if include_ndvi:
        fixed_terms.append("dyad_daily_ndvi_mean_z")
        fixed_arrays.append(model_df["dyad_daily_ndvi_mean_z"].to_numpy(float))

    random_arrays = []
    random_terms = []
    if hierarchical:
        for col, prefix in [
            ("group_a", "group1"),
            ("group_b", "group2"),
            ("pair_key", "dyad"),
        ]:
            levels = sorted(model_df[col].astype(str).dropna().unique())
            for level in levels:
                random_terms.append(f"{prefix}[{level}]")
                random_arrays.append(model_df[col].astype(str).eq(level).astype(float).to_numpy())

    term_names = fixed_terms + random_terms
    X = np.column_stack(fixed_arrays + random_arrays)
    y = model_df["any_interaction"].to_numpy(float)
    base_probability = float(np.clip(np.nanmean(y), 1e-4, 1 - 1e-4))
    prior_mean_by_term = {
        "Intercept": float(np.log(base_probability / (1 - base_probability))),
        "log1p_range_mean_centroid_distance_m_z": -0.35,
        "log1p_prior_pair_interaction_days_z": 0.20,
        "log1p_group_size_total_z": 0.10,
        "log1p_group_size_abs_diff_z": -0.05,
        "dyad_daily_ndvi_mean_z": 0.05,
    }
    prior_mean = np.array([prior_mean_by_term.get(term, 0.0) for term in term_names], dtype=float)
    prior_sd = np.array(
        [5.0] + [1.5] * (len(fixed_terms) - 1) + [0.75] * len(random_terms),
        dtype=float,
    )
    return X, y, term_names, fixed_terms, prior_mean, prior_sd


def fit_meeting_probability_model(
    weekly_rows: pd.DataFrame,
    include_ndvi: bool,
    hierarchical: bool,
) -> tuple[object, str, str, pd.DataFrame]:
    formula = build_meeting_probability_formula(include_ndvi, hierarchical)
    cols = [
        "pair_key",
        "group_a",
        "group_b",
        "any_interaction",
        "log1p_range_mean_centroid_distance_m_z",
        "log1p_prior_pair_interaction_days_z",
        "log1p_group_size_total_z",
        "log1p_group_size_abs_diff_z",
    ]
    if include_ndvi:
        cols.append("dyad_daily_ndvi_mean_z")
    model_df = weekly_rows[cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if len(model_df) < MIN_MODEL_ROWS or model_df["pair_key"].nunique() < 2:
        raise ValueError("Not enough dyad intervals to fit meeting-probability model.")

    X, y, term_names, fixed_terms, prior_mean, prior_sd = make_bayesian_meeting_design(
        model_df,
        include_ndvi,
        hierarchical,
    )
    prior_precision = np.diag(1.0 / (prior_sd**2))
    beta = prior_mean.copy()
    for _ in range(100):
        eta = X @ beta
        p = sigmoid(eta)
        weights = np.clip(p * (1 - p), 1e-6, np.inf)
        hessian = X.T @ (weights[:, None] * X) + prior_precision
        score = X.T @ (y - p) - prior_precision @ (beta - prior_mean)
        step = np.linalg.solve(hessian, score)
        beta = beta + step
        if float(np.max(np.abs(step))) < 1e-7:
            break
    eta = X @ beta
    p = sigmoid(eta)
    weights = np.clip(p * (1 - p), 1e-6, np.inf)
    hessian = X.T @ (weights[:, None] * X) + prior_precision
    posterior_cov = np.linalg.inv(hessian)
    rng = np.random.default_rng(20260710 if hierarchical else 20260711)
    posterior_draws = rng.multivariate_normal(beta, posterior_cov, size=BAYESIAN_DRAWS)
    result = BayesianLogisticMixedResult(
        term_names=term_names,
        fixed_terms=fixed_terms,
        posterior_draws=posterior_draws,
        formula=formula,
        prior_mean=prior_mean,
        prior_sd=prior_sd,
    )
    if hierarchical:
        note = (
            "Approximate Bayesian logistic mixed model with weak informed Normal fixed-effect priors "
            "and Gaussian random-intercept priors for group 1, group 2, and dyad"
        )
    else:
        note = "Approximate Bayesian logistic model with weak informed Normal fixed-effect priors"
    return result, formula, note, model_df


def make_bayesian_event_design(
    model_df: pd.DataFrame,
    include_ndvi: bool,
    response_col: str,
    model_kind: str,
    distance_window: str = "before_after",
    extra_terms: list[str] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str], np.ndarray, np.ndarray]:
    distance_term = EVENT_DISTANCE_WINDOWS[distance_window]["term"]
    fixed_terms = [
        "Intercept",
        distance_term,
        "log1p_prior_meetings_z",
        "log1p_group_size_total_1wk_before_after_z",
        "log1p_group_size_abs_diff_1wk_before_after_z",
    ]
    fixed_arrays = [
        np.ones(len(model_df)),
        model_df[distance_term].to_numpy(float),
        model_df["log1p_prior_meetings_z"].to_numpy(float),
        model_df["log1p_group_size_total_1wk_before_after_z"].to_numpy(float),
        model_df["log1p_group_size_abs_diff_1wk_before_after_z"].to_numpy(float),
    ]
    if include_ndvi:
        fixed_terms.append("dyad_ndvi_mean_1wk_before_after_z")
        fixed_arrays.append(model_df["dyad_ndvi_mean_1wk_before_after_z"].to_numpy(float))
    if extra_terms:
        for term in extra_terms:
            fixed_terms.append(term)
            fixed_arrays.append(model_df[term].to_numpy(float))

    random_arrays = []
    random_terms = []
    for col, prefix in [
        ("group_a", "group1"),
        ("group_b", "group2"),
        ("pair_key", "dyad"),
    ]:
        levels = sorted(model_df[col].astype(str).dropna().unique())
        for level in levels:
            random_terms.append(f"{prefix}[{level}]")
            random_arrays.append(model_df[col].astype(str).eq(level).astype(float).to_numpy())

    term_names = fixed_terms + random_terms
    X = np.column_stack(fixed_arrays + random_arrays)
    y = model_df[response_col].to_numpy(float)
    prior_mean_by_term = {
        "Intercept": float(np.nanmean(y)),
        distance_term: -0.15,
        "log1p_prior_meetings_z": 0.10,
        "log1p_group_size_total_1wk_before_after_z": 0.10,
        "log1p_group_size_abs_diff_1wk_before_after_z": -0.05,
        "dyad_ndvi_mean_1wk_before_after_z": 0.05,
        "log1p_event_gps_movement_mean_m_per_h_z": 0.10,
    }
    if model_kind == "integration":
        prior_mean_by_term.update(
            {
                distance_term: -0.25,
                "log1p_prior_meetings_z": 0.15,
            }
        )
    prior_mean = np.array([prior_mean_by_term.get(term, 0.0) for term in term_names], dtype=float)
    prior_sd = np.array(
        [5.0] + [1.5] * (len(fixed_terms) - 1) + [0.75] * len(random_terms),
        dtype=float,
    )
    return X, y, term_names, fixed_terms, prior_mean, prior_sd


def fit_event_bayesian_model(
    event_rows: pd.DataFrame,
    response_col: str,
    include_ndvi: bool,
    model_kind: str,
    distance_window: str = "before_after",
    extra_terms: list[str] | None = None,
) -> tuple[object, str, str, pd.DataFrame]:
    distance_term = EVENT_DISTANCE_WINDOWS[distance_window]["term"]
    formula = build_event_bayesian_formula(
        response_col,
        include_ndvi,
        distance_window=distance_window,
        extra_terms=extra_terms,
    )
    cols = [
        "pair_key",
        "group_a",
        "group_b",
        response_col,
        distance_term,
        "log1p_prior_meetings_z",
        "log1p_group_size_total_1wk_before_after_z",
        "log1p_group_size_abs_diff_1wk_before_after_z",
    ]
    if include_ndvi:
        cols.append("dyad_ndvi_mean_1wk_before_after_z")
    if extra_terms:
        cols.extend(extra_terms)
    model_df = event_rows[cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if len(model_df) < MIN_MODEL_ROWS or model_df["pair_key"].nunique() < 2:
        raise ValueError(f"Not enough interaction episodes to fit {model_kind} model.")

    X, y, term_names, fixed_terms, prior_mean, prior_sd = make_bayesian_event_design(
        model_df,
        include_ndvi,
        response_col,
        model_kind,
        distance_window=distance_window,
        extra_terms=extra_terms,
    )
    prior_precision = np.diag(1.0 / (prior_sd**2))
    posterior_precision = prior_precision + X.T @ X
    posterior_cov_unit = np.linalg.inv(posterior_precision)
    posterior_mean = posterior_cov_unit @ (prior_precision @ prior_mean + X.T @ y)
    prior_shape = BAYESIAN_PRIOR_SHAPE
    prior_scale = BAYESIAN_PRIOR_SCALE
    post_shape = prior_shape + len(y) / 2.0
    post_scale = prior_scale + 0.5 * (
        y @ y
        + prior_mean @ prior_precision @ prior_mean
        - posterior_mean @ posterior_precision @ posterior_mean
    )
    post_scale = float(max(post_scale, 1e-9))
    rng = np.random.default_rng(20260709)
    sigma2_draws = post_scale / rng.gamma(post_shape, 1.0, size=BAYESIAN_DRAWS)
    chol = np.linalg.cholesky(posterior_cov_unit)
    standard = rng.normal(size=(BAYESIAN_DRAWS, len(term_names)))
    posterior_draws = posterior_mean + (standard @ chol.T) * np.sqrt(sigma2_draws)[:, None]
    result = BayesianLinearMixedResult(
        term_names=term_names,
        fixed_terms=fixed_terms,
        posterior_draws=posterior_draws,
        sigma2_draws=sigma2_draws,
        formula=formula,
        prior_mean=prior_mean,
        prior_sd=prior_sd,
    )
    note = (
        "Bayesian Gaussian linear mixed model with weak informed Normal fixed-effect priors, "
        "Gaussian random-intercept priors for group 1, group 2, and dyad, "
        "and inverse-gamma residual variance prior"
    )
    return result, formula, note, model_df


def fit_event_duration_model(event_rows: pd.DataFrame, include_ndvi: bool) -> tuple[object, str, str, pd.DataFrame]:
    return fit_event_bayesian_model(
        event_rows,
        response_col="log1p_duration_hours",
        include_ndvi=include_ndvi,
        model_kind="duration",
        distance_window="before_after",
    )


def fit_event_integration_model(event_rows: pd.DataFrame, include_ndvi: bool) -> tuple[object, str, str, pd.DataFrame]:
    return fit_event_bayesian_model(
        event_rows,
        response_col="logit_integration_5m_fraction",
        include_ndvi=include_ndvi,
        model_kind="integration",
        distance_window="before_after",
    )


def meeting_probability_summary_row(
    model_rows_all: pd.DataFrame,
    model_df: pd.DataFrame,
    formula: str,
    note: str,
    model_name: str,
    model_type: str,
    time_unit: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model": model_name,
                "model_type": model_type,
                "time_unit": time_unit,
                "response": "any_interaction",
                "n_dyad_intervals": int(len(model_df)),
                "n_dyads": int(model_df["pair_key"].nunique()),
                "n_meeting_intervals": int(model_df["any_interaction"].sum()),
                "meeting_probability": float(model_df["any_interaction"].mean()),
                "median_prior_meetings": float(model_rows_all["prior_pair_interaction_days"].median()),
                "formula": formula,
                "fit_note": note,
            }
        ]
    )


def meeting_probability_prediction_grid(
    weekly_rows: pd.DataFrame,
    result: object,
    include_ndvi: bool,
    predictor: str,
    model_name: str,
    model_type: str,
) -> pd.DataFrame:
    reference_cols = [
        "range_mean_centroid_distance_m",
        "log1p_range_mean_centroid_distance_m",
        "prior_pair_interaction_days",
        "log1p_prior_pair_interaction_days",
        "log1p_group_size_total_z",
        "log1p_group_size_abs_diff_z",
    ]
    if include_ndvi:
        reference_cols.append("dyad_daily_ndvi_mean")
    reference = weekly_rows.replace([np.inf, -np.inf], np.nan).dropna(subset=reference_cols).copy()
    n_grid = 140
    grid = pd.DataFrame(
        {
            "range_mean_centroid_distance_m": np.repeat(
                reference["range_mean_centroid_distance_m"].median(), n_grid
            ),
            "prior_pair_interaction_days": np.repeat(reference["prior_pair_interaction_days"].median(), n_grid),
            "log1p_group_size_total_z": 0.0,
            "log1p_group_size_abs_diff_z": 0.0,
        }
    )
    if include_ndvi:
        grid["dyad_daily_ndvi_mean"] = reference["dyad_daily_ndvi_mean"].mean()

    predictor_columns = {
        "centroid_distance": (
            "range_mean_centroid_distance_m",
            "Mean centroid distance for dyad range (m)",
            True,
        ),
        "history": ("prior_pair_interaction_days", "Prior meeting weeks for dyad", False),
        "ndvi_mean": ("dyad_daily_ndvi_mean", "Dyad mean NDVI", False),
    }
    raw_col, label, log_x = predictor_columns[predictor]
    if predictor in {"centroid_distance", "history"}:
        lo = 0.0 if predictor == "history" else reference[raw_col].quantile(0.01)
        hi = max(float(reference[raw_col].quantile(0.99)), 1.0)
        values = np.linspace(lo, hi, n_grid)
    else:
        if not include_ndvi:
            return pd.DataFrame()
        values = np.linspace(reference[raw_col].quantile(0.01), reference[raw_col].quantile(0.99), n_grid)
    grid[raw_col] = values

    grid["log1p_range_mean_centroid_distance_m"] = np.log1p(
        grid["range_mean_centroid_distance_m"].clip(lower=0)
    )
    grid["log1p_prior_pair_interaction_days"] = np.log1p(grid["prior_pair_interaction_days"])
    for raw in [
        "log1p_range_mean_centroid_distance_m",
        "log1p_prior_pair_interaction_days",
        "dyad_daily_ndvi_mean",
    ]:
        if raw not in grid:
            continue
        mean = weekly_rows[raw].mean()
        sd = weekly_rows[raw].std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        grid[f"{raw}_z"] = (grid[raw] - mean) / sd

    fixed_arrays = [
        np.ones(len(grid)),
        grid["log1p_range_mean_centroid_distance_m_z"].to_numpy(float),
        grid["log1p_prior_pair_interaction_days_z"].to_numpy(float),
        grid["log1p_group_size_total_z"].to_numpy(float),
        grid["log1p_group_size_abs_diff_z"].to_numpy(float),
    ]
    if include_ndvi:
        fixed_arrays.append(grid["dyad_daily_ndvi_mean_z"].to_numpy(float))
    fixed_design = np.column_stack(fixed_arrays)
    pred_draws = sigmoid(result.fixed_draws() @ fixed_design.T)
    return pd.DataFrame(
        {
            "predictor": predictor,
            "predictor_value": values,
            "predictor_label": label,
            "x_log_scale": log_x,
            "predicted_probability": pred_draws.mean(axis=0),
            "ci_low": np.quantile(pred_draws, 0.025, axis=0),
            "ci_high": np.quantile(pred_draws, 0.975, axis=0),
            "model": model_name,
            "model_type": model_type,
        }
    )


def plot_meeting_probability_predictions(
    weekly_rows: pd.DataFrame,
    predictions: pd.DataFrame,
    out_dir: Path,
    filename_prefix: str = "meeting_probability",
    y_label: str = "Weekly probability of meeting",
) -> list[Path]:
    paths = []
    plot_specs = [
        (
            "centroid_distance",
            "range_mean_centroid_distance_m",
            "Mean centroid distance for dyad range (m, log scale)",
            True,
            f"{filename_prefix}_centroid_distance.png",
        ),
        ("history", "prior_pair_interaction_days", "Prior meeting intervals for dyad", False, f"{filename_prefix}_history.png"),
        ("ndvi_mean", "dyad_daily_ndvi_mean", "Dyad mean NDVI", False, f"{filename_prefix}_ndvi_mean.png"),
    ]
    rng = np.random.default_rng(20260712)
    for predictor, x_col, xlabel, x_log, filename in plot_specs:
        pred = predictions[predictions["predictor"].eq(predictor)]
        if pred.empty or x_col not in weekly_rows.columns:
            continue
        plot_df = weekly_rows.dropna(subset=[x_col, "any_interaction"]).copy()
        fig, ax = plt.subplots(figsize=(9.5, 5.8))
        ax.scatter(
            plot_df[x_col],
            plot_df["any_interaction"] + rng.normal(0, 0.015, size=len(plot_df)),
            s=14,
            alpha=0.18,
            color="#f28e2b",
            edgecolor="none",
        )
        for model_type, sub in pred.groupby("model_type", observed=True):
            sub = sub.sort_values("predictor_value")
            line_style = "-" if model_type == "hierarchical" else "--"
            ax.plot(
                sub["predictor_value"],
                sub["predicted_probability"],
                linewidth=2.3,
                linestyle=line_style,
                label=model_type,
            )
            ax.fill_between(sub["predictor_value"], sub["ci_low"], sub["ci_high"], alpha=0.13)
        if x_log:
            ax.set_xscale("log")
        ax.set_ylim(-0.08, 1.08)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(y_label)
        ax.set_title("Probability of meeting")
        ax.grid(True, alpha=0.22)
        ax.legend(frameon=False)
        fig.tight_layout()
        path = out_dir / filename
        fig.savefig(path, dpi=220)
        plt.close(fig)
        paths.append(path)
    return paths


def run_meeting_probability_analyses(
    weekly_rows: pd.DataFrame,
    out_dir: Path,
    include_ndvi: bool,
    time_unit: str = "weekly",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[Path]]:
    outputs = []
    for hierarchical, model_name, model_type in [
        (False, f"meeting_probability_{time_unit}_separate", "separate"),
        (True, f"meeting_probability_{time_unit}_hierarchical", "hierarchical"),
    ]:
        result, formula, note, model_df = fit_meeting_probability_model(
            weekly_rows,
            include_ndvi=include_ndvi,
            hierarchical=hierarchical,
        )
        summary = meeting_probability_summary_row(
            weekly_rows,
            model_df,
            formula,
            note,
            model_name=model_name,
            model_type=model_type,
            time_unit=time_unit,
        )
        coef = coefficient_table(result, model_name, "any_interaction")
        coef["model_type"] = model_type
        coef["time_unit"] = time_unit
        predictions = pd.concat(
            [
                grid
                for grid in (
                    meeting_probability_prediction_grid(
                        weekly_rows,
                        result,
                        include_ndvi,
                        predictor,
                        model_name=model_name,
                        model_type=model_type,
                    )
                    for predictor in ["centroid_distance", "history", "ndvi_mean"]
                )
                if not grid.empty
            ],
            ignore_index=True,
        )
        predictions["time_unit"] = time_unit
        outputs.append((summary, coef, predictions))
    summary = pd.concat([item[0] for item in outputs], ignore_index=True)
    coef = pd.concat([item[1] for item in outputs], ignore_index=True)
    predictions = pd.concat([item[2] for item in outputs], ignore_index=True)
    plots = plot_meeting_probability_predictions(
        weekly_rows,
        predictions,
        out_dir,
        filename_prefix=f"meeting_probability_{time_unit}",
        y_label=f"{time_unit.capitalize()} probability of meeting",
    )
    return summary, coef, predictions, plots


def event_duration_summary_row(
    event_rows: pd.DataFrame,
    model_df: pd.DataFrame,
    formula: str,
    note: str,
    model_name: str = "event_duration",
    merge_size_class: str = "all",
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model": model_name,
                "merge_size_class": merge_size_class,
                "response": "log1p_duration_hours",
                "n_events": int(len(model_df)),
                "n_dyads": int(model_df["pair_key"].nunique()),
                "median_duration_hours": float(event_rows["duration_hours"].median()),
                "max_duration_hours": float(event_rows["duration_hours"].max()),
                "max_duration_days": float(event_rows["duration_days"].max()),
                "median_active_contact_hours": float(event_rows["active_contact_hours"].median()),
                "median_prior_meetings": float(event_rows["prior_meetings"].median()),
                "formula": formula,
                "fit_note": note,
            }
        ]
    )


def event_integration_summary_row(
    event_rows: pd.DataFrame,
    model_df: pd.DataFrame,
    formula: str,
    note: str,
    model_name: str = "event_integration_5m",
    merge_size_class: str = "all",
    ndvi_used: bool = False,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model": model_name,
                "merge_size_class": merge_size_class,
                "response": "logit_integration_5m_fraction",
                "n_events": int(len(model_df)),
                "n_dyads": int(model_df["pair_key"].nunique()),
                "median_integration_5m_fraction": float(event_rows["integration_5m_fraction"].median()),
                "mean_integration_5m_fraction": float(event_rows["integration_5m_fraction"].mean()),
                "median_cross_group_5m_edges": float(event_rows["total_cross_group_5m_edges"].median()),
                "median_candidate_5m_edges": float(event_rows["total_candidate_5m_edges"].median()),
                "median_duration_hours": float(event_rows["duration_hours"].median()),
                "median_prior_meetings": float(event_rows["prior_meetings"].median()),
                "ndvi_used": ndvi_used,
                "formula": formula,
                "fit_note": note,
            }
        ]
    )


def event_duration_prediction_grid(
    event_rows: pd.DataFrame,
    result: object,
    include_ndvi: bool,
    predictor: str,
    model_name: str = "event_duration",
    merge_size_class: str = "all",
    response_kind: str = "duration",
) -> pd.DataFrame:
    reference_cols = [
        "centroid_distance_median_1wk_before_after_m",
        "log1p_centroid_distance_median_1wk_before_after_m",
        "prior_meetings",
        "log1p_prior_meetings",
        "log1p_group_size_total_1wk_before_after_z",
        "log1p_group_size_abs_diff_1wk_before_after_z",
        "pair_key",
    ]
    if include_ndvi:
        reference_cols.append("dyad_ndvi_mean_1wk_before_after")
    reference = event_rows.replace([np.inf, -np.inf], np.nan).dropna(subset=reference_cols).copy()
    n_grid = 140
    grid = pd.DataFrame(
        {
            "centroid_distance_median_1wk_before_after_m": np.repeat(
                reference["centroid_distance_median_1wk_before_after_m"].median(), n_grid
            ),
            "prior_meetings": np.repeat(reference["prior_meetings"].median(), n_grid),
            "log1p_group_size_total_1wk_before_after_z": 0.0,
            "log1p_group_size_abs_diff_1wk_before_after_z": 0.0,
            "pair_key": reference["pair_key"].mode().iloc[0],
        }
    )
    if include_ndvi:
        grid["dyad_ndvi_mean_1wk_before_after"] = event_rows["dyad_ndvi_mean_1wk_before_after"].mean()

    predictor_columns = {
        "centroid_distance": (
            "centroid_distance_median_1wk_before_after_m",
            "Median centroid distance in +/- 1 week event window (m)",
            True,
        ),
        "history": ("prior_meetings", "Prior meetings for dyad", False),
        "ndvi_mean": ("dyad_ndvi_mean_1wk_before_after", "Dyad mean NDVI in +/- 1 week event window", False),
    }
    raw_col, label, log_x = predictor_columns[predictor]
    if predictor in {"centroid_distance", "history"}:
        lo = 0.0 if predictor == "history" else reference[raw_col].quantile(0.01)
        hi = max(float(reference[raw_col].quantile(0.99)), 1.0)
        values = np.linspace(lo, hi, n_grid)
    else:
        if not include_ndvi:
            return pd.DataFrame()
        values = np.linspace(reference[raw_col].quantile(0.01), reference[raw_col].quantile(0.99), n_grid)
    grid[raw_col] = values

    grid["log1p_centroid_distance_median_1wk_before_after_m"] = np.log1p(
        grid["centroid_distance_median_1wk_before_after_m"].clip(lower=0)
    )
    grid["log1p_prior_meetings"] = np.log1p(grid["prior_meetings"])
    for raw in [
        "log1p_centroid_distance_median_1wk_before_after_m",
        "log1p_prior_meetings",
        "dyad_ndvi_mean_1wk_before_after",
    ]:
        if raw not in grid:
            continue
        mean = event_rows[raw].mean()
        sd = event_rows[raw].std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        grid[f"{raw}_z"] = (grid[raw] - mean) / sd

    fixed_arrays = [
        np.ones(len(grid)),
        grid["log1p_centroid_distance_median_1wk_before_after_m_z"].to_numpy(float),
        grid["log1p_prior_meetings_z"].to_numpy(float),
        grid["log1p_group_size_total_1wk_before_after_z"].to_numpy(float),
        grid["log1p_group_size_abs_diff_1wk_before_after_z"].to_numpy(float),
    ]
    if include_ndvi:
        fixed_arrays.append(grid["dyad_ndvi_mean_1wk_before_after_z"].to_numpy(float))
    fixed_design = np.column_stack(fixed_arrays)
    pred_draws_linear = result.fixed_draws() @ fixed_design.T
    if response_kind == "integration":
        pred_draws = 1 / (1 + np.exp(-pred_draws_linear))
    else:
        pred_draws = np.expm1(pred_draws_linear)
    pred = pred_draws.mean(axis=0)
    low = np.quantile(pred_draws, 0.025, axis=0)
    high = np.quantile(pred_draws, 0.975, axis=0)
    out = pd.DataFrame(
        {
            "predictor": predictor,
            "predictor_value": values,
            "predictor_label": label,
            "x_log_scale": log_x,
            "predicted": pred,
            "ci_low": low,
            "ci_high": high,
            "model": model_name,
            "merge_size_class": merge_size_class,
            "response_kind": response_kind,
        }
    )
    if response_kind == "integration":
        out["predicted_integration_5m_fraction"] = out["predicted"].clip(0, 1)
        out["ci_low"] = out["ci_low"].clip(0, 1)
        out["ci_high"] = out["ci_high"].clip(0, 1)
    else:
        out["predicted_duration_hours"] = out["predicted"]
        out["ci_low"] = out["ci_low"].clip(EVENT_BIN_WIDTH_MINUTES / 60.0, np.inf)
        out["ci_high"] = out["ci_high"].clip(EVENT_BIN_WIDTH_MINUTES / 60.0, np.inf)
    return out


def plot_event_duration_predictions(
    event_rows: pd.DataFrame,
    predictions: pd.DataFrame,
    out_dir: Path,
    filename_prefix: str = "event_duration",
    title_prefix: str = "All merge events",
    response_kind: str = "duration",
) -> list[Path]:
    paths = []
    plot_specs = [
        (
            "centroid_distance",
            "centroid_distance_median_1wk_before_after_m",
            "Median centroid distance in +/- 1 week event window (m, log scale)",
            True,
            f"{filename_prefix}_centroid_distance.png",
        ),
        ("history", "prior_meetings", "Prior meetings for dyad", False, f"{filename_prefix}_history.png"),
        (
            "ndvi_mean",
            "dyad_ndvi_mean_1wk_before_after",
            "Dyad mean NDVI in +/- 1 week event window",
            False,
            f"{filename_prefix}_ndvi_mean.png",
        ),
    ]
    for predictor, x_col, xlabel, x_log, filename in plot_specs:
        if x_col not in event_rows.columns:
            continue
        pred = predictions[predictions["predictor"].eq(predictor)]
        if pred.empty:
            continue
        y_col = "integration_5m_fraction" if response_kind == "integration" else "duration_hours"
        plot_df = event_rows.dropna(subset=[x_col, y_col]).copy()
        fig, ax = plt.subplots(figsize=(9.5, 5.8))
        sizes = np.clip(np.sqrt(plot_df["positive_2min_bins"]) * 8, 16, 110)
        color = "#59a14f" if response_kind == "integration" else "#4c78a8"
        ax.scatter(
            plot_df[x_col],
            plot_df[y_col],
            s=sizes,
            alpha=0.35,
            color=color,
            edgecolor="white",
            linewidth=0.35,
        )
        pred = pred.sort_values("predictor_value")
        pred_col = (
            "predicted_integration_5m_fraction"
            if response_kind == "integration"
            else "predicted_duration_hours"
        )
        ax.plot(pred["predictor_value"], pred[pred_col], color="#e15759", linewidth=2.3)
        ax.fill_between(pred["predictor_value"], pred["ci_low"], pred["ci_high"], color="#e15759", alpha=0.16)
        if x_log:
            ax.set_xscale("log")
        if response_kind == "integration":
            ax.set_ylim(-0.04, 1.04)
            ax.set_ylabel("Event 5 m cross-group integration fraction")
        else:
            ax.set_yscale("log")
            ax.set_ylabel("Real event duration (hours, log scale)")
        ax.set_xlabel(xlabel)
        ax.set_title(title_prefix)
        ax.grid(True, alpha=0.22)
        fig.tight_layout()
        path = out_dir / filename
        fig.savefig(path, dpi=220)
        plt.close(fig)
        paths.append(path)
    return paths


def run_event_duration_analysis(
    event_rows: pd.DataFrame,
    out_dir: Path,
    merge_size_class: str,
    filename_prefix: str,
    title_prefix: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[Path], bool]:
    subset = (
        event_rows.copy()
        if merge_size_class == "all"
        else event_rows[event_rows["merge_size_class"].eq(merge_size_class)].copy()
    )
    include_event_ndvi = bool(
        subset["dyad_ndvi_mean_1wk_before_after_z"].notna().sum()
        >= MIN_MODEL_ROWS
    )
    result, formula, note, model_df = fit_event_duration_model(
        subset,
        include_ndvi=include_event_ndvi,
    )
    model_name = "event_duration" if merge_size_class == "all" else f"event_duration_{merge_size_class}"
    summary = event_duration_summary_row(
        subset,
        model_df,
        formula,
        note,
        model_name=model_name,
        merge_size_class=merge_size_class,
    )
    coef = coefficient_table(result, model_name, "log1p_duration_hours")
    coef["merge_size_class"] = merge_size_class
    predictions = pd.concat(
        [
            grid
            for grid in (
                event_duration_prediction_grid(
                    subset,
                    result,
                    include_event_ndvi,
                    predictor,
                    model_name=model_name,
                    merge_size_class=merge_size_class,
                    response_kind="duration",
                )
                for predictor in ["centroid_distance", "history", "ndvi_mean"]
            )
            if not grid.empty
        ],
        ignore_index=True,
    )
    plots = plot_event_duration_predictions(
        subset,
        predictions,
        out_dir,
        filename_prefix=filename_prefix,
        title_prefix=title_prefix,
    )
    return summary, coef, predictions, plots, include_event_ndvi


def run_event_integration_analysis(
    event_rows: pd.DataFrame,
    out_dir: Path,
    filename_prefix: str = "event_integration_5m",
    title_prefix: str = "All merge events: 5 m integration",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[Path], bool]:
    subset = event_rows[event_rows["integration_5m_fraction"].notna()].copy()
    include_event_ndvi = bool(
        subset["dyad_ndvi_mean_1wk_before_after_z"].notna().all()
        and len(subset) >= MIN_MODEL_ROWS
    )
    result, formula, note, model_df = fit_event_integration_model(
        subset,
        include_ndvi=include_event_ndvi,
    )
    model_name = "event_integration_5m"
    summary = event_integration_summary_row(
        subset,
        model_df,
        formula,
        note,
        model_name=model_name,
        merge_size_class="all",
        ndvi_used=include_event_ndvi,
    )
    coef = coefficient_table(result, model_name, "logit_integration_5m_fraction")
    coef["merge_size_class"] = "all"
    predictions = pd.concat(
        [
            grid
            for grid in (
                event_duration_prediction_grid(
                    subset,
                    result,
                    include_event_ndvi,
                    predictor,
                    model_name=model_name,
                    merge_size_class="all",
                    response_kind="integration",
                )
                for predictor in ["centroid_distance", "history", "ndvi_mean"]
            )
            if not grid.empty
        ],
        ignore_index=True,
    )
    plots = plot_event_duration_predictions(
        subset,
        predictions,
        out_dir,
        filename_prefix=filename_prefix,
        title_prefix=title_prefix,
        response_kind="integration",
    )
    return summary, coef, predictions, plots, include_event_ndvi


INTEGRATION_METRIC_SPECS = [
    {
        "metric": "integration_5m_fraction",
        "response_col": "logit_integration_5m_fraction",
        "label": "5 m cross-group link fraction",
        "interpretation": "Higher values mean a larger fraction of 5 m candidate links are between the two groups.",
    },
    {
        "metric": "mean_cross_edge_fraction",
        "response_col": "logit_mean_cross_edge_fraction",
        "label": "Mean cross-edge fraction",
        "interpretation": "Higher values mean more observed proximity edges cross between groups during the event.",
    },
    {
        "metric": "mean_composition_entropy_norm",
        "response_col": "logit_mean_composition_entropy_norm",
        "label": "Composition entropy",
        "interpretation": "Higher values mean clusters are more compositionally mixed rather than dominated by one group.",
    },
    {
        "metric": "mean_pair_balance",
        "response_col": "logit_mean_pair_balance",
        "label": "Pair balance",
        "interpretation": "Higher values mean the two groups contribute more evenly to mixed clusters.",
    },
    {
        "metric": "mean_edge_modularity_q",
        "response_col": "mean_edge_modularity_q",
        "label": "Edge modularity Q",
        "interpretation": "Higher values mean stronger separation by group, so negative effects on this metric indicate more integration.",
    },
]


def integration_metric_summary_row(
    event_rows: pd.DataFrame,
    model_df: pd.DataFrame,
    formula: str,
    note: str,
    spec: dict[str, str],
    ndvi_used: bool,
    model_variant: str,
) -> pd.DataFrame:
    metric = spec["metric"]
    return pd.DataFrame(
        [
            {
                "metric": metric,
                "metric_label": spec["label"],
                "model_variant": model_variant,
                "response": spec["response_col"],
                "n_events": int(len(model_df)),
                "n_dyads": int(model_df["pair_key"].nunique()),
                "mean_metric": float(event_rows[metric].mean()),
                "median_metric": float(event_rows[metric].median()),
                "min_metric": float(event_rows[metric].min()),
                "max_metric": float(event_rows[metric].max()),
                "ndvi_used": ndvi_used,
                "n_events_with_ndvi": int(event_rows["dyad_ndvi_mean_1wk_before_after_z"].notna().sum()),
                "formula": formula,
                "fit_note": note,
                "interpretation": spec["interpretation"],
            }
        ]
    )


def plot_integration_metric_distributions(event_rows: pd.DataFrame, out_dir: Path) -> Path:
    specs = [spec for spec in INTEGRATION_METRIC_SPECS if spec["metric"] in event_rows.columns]
    n_cols = 2
    n_rows = int(np.ceil(len(specs) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11.0, max(4.2, 3.2 * n_rows)), squeeze=False)
    colors = ["#59a14f", "#4c78a8", "#f28e2b", "#af7aa1", "#e15759"]
    for ax, spec, color in zip(axes.ravel(), specs, colors):
        values = pd.to_numeric(event_rows[spec["metric"]], errors="coerce").dropna()
        ax.hist(values, bins=min(18, max(6, int(np.sqrt(len(values))))), color=color, alpha=0.78, edgecolor="white")
        ax.axvline(values.median(), color="#222222", linewidth=1.2, linestyle="--")
        ax.set_title(spec["label"])
        ax.set_xlabel("Metric value")
        ax.set_ylabel("Events")
        ax.grid(True, axis="y", alpha=0.18)
    for ax in axes.ravel()[len(specs):]:
        ax.axis("off")
    fig.tight_layout()
    path = out_dir / "integration_metric_distributions.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_integration_metric_relationships(event_rows: pd.DataFrame, out_dir: Path) -> Path:
    specs = [spec for spec in INTEGRATION_METRIC_SPECS if spec["metric"] in event_rows.columns]
    base = "integration_5m_fraction"
    compare_specs = [spec for spec in specs if spec["metric"] != base]
    n_cols = 2
    n_rows = int(np.ceil(len(compare_specs) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11.0, max(4.2, 3.4 * n_rows)), squeeze=False)
    for ax, spec in zip(axes.ravel(), compare_specs):
        plot_df = event_rows[[base, spec["metric"], "duration_hours"]].replace([np.inf, -np.inf], np.nan).dropna()
        sizes = np.clip(np.sqrt(plot_df["duration_hours"]) * 18, 20, 130)
        ax.scatter(
            plot_df[base],
            plot_df[spec["metric"]],
            s=sizes,
            alpha=0.42,
            color="#4c78a8",
            edgecolor="white",
            linewidth=0.35,
        )
        if len(plot_df) >= 3:
            corr = plot_df[base].corr(plot_df[spec["metric"]], method="spearman")
            ax.text(
                0.03,
                0.94,
                f"Spearman r = {corr:.2f}",
                transform=ax.transAxes,
                va="top",
                fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.9},
            )
        ax.set_xlabel("5 m cross-group link fraction")
        ax.set_ylabel(spec["label"])
        ax.grid(True, alpha=0.2)
    for ax in axes.ravel()[len(compare_specs):]:
        ax.axis("off")
    fig.tight_layout()
    path = out_dir / "integration_metric_relationships.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_integration_metric_effect_comparison(coef: pd.DataFrame, out_dir: Path) -> Path:
    fixed = fixed_effects_for_plot(coef, EVENT_EFFECT_LABELS)
    if fixed.empty:
        return out_dir / "integration_metric_effect_comparison.png"
    order = list(EVENT_EFFECT_LABELS.values())
    metric_order = [spec["label"] for spec in INTEGRATION_METRIC_SPECS]
    fixed["effect_label"] = pd.Categorical(fixed["effect_label"], categories=order, ordered=True)
    fixed["metric_label"] = pd.Categorical(fixed["metric_label"], categories=metric_order, ordered=True)
    effects = [effect for effect in order if effect in set(fixed["effect_label"].astype(str))]
    metrics = [metric for metric in metric_order if metric in set(fixed["metric_label"].astype(str))]
    fig, axes = plt.subplots(1, len(effects), figsize=(4.4 * len(effects), max(4.5, 0.48 * len(metrics) + 2.0)), squeeze=False)
    color = "#59a14f"
    for ax, effect in zip(axes.ravel(), effects):
        sub = fixed[fixed["effect_label"].astype(str).eq(effect)].sort_values("metric_label")
        y = np.arange(len(sub))
        ax.axvline(0, color="#333333", linewidth=1.0, alpha=0.65)
        ax.errorbar(
            sub["estimate"],
            y,
            xerr=[
                sub["estimate"] - sub["ci_low"],
                sub["ci_high"] - sub["estimate"],
            ],
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=1.8,
            capsize=3.5,
            markersize=6,
            alpha=0.92,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(sub["metric_label"] if effect == effects[0] else [])
        ax.invert_yaxis()
        ax.set_title(effect)
        ax.set_xlabel("Std. coefficient")
        ax.grid(True, axis="x", alpha=0.22)
    fig.suptitle("Integration metric effect-size comparison", y=1.02)
    fig.tight_layout()
    path = out_dir / "integration_metric_effect_comparison.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_integration_gps_control_effect_comparison(coef: pd.DataFrame, out_dir: Path) -> Path:
    fixed = fixed_effects_for_plot(coef, EVENT_EFFECT_LABELS)
    fixed = fixed[fixed["model_variant"].isin(["baseline", "gps_control"])].copy()
    if fixed.empty:
        return out_dir / "integration_metric_gps_control_effect_comparison.png"
    effect_order = [
        "Centroid distance",
        "Prior meetings",
        "Total group size",
        "Group size difference",
        "GPS movement",
    ]
    metric_order = [spec["label"] for spec in INTEGRATION_METRIC_SPECS]
    fixed["effect_label"] = pd.Categorical(fixed["effect_label"], categories=effect_order, ordered=True)
    fixed["metric_label"] = pd.Categorical(fixed["metric_label"], categories=metric_order, ordered=True)
    metrics = [metric for metric in metric_order if metric in set(fixed["metric_label"].astype(str))]
    effects = [effect for effect in effect_order if effect in set(fixed["effect_label"].astype(str))]
    fig, axes = plt.subplots(
        1,
        len(effects),
        figsize=(4.2 * len(effects), max(4.8, 0.5 * len(metrics) + 2.0)),
        squeeze=False,
    )
    colors = {"baseline": "#9aa0a6", "gps_control": "#f28e2b"}
    offsets = {"baseline": -0.11, "gps_control": 0.11}
    labels = {"baseline": "Baseline", "gps_control": "+ GPS movement"}
    for ax, effect in zip(axes.ravel(), effects):
        sub = fixed[fixed["effect_label"].astype(str).eq(effect)].sort_values(["metric_label", "model_variant"])
        ax.axvline(0, color="#333333", linewidth=1.0, alpha=0.65)
        for variant, variant_df in sub.groupby("model_variant", observed=True):
            positions = np.array([metrics.index(label) for label in variant_df["metric_label"].astype(str)]) + offsets[variant]
            ax.errorbar(
                variant_df["estimate"],
                positions,
                xerr=[
                    variant_df["estimate"] - variant_df["ci_low"],
                    variant_df["ci_high"] - variant_df["estimate"],
                ],
                fmt="o",
                color=colors[variant],
                ecolor=colors[variant],
                elinewidth=1.8,
                capsize=3,
                markersize=5.5,
                alpha=0.92,
                label=labels[variant],
            )
        ax.set_yticks(np.arange(len(metrics)))
        ax.set_yticklabels(metrics if effect == effects[0] else [])
        ax.invert_yaxis()
        ax.set_title(effect)
        ax.set_xlabel("Std. coefficient")
        ax.grid(True, axis="x", alpha=0.22)
    axes.ravel()[-1].legend(frameon=False, loc="lower right")
    fig.suptitle("Integration effects before and after controlling for GPS movement", y=1.02)
    fig.tight_layout()
    path = out_dir / "integration_metric_gps_control_effect_comparison.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def summarize_event_activity(
    event_rows: pd.DataFrame,
    proximity_path: Path,
    vedba_dir: Path,
) -> pd.DataFrame:
    if event_rows.empty or not vedba_dir.exists():
        return pd.DataFrame()
    events = event_rows[
        [
            "interaction_episode_id",
            "pair_key",
            "group_a",
            "group_b",
            "episode_start",
            "episode_end_exclusive",
            "duration_hours",
            "integration_5m_fraction",
        ]
    ].copy()
    events["episode_start"] = pd.to_datetime(events["episode_start"])
    events["episode_end_exclusive"] = pd.to_datetime(events["episode_end_exclusive"])
    min_start = events["episode_start"].min()
    max_end = events["episode_end_exclusive"].max()
    proximity = pd.read_parquet(
        proximity_path,
        columns=["animal_id", "timestamp", "dynamic_social_unit"],
    )
    proximity["timestamp"] = pd.to_datetime(proximity["timestamp"])
    proximity = proximity[
        proximity["timestamp"].between(min_start, max_end)
        & proximity["animal_id"].notna()
        & proximity["dynamic_social_unit"].notna()
    ].copy()
    proximity["animal_id"] = proximity["animal_id"].astype(str)
    proximity["dynamic_social_unit"] = proximity["dynamic_social_unit"].astype(str)

    vedba_cache: dict[str, pd.DataFrame | None] = {}

    def load_vedba(animal_id: str) -> pd.DataFrame | None:
        if animal_id not in vedba_cache:
            path = vedba_dir / f"{animal_id}.parquet"
            if not path.exists():
                vedba_cache[animal_id] = None
            else:
                data = pd.read_parquet(path, columns=["timestamp", "vedba", "logvedba"])
                data["timestamp"] = pd.to_datetime(data["timestamp"])
                vedba_cache[animal_id] = data
        return vedba_cache[animal_id]

    def summarize_group(event: pd.Series, group_col: str) -> dict[str, object]:
        group_name = str(event[group_col])
        present = proximity[
            proximity["dynamic_social_unit"].eq(group_name)
            & proximity["timestamp"].between(event["episode_start"], event["episode_end_exclusive"])
        ]
        animal_ids = sorted(present["animal_id"].unique())
        animal_summaries = []
        for animal_id in animal_ids:
            data = load_vedba(animal_id)
            if data is None:
                continue
            chunk = data[
                data["timestamp"].between(event["episode_start"], event["episode_end_exclusive"])
                & data["vedba"].notna()
            ]
            if chunk.empty:
                continue
            animal_summaries.append(
                {
                    "animal_id": animal_id,
                    "n_vedba_minutes": int(len(chunk)),
                    "mean_vedba": float(chunk["vedba"].mean()),
                    "median_vedba": float(chunk["vedba"].median()),
                    "sd_vedba": float(chunk["vedba"].std(ddof=0)),
                    "max_vedba": float(chunk["vedba"].max()),
                    "mean_logvedba": float(chunk["logvedba"].mean()),
                }
            )
        if not animal_summaries:
            return {
                f"{group_col}_activity_group": group_name,
                f"{group_col}_activity_animals_expected": int(len(animal_ids)),
                f"{group_col}_activity_animals_with_vedba": 0,
                f"{group_col}_activity_minutes": 0,
                f"{group_col}_vedba_mean": np.nan,
                f"{group_col}_vedba_median": np.nan,
                f"{group_col}_vedba_sd_between_animals": np.nan,
                f"{group_col}_vedba_cv_between_animals": np.nan,
                f"{group_col}_vedba_mean_within_animal_sd": np.nan,
                f"{group_col}_vedba_max": np.nan,
                f"{group_col}_logvedba_mean": np.nan,
                f"{group_col}_activity_animals": "",
            }
        animal_df = pd.DataFrame(animal_summaries)
        mean_activity = float(animal_df["mean_vedba"].mean())
        sd_between = float(animal_df["mean_vedba"].std(ddof=0))
        return {
            f"{group_col}_activity_group": group_name,
            f"{group_col}_activity_animals_expected": int(len(animal_ids)),
            f"{group_col}_activity_animals_with_vedba": int(len(animal_df)),
            f"{group_col}_activity_minutes": int(animal_df["n_vedba_minutes"].sum()),
            f"{group_col}_vedba_mean": mean_activity,
            f"{group_col}_vedba_median": float(animal_df["median_vedba"].mean()),
            f"{group_col}_vedba_sd_between_animals": sd_between,
            f"{group_col}_vedba_cv_between_animals": sd_between / mean_activity if mean_activity > 0 else np.nan,
            f"{group_col}_vedba_mean_within_animal_sd": float(animal_df["sd_vedba"].mean()),
            f"{group_col}_vedba_max": float(animal_df["max_vedba"].max()),
            f"{group_col}_logvedba_mean": float(animal_df["mean_logvedba"].mean()),
            f"{group_col}_activity_animals": ", ".join(animal_df["animal_id"]),
        }

    rows = []
    for _, event in events.iterrows():
        left = summarize_group(event, "group_a")
        right = summarize_group(event, "group_b")
        combined_means = [
            left.get("group_a_vedba_mean", np.nan),
            right.get("group_b_vedba_mean", np.nan),
        ]
        combined_sds = [
            left.get("group_a_vedba_sd_between_animals", np.nan),
            right.get("group_b_vedba_sd_between_animals", np.nan),
        ]
        rows.append(
            {
                "interaction_episode_id": event["interaction_episode_id"],
                "pair_key": event["pair_key"],
                "episode_start": event["episode_start"],
                "episode_end_exclusive": event["episode_end_exclusive"],
                "duration_hours": event["duration_hours"],
                "integration_5m_fraction": event["integration_5m_fraction"],
                **left,
                **right,
                "event_vedba_mean": float(np.nanmean(combined_means)),
                "event_vedba_abs_group_difference": float(
                    abs(left.get("group_a_vedba_mean", np.nan) - right.get("group_b_vedba_mean", np.nan))
                ),
                "event_vedba_mean_between_animal_sd": float(np.nanmean(combined_sds)),
            }
        )
    out = pd.DataFrame(rows)
    for col in [
        "event_vedba_mean",
        "event_vedba_abs_group_difference",
        "event_vedba_mean_between_animal_sd",
    ]:
        out[f"log1p_{col}"] = np.log1p(pd.to_numeric(out[col], errors="coerce"))
    return out


def summarize_event_gps_movement(event_rows: pd.DataFrame, proximity_path: Path) -> pd.DataFrame:
    if event_rows.empty:
        return pd.DataFrame()
    events = event_rows[
        [
            "interaction_episode_id",
            "pair_key",
            "group_a",
            "group_b",
            "episode_start",
            "episode_end_exclusive",
            "duration_hours",
            "integration_5m_fraction",
        ]
    ].copy()
    events["episode_start"] = pd.to_datetime(events["episode_start"])
    events["episode_end_exclusive"] = pd.to_datetime(events["episode_end_exclusive"])
    min_start = events["episode_start"].min()
    max_end = events["episode_end_exclusive"].max()
    proximity = pd.read_parquet(
        proximity_path,
        columns=["animal_id", "timestamp", "median_latitude", "median_longitude", "dynamic_social_unit"],
    )
    proximity["timestamp"] = pd.to_datetime(proximity["timestamp"])
    proximity = proximity[
        proximity["timestamp"].between(min_start, max_end)
        & proximity["animal_id"].notna()
        & proximity["dynamic_social_unit"].notna()
        & proximity["median_latitude"].notna()
        & proximity["median_longitude"].notna()
    ].copy()
    proximity["animal_id"] = proximity["animal_id"].astype(str)
    proximity["dynamic_social_unit"] = proximity["dynamic_social_unit"].astype(str)

    def summarize_group(event: pd.Series, group_col: str) -> dict[str, object]:
        group_name = str(event[group_col])
        present = proximity[
            proximity["dynamic_social_unit"].eq(group_name)
            & proximity["timestamp"].between(event["episode_start"], event["episode_end_exclusive"])
        ].copy()
        animal_summaries = []
        for animal_id, animal in present.sort_values("timestamp").groupby("animal_id", observed=True):
            animal = animal.drop_duplicates("timestamp").sort_values("timestamp")
            if len(animal) < 2:
                continue
            lat_prev = animal["median_latitude"].shift()
            lon_prev = animal["median_longitude"].shift()
            step_m = haversine_m(lat_prev, lon_prev, animal["median_latitude"], animal["median_longitude"])
            total_m = float(step_m.iloc[1:].sum())
            elapsed_h = float((animal["timestamp"].max() - animal["timestamp"].min()).total_seconds() / 3600.0)
            animal_summaries.append(
                {
                    "animal_id": str(animal_id),
                    "n_gps_fixes": int(len(animal)),
                    "gps_path_m": total_m,
                    "gps_m_per_h": total_m / elapsed_h if elapsed_h > 0 else np.nan,
                }
            )
        if not animal_summaries:
            return {
                f"{group_col}_gps_group": group_name,
                f"{group_col}_gps_animals_with_movement": 0,
                f"{group_col}_gps_movement_mean_m_per_h": np.nan,
                f"{group_col}_gps_movement_sd_between_animals": np.nan,
                f"{group_col}_gps_path_mean_m": np.nan,
                f"{group_col}_gps_animals": "",
            }
        animal_df = pd.DataFrame(animal_summaries)
        return {
            f"{group_col}_gps_group": group_name,
            f"{group_col}_gps_animals_with_movement": int(len(animal_df)),
            f"{group_col}_gps_movement_mean_m_per_h": float(animal_df["gps_m_per_h"].mean()),
            f"{group_col}_gps_movement_sd_between_animals": float(animal_df["gps_m_per_h"].std(ddof=0)),
            f"{group_col}_gps_path_mean_m": float(animal_df["gps_path_m"].mean()),
            f"{group_col}_gps_animals": ", ".join(animal_df["animal_id"]),
        }

    rows = []
    for _, event in events.iterrows():
        left = summarize_group(event, "group_a")
        right = summarize_group(event, "group_b")
        combined_means = [
            left.get("group_a_gps_movement_mean_m_per_h", np.nan),
            right.get("group_b_gps_movement_mean_m_per_h", np.nan),
        ]
        combined_sds = [
            left.get("group_a_gps_movement_sd_between_animals", np.nan),
            right.get("group_b_gps_movement_sd_between_animals", np.nan),
        ]
        rows.append(
            {
                "interaction_episode_id": event["interaction_episode_id"],
                "pair_key": event["pair_key"],
                "episode_start": event["episode_start"],
                "episode_end_exclusive": event["episode_end_exclusive"],
                "duration_hours": event["duration_hours"],
                "integration_5m_fraction": event["integration_5m_fraction"],
                **left,
                **right,
                "event_gps_movement_mean_m_per_h": float(np.nanmean(combined_means)),
                "event_gps_movement_abs_group_difference_m_per_h": float(
                    abs(
                        left.get("group_a_gps_movement_mean_m_per_h", np.nan)
                        - right.get("group_b_gps_movement_mean_m_per_h", np.nan)
                    )
                ),
                "event_gps_movement_mean_between_animal_sd_m_per_h": float(np.nanmean(combined_sds)),
            }
        )
    out = pd.DataFrame(rows)
    for col in [
        "event_gps_movement_mean_m_per_h",
        "event_gps_movement_abs_group_difference_m_per_h",
        "event_gps_movement_mean_between_animal_sd_m_per_h",
    ]:
        out[f"log1p_{col}"] = np.log1p(pd.to_numeric(out[col], errors="coerce"))
    return out


def plot_event_activity_summary(activity: pd.DataFrame, out_dir: Path) -> list[Path]:
    if activity.empty:
        return []
    paths = []
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4))
    specs = [
        ("event_vedba_mean", "Event mean VeDBA"),
        ("event_vedba_abs_group_difference", "Absolute group difference in mean VeDBA"),
        ("event_vedba_mean_between_animal_sd", "Mean between-animal SD within groups"),
    ]
    for ax, (col, label) in zip(axes, specs):
        values = pd.to_numeric(activity[col], errors="coerce").dropna()
        ax.hist(values, bins=min(18, max(6, int(np.sqrt(len(values))))), color="#4c78a8", alpha=0.78, edgecolor="white")
        ax.axvline(values.median(), color="#222222", linewidth=1.2, linestyle="--")
        ax.set_xlabel(label)
        ax.set_ylabel("Events")
        ax.grid(True, axis="y", alpha=0.18)
    fig.tight_layout()
    path = out_dir / "event_activity_vedba_distributions.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)

    plot_df = activity.dropna(subset=["event_vedba_mean", "event_vedba_mean_between_animal_sd", "integration_5m_fraction"]).copy()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    for ax, y_col, ylabel in [
        (axes[0], "event_vedba_mean", "Event mean VeDBA"),
        (axes[1], "event_vedba_mean_between_animal_sd", "Activity variation among animals"),
    ]:
        sizes = np.clip(np.sqrt(plot_df["duration_hours"]) * 18, 20, 130)
        ax.scatter(
            plot_df["integration_5m_fraction"],
            plot_df[y_col],
            s=sizes,
            alpha=0.45,
            color="#59a14f",
            edgecolor="white",
            linewidth=0.35,
        )
        corr = plot_df["integration_5m_fraction"].corr(plot_df[y_col], method="spearman")
        ax.text(
            0.03,
            0.94,
            f"Spearman r = {corr:.2f}",
            transform=ax.transAxes,
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.9},
        )
        ax.set_xlabel("5 m integration fraction")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.2)
    fig.tight_layout()
    path = out_dir / "event_activity_vs_integration.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_event_gps_movement_summary(gps_activity: pd.DataFrame, out_dir: Path) -> list[Path]:
    if gps_activity.empty:
        return []
    paths = []
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4))
    specs = [
        ("event_gps_movement_mean_m_per_h", "Event GPS movement (m/h)"),
        ("event_gps_movement_abs_group_difference_m_per_h", "Absolute group difference in movement"),
        ("event_gps_movement_mean_between_animal_sd_m_per_h", "Movement variation among animals"),
    ]
    for ax, (col, label) in zip(axes, specs):
        values = pd.to_numeric(gps_activity[col], errors="coerce").dropna()
        ax.hist(values, bins=min(18, max(6, int(np.sqrt(len(values))))), color="#f28e2b", alpha=0.78, edgecolor="white")
        ax.axvline(values.median(), color="#222222", linewidth=1.2, linestyle="--")
        ax.set_xlabel(label)
        ax.set_ylabel("Events")
        ax.grid(True, axis="y", alpha=0.18)
    fig.tight_layout()
    path = out_dir / "event_gps_movement_distributions.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)

    plot_df = gps_activity.dropna(
        subset=[
            "event_gps_movement_mean_m_per_h",
            "event_gps_movement_mean_between_animal_sd_m_per_h",
            "integration_5m_fraction",
        ]
    ).copy()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    for ax, y_col, ylabel in [
        (axes[0], "event_gps_movement_mean_m_per_h", "Event GPS movement (m/h)"),
        (axes[1], "event_gps_movement_mean_between_animal_sd_m_per_h", "GPS movement variation among animals"),
    ]:
        sizes = np.clip(np.sqrt(plot_df["duration_hours"]) * 18, 20, 130)
        ax.scatter(
            plot_df["integration_5m_fraction"],
            plot_df[y_col],
            s=sizes,
            alpha=0.45,
            color="#f28e2b",
            edgecolor="white",
            linewidth=0.35,
        )
        corr = plot_df["integration_5m_fraction"].corr(plot_df[y_col], method="spearman")
        ax.text(
            0.03,
            0.94,
            f"Spearman r = {corr:.2f}",
            transform=ax.transAxes,
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.9},
        )
        ax.set_xlabel("5 m integration fraction")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.2)
    fig.tight_layout()
    path = out_dir / "event_gps_movement_vs_integration.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)
    return paths


def run_integration_metric_analyses(
    event_rows: pd.DataFrame,
    out_dir: Path,
    gps_movement_summary: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    summaries = []
    coefficients = []
    model_rows = event_rows.copy()
    if gps_movement_summary is not None and not gps_movement_summary.empty:
        gps_cols = [
            "interaction_episode_id",
            "event_gps_movement_mean_m_per_h",
        ]
        model_rows = model_rows.merge(
            gps_movement_summary[gps_cols],
            on="interaction_episode_id",
            how="left",
        )
        model_rows["log1p_event_gps_movement_mean_m_per_h"] = np.log1p(
            pd.to_numeric(model_rows["event_gps_movement_mean_m_per_h"], errors="coerce")
        )
        model_rows["log1p_event_gps_movement_mean_m_per_h_z"] = zscore(
            model_rows["log1p_event_gps_movement_mean_m_per_h"]
        )
    for spec in INTEGRATION_METRIC_SPECS:
        if spec["metric"] not in model_rows.columns or spec["response_col"] not in model_rows.columns:
            continue
        subset = model_rows[model_rows[spec["metric"]].notna()].copy()
        variants = [("baseline", [])]
        if "log1p_event_gps_movement_mean_m_per_h_z" in subset.columns:
            variants.append(("gps_control", ["log1p_event_gps_movement_mean_m_per_h_z"]))
        for model_variant, extra_terms in variants:
            include_event_ndvi = bool(
                subset["dyad_ndvi_mean_1wk_before_after_z"].notna().all()
                and len(subset) >= MIN_MODEL_ROWS
            )
            result, formula, note, model_df = fit_event_bayesian_model(
                subset,
                response_col=spec["response_col"],
                include_ndvi=include_event_ndvi,
                model_kind="integration",
                distance_window="before_after",
                extra_terms=extra_terms,
            )
            model_name = f"integration_metric_{spec['metric']}_{model_variant}"
            summary = integration_metric_summary_row(
                subset,
                model_df,
                formula,
                note,
                spec,
                include_event_ndvi,
                model_variant,
            )
            coef = coefficient_table(result, model_name, spec["response_col"])
            coef["metric"] = spec["metric"]
            coef["metric_label"] = spec["label"]
            coef["model_variant"] = model_variant
            summaries.append(summary)
            coefficients.append(coef)
    summary_df = pd.concat(summaries, ignore_index=True)
    coef_df = pd.concat(coefficients, ignore_index=True)
    plots = [
        plot_integration_metric_distributions(model_rows, out_dir),
        plot_integration_metric_relationships(model_rows, out_dir),
        plot_integration_metric_effect_comparison(coef_df[coef_df["model_variant"].eq("baseline")], out_dir),
        plot_integration_gps_control_effect_comparison(coef_df, out_dir),
    ]
    return summary_df, coef_df, plots


def write_integration_metrics_dashboard(
    out_dir: Path,
    summary: pd.DataFrame,
    coef: pd.DataFrame,
    plots: list[Path],
    activity_summary: pd.DataFrame | None = None,
    activity_plots: list[Path] | None = None,
    gps_movement_summary: pd.DataFrame | None = None,
    gps_movement_plots: list[Path] | None = None,
) -> Path:
    metric_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "metric_label",
                "model_variant",
                "n_events",
                "n_dyads",
                "median_metric",
                "mean_metric",
                "min_metric",
                "max_metric",
                "ndvi_used",
                "n_events_with_ndvi",
                "interpretation",
            ]
        )
        + "</tr>"
        for _, row in summary.round(
            {
                "median_metric": 3,
                "mean_metric": 3,
                "min_metric": 3,
                "max_metric": 3,
            }
        ).iterrows()
    )
    effect_rows = fixed_effects_for_plot(coef, EVENT_EFFECT_LABELS)
    effect_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "metric_label",
                "model_variant",
                "effect_label",
                "estimate",
                "ci_low",
                "ci_high",
            ]
        )
        + "</tr>"
        for _, row in effect_rows.round({"estimate": 3, "ci_low": 3, "ci_high": 3}).iterrows()
    )
    images = "\n".join(f"<img src='{path.name}' alt='{path.stem}'>" for path in plots)
    activity_summary = pd.DataFrame() if activity_summary is None else activity_summary
    activity_plots = [] if activity_plots is None else activity_plots
    gps_movement_summary = pd.DataFrame() if gps_movement_summary is None else gps_movement_summary
    gps_movement_plots = [] if gps_movement_plots is None else gps_movement_plots
    activity_images = "\n".join(f"<img src='{path.name}' alt='{path.stem}'>" for path in activity_plots)
    gps_movement_images = "\n".join(f"<img src='{path.name}' alt='{path.stem}'>" for path in gps_movement_plots)
    activity_table = ""
    activity_note = "<p class='note'>No VeDBA activity summary was available.</p>"
    if not activity_summary.empty:
        activity_note = (
            f"<p class='note'>VeDBA activity was summarized for "
            f"{activity_summary['event_vedba_mean'].notna().sum()} events. "
            "Group-level activity is the average of per-animal mean VeDBA values during the event window. "
            "Activity variation is measured as the between-animal standard deviation of per-animal mean VeDBA within each group.</p>"
        )
        activity_table = "".join(
            "<tr>"
            + "".join(
                f"<td>{html.escape(str(row[col]))}</td>"
                for col in [
                    "interaction_episode_id",
                    "pair_key",
                    "duration_hours",
                    "group_a_vedba_mean",
                    "group_b_vedba_mean",
                    "event_vedba_mean",
                    "event_vedba_abs_group_difference",
                    "event_vedba_mean_between_animal_sd",
                    "group_a_activity_animals_with_vedba",
                    "group_b_activity_animals_with_vedba",
                ]
            )
            + "</tr>"
            for _, row in activity_summary.sort_values("event_vedba_mean", ascending=False)
            .head(20)
            .round(
                {
                    "duration_hours": 2,
                    "group_a_vedba_mean": 2,
                    "group_b_vedba_mean": 2,
                    "event_vedba_mean": 2,
                    "event_vedba_abs_group_difference": 2,
                    "event_vedba_mean_between_animal_sd": 2,
                }
            )
            .iterrows()
        )
    gps_movement_table = ""
    gps_movement_note = "<p class='note'>No GPS movement summary was available.</p>"
    if not gps_movement_summary.empty:
        gps_movement_note = (
            f"<p class='note'>GPS movement was summarized for "
            f"{gps_movement_summary['event_gps_movement_mean_m_per_h'].notna().sum()} events. "
            "Movement is estimated as per-animal path length inside the event window divided by the elapsed time covered by that animal's fixes. "
            "Group-level movement is the average of animal movement rates within the group.</p>"
        )
        gps_movement_table = "".join(
            "<tr>"
            + "".join(
                f"<td>{html.escape(str(row[col]))}</td>"
                for col in [
                    "interaction_episode_id",
                    "pair_key",
                    "duration_hours",
                    "group_a_gps_movement_mean_m_per_h",
                    "group_b_gps_movement_mean_m_per_h",
                    "event_gps_movement_mean_m_per_h",
                    "event_gps_movement_abs_group_difference_m_per_h",
                    "event_gps_movement_mean_between_animal_sd_m_per_h",
                    "group_a_gps_animals_with_movement",
                    "group_b_gps_animals_with_movement",
                ]
            )
            + "</tr>"
            for _, row in gps_movement_summary.sort_values("event_gps_movement_mean_m_per_h", ascending=False)
            .head(20)
            .round(
                {
                    "duration_hours": 2,
                    "group_a_gps_movement_mean_m_per_h": 1,
                    "group_b_gps_movement_mean_m_per_h": 1,
                    "event_gps_movement_mean_m_per_h": 1,
                    "event_gps_movement_abs_group_difference_m_per_h": 1,
                    "event_gps_movement_mean_between_animal_sd_m_per_h": 1,
                }
            )
            .iterrows()
        )
    out = out_dir / "integration_metrics_dashboard.html"
    out.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Integration metrics</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;color:#222;background:#fff}}
.wrap{{padding:24px 30px 44px;max-width:1420px}}
.note{{color:#555;max-width:1120px;line-height:1.45}}
table{{border-collapse:collapse;font-size:13px;margin:16px 0 26px}}
td,th{{border-bottom:1px solid #e5e5e5;padding:6px 9px;text-align:left;vertical-align:top}}
img{{max-width:100%;display:block;margin:20px 0 34px;border:1px solid #ddd;border-radius:6px}}
code{{background:#f4f4f4;padding:1px 4px;border-radius:3px}}
</style>
</head>
<body><div class="wrap">
<h1>Integration metrics only</h1>
<p class="note">This page compares event-level metrics that can represent integration after two groups meet. Each metric is modeled with the same hierarchical event structure: centroid distance in the week-before through week-after event window, prior meetings, Total group size, group-size difference, plus random intercepts for group 1, group 2, and dyad. The GPS-control model adds event-level GPS movement as a general activity control. NDVI is included only when it is available for all events in that metric; otherwise it is omitted so the model keeps all 108 interaction events where possible. Higher modularity means separation, so effects on modularity have the opposite biological direction from the other metrics.</p>
<p><a href="integration_metric_summary.csv">summary CSV</a> |
<a href="integration_metric_coefficients.csv">coefficients CSV</a> |
<a href="daily_interaction_hurdle_dashboard.html">full model dashboard</a></p>
<h2>Metric definitions</h2>
<table>
<tr><th>metric</th><th>model</th><th>events</th><th>dyads</th><th>median</th><th>mean</th><th>min</th><th>max</th><th>NDVI used</th><th>events with NDVI</th><th>meaning</th></tr>
{metric_table}
</table>
<h2>Plots</h2>
{images}
<h2>Event activity</h2>
{activity_note}
<p><a href="event_activity_vedba_summary.csv">event activity VeDBA summary CSV</a></p>
{activity_images}
<table>
<tr><th>event</th><th>dyad</th><th>duration h</th><th>group A mean</th><th>group B mean</th><th>event mean</th><th>group difference</th><th>activity variation</th><th>group A animals</th><th>group B animals</th></tr>
{activity_table}
</table>
<h2>GPS movement control</h2>
{gps_movement_note}
<p><a href="event_gps_movement_summary.csv">event GPS movement summary CSV</a></p>
{gps_movement_images}
<table>
<tr><th>event</th><th>dyad</th><th>duration h</th><th>group A m/h</th><th>group B m/h</th><th>event m/h</th><th>group difference</th><th>movement variation</th><th>group A animals</th><th>group B animals</th></tr>
{gps_movement_table}
</table>
<h2>Effect table</h2>
<table>
<tr><th>metric</th><th>model</th><th>effect</th><th>estimate</th><th>CI low</th><th>CI high</th></tr>
{effect_table}
</table>
</div></body></html>
""",
        encoding="utf-8",
    )
    return out


def build_event_distance_window_example(event_rows: pd.DataFrame, n_events: int = 8) -> pd.DataFrame:
    cols = [
        "pair_key",
        "group_a",
        "group_b",
        "episode_start",
        "episode_end",
        "duration_hours",
        "active_contact_hours",
        "prior_meetings",
        "centroid_distance_median_1wk_before_m",
        "centroid_distance_median_1wk_before_after_m",
        "centroid_distance_median_1wk_after_m",
        "group_size_total_1wk_before_after",
        "group_size_abs_diff_1wk_before_after",
        "dyad_ndvi_mean_1wk_before_after",
        "integration_5m_fraction",
    ]
    available = [col for col in cols if col in event_rows.columns]
    example = (
        event_rows.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["centroid_distance_median_1wk_before_after_m"])
        .nlargest(n_events, "centroid_distance_median_1wk_before_after_m")[available]
        .copy()
    )
    if example.empty:
        return example
    example["before_minus_after_m"] = (
        example["centroid_distance_median_1wk_before_m"]
        - example["centroid_distance_median_1wk_after_m"]
    )
    example["before_after_minus_before_m"] = (
        example["centroid_distance_median_1wk_before_after_m"]
        - example["centroid_distance_median_1wk_before_m"]
    )
    return example


def plot_event_distance_window_comparison(comparison: pd.DataFrame, out_dir: Path) -> Path:
    plot_df = comparison.dropna(subset=["distance_estimate", "distance_ci_low", "distance_ci_high"]).copy()
    plot_df["response_label"] = plot_df["response"].map(
        {
            "log1p_duration_hours": "Duration",
            "logit_integration_5m_fraction": "5 m integration",
        }
    )
    order = ["before", "before_after", "after"]
    labels = {"before": "Before", "before_after": "Before + after", "after": "After"}
    responses = [r for r in ["Duration", "5 m integration"] if r in set(plot_df["response_label"])]
    fig, axes = plt.subplots(1, max(1, len(responses)), figsize=(5.4 * max(1, len(responses)), 4.8), squeeze=False)
    colors = {"Duration": "#4c78a8", "5 m integration": "#59a14f"}
    for ax, response_label in zip(axes.ravel(), responses):
        sub = plot_df[plot_df["response_label"].eq(response_label)].copy()
        sub["distance_window"] = pd.Categorical(sub["distance_window"], categories=order, ordered=True)
        sub = sub.sort_values("distance_window")
        y = np.arange(len(sub))
        ax.axvline(0, color="#333333", linewidth=1.0, alpha=0.65)
        ax.errorbar(
            sub["distance_estimate"],
            y,
            xerr=[
                sub["distance_estimate"] - sub["distance_ci_low"],
                sub["distance_ci_high"] - sub["distance_estimate"],
            ],
            fmt="o",
            color=colors[response_label],
            ecolor=colors[response_label],
            elinewidth=2.0,
            capsize=4,
        )
        ax.set_yticks(y)
        ax.set_yticklabels([labels[value] for value in sub["distance_window"].astype(str)])
        ax.invert_yaxis()
        ax.set_xlabel("Standardized centroid-distance coefficient")
        ax.set_title(response_label)
        ax.grid(True, axis="x", alpha=0.22)
    fig.tight_layout()
    path = out_dir / "event_distance_window_effect_comparison.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def run_event_distance_window_comparison(
    event_rows: pd.DataFrame,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    rows = []
    coefs = []
    for response_col, model_kind in [
        ("log1p_duration_hours", "duration"),
        ("logit_integration_5m_fraction", "integration"),
    ]:
        subset = event_rows.copy()
        if model_kind == "integration":
            subset = subset[subset["integration_5m_fraction"].notna()].copy()
        include_event_ndvi = bool(subset["dyad_ndvi_mean_1wk_before_after_z"].notna().sum() >= MIN_MODEL_ROWS)
        for distance_window in ["before", "before_after", "after"]:
            result, formula, note, model_df = fit_event_bayesian_model(
                subset,
                response_col=response_col,
                include_ndvi=include_event_ndvi,
                model_kind=model_kind,
                distance_window=distance_window,
            )
            model_name = f"event_{model_kind}_{distance_window}"
            coef = coefficient_table(result, model_name, response_col)
            coef["distance_window"] = distance_window
            coefs.append(coef)
            distance_term = EVENT_DISTANCE_WINDOWS[distance_window]["term"]
            coef_index = coef.set_index("term")
            rows.append(
                {
                    "model": model_name,
                    "response": response_col,
                    "distance_window": distance_window,
                    "n_events": int(len(model_df)),
                    "n_dyads": int(model_df["pair_key"].nunique()),
                    "distance_estimate": float(coef_index.loc[distance_term, "estimate"]),
                    "distance_ci_low": float(coef_index.loc[distance_term, "ci_low"]),
                    "distance_ci_high": float(coef_index.loc[distance_term, "ci_high"]),
                    "total_size_estimate": float(coef_index.loc["log1p_group_size_total_1wk_before_after_z", "estimate"]),
                    "size_difference_estimate": float(coef_index.loc["log1p_group_size_abs_diff_1wk_before_after_z", "estimate"]),
                    "formula": formula,
                    "fit_note": note,
                }
            )
    comparison = pd.DataFrame(rows)
    plots = [plot_event_distance_window_comparison(comparison, out_dir)]
    return comparison, pd.concat(coefs, ignore_index=True), plots


def build_example_rows(rows: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "period_start",
        "pair_key",
        "daily_centroid_distance_m",
        "range_mean_centroid_distance_m",
        "range_median_centroid_distance_m",
        "range_overlap_days",
        "prior_pair_interaction_days",
        "any_interaction",
        "positive_duration_hours",
        "positive_2min_bins",
        "eligible_2min_bins",
        "total_cross_edges",
        "total_candidate_edges",
        "integration_5m_fraction",
        "edge_interaction_probability",
        "dyad_daily_ndvi_mean",
        "dyad_daily_ndvi_abs_diff",
    ]
    chunks = []
    positive = rows[rows["any_interaction"].eq(1)].nlargest(6, "positive_duration_hours").copy()
    positive["example_type"] = "positive dyad-week, longest durations"
    chunks.append(positive)

    eligible_zero = rows[
        rows["any_interaction"].eq(0) & rows["eligible_2min_bins"].gt(0)
    ].nsmallest(6, "range_mean_centroid_distance_m").copy()
    eligible_zero["example_type"] = "eligible 2-min rows, no cross-edge"
    chunks.append(eligible_zero)

    candidate_zero = rows[
        rows["any_interaction"].eq(0) & rows["eligible_2min_bins"].eq(0)
    ].nsmallest(6, "range_mean_centroid_distance_m").copy()
    candidate_zero["example_type"] = "centroid candidate only, no eligible 2-min rows"
    chunks.append(candidate_zero)

    out = pd.concat(chunks, ignore_index=True)
    return out[["example_type", *cols]].copy()


def plot_example_data(rows: pd.DataFrame, out_dir: Path) -> Path:
    plot_df = rows.copy()
    plot_df["plot_duration_hours"] = plot_df["positive_duration_hours"].where(
        plot_df["any_interaction"].eq(1), 0.0
    )
    rng = np.random.default_rng(7)
    sample_zero = plot_df[plot_df["any_interaction"].eq(0)].sample(
        min(2500, int(plot_df["any_interaction"].eq(0).sum())),
        random_state=7,
    )
    positive = plot_df[plot_df["any_interaction"].eq(1)]
    sample = pd.concat([sample_zero, positive], ignore_index=True)
    sample["jittered_duration"] = sample["plot_duration_hours"] + rng.normal(0, 0.035, len(sample))
    sample.loc[sample["any_interaction"].eq(0), "jittered_duration"] = rng.normal(
        0, 0.035, int(sample["any_interaction"].eq(0).sum())
    )

    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    zeros = sample[sample["any_interaction"].eq(0)]
    positives = sample[sample["any_interaction"].eq(1)]
    ax.scatter(
        zeros["range_mean_centroid_distance_m"],
        zeros["jittered_duration"],
        s=12,
        alpha=0.18,
        color="#9aa0a6",
        edgecolor="none",
        label="no interaction that week",
    )
    ax.scatter(
        positives["range_mean_centroid_distance_m"],
        positives["positive_duration_hours"],
        s=38,
        alpha=0.74,
        color="#2f7ed8",
        edgecolor="white",
        linewidth=0.45,
        label="positive dyad-week duration",
    )
    ax.axhline(0, color="#333333", linewidth=0.8, alpha=0.45)
    ax.set_xscale("log")
    ax.set_xlabel("Average daytime centroid distance for group pair (m, log scale)")
    ax.set_ylabel("Interaction input: 0/1 occurrence and positive duration (h)")
    ax.set_title("Example model input rows: distance paired with interaction outcome")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.22)
    fig.tight_layout()
    path = out_dir / "daily_interaction_examples_distance.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def lonlat_to_local_m(lon: pd.Series, lat: pd.Series, ref_lon: float, ref_lat: float) -> tuple[pd.Series, pd.Series]:
    x = (lon.astype(float) - ref_lon) * 111_320.0 * np.cos(np.radians(ref_lat))
    y = (lat.astype(float) - ref_lat) * 110_540.0
    return x, y


def safe_slug(value: object) -> str:
    text = str(value).strip().lower()
    slug = "".join(char if char.isalnum() else "_" for char in text)
    return "_".join(part for part in slug.split("_") if part)


def plot_specific_event(
    rows: pd.DataFrame,
    raw_path: Path,
    proximity_path: Path,
    out_dir: Path,
    event: pd.Series | None = None,
    output_prefix: str = "specific_interaction_event",
) -> tuple[Path, pd.DataFrame]:
    if event is None:
        event = rows[rows["any_interaction"].eq(1)].nlargest(1, "positive_duration_hours").iloc[0]
    event_day = pd.to_datetime(event["period_start"]).floor("D")
    pair_key = str(event["pair_key"])
    group_a = str(event["group_a"])
    group_b = str(event["group_b"])

    raw = pd.read_csv(raw_path, parse_dates=["timestamp", "bin_2min"])
    metric = raw[
        raw["pair_key"].eq(pair_key)
        & raw["bin_2min"].dt.floor("D").eq(event_day)
        & in_daytime_window(raw["bin_2min"])
    ].copy()
    if metric.empty:
        raise ValueError(f"No 2-min metrics found for {pair_key} on {event_day.date()}.")
    metric["positive_bin"] = metric["cross_edges"].gt(0).astype(int)
    metric["cross_edge_fraction"] = np.where(
        metric["total_edges"].gt(0),
        metric["cross_edges"] / metric["total_edges"],
        np.nan,
    )

    hourly = (
        metric.groupby("timestamp", observed=True)
        .agg(
            bins=("bin_2min", "size"),
            positive_bins=("positive_bin", "sum"),
            cross_edges=("cross_edges", "sum"),
            total_edges=("total_edges", "sum"),
            mean_cross_edge_fraction=("cross_edge_fraction", "mean"),
            mean_pair_n=("pair_n", "mean"),
            mean_cluster_size_total=("cluster_size_total", "mean"),
            cluster_groups=("cluster_groups", lambda s: " | ".join(sorted(set(map(str, s))))),
        )
        .reset_index()
    )
    hourly["positive_duration_hours"] = hourly["positive_bins"] * (2.0 / 60.0)
    available_hours = hourly["timestamp"].sort_values().tolist()
    candidate_snapshot_hours = [
        available_hours[0],
        available_hours[len(available_hours) // 2],
        available_hours[-1],
    ]
    snapshot_hours = []
    for ts in candidate_snapshot_hours:
        if ts not in snapshot_hours:
            snapshot_hours.append(ts)
    snapshot_hours = sorted(snapshot_hours)

    prox_cols = [
        "animal_id",
        "timestamp",
        "median_longitude",
        "median_latitude",
        "dbscan_cluster",
        "dynamic_social_unit",
    ]
    prox = pd.read_parquet(proximity_path, columns=prox_cols)
    prox["timestamp"] = pd.to_datetime(prox["timestamp"])
    prox_day = prox[prox["timestamp"].isin(snapshot_hours)].dropna(
        subset=["median_longitude", "median_latitude", "dbscan_cluster", "dynamic_social_unit"]
    ).copy()

    snapshot_rows = []
    for ts in snapshot_hours:
        hour_rows = prox_day[prox_day["timestamp"].eq(ts)].copy()
        pair_clusters = hour_rows.loc[
            hour_rows["dynamic_social_unit"].isin([group_a, group_b]), "dbscan_cluster"
        ].dropna().unique()
        cluster_rows = hour_rows[hour_rows["dbscan_cluster"].isin(pair_clusters)].copy()
        snapshot_rows.append(cluster_rows)
    spatial = pd.concat(snapshot_rows, ignore_index=True) if snapshot_rows else pd.DataFrame()
    if spatial.empty:
        raise ValueError(f"No hourly spatial rows found for {pair_key} on {event_day.date()}.")

    ref_lon = float(spatial["median_longitude"].mean())
    ref_lat = float(spatial["median_latitude"].mean())
    spatial["x_m"], spatial["y_m"] = lonlat_to_local_m(
        spatial["median_longitude"], spatial["median_latitude"], ref_lon, ref_lat
    )

    colors = {
        group_a: "#e15759",
        group_b: "#4e79a7",
        "Lilac": "#b07aa1",
    }
    fallback_colors = ["#59a14f", "#f28e2b", "#76b7b2", "#edc948", "#bab0ac"]

    fig = plt.figure(figsize=(13.6, 8.2))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.05, 1.0], hspace=0.34, wspace=0.24)
    fig.suptitle(
        f"Specific interaction event: {pair_key} on {event_day.date()}",
        fontsize=16,
        y=0.98,
    )

    for i, ts in enumerate(snapshot_hours):
        ax = fig.add_subplot(grid[0, i])
        sub = spatial[spatial["timestamp"].eq(ts)]
        for j, (unit, unit_rows) in enumerate(sub.groupby("dynamic_social_unit", observed=True)):
            color = colors.get(str(unit), fallback_colors[j % len(fallback_colors)])
            ax.scatter(
                unit_rows["x_m"],
                unit_rows["y_m"],
                s=46,
                color=color,
                alpha=0.86,
                edgecolor="white",
                linewidth=0.55,
                label=str(unit),
            )
        centroids = (
            sub[sub["dynamic_social_unit"].isin([group_a, group_b])]
            .groupby("dynamic_social_unit", observed=True)
            .agg(x_m=("x_m", "mean"), y_m=("y_m", "mean"), n=("animal_id", "nunique"))
            .reset_index()
        )
        if len(centroids) == 2:
            ax.plot(centroids["x_m"], centroids["y_m"], color="#222222", linewidth=1.4, alpha=0.75)
            for _, centroid in centroids.iterrows():
                ax.scatter(
                    centroid["x_m"],
                    centroid["y_m"],
                    marker="X",
                    s=120,
                    color=colors.get(str(centroid["dynamic_social_unit"]), "#333333"),
                    edgecolor="#111111",
                    linewidth=0.7,
                    zorder=5,
                )
        cluster_ids = ", ".join(map(str, sorted(sub["dbscan_cluster"].dropna().astype(int).unique())))
        ax.set_title(f"{pd.to_datetime(ts).strftime('%H:%M')} cluster {cluster_ids}")
        ax.set_xlabel("east-west position (m)")
        if i == 0:
            ax.set_ylabel("north-south position (m)")
        ax.axhline(0, color="#cccccc", linewidth=0.7, zorder=0)
        ax.axvline(0, color="#cccccc", linewidth=0.7, zorder=0)
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(True, alpha=0.18)
        if i == len(snapshot_hours) - 1:
            ax.legend(frameon=False, fontsize=8, loc="best")

    ax1 = fig.add_subplot(grid[1, :])
    ax1.plot(metric["bin_2min"], metric["cross_edges"], color="#4e79a7", linewidth=1.3, label="cross-group 5 m edges")
    ax1.fill_between(
        metric["bin_2min"],
        0,
        metric["positive_bin"] * max(float(metric["cross_edges"].max()), 1.0),
        color="#4e79a7",
        alpha=0.08,
        step="mid",
        label="positive 2-min bin",
    )
    ax1.set_ylabel("cross edges per 2-min bin")
    ax1.set_xlabel("time")
    ax1.grid(True, alpha=0.2)
    ax2 = ax1.twinx()
    ax2.plot(metric["bin_2min"], metric["cross_edge_fraction"], color="#f28e2b", linewidth=1.2, alpha=0.85, label="cross-edge fraction")
    ax2.set_ylabel("cross-edge fraction")
    ax2.set_ylim(-0.02, min(1.02, max(0.2, float(metric["cross_edge_fraction"].max()) + 0.12)))

    summary_text = (
        f"average range distance: {event['range_mean_centroid_distance_m']:.1f} m\n"
        f"event-day centroid distance: {event['daily_centroid_distance_m']:.1f} m\n"
        f"positive duration: {event['positive_duration_hours']:.2f} h "
        f"({int(event['positive_2min_bins'])}/{int(event['eligible_2min_bins'])} eligible 2-min bins)\n"
        f"cross edges: {int(event['total_cross_edges'])}/{int(event['total_candidate_edges'])} "
        f"({event['edge_interaction_probability']:.3f})\n"
        f"hourly cluster composition: {hourly['cluster_groups'].mode().iloc[0]}"
    )
    ax1.text(
        0.01,
        0.98,
        summary_text,
        transform=ax1.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.92, "boxstyle": "round,pad=0.45"},
    )
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, frameon=False, loc="upper right")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    path = out_dir / f"{output_prefix}_{safe_slug(pair_key)}_{event_day.date()}.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)

    detail = pd.DataFrame(
        [
            {
                "period_start": event_day.date().isoformat(),
                "pair_key": pair_key,
                "group_a": group_a,
                "group_b": group_b,
                "daily_centroid_distance_m": event["daily_centroid_distance_m"],
                "range_mean_centroid_distance_m": event["range_mean_centroid_distance_m"],
                "range_median_centroid_distance_m": event["range_median_centroid_distance_m"],
                "range_overlap_days": event["range_overlap_days"],
                "prior_pair_interaction_days": event["prior_pair_interaction_days"],
                "prior_pair_interaction_duration_hours": event["prior_pair_interaction_duration_hours"],
                "positive_duration_hours": event["positive_duration_hours"],
                "positive_2min_bins": event["positive_2min_bins"],
                "eligible_2min_bins": event["eligible_2min_bins"],
                "total_cross_edges": event["total_cross_edges"],
                "total_candidate_edges": event["total_candidate_edges"],
                "edge_interaction_probability": event["edge_interaction_probability"],
                "integration_5m_fraction": event["integration_5m_fraction"],
                "first_2min_bin": metric["bin_2min"].min(),
                "last_2min_bin": metric["bin_2min"].max(),
                "snapshot_hours": ", ".join(pd.to_datetime(snapshot_hours).strftime("%H:%M")),
                "cluster_groups": hourly["cluster_groups"].mode().iloc[0],
            }
        ]
    )
    return path, detail


def plot_long_distance_events(
    rows: pd.DataFrame,
    raw_path: Path,
    proximity_path: Path,
    out_dir: Path,
    n_events: int = 5,
) -> tuple[list[Path], pd.DataFrame]:
    events = rows[rows["any_interaction"].eq(1)].nlargest(n_events, "range_mean_centroid_distance_m")
    paths: list[Path] = []
    details = []
    for rank, (_, event) in enumerate(events.iterrows(), start=1):
        path, detail = plot_specific_event(
            rows,
            raw_path,
            proximity_path,
            out_dir,
            event=event,
            output_prefix=f"long_centroid_distance_event_{rank:02d}",
        )
        detail.insert(0, "rank", rank)
        paths.append(path)
        details.append(detail)
    return paths, pd.concat(details, ignore_index=True) if details else pd.DataFrame()


def write_dashboard(
    out_dir: Path,
    summary: pd.DataFrame,
    meeting_probability_summary: pd.DataFrame,
    event_duration_summary: pd.DataFrame,
    event_integration_summary: pd.DataFrame,
    image_paths: list[Path],
    meeting_probability_plots: list[Path],
    event_duration_plots: list[Path],
    event_integration_plots: list[Path],
    event_effect_plots: list[Path],
    example_rows: pd.DataFrame,
    example_plot: Path,
    event_plot: Path,
    event_detail: pd.DataFrame,
    long_event_plots: list[Path],
    long_event_detail: pd.DataFrame,
    event_distance_window_comparison: pd.DataFrame,
    event_distance_window_example: pd.DataFrame,
    model_rows: pd.DataFrame,
    episodes: pd.DataFrame,
    include_ndvi: bool,
) -> Path:
    table_rows = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "model",
                "response",
                "n_dyad_days",
                "n_dyads",
                "n_positive_dyad_days",
                "overall_daily_interaction_probability",
                "overall_5m_integration_fraction",
                "n_eligible_5m_dyad_days",
                "median_positive_duration_hours",
                "fit_note",
            ]
        )
        + "</tr>"
        for _, row in summary.round(
            {
                "overall_daily_interaction_probability": 4,
                "overall_5m_integration_fraction": 4,
                "median_positive_duration_hours": 3,
            }
        ).iterrows()
    )
    images = "\n".join(f"<img src='{path.name}' alt='{path.stem}'>" for path in image_paths)
    meeting_probability_images = "\n".join(
        f"<img src='{path.name}' alt='{path.stem}'>" for path in meeting_probability_plots
    )
    event_duration_images = "\n".join(
        f"<img src='{path.name}' alt='{path.stem}'>" for path in event_duration_plots
    )
    event_integration_images = "\n".join(
        f"<img src='{path.name}' alt='{path.stem}'>" for path in event_integration_plots
    )
    event_effect_images = "\n".join(
        f"<img src='{path.name}' alt='{path.stem}'>" for path in event_effect_plots
    )
    event_duration_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "merge_size_class",
                "n_events",
                "n_dyads",
                "median_duration_hours",
                "max_duration_hours",
                "max_duration_days",
                "median_active_contact_hours",
                "median_prior_meetings",
                "fit_note",
            ]
        )
        + "</tr>"
        for _, row in event_duration_summary.round(
            {
                "median_duration_hours": 3,
                "max_duration_hours": 2,
                "max_duration_days": 2,
                "median_active_contact_hours": 3,
                "median_prior_meetings": 1,
            }
        ).iterrows()
    )
    meeting_probability_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "time_unit",
                "model_type",
                "n_dyad_intervals",
                "n_dyads",
                "n_meeting_intervals",
                "meeting_probability",
                "median_prior_meetings",
                "fit_note",
            ]
        )
        + "</tr>"
        for _, row in meeting_probability_summary.round(
            {
                "meeting_probability": 4,
                "median_prior_meetings": 1,
            }
        ).iterrows()
    )
    event_integration_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "merge_size_class",
                "n_events",
                "n_dyads",
                "median_integration_5m_fraction",
                "mean_integration_5m_fraction",
                "median_cross_group_5m_edges",
                "median_candidate_5m_edges",
                "median_duration_hours",
                "median_prior_meetings",
                "fit_note",
            ]
        )
        + "</tr>"
        for _, row in event_integration_summary.round(
            {
                "median_integration_5m_fraction": 3,
                "mean_integration_5m_fraction": 3,
                "median_cross_group_5m_edges": 1,
                "median_candidate_5m_edges": 1,
                "median_duration_hours": 3,
                "median_prior_meetings": 1,
            }
        ).iterrows()
    )
    example_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "example_type",
                "period_start",
                "pair_key",
                "range_mean_centroid_distance_m",
                "daily_centroid_distance_m",
                "prior_pair_interaction_days",
                "any_interaction",
                "positive_duration_hours",
                "positive_2min_bins",
                "eligible_2min_bins",
                "total_cross_edges",
                "total_candidate_edges",
            ]
        )
        + "</tr>"
        for _, row in example_rows.round(
            {
                "range_mean_centroid_distance_m": 2,
                "daily_centroid_distance_m": 2,
                "positive_duration_hours": 2,
                "edge_interaction_probability": 4,
            }
        ).iterrows()
    )
    event = event_detail.iloc[0]
    long_event_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "rank",
                "period_start",
                "pair_key",
                "range_mean_centroid_distance_m",
                "daily_centroid_distance_m",
                "prior_pair_interaction_days",
                "positive_duration_hours",
                "positive_2min_bins",
                "eligible_2min_bins",
                "integration_5m_fraction",
                "edge_interaction_probability",
                "cluster_groups",
            ]
        )
        + "</tr>"
        for _, row in long_event_detail.round(
            {
                "range_mean_centroid_distance_m": 1,
                "daily_centroid_distance_m": 1,
                "integration_5m_fraction": 3,
                "positive_duration_hours": 2,
                "edge_interaction_probability": 3,
            }
        ).iterrows()
    )
    long_event_images = "\n".join(
        f"<img src='{path.name}' alt='{path.stem}'>" for path in long_event_plots
    )
    distance_window_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "response",
                "distance_window",
                "n_events",
                "n_dyads",
                "distance_estimate",
                "distance_ci_low",
                "distance_ci_high",
                "total_size_estimate",
                "size_difference_estimate",
            ]
        )
        + "</tr>"
        for _, row in event_distance_window_comparison.round(
            {
                "distance_estimate": 3,
                "distance_ci_low": 3,
                "distance_ci_high": 3,
                "total_size_estimate": 3,
                "size_difference_estimate": 3,
            }
        ).iterrows()
    )
    distance_window_example_table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "pair_key",
                "episode_start",
                "duration_hours",
                "centroid_distance_median_1wk_before_m",
                "centroid_distance_median_1wk_before_after_m",
                "centroid_distance_median_1wk_after_m",
                "before_minus_after_m",
                "group_size_total_1wk_before_after",
                "group_size_abs_diff_1wk_before_after",
            ]
            if col in event_distance_window_example.columns
        )
        + "</tr>"
        for _, row in event_distance_window_example.round(
            {
                "duration_hours": 2,
                "centroid_distance_median_1wk_before_m": 1,
                "centroid_distance_median_1wk_before_after_m": 1,
                "centroid_distance_median_1wk_after_m": 1,
                "before_minus_after_m": 1,
                "group_size_total_1wk_before_after": 1,
                "group_size_abs_diff_1wk_before_after": 1,
            }
        ).iterrows()
    )
    episode_note = "No positive episodes were detected."
    if not episodes.empty:
        episode_note = (
            f"{len(episodes):,} positive interaction episodes; "
            f"median episode duration {episodes['duration_hours'].median():.2f} h; "
            f"max episode duration {episodes['duration_hours'].max():.2f} h."
        )
    ndvi_note = "NDVI terms included." if include_ndvi else "NDVI terms not included."
    out = out_dir / "daily_interaction_hurdle_dashboard.html"
    out.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Bayesian meeting models</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;color:#222;background:#fff}}
.wrap{{padding:24px 30px 44px;max-width:1420px}}
.note{{color:#555;max-width:1120px;line-height:1.45}}
table{{border-collapse:collapse;font-size:13px;margin:16px 0 26px}}
td,th{{border-bottom:1px solid #e5e5e5;padding:6px 9px;text-align:left;vertical-align:top}}
img{{max-width:100%;display:block;margin:20px 0 34px;border:1px solid #ddd;border-radius:6px}}
code{{background:#f4f4f4;padding:1px 4px;border-radius:3px}}
</style>
</head>
<body><div class="wrap">
<h1>Bayesian meeting models</h1>
<p class="note">The meeting-probability models use all dyad intervals, including intervals with no detected meeting. The duration and integration models use only dyads that interacted, and all merge sizes are fit together. Each meeting-event row has a real duration measured as elapsed time from the first positive 2-minute bin to the last positive 2-minute bin plus one bin width. Positive bins separated by up to {EVENT_MAX_GAP_HOURS:.0f} hours are treated as the same meeting, so an event can last from minutes to multiple days when interaction evidence continues across daytime windows. The duration model is a linear Bayesian Gaussian mixed model on <code>log1p(duration_hours)</code>. The integration model is a linear Bayesian Gaussian mixed model on <code>logit(total cross-group 5 m links / total candidate 5 m links)</code> for the same event episodes. Fixed effects are centroid distance, previous meetings, Total group size, absolute group-size difference, and mean NDVI. Random intercepts in the hierarchical versions are included for group 1, group 2, and dyad. Weather and NDVI difference/change are not included.</p>
<p><a href="event_duration_model_rows.csv">event-duration rows CSV</a> |
<a href="meeting_probability_predictions.csv">meeting-probability predictions CSV</a> |
<a href="meeting_probability_coefficients.csv">meeting-probability coefficients CSV</a> |
<a href="meeting_probability_summary.csv">meeting-probability summary CSV</a> |
<a href="event_duration_predictions.csv">event-duration predictions CSV</a> |
<a href="event_duration_coefficients.csv">event-duration coefficients CSV</a> |
<a href="event_duration_summary.csv">event-duration summary CSV</a> |
<a href="event_integration_predictions.csv">event-integration predictions CSV</a> |
<a href="event_integration_coefficients.csv">event-integration coefficients CSV</a> |
<a href="event_integration_summary.csv">event-integration summary CSV</a> |
<a href="event_distance_window_comparison.csv">distance-window comparison CSV</a> |
<a href="event_distance_window_example.csv">distance-window example CSV</a> |
<a href="integration_metrics_dashboard.html">integration metrics only</a></p>
<h2>Effect sizes</h2>
<p class="note">These forest plots summarize the standardized fixed-effect coefficients. Points are posterior means and horizontal bars are 95% posterior intervals. The three-part hierarchical comparison uses the hierarchical meeting-probability model, then duration and integration conditional on a detected meeting.</p>
{event_effect_images}
<h2>Probability of meeting</h2>
<p class="note">The probability step is fit twice: once on dyad-days and once on dyad-weeks, so we can see whether the time unit changes the story. The separate version has fixed effects only. The hierarchical version adds random intercepts for group 1, group 2, and dyad.</p>
<table>
<tr><th>time unit</th><th>model type</th><th>dyad intervals</th><th>dyads</th><th>meeting intervals</th><th>meeting probability</th><th>median prior meetings</th><th>fit note</th></tr>
{meeting_probability_table}
</table>
{meeting_probability_images}
<h2>Hierarchical three-part structure</h2>
<p class="note">The three-part structure follows the biology of the process. First, two groups either meet or they do not, so the probability model uses all dyad intervals and asks what makes a meeting happen. Second, duration is only meaningful after a meeting exists, so the duration model is conditional on detected meeting events. Third, 5 m integration is a property of how mixed the animals were during that event, so it is modeled conditional on the same events. The shared group and dyad random intercepts separate general group/dyad tendencies from the effects of distance, history, Total group size, size imbalance, and NDVI.</p>
<table>
<tr><th>step</th><th>row unit</th><th>response</th><th>meaning</th></tr>
<tr><td>1</td><td>dyad-week</td><td><code>Pr(any meeting)</code></td><td>Whether the two groups meet at least once in that week.</td></tr>
<tr><td>2</td><td>meeting event</td><td><code>log1p(duration_hours)</code></td><td>How long the detected meeting episode lasts, conditional on a meeting.</td></tr>
<tr><td>3</td><td>meeting event</td><td><code>logit(5 m integration)</code></td><td>How mixed the event is at 5 m, conditional on a meeting.</td></tr>
</table>
<h2>Duration</h2>
<table>
<tr><th>merge scale</th><th>events</th><th>dyads</th><th>median real duration h</th><th>max real duration h</th><th>max real duration d</th><th>median active contact h</th><th>median prior meetings</th><th>fit note</th></tr>
{event_duration_table}
</table>
{event_duration_images}
<h2>5 m integration</h2>
<table>
<tr><th>merge scale</th><th>events</th><th>dyads</th><th>median integration</th><th>mean integration</th><th>median cross-group 5 m links</th><th>median candidate 5 m links</th><th>median real duration h</th><th>median prior meetings</th><th>fit note</th></tr>
{event_integration_table}
</table>
{event_integration_images}
<h2>Distance-window comparison</h2>
<p class="note">This checks whether the centroid-distance predictor behaves differently when it is measured only in the week before the event, in the combined week-before through week-after window, or only in the week after the event. The main duration and integration plots above use the before-plus-after window.</p>
<table>
<tr><th>response</th><th>distance window</th><th>events</th><th>dyads</th><th>distance estimate</th><th>CI low</th><th>CI high</th><th>total size estimate</th><th>size-difference estimate</th></tr>
{distance_window_table}
</table>
<h3>Example event rows</h3>
<table>
<tr><th>dyad</th><th>episode start</th><th>duration h</th><th>before distance m</th><th>before+after distance m</th><th>after distance m</th><th>before-after difference m</th><th>total group size</th><th>size difference</th></tr>
{distance_window_example_table}
</table>
</div></body></html>
""",
        encoding="utf-8",
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit weekly two-part intergroup interaction hurdle models.")
    parser.add_argument("--daily-dir", type=Path, default=DEFAULT_DAILY_DIR)
    parser.add_argument("--raw-interactions", type=Path, default=DEFAULT_RAW_INTERACTIONS)
    parser.add_argument("--proximity-status", type=Path, default=DEFAULT_PROXIMITY_STATUS)
    parser.add_argument("--vedba-dir", type=Path, default=DEFAULT_VEDBA_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--size-source",
        choices=["collars", "demographic"],
        default="collars",
        help=(
            "Which size variable enters the models. 'collars' (default) sums the observed "
            "collar counts and reproduces every run before 2026-09-03. 'demographic' uses "
            "independent group sizes from EAS and additionally controls for collar count, "
            "so a size effect is not just sampling effort. Collar coverage ranges 3%%-42%% "
            "across groups, so the two are not interchangeable."
        ),
    )
    args = parser.parse_args()

    global SIZE_SOURCE
    SIZE_SOURCE = args.size_source
    print(f"size covariate source: {SIZE_SOURCE}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    daily_rows, episodes = build_model_rows(args.daily_dir, args.raw_interactions, args.proximity_status)
    rows = aggregate_to_weekly_rows(daily_rows)
    include_ndvi = bool(
        rows[["dyad_daily_ndvi_mean_z", "dyad_daily_ndvi_abs_diff_z"]]
        .notna()
        .all(axis=1)
        .sum()
        >= MIN_MODEL_ROWS
    )
    (
        weekly_meeting_probability_summary,
        weekly_meeting_probability_coef,
        weekly_meeting_probability_predictions,
        weekly_meeting_probability_plots,
    ) = run_meeting_probability_analyses(
        rows,
        args.out_dir,
        include_ndvi=include_ndvi,
        time_unit="weekly",
    )
    include_daily_ndvi = bool(
        daily_rows[["dyad_daily_ndvi_mean_z", "dyad_daily_ndvi_abs_diff_z"]]
        .notna()
        .all(axis=1)
        .sum()
        >= MIN_MODEL_ROWS
    )
    (
        daily_meeting_probability_summary,
        daily_meeting_probability_coef,
        daily_meeting_probability_predictions,
        daily_meeting_probability_plots,
    ) = run_meeting_probability_analyses(
        daily_rows,
        args.out_dir,
        include_ndvi=include_daily_ndvi,
        time_unit="daily",
    )
    meeting_probability_summary = pd.concat(
        [weekly_meeting_probability_summary, daily_meeting_probability_summary],
        ignore_index=True,
    )
    meeting_probability_coef = pd.concat(
        [weekly_meeting_probability_coef, daily_meeting_probability_coef],
        ignore_index=True,
    )
    meeting_probability_predictions = pd.concat(
        [weekly_meeting_probability_predictions, daily_meeting_probability_predictions],
        ignore_index=True,
    )
    meeting_probability_plots = weekly_meeting_probability_plots + daily_meeting_probability_plots
    event_duration_rows = build_event_duration_rows(episodes, daily_rows)
    (
        event_duration_summary,
        event_duration_coef,
        event_duration_predictions,
        event_duration_plots,
        event_duration_ndvi_used,
    ) = run_event_duration_analysis(
        event_duration_rows,
        args.out_dir,
        "all",
        "event_duration",
        "All merge events: duration",
    )
    (
        event_integration_summary,
        event_integration_coef,
        event_integration_predictions,
        event_integration_plots,
        event_integration_ndvi_used,
    ) = run_event_integration_analysis(
        event_duration_rows,
        args.out_dir,
    )
    (
        event_distance_window_comparison,
        event_distance_window_coef,
        event_distance_window_plots,
    ) = run_event_distance_window_comparison(event_duration_rows, args.out_dir)
    event_distance_window_example = build_event_distance_window_example(event_duration_rows)
    event_gps_movement_summary = summarize_event_gps_movement(
        event_duration_rows,
        args.proximity_status,
    )
    event_gps_movement_plots = plot_event_gps_movement_summary(event_gps_movement_summary, args.out_dir)
    (
        integration_metric_summary,
        integration_metric_coef,
        integration_metric_plots,
    ) = run_integration_metric_analyses(
        event_duration_rows,
        args.out_dir,
        gps_movement_summary=event_gps_movement_summary,
    )
    event_activity_summary = summarize_event_activity(
        event_duration_rows,
        args.proximity_status,
        args.vedba_dir,
    )
    event_activity_plots = plot_event_activity_summary(event_activity_summary, args.out_dir)
    integration_metrics_dashboard = write_integration_metrics_dashboard(
        args.out_dir,
        integration_metric_summary,
        integration_metric_coef,
        integration_metric_plots,
        activity_summary=event_activity_summary,
        activity_plots=event_activity_plots,
        gps_movement_summary=event_gps_movement_summary,
        gps_movement_plots=event_gps_movement_plots,
    )

    fit_specs = [
        ("probability", "any_interaction", rows, Binomial()),
        ("duration_positive_weeks", "log1p_positive_duration_hours", rows[rows["any_interaction"].eq(1)].copy(), Gaussian()),
        ("integration_5m", "integration_5m_fraction", rows[rows["total_candidate_edges"].gt(0)].copy(), Binomial()),
    ]
    summaries = []
    coefficients = []
    predictions = []
    for model_name, response, model_rows, family in fit_specs:
        result, formula, note = fit_gee(model_rows, response, include_ndvi=include_ndvi, family=family)
        summaries.append(summary_row(rows, model_name, response, formula, note, result))
        coefficients.append(coefficient_table(result, model_name, response))
        for predictor in ["distance", "history", "ndvi_mean", "ndvi_abs_diff"]:
            grid = prediction_grid(
                rows,
                result,
                model_name,
                response,
                include_ndvi=include_ndvi,
                predictor=predictor,
            )
            if not grid.empty:
                predictions.append(grid)

    summary = pd.DataFrame(summaries)
    coef = pd.concat(coefficients, ignore_index=True)
    pred = pd.concat(predictions, ignore_index=True)

    rows.to_csv(args.out_dir / "daily_interaction_hurdle_model_rows.csv", index=False)
    episodes.to_csv(args.out_dir / "daily_interaction_positive_episodes.csv", index=False)
    summary.to_csv(args.out_dir / "daily_interaction_hurdle_summary.csv", index=False)
    coef.to_csv(args.out_dir / "daily_interaction_hurdle_coefficients.csv", index=False)
    pred.to_csv(args.out_dir / "daily_interaction_hurdle_predictions.csv", index=False)
    daily_rows.to_csv(args.out_dir / "daily_interaction_hurdle_daily_event_rows.csv", index=False)
    meeting_probability_summary.to_csv(args.out_dir / "meeting_probability_summary.csv", index=False)
    meeting_probability_coef.to_csv(args.out_dir / "meeting_probability_coefficients.csv", index=False)
    meeting_probability_predictions.to_csv(args.out_dir / "meeting_probability_predictions.csv", index=False)
    event_duration_rows.to_csv(args.out_dir / "event_duration_model_rows.csv", index=False)
    event_duration_summary.to_csv(args.out_dir / "event_duration_summary.csv", index=False)
    event_duration_coef.to_csv(args.out_dir / "event_duration_coefficients.csv", index=False)
    event_duration_predictions.to_csv(args.out_dir / "event_duration_predictions.csv", index=False)
    event_integration_summary.to_csv(args.out_dir / "event_integration_summary.csv", index=False)
    event_integration_coef.to_csv(args.out_dir / "event_integration_coefficients.csv", index=False)
    event_integration_predictions.to_csv(args.out_dir / "event_integration_predictions.csv", index=False)
    event_distance_window_comparison.to_csv(args.out_dir / "event_distance_window_comparison.csv", index=False)
    event_distance_window_coef.to_csv(args.out_dir / "event_distance_window_coefficients.csv", index=False)
    event_distance_window_example.to_csv(args.out_dir / "event_distance_window_example.csv", index=False)
    integration_metric_summary.to_csv(args.out_dir / "integration_metric_summary.csv", index=False)
    integration_metric_coef.to_csv(args.out_dir / "integration_metric_coefficients.csv", index=False)
    event_activity_summary.to_csv(args.out_dir / "event_activity_vedba_summary.csv", index=False)
    event_gps_movement_summary.to_csv(args.out_dir / "event_gps_movement_summary.csv", index=False)
    weekly_hierarchical_meeting_coef = meeting_probability_coef[
        meeting_probability_coef["model_type"].eq("hierarchical")
        & meeting_probability_coef["time_unit"].eq("weekly")
    ]
    weekly_separate_meeting_coef = meeting_probability_coef[
        meeting_probability_coef["model_type"].eq("separate")
        & meeting_probability_coef["time_unit"].eq("weekly")
    ]
    event_effect_plots = [
        plot_event_effect_sizes(
            weekly_separate_meeting_coef,
            args.out_dir,
            "meeting_probability_separate_effect_sizes.png",
            "Meeting-probability effect sizes: separate model",
            "Standardized coefficient on logit(meeting probability)",
            "#f28e2b",
            label_map=MEETING_EFFECT_LABELS,
        ),
        plot_event_effect_sizes(
            weekly_hierarchical_meeting_coef,
            args.out_dir,
            "meeting_probability_hierarchical_effect_sizes.png",
            "Meeting-probability effect sizes: hierarchical model",
            "Standardized coefficient on logit(meeting probability)",
            "#f28e2b",
            label_map=MEETING_EFFECT_LABELS,
        ),
        plot_event_effect_sizes(
            event_duration_coef,
            args.out_dir,
            "event_duration_effect_sizes.png",
            "Duration model effect sizes",
            "Standardized coefficient on log1p(duration hours)",
            "#4c78a8",
        ),
        plot_event_effect_sizes(
            event_integration_coef,
            args.out_dir,
            "event_integration_5m_effect_sizes.png",
            "5 m integration model effect sizes",
            "Standardized coefficient on logit(integration fraction)",
            "#59a14f",
        ),
        plot_event_effect_size_comparison(
            event_duration_coef,
            event_integration_coef,
            args.out_dir,
        ),
        plot_three_part_effect_size_comparison(
            weekly_hierarchical_meeting_coef,
            event_duration_coef,
            event_integration_coef,
            args.out_dir,
        ),
        *event_distance_window_plots,
    ]
    images = plot_predictions(rows, pred, args.out_dir)
    example_rows = build_example_rows(rows)
    example_rows.to_csv(args.out_dir / "daily_interaction_example_rows.csv", index=False)
    example_plot = plot_example_data(rows, args.out_dir)
    event_plot, event_detail = plot_specific_event(daily_rows, args.raw_interactions, args.proximity_status, args.out_dir)
    event_detail.to_csv(args.out_dir / "specific_interaction_event_detail.csv", index=False)
    long_event_plots, long_event_detail = plot_long_distance_events(
        daily_rows,
        args.raw_interactions,
        args.proximity_status,
        args.out_dir,
    )
    long_event_detail.to_csv(args.out_dir / "long_centroid_distance_event_detail.csv", index=False)
    dashboard = write_dashboard(
        args.out_dir,
        summary,
        meeting_probability_summary,
        event_duration_summary,
        event_integration_summary,
        images,
        meeting_probability_plots,
        event_duration_plots,
        event_integration_plots,
        event_effect_plots,
        example_rows,
        example_plot,
        event_plot,
        event_detail,
        long_event_plots,
        long_event_detail,
        event_distance_window_comparison,
        event_distance_window_example,
        rows,
        episodes,
        include_ndvi,
    )
    metadata = {
        "daily_dir": str(args.daily_dir),
        "raw_interactions": str(args.raw_interactions),
        "proximity_status": str(args.proximity_status),
        "vedba_dir": str(args.vedba_dir),
        "ndvi_used": include_ndvi,
        "event_duration_ndvi_used": event_duration_ndvi_used,
        "event_integration_ndvi_used": event_integration_ndvi_used,
        "meeting_probability_weekly_ndvi_used": include_ndvi,
        "meeting_probability_daily_ndvi_used": include_daily_ndvi,
        "modeling_note": (
            "Python first-pass hurdle model. Occurrence and positive duration are fit separately "
            "with dyad-clustered GEE. This is an approximation to a hierarchical hurdle model. "
            "The fitted row unit is a group-pair week aggregated from daytime dyad-days. "
            f"Centroids and interactions are restricted to {DAY_START_HOUR:02d}:00-{DAY_END_HOUR:02d}:00; "
            "the sparse overnight period is treated as held constant rather than used to update locations."
        ),
        "event_duration_modeling_note": (
            "Event-duration model uses only dyads that interacted. Positive bins separated by up to "
            f"{EVENT_MAX_GAP_HOURS:.1f} hours are merged into one meeting episode. The duration response "
            "is elapsed time from first positive bin to last positive bin plus one 2-minute bin width. "
            "The fitted duration model is a linear Bayesian Gaussian mixed model on log1p(duration_hours). "
            "All merge sizes are pooled in one model. "
            "Fixed predictors include log median centroid distance in the one-week-before through "
            "one-week-after event window, mean NDVI in that window, total group size, absolute group-size difference, "
            "and prior meeting count. "
            "Random intercepts are included for group_a, group_b, and pair_key. Weather and NDVI difference/change "
            "are not included in the duration model. Priors are weakly informed with expected negative distance "
            "effects and positive prior-meeting effects."
        ),
        "event_integration_modeling_note": (
            "Event-integration model uses the same pooled meeting episodes, restricted to events with candidate "
            "5 m links. The response is logit(total cross-group 5 m links / total candidate 5 m links). "
            "Fixed predictors and random intercepts match the event-duration model, with weakly informed priors."
        ),
        "three_part_bayesian_modeling_note": (
            "Meeting probability is fit on dyad-days and dyad-weeks with both a fixed-effect separate Bayesian logistic model "
            "and a hierarchical Bayesian logistic model. The hierarchical three-part view uses the hierarchical "
            "meeting-probability model, followed by event-duration and event-integration models conditional on "
            "a detected meeting event."
        ),
        "daytime_window": {
            "start_hour": DAY_START_HOUR,
            "end_hour": DAY_END_HOUR,
            "overnight_assumption": "location held constant from 16:00 to 03:00; no sparse night fixes used for centroid distance or interaction evidence",
        },
        "equivalent_hierarchical_formula": {
            "occurrence": "interval_any_interaction ~ s(log1p(range_mean_centroid_distance_m)) + prior positive intervals + s(ndvi_mean) + s(ndvi_abs_diff) + total group size + group size difference + (1|group_a) + (1|group_b) + (1|pair_key)",
            "duration": "log1p(weekly_duration_hours) ~ s(log1p(range_mean_centroid_distance_m)) + prior positive weeks + s(ndvi_mean) + s(ndvi_abs_diff) + total group size + group size difference + (1|group_a) + (1|group_b) + (1|pair_key), conditional on weekly_any_interaction == 1",
            "integration_5m": "weekly_cross_group_5m_edges / weekly_candidate_5m_edges ~ s(log1p(range_mean_centroid_distance_m)) + prior positive weeks + s(ndvi_mean) + s(ndvi_abs_diff) + total group size + group size difference + (1|group_a) + (1|group_b) + (1|pair_key)",
            "event_duration_bayesian": "log1p(duration_hours) ~ log1p(median_centroid_distance_1wk_before_after) + log1p(prior_meetings) + mean_ndvi_1wk_before_after + total group size + group size difference + (1|group_a) + (1|group_b) + (1|pair_key), all merge sizes pooled",
            "event_integration_bayesian": "logit(total_cross_group_5m_edges / total_candidate_5m_edges) ~ log1p(median_centroid_distance_1wk_before_after) + log1p(prior_meetings) + mean_ndvi_1wk_before_after + total group size + group size difference + (1|group_a) + (1|group_b) + (1|pair_key), all merge sizes pooled",
            "meeting_probability_bayesian_separate": "interval_any_interaction ~ log1p(range_mean_centroid_distance_m) + log1p(prior_meetings) + mean_ndvi + total group size + group size difference",
            "meeting_probability_bayesian_hierarchical": "interval_any_interaction ~ log1p(range_mean_centroid_distance_m) + log1p(prior_meetings) + mean_ndvi + total group size + group size difference + (1|group_a) + (1|group_b) + (1|pair_key)",
        },
        "outputs": {
            "dashboard": str(dashboard),
            "model_rows": str(args.out_dir / "daily_interaction_hurdle_model_rows.csv"),
            "daily_event_rows": str(args.out_dir / "daily_interaction_hurdle_daily_event_rows.csv"),
            "episodes": str(args.out_dir / "daily_interaction_positive_episodes.csv"),
            "example_rows": str(args.out_dir / "daily_interaction_example_rows.csv"),
            "specific_event_detail": str(args.out_dir / "specific_interaction_event_detail.csv"),
            "long_centroid_distance_event_detail": str(args.out_dir / "long_centroid_distance_event_detail.csv"),
            "summary": str(args.out_dir / "daily_interaction_hurdle_summary.csv"),
            "coefficients": str(args.out_dir / "daily_interaction_hurdle_coefficients.csv"),
            "predictions": str(args.out_dir / "daily_interaction_hurdle_predictions.csv"),
            "meeting_probability_summary": str(args.out_dir / "meeting_probability_summary.csv"),
            "meeting_probability_coefficients": str(args.out_dir / "meeting_probability_coefficients.csv"),
            "meeting_probability_predictions": str(args.out_dir / "meeting_probability_predictions.csv"),
            "event_duration_rows": str(args.out_dir / "event_duration_model_rows.csv"),
            "event_duration_summary": str(args.out_dir / "event_duration_summary.csv"),
            "event_duration_coefficients": str(args.out_dir / "event_duration_coefficients.csv"),
            "event_duration_predictions": str(args.out_dir / "event_duration_predictions.csv"),
            "event_integration_summary": str(args.out_dir / "event_integration_summary.csv"),
            "event_integration_coefficients": str(args.out_dir / "event_integration_coefficients.csv"),
            "event_integration_predictions": str(args.out_dir / "event_integration_predictions.csv"),
            "event_distance_window_comparison": str(args.out_dir / "event_distance_window_comparison.csv"),
            "event_distance_window_coefficients": str(args.out_dir / "event_distance_window_coefficients.csv"),
            "event_distance_window_example": str(args.out_dir / "event_distance_window_example.csv"),
            "integration_metrics_dashboard": str(integration_metrics_dashboard),
            "integration_metric_summary": str(args.out_dir / "integration_metric_summary.csv"),
            "integration_metric_coefficients": str(args.out_dir / "integration_metric_coefficients.csv"),
            "event_activity_vedba_summary": str(args.out_dir / "event_activity_vedba_summary.csv"),
            "event_gps_movement_summary": str(args.out_dir / "event_gps_movement_summary.csv"),
        },
    }
    (args.out_dir / "daily_interaction_hurdle_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print(f"Wrote dashboard: {dashboard.resolve()}")
    print(f"Wrote model rows: {(args.out_dir / 'daily_interaction_hurdle_model_rows.csv').resolve()}")


if __name__ == "__main__":
    main()
