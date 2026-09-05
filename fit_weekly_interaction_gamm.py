from __future__ import annotations

import argparse
import html
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.tools.sm_exceptions import DomainWarning


BASE = Path(r"C:\Users\rharel\Documents")
DEFAULT_CANONICAL = (
    BASE
    / "New project"
    / "outputs"
    / "canonical_robust_hourly_membership_local_2h_support"
    / "canonical_hourly_membership.parquet"
)
DEFAULT_INTERACTIONS = (
    BASE
    / "group_mebership"
    / "outputs"
    / "bigmerge_5m_long_events_2min_days"
    / "bigmerge_5m_hourly_long_event_assignments.csv"
)
DEFAULT_NDVI = BASE / "Github" / "detect_group_membership" / "indices_sentinel2.csv"
DEFAULT_OUT_DIR = BASE / "group_mebership" / "outputs" / "weekly_interaction_gamm"

EARTH_RADIUS_M = 6_371_000.0
SPLINE_DF_DISTANCE = 5
SPLINE_DF_NDVI = 4
MIN_MODEL_ROWS = 20
MIN_TOTAL_EDGES = 1
MIN_GROUP_CONTROL_ROWS = 5
N_DISTANCE_ANNOTATIONS = 3


def normalize_group_name(value: object) -> str:
    return "".join(str(value).split()).lower()


def period_start(series: pd.Series, period: str) -> pd.Series:
    values = pd.to_datetime(series)
    if period == "daily":
        return values.dt.floor("D")
    return values.dt.to_period("W-SUN").dt.start_time


def haversine_m(lat1: pd.Series, lon1: pd.Series, lat2: pd.Series, lon2: pd.Series) -> pd.Series:
    lat1_rad = np.radians(lat1.astype(float))
    lon1_rad = np.radians(lon1.astype(float))
    lat2_rad = np.radians(lat2.astype(float))
    lon2_rad = np.radians(lon2.astype(float))
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(a))


def load_weekly_group_centroids(canonical_path: Path, period: str = "weekly") -> pd.DataFrame:
    columns = [
        "animal_id",
        "window_start",
        "mean_latitude",
        "mean_longitude",
        "dynamic_social_unit",
    ]
    data = pd.read_parquet(canonical_path, columns=columns)
    data["window_start"] = pd.to_datetime(data["window_start"])
    data = data.dropna(
        subset=["window_start", "mean_latitude", "mean_longitude", "dynamic_social_unit"]
    ).copy()
    data["period_start"] = period_start(data["window_start"], period)
    data["dynamic_social_unit"] = data["dynamic_social_unit"].astype(str)
    centroids = (
        data.groupby(["period_start", "dynamic_social_unit"], observed=True)
        .agg(
            centroid_latitude=("mean_latitude", "mean"),
            centroid_longitude=("mean_longitude", "mean"),
            n_position_rows=("animal_id", "size"),
            n_animals=("animal_id", "nunique"),
            first_timestamp=("window_start", "min"),
            last_timestamp=("window_start", "max"),
        )
        .reset_index()
    )
    centroids["group_key"] = centroids["dynamic_social_unit"].map(normalize_group_name)
    return centroids


def load_weekly_interactions(interactions_path: Path, period: str = "weekly") -> pd.DataFrame:
    rows = pd.read_csv(interactions_path, parse_dates=["timestamp"])
    required = ["timestamp", "group_a", "group_b", "pair_key", "cross_edges", "total_edges"]
    missing = [col for col in required if col not in rows.columns]
    if missing:
        raise ValueError(f"Interaction file is missing required columns: {missing}")

    rows["period_start"] = period_start(rows["timestamp"], period)
    rows["group_a"] = rows["group_a"].astype(str)
    rows["group_b"] = rows["group_b"].astype(str)
    rows["pair_key"] = rows["pair_key"].astype(str)
    rows["cross_edges"] = pd.to_numeric(rows["cross_edges"], errors="coerce").fillna(0)
    rows["total_edges"] = pd.to_numeric(rows["total_edges"], errors="coerce").fillna(0)
    rows["non_cross_edges"] = (rows["total_edges"] - rows["cross_edges"]).clip(lower=0)

    agg_spec: dict[str, tuple[str, str]] = {
        "interaction_edges": ("cross_edges", "sum"),
        "non_interaction_edges": ("non_cross_edges", "sum"),
        "total_edges": ("total_edges", "sum"),
        "n_interaction_rows": ("pair_key", "size"),
        "first_timestamp": ("timestamp", "min"),
        "last_timestamp": ("timestamp", "max"),
    }
    optional_means = [
        "pair_centroid_distance_m",
        "pair_n",
        "cluster_size_total",
        "group_a_fraction_present",
        "group_b_fraction_present",
        "min_group_fraction_present",
        "composition_entropy_norm",
        "pair_balance",
    ]
    for col in optional_means:
        if col in rows.columns:
            agg_spec[f"mean_{col}"] = (col, "mean")

    weekly = (
        rows.groupby(["period_start", "group_a", "group_b", "pair_key"], observed=True)
        .agg(**agg_spec)
        .reset_index()
    )
    weekly = weekly[weekly["total_edges"].ge(MIN_TOTAL_EDGES)].copy()
    weekly["interaction_probability"] = weekly["interaction_edges"] / weekly["total_edges"]
    weekly["group_a_key"] = weekly["group_a"].map(normalize_group_name)
    weekly["group_b_key"] = weekly["group_b"].map(normalize_group_name)
    return weekly


def add_weekly_centroid_distance(weekly: pd.DataFrame, centroids: pd.DataFrame) -> pd.DataFrame:
    centroid_cols = [
        "period_start",
        "group_key",
        "dynamic_social_unit",
        "centroid_latitude",
        "centroid_longitude",
        "n_position_rows",
        "n_animals",
    ]
    left = centroids[centroid_cols].rename(
        columns={
            "group_key": "group_a_key",
            "dynamic_social_unit": "group_a_centroid_label",
            "centroid_latitude": "group_a_centroid_latitude",
            "centroid_longitude": "group_a_centroid_longitude",
            "n_position_rows": "group_a_centroid_position_rows",
            "n_animals": "group_a_centroid_animals",
        }
    )
    right = centroids[centroid_cols].rename(
        columns={
            "group_key": "group_b_key",
            "dynamic_social_unit": "group_b_centroid_label",
            "centroid_latitude": "group_b_centroid_latitude",
            "centroid_longitude": "group_b_centroid_longitude",
            "n_position_rows": "group_b_centroid_position_rows",
            "n_animals": "group_b_centroid_animals",
        }
    )
    out = weekly.merge(left, on=["period_start", "group_a_key"], how="left")
    out = out.merge(right, on=["period_start", "group_b_key"], how="left")
    have_centroids = out[
        [
            "group_a_centroid_latitude",
            "group_a_centroid_longitude",
            "group_b_centroid_latitude",
            "group_b_centroid_longitude",
        ]
    ].notna().all(axis=1)
    out["weekly_centroid_distance_m"] = np.nan
    out.loc[have_centroids, "weekly_centroid_distance_m"] = haversine_m(
        out.loc[have_centroids, "group_a_centroid_latitude"],
        out.loc[have_centroids, "group_a_centroid_longitude"],
        out.loc[have_centroids, "group_b_centroid_latitude"],
        out.loc[have_centroids, "group_b_centroid_longitude"],
    )
    return out


def infer_column(columns: list[str], candidates: list[str]) -> str | None:
    lower = {col.lower(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def load_weekly_ndvi(ndvi_path: Path, period: str = "weekly") -> pd.DataFrame:
    if not ndvi_path.exists():
        return pd.DataFrame()
    ndvi = pd.read_csv(ndvi_path)
    group_col = infer_column(ndvi.columns.tolist(), ["group_id", "dynamic_social_unit", "group", "group_name"])
    date_col = infer_column(ndvi.columns.tolist(), ["date", "image_date", "period_start", "week_start"])
    value_col = infer_column(ndvi.columns.tolist(), ["NDVI", "ndvi", "mean_ndvi", "ndvi_mean"])
    if not group_col or not date_col or not value_col:
        raise ValueError(
            "NDVI table must include group, date/week, and NDVI columns. "
            f"Found columns: {list(ndvi.columns)}"
        )

    ndvi = ndvi.rename(columns={group_col: "group_name", date_col: "ndvi_date", value_col: "ndvi"})
    ndvi["ndvi_date"] = pd.to_datetime(ndvi["ndvi_date"], errors="coerce")
    ndvi["ndvi"] = pd.to_numeric(ndvi["ndvi"], errors="coerce")
    ndvi = ndvi.dropna(subset=["group_name", "ndvi_date", "ndvi"]).copy()
    ndvi["period_start"] = period_start(ndvi["ndvi_date"], period)
    ndvi["group_key"] = ndvi["group_name"].map(normalize_group_name)
    weekly = (
        ndvi.groupby(["period_start", "group_key"], observed=True)
        .agg(
            weekly_ndvi_mean=("ndvi", "mean"),
            weekly_ndvi_median=("ndvi", "median"),
            weekly_ndvi_n=("ndvi", "size"),
            weekly_ndvi_first_date=("ndvi_date", "min"),
            weekly_ndvi_last_date=("ndvi_date", "max"),
        )
        .reset_index()
    )
    return weekly


def add_weekly_ndvi(weekly: pd.DataFrame, ndvi: pd.DataFrame) -> pd.DataFrame:
    if ndvi.empty:
        return weekly.copy()
    left = ndvi.rename(
        columns={
            "group_key": "group_a_key",
            "weekly_ndvi_mean": "group_a_weekly_ndvi_mean",
            "weekly_ndvi_median": "group_a_weekly_ndvi_median",
            "weekly_ndvi_n": "group_a_weekly_ndvi_n",
            "weekly_ndvi_first_date": "group_a_weekly_ndvi_first_date",
            "weekly_ndvi_last_date": "group_a_weekly_ndvi_last_date",
        }
    )
    right = ndvi.rename(
        columns={
            "group_key": "group_b_key",
            "weekly_ndvi_mean": "group_b_weekly_ndvi_mean",
            "weekly_ndvi_median": "group_b_weekly_ndvi_median",
            "weekly_ndvi_n": "group_b_weekly_ndvi_n",
            "weekly_ndvi_first_date": "group_b_weekly_ndvi_first_date",
            "weekly_ndvi_last_date": "group_b_weekly_ndvi_last_date",
        }
    )
    out = weekly.merge(left, on=["period_start", "group_a_key"], how="left")
    out = out.merge(right, on=["period_start", "group_b_key"], how="left")
    out["dyad_weekly_ndvi_mean"] = out[
        ["group_a_weekly_ndvi_mean", "group_b_weekly_ndvi_mean"]
    ].mean(axis=1)
    out["dyad_weekly_ndvi_abs_diff"] = (
        out["group_a_weekly_ndvi_mean"] - out["group_b_weekly_ndvi_mean"]
    ).abs()
    return out


def zscore(series: pd.Series) -> pd.Series:
    values = series.astype(float)
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index)
    return (values - values.mean()) / sd


def prepare_model_df(rows: pd.DataFrame, include_ndvi: bool) -> pd.DataFrame:
    cols = [
        "period_start",
        "pair_key",
        "group_a",
        "group_b",
        "interaction_edges",
        "non_interaction_edges",
        "total_edges",
        "interaction_probability",
        "weekly_centroid_distance_m",
    ]
    if include_ndvi:
        cols += ["dyad_weekly_ndvi_mean", "dyad_weekly_ndvi_abs_diff"]
    model_df = rows[cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    model_df = model_df[model_df["total_edges"].gt(0)]
    model_df["log1p_weekly_centroid_distance_m"] = np.log1p(
        model_df["weekly_centroid_distance_m"].clip(lower=0)
    )
    model_df["log1p_weekly_centroid_distance_m_z"] = zscore(
        model_df["log1p_weekly_centroid_distance_m"]
    )
    for col in ["group_a", "group_b"]:
        counts = model_df[col].value_counts()
        keep = set(counts[counts.ge(MIN_GROUP_CONTROL_ROWS)].index)
        model_df[f"{col}_control"] = np.where(
            model_df[col].isin(keep),
            model_df[col],
            f"other_{col}",
        )
    if include_ndvi:
        model_df["dyad_weekly_ndvi_mean_z"] = zscore(model_df["dyad_weekly_ndvi_mean"])
        model_df["dyad_weekly_ndvi_abs_diff_z"] = zscore(model_df["dyad_weekly_ndvi_abs_diff"])
    return model_df.dropna().copy()


def build_formula(include_ndvi: bool) -> str:
    terms = [
        f"bs(log1p_weekly_centroid_distance_m_z, df={SPLINE_DF_DISTANCE}, degree=3, include_intercept=False)",
    ]
    if include_ndvi:
        terms.extend(
            [
                f"bs(dyad_weekly_ndvi_mean_z, df={SPLINE_DF_NDVI}, degree=3, include_intercept=False)",
                f"bs(dyad_weekly_ndvi_abs_diff_z, df={SPLINE_DF_NDVI}, degree=3, include_intercept=False)",
            ]
        )
    terms.extend(["C(group_a_control)", "C(group_b_control)"])
    return "interaction_probability ~ " + " + ".join(terms)


def fit_binomial_spline(model_df: pd.DataFrame, include_ndvi: bool) -> tuple[object, str, str]:
    formula = build_formula(include_ndvi)
    if len(model_df) < MIN_MODEL_ROWS or model_df["pair_key"].nunique() < 2:
        raise ValueError("Not enough weekly dyad rows to fit model.")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DomainWarning)
        try:
            model = GEE.from_formula(
                formula,
                groups="pair_key",
                data=model_df,
                family=Binomial(),
                weights=model_df["total_edges"],
            )
            result = model.fit(maxiter=200)
            return result, formula, "GEE binomial spline; dyad-clustered working correlation"
        except Exception as exc:
            glm = sm.GLM.from_formula(
                formula,
                data=model_df,
                family=Binomial(),
                freq_weights=model_df["total_edges"],
            )
            result = glm.fit(maxiter=200)
            return result, formula, f"GLM binomial spline fallback after GEE error: {exc}"


def coefficient_table(result: object, model_name: str) -> pd.DataFrame:
    params = result.params
    bse = result.bse
    out = pd.DataFrame(
        {
            "model": model_name,
            "term": params.index,
            "estimate": params.to_numpy(),
            "se": bse.reindex(params.index).to_numpy(),
        }
    )
    out["z"] = out["estimate"] / out["se"]
    out["p_value"] = 2 * (1 - norm.cdf(out["z"].abs()))
    return out


def model_summary(result: object, model_df: pd.DataFrame, model_name: str, formula: str, fit_note: str) -> dict[str, object]:
    return {
        "model": model_name,
        "fit_note": fit_note,
        "formula": formula,
        "n_period_dyad_rows": int(len(model_df)),
        "n_dyads": int(model_df["pair_key"].nunique()),
        "n_group_a": int(model_df["group_a"].nunique()),
        "n_group_b": int(model_df["group_b"].nunique()),
        "total_edges": float(model_df["total_edges"].sum()),
        "interaction_edges": float(model_df["interaction_edges"].sum()),
        "overall_interaction_probability": float(
            model_df["interaction_edges"].sum() / model_df["total_edges"].sum()
        ),
        "distance_min_m": float(model_df["weekly_centroid_distance_m"].min()),
        "distance_median_m": float(model_df["weekly_centroid_distance_m"].median()),
        "distance_max_m": float(model_df["weekly_centroid_distance_m"].max()),
        "aic": float(getattr(result, "aic", np.nan)) if np.isfinite(getattr(result, "aic", np.nan)) else np.nan,
    }


def prediction_grid(model_df: pd.DataFrame, result: object, include_ndvi: bool, model_name: str) -> pd.DataFrame:
    distance_grid = np.linspace(
        model_df["weekly_centroid_distance_m"].quantile(0.02),
        model_df["weekly_centroid_distance_m"].quantile(0.98),
        120,
    )
    group_a = model_df["group_a"].mode().iloc[0]
    group_b = model_df["group_b"].mode().iloc[0]
    group_a_control = model_df.loc[model_df["group_a"].eq(group_a), "group_a_control"].mode().iloc[0]
    group_b_control = model_df.loc[model_df["group_b"].eq(group_b), "group_b_control"].mode().iloc[0]
    grid = pd.DataFrame(
        {
            "weekly_centroid_distance_m": distance_grid,
            "log1p_weekly_centroid_distance_m": np.log1p(distance_grid),
            "group_a": group_a,
            "group_b": group_b,
            "group_a_control": group_a_control,
            "group_b_control": group_b_control,
            "pair_key": f"{group_a} - {group_b}",
        }
    )
    dist_mean = model_df["log1p_weekly_centroid_distance_m"].mean()
    dist_sd = model_df["log1p_weekly_centroid_distance_m"].std(ddof=0)
    grid["log1p_weekly_centroid_distance_m_z"] = (
        grid["log1p_weekly_centroid_distance_m"] - dist_mean
    ) / dist_sd

    if include_ndvi:
        for col in ["dyad_weekly_ndvi_mean", "dyad_weekly_ndvi_abs_diff"]:
            zcol = f"{col}_z"
            grid[col] = model_df[col].median()
            grid[zcol] = (grid[col] - model_df[col].mean()) / model_df[col].std(ddof=0)

    design = patsy.build_design_matrices([result.model.data.design_info], grid)[0]
    exog = np.asarray(design)
    beta = result.params.to_numpy()
    cov = result.cov_params()
    if isinstance(cov, pd.DataFrame):
        cov = cov.loc[result.params.index, result.params.index].to_numpy()
    linear = exog @ beta
    linear_se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", exog, cov, exog), 0))
    grid["predicted_probability"] = 1 / (1 + np.exp(-linear))
    grid["ci_low"] = 1 / (1 + np.exp(-(linear - 1.96 * linear_se)))
    grid["ci_high"] = 1 / (1 + np.exp(-(linear + 1.96 * linear_se)))
    grid["model"] = model_name
    return grid


def plot_predictions(model_df: pd.DataFrame, pred: pd.DataFrame, out_dir: Path, model_name: str) -> Path:
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    plot_df = model_df.copy()
    if len(plot_df) > 5000:
        plot_df = plot_df.sample(5000, random_state=42)
    sizes = np.clip(np.sqrt(plot_df["total_edges"]) * 7, 12, 95)
    period_label = model_df.attrs.get("period_label", "weekly")
    period_title = str(period_label).title()
    ax.scatter(
        plot_df["weekly_centroid_distance_m"],
        plot_df["interaction_probability"],
        s=sizes,
        color="#4C78A8",
        alpha=0.2,
        linewidth=0,
        label=f"{period_label} dyad rows",
    )
    ax.plot(
        pred["weekly_centroid_distance_m"],
        pred["predicted_probability"],
        color="#D95F02",
        linewidth=2.6,
        label="model prediction",
    )
    ax.fill_between(
        pred["weekly_centroid_distance_m"],
        pred["ci_low"],
        pred["ci_high"],
        color="#D95F02",
        alpha=0.18,
        linewidth=0,
    )
    annotate_far_right_points(ax, model_df)
    ax.set_xlabel(f"{period_title} distance between group centroids (m)")
    ax.set_ylabel("Probability of interaction: summed cross edges / summed total edges")
    label = model_df.attrs.get("analysis_label", "")
    title = f"{period_title} intergroup interaction probability"
    if label:
        title += f"\n{label}"
    ax.set_title(
        title,
        fontsize=14,
        fontweight="bold",
    )
    ax.text(
        0.01,
        0.98,
        f"Binomial spline with group controls; point size scales with {period_label} total edges",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color="#555",
        fontsize=9.5,
    )
    ax.grid(True, color="#e6e6e6", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="best", frameon=True)
    fig.tight_layout()
    out = out_dir / f"{model_name}_prediction_distance.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def far_right_dyad_weeks(rows: pd.DataFrame, n: int = N_DISTANCE_ANNOTATIONS) -> pd.DataFrame:
    cols = [
        "period_start",
        "pair_key",
        "interaction_edges",
        "total_edges",
        "interaction_probability",
        "weekly_centroid_distance_m",
    ]
    optional = ["dyad_weekly_ndvi_mean", "dyad_weekly_ndvi_abs_diff"]
    cols.extend([col for col in optional if col in rows.columns])
    return (
        rows.dropna(subset=["weekly_centroid_distance_m"])
        .sort_values("weekly_centroid_distance_m", ascending=False)
        .loc[:, cols]
        .head(n)
        .copy()
    )


def annotate_far_right_points(ax: plt.Axes, model_df: pd.DataFrame) -> None:
    annotations = far_right_dyad_weeks(model_df)
    if annotations.empty:
        return

    offsets = [(-128, 38), (-118, -34), (-132, 18), (-122, -48), (-116, 54)]
    for i, (_, row) in enumerate(annotations.iterrows()):
        period = pd.to_datetime(row["period_start"]).date().isoformat()
        label = (
            f"{row['pair_key']}\n"
            f"{period}; {int(row['interaction_edges'])}/{int(row['total_edges'])}"
        )
        ax.annotate(
            label,
            xy=(row["weekly_centroid_distance_m"], row["interaction_probability"]),
            xytext=offsets[i % len(offsets)],
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize=8.2,
            color="#333333",
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#B8B8B8",
                "alpha": 0.92,
            },
            arrowprops={
                "arrowstyle": "->",
                "color": "#777777",
                "linewidth": 0.8,
                "shrinkA": 2,
                "shrinkB": 3,
            },
        )


def write_html(
    out_dir: Path,
    summary: pd.DataFrame,
    image_paths: list[Path],
    ndvi_used: bool,
    analysis_label: str,
    period: str,
) -> Path:
    out = out_dir / "weekly_interaction_gamm_dashboard.html"
    table_rows = "\n".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "model",
                "n_period_dyad_rows",
                "n_dyads",
                "total_edges",
                "overall_interaction_probability",
                "distance_median_m",
                "fit_note",
            ]
        )
        + "</tr>"
        for _, row in summary.round(
            {
                "total_edges": 0,
                "overall_interaction_probability": 4,
                "distance_median_m": 1,
            }
        ).iterrows()
    )
    images = "\n".join(
        f"<img src='{html.escape(path.name)}' alt='{html.escape(path.stem)}'>" for path in image_paths
    )
    period_label = "daily" if period == "daily" else "weekly"
    period_title = period_label.title()
    ndvi_note = (
        f"NDVI terms included where both groups had {period_label} NDVI."
        if ndvi_used
        else f"NDVI was not included because no usable {period_label} dyad NDVI was available."
    )
    label_note = f"<p class=\"note\"><b>Analysis label:</b> {html.escape(analysis_label)}</p>" if analysis_label else ""
    out.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Weekly interaction GAMM</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;background:#fff;color:#222}}
.wrap{{padding:24px 30px 44px;max-width:1320px}}
.note{{line-height:1.45;color:#555;max-width:980px}}
img{{max-width:100%;display:block;margin:18px 0 34px;border:1px solid #ddd;border-radius:6px}}
table{{border-collapse:collapse;margin:16px 0 32px;font-size:13px}}
th,td{{border:1px solid #ddd;padding:6px 8px;text-align:left;vertical-align:top}}
th{{background:#f2f2f2}}
</style></head><body><div class="wrap">
<h1>{period_title} intergroup interaction model</h1>
{label_note}
<p class="note">Response is a binomial probability built from summed {period_label} interactions:
cross-group edges divided by total edges for each group pair and {period_label.rstrip("ly")}. Predictors are {period_label}
group-centroid distance, group controls, and NDVI when available. {html.escape(ndvi_note)}</p>
<table>
<tr><th>model</th><th>rows</th><th>dyads</th><th>total edges</th><th>overall probability</th><th>median distance m</th><th>fit note</th></tr>
{table_rows}
</table>
{images}
</div></body></html>""",
        encoding="utf-8",
    )
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit weekly intergroup interaction GAMM-like binomial spline models."
    )
    parser.add_argument("--canonical", type=Path, default=DEFAULT_CANONICAL)
    parser.add_argument("--interactions", type=Path, default=DEFAULT_INTERACTIONS)
    parser.add_argument("--ndvi", type=Path, default=DEFAULT_NDVI)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--period", choices=["weekly", "daily"], default="weekly")
    parser.add_argument(
        "--no-ndvi",
        action="store_true",
        help="Skip NDVI join and fit only the centroid-distance model.",
    )
    parser.add_argument(
        "--analysis-label",
        default="",
        help="Short label to show in output dashboards and plot titles.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    centroids = load_weekly_group_centroids(args.canonical, period=args.period)
    interactions = load_weekly_interactions(args.interactions, period=args.period)
    rows = add_weekly_centroid_distance(interactions, centroids)

    ndvi = pd.DataFrame()
    if not args.no_ndvi:
        ndvi = load_weekly_ndvi(args.ndvi, period=args.period)
        rows = add_weekly_ndvi(rows, ndvi)

    centroids.to_csv(args.out_dir / "weekly_group_centroids.csv", index=False)
    if not ndvi.empty:
        ndvi.to_csv(args.out_dir / "weekly_group_ndvi.csv", index=False)
    rows.to_csv(args.out_dir / "weekly_interaction_model_input_rows.csv", index=False)
    far_right_dyad_weeks(rows, n=12).to_csv(
        args.out_dir / "far_right_dyad_weeks_by_centroid_distance.csv",
        index=False,
    )

    model_outputs: list[pd.DataFrame] = []
    coefficient_outputs: list[pd.DataFrame] = []
    prediction_outputs: list[pd.DataFrame] = []
    image_paths: list[Path] = []

    model_specs = [("distance_group_controls", False)]
    ndvi_candidate = prepare_model_df(rows, include_ndvi=True) if not ndvi.empty else pd.DataFrame()
    ndvi_used = len(ndvi_candidate) >= MIN_MODEL_ROWS and ndvi_candidate["pair_key"].nunique() >= 2
    if ndvi_used:
        model_specs.append(("distance_ndvi_group_controls", True))

    for model_name, include_ndvi in model_specs:
        model_df = prepare_model_df(rows, include_ndvi=include_ndvi)
        model_df.attrs["analysis_label"] = args.analysis_label
        model_df.attrs["period_label"] = args.period
        result, formula, fit_note = fit_binomial_spline(model_df, include_ndvi=include_ndvi)
        summary_row = model_summary(result, model_df, model_name, formula, fit_note)
        model_outputs.append(pd.DataFrame([summary_row]))
        coefficient_outputs.append(coefficient_table(result, model_name))
        pred = prediction_grid(model_df, result, include_ndvi, model_name)
        prediction_outputs.append(pred)
        image_paths.append(plot_predictions(model_df, pred, args.out_dir, model_name))
        model_df.to_csv(args.out_dir / f"{model_name}_model_rows.csv", index=False)

    summary = pd.concat(model_outputs, ignore_index=True)
    coefficients = pd.concat(coefficient_outputs, ignore_index=True)
    predictions = pd.concat(prediction_outputs, ignore_index=True)

    summary.to_csv(args.out_dir / "weekly_interaction_gamm_model_summary.csv", index=False)
    coefficients.to_csv(args.out_dir / "weekly_interaction_gamm_coefficients.csv", index=False)
    predictions.to_csv(args.out_dir / "weekly_interaction_gamm_predictions.csv", index=False)
    dashboard = write_html(args.out_dir, summary, image_paths, ndvi_used, args.analysis_label, args.period)

    metadata = {
        "canonical": str(args.canonical),
        "interactions": str(args.interactions),
        "ndvi": str(args.ndvi) if not args.no_ndvi else None,
        "ndvi_used": ndvi_used,
        "analysis_label": args.analysis_label,
        "period": args.period,
        "modeling_note": (
            "Python fit uses binomial spline terms with group fixed effects and dyad-clustered GEE. "
            "Rscript/mgcv is not available in this environment."
        ),
        "equivalent_mgcv_formula": (
            "mgcv::bam(cbind(interaction_edges, non_interaction_edges) ~ "
            "s(log1p(weekly_centroid_distance_m), k=5) + "
            "s(dyad_weekly_ndvi_mean, k=4) + "
            "s(dyad_weekly_ndvi_abs_diff, k=4) + "
            "s(group_a, bs='re') + s(group_b, bs='re') + s(pair_key, bs='re'), "
            "family=binomial(), method='fREML', data=weekly_rows)"
        ),
        "python_group_control_note": (
            f"Python model collapses group_a/group_b levels with fewer than {MIN_GROUP_CONTROL_ROWS} "
            "weekly rows into side-specific other_group controls to reduce sparse fixed-effect instability."
        ),
        "outputs": {
            "dashboard": str(dashboard),
            "summary": str(args.out_dir / "weekly_interaction_gamm_model_summary.csv"),
            "coefficients": str(args.out_dir / "weekly_interaction_gamm_coefficients.csv"),
            "predictions": str(args.out_dir / "weekly_interaction_gamm_predictions.csv"),
            "model_input": str(args.out_dir / "weekly_interaction_model_input_rows.csv"),
            "far_right_dyad_weeks": str(args.out_dir / "far_right_dyad_weeks_by_centroid_distance.csv"),
        },
    }
    (args.out_dir / "weekly_interaction_gamm_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print(f"Wrote dashboard: {dashboard.resolve()}")
    print(f"Wrote summary: {(args.out_dir / 'weekly_interaction_gamm_model_summary.csv').resolve()}")
    for path in image_paths:
        print(path.resolve())


if __name__ == "__main__":
    main()
