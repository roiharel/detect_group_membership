from __future__ import annotations

import json
import warnings
from itertools import count
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from patsy import build_design_matrices, dmatrix
from statsmodels.regression.mixed_linear_model import MixedLM


BASE = Path(r"C:\Users\rharel\Documents")
CANONICAL = (
    BASE
    / "New project"
    / "outputs"
    / "canonical_robust_hourly_membership_local_2h_support_conservative_dynamic"
    / "canonical_hourly_membership.parquet"
)
OUT_DIR = BASE / "group_mebership" / "outputs" / "canonical_integration_gamm_conservative_14d"

EARTH_RADIUS_M = 6_371_000.0
SUSTAINED = "sustained_non_origin_association"
ISOLATED = "ISOLATED"

MAX_GAP_HOURS = 36
MIN_SPLIT_TRACE_HOURS = 6
MIN_SPLIT_TRACE_TIMESTAMPS = 6
MIN_CONTEXT_ANIMALS = 2
SPLINE_DF = 5
MAX_MODEL_DAYS = 60
RAW_PLOT_SAMPLE = 25_000

TRACE_COLORS = {
    "dispersal_join": "#8d5a9e",
    "within_group_split_off": "#009e73",
}
RESPONSE_SPECS = {
    "centroid_centrality_percentile": {
        "label": "Centroid centrality percentile",
        "interpretation": "higher = more integrated / central",
        "ylim": (-3, 103),
    },
    "normalized_centroid_distance": {
        "label": "Normalized centroid distance",
        "interpretation": "lower = more integrated / central",
        "ylim": None,
    },
    "log1p_normalized_centroid_distance": {
        "label": "log1p normalized centroid distance",
        "interpretation": "lower = more integrated / central; log scale dampens extreme split distances",
        "ylim": None,
    },
}


def haversine_m(
    lat1: float | np.ndarray,
    lon1: float | np.ndarray,
    lat2: np.ndarray,
    lon2: np.ndarray,
) -> np.ndarray:
    lat1_rad = np.radians(lat1)
    lon1_rad = np.radians(lon1)
    lat2_rad = np.radians(lat2.astype(float))
    lon2_rad = np.radians(lon2.astype(float))
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(a))


def load_canonical() -> pd.DataFrame:
    columns = [
        "animal_id",
        "window_start",
        "sex",
        "age",
        "mean_latitude",
        "mean_longitude",
        "temp_group_id",
        "temp_group_type",
        "temp_group_size",
        "assigned_social_unit",
        "dynamic_social_unit",
        "origin_group",
        "dynamic_assignment",
        "dynamic_target_group",
        "dynamic_segment_id",
        "dynamic_contact_context",
        "is_carried_night",
        "is_local_2h_supported",
        "position_support_confidence",
        "membership_confidence",
    ]
    data = pd.read_parquet(CANONICAL, columns=columns)
    data["window_start"] = pd.to_datetime(data["window_start"])
    data = data.dropna(
        subset=[
            "animal_id",
            "window_start",
            "mean_latitude",
            "mean_longitude",
            "temp_group_id",
            "dynamic_social_unit",
        ]
    ).copy()
    data["animal_id"] = data["animal_id"].astype(str)
    data["dynamic_social_unit"] = data["dynamic_social_unit"].astype(str)
    data["temp_group_id"] = data["temp_group_id"].astype(str)
    return data


def add_split_component_flags(data: pd.DataFrame) -> pd.DataFrame:
    components = (
        data.groupby(["window_start", "dynamic_social_unit", "temp_group_id"], observed=True)
        .agg(component_size=("animal_id", "nunique"))
        .reset_index()
        .sort_values(
            ["window_start", "dynamic_social_unit", "component_size", "temp_group_id"],
            ascending=[True, True, False, True],
        )
    )
    components["component_rank"] = (
        components.groupby(["window_start", "dynamic_social_unit"], observed=True).cumcount() + 1
    )
    n_components = (
        components.groupby(["window_start", "dynamic_social_unit"], observed=True)["temp_group_id"]
        .nunique()
        .rename("n_dynamic_unit_components")
        .reset_index()
    )
    components = components.merge(n_components, on=["window_start", "dynamic_social_unit"], how="left")
    components["is_main_dynamic_unit_component"] = components["component_rank"].eq(1)
    out = data.merge(
        components,
        on=["window_start", "dynamic_social_unit", "temp_group_id"],
        how="left",
    )
    out["is_split_off_component"] = (
        out["n_dynamic_unit_components"].fillna(1).gt(1)
        & ~out["is_main_dynamic_unit_component"].fillna(True)
    )
    return out


def assign_contiguous_trace_ids(
    rows: pd.DataFrame,
    grouping_cols: list[str],
    prefix: str,
    min_hours: int,
    min_timestamps: int,
) -> pd.DataFrame:
    if rows.empty:
        return rows.copy()
    out = rows.sort_values(grouping_cols + ["window_start"]).copy()
    trace_counter = count(1)
    chunks = []
    for _, group in out.groupby(grouping_cols, sort=True, observed=True):
        group = group.sort_values("window_start").copy()
        gap_hours = group["window_start"].diff().dt.total_seconds().div(3600)
        group["_trace_run"] = gap_hours.gt(MAX_GAP_HOURS).fillna(False).cumsum()
        for _, run in group.groupby("_trace_run", sort=True):
            duration_hours = (
                (run["window_start"].max() - run["window_start"].min()).total_seconds() / 3600
            ) + 1
            n_timestamps = run["window_start"].nunique()
            if duration_hours < min_hours or n_timestamps < min_timestamps:
                continue
            run = run.copy()
            run["trace_id"] = f"{prefix}_{next(trace_counter):05d}"
            run["trace_start_time"] = run["window_start"].min()
            run["trace_end_time"] = run["window_start"].max()
            run["trace_duration_hours"] = duration_hours
            run["trace_n_timestamps"] = n_timestamps
            chunks.append(run.drop(columns="_trace_run"))
    if not chunks:
        return pd.DataFrame(columns=list(rows.columns) + ["trace_id"])
    return pd.concat(chunks, ignore_index=True)


def build_dispersal_traces(data: pd.DataFrame) -> pd.DataFrame:
    rows = data[data["dynamic_assignment"].eq(SUSTAINED)].copy()
    if rows.empty:
        return rows
    rows["trace_type"] = "dispersal_join"
    rows["context_label"] = rows["dynamic_social_unit"]
    rows["trace_grouping_key"] = rows["dynamic_segment_id"].fillna(
        rows["animal_id"] + "__" + rows["dynamic_social_unit"]
    )
    return assign_contiguous_trace_ids(
        rows,
        grouping_cols=["animal_id", "trace_grouping_key", "context_label"],
        prefix="DISP",
        min_hours=24,
        min_timestamps=6,
    )


def build_split_traces(data: pd.DataFrame) -> pd.DataFrame:
    rows = data[data["is_split_off_component"]].copy()
    if rows.empty:
        return rows
    rows["trace_type"] = "within_group_split_off"
    rows["context_label"] = rows["dynamic_social_unit"]
    return assign_contiguous_trace_ids(
        rows,
        grouping_cols=["animal_id", "dynamic_social_unit"],
        prefix="SPLIT",
        min_hours=MIN_SPLIT_TRACE_HOURS,
        min_timestamps=MIN_SPLIT_TRACE_TIMESTAMPS,
    )


def context_for_row(row: pd.Series, same_time: pd.DataFrame) -> pd.DataFrame:
    if row["trace_type"] == "dispersal_join":
        context = same_time[
            same_time["temp_group_id"].eq(row["temp_group_id"])
            & same_time["dynamic_social_unit"].eq(row["context_label"])
            & same_time["animal_id"].ne(row["animal_id"])
        ].copy()
    else:
        context = same_time[
            same_time["dynamic_social_unit"].eq(row["dynamic_social_unit"])
            & same_time["animal_id"].ne(row["animal_id"])
        ].copy()
    return context.dropna(subset=["mean_latitude", "mean_longitude"])


def centrality_metrics_for_row(row: pd.Series, context: pd.DataFrame) -> dict[str, object] | None:
    if len(context) < MIN_CONTEXT_ANIMALS:
        return None
    lat = context["mean_latitude"].astype(float).to_numpy()
    lon = context["mean_longitude"].astype(float).to_numpy()
    centroid_lat = float(np.mean(lat))
    centroid_lon = float(np.mean(lon))
    member_to_centroid = haversine_m(centroid_lat, centroid_lon, lat, lon)
    focal_to_centroid = float(
        haversine_m(
            float(row["mean_latitude"]),
            float(row["mean_longitude"]),
            np.array([centroid_lat]),
            np.array([centroid_lon]),
        )[0]
    )
    focal_to_members = haversine_m(
        float(row["mean_latitude"]),
        float(row["mean_longitude"]),
        lat,
        lon,
    )
    median_radius = float(np.median(member_to_centroid)) if len(member_to_centroid) else np.nan
    normalized = focal_to_centroid / median_radius if median_radius > 0 else np.nan
    nearest_idx = int(np.argmin(focal_to_members))

    return {
        "n_context_animals": int(len(context)),
        "context_animals": ",".join(sorted(context["animal_id"].astype(str).unique())),
        "focal_to_centroid_m": focal_to_centroid,
        "context_median_radius_m": median_radius,
        "normalized_centroid_distance": normalized,
        "centroid_centrality_percentile": float(100 * np.mean(member_to_centroid >= focal_to_centroid)),
        "nearest_member_id": str(context["animal_id"].astype(str).to_numpy()[nearest_idx]),
        "nearest_member_distance_m": float(focal_to_members[nearest_idx]),
    }


def compute_trace_metrics(data: pd.DataFrame, traces: pd.DataFrame) -> pd.DataFrame:
    if traces.empty:
        return pd.DataFrame()
    by_time = {timestamp: group.copy() for timestamp, group in data.groupby("window_start", observed=True)}
    metric_rows = []
    trace_columns = [
        "trace_id",
        "trace_type",
        "trace_start_time",
        "trace_end_time",
        "trace_duration_hours",
        "trace_n_timestamps",
        "context_label",
        "animal_id",
        "sex",
        "age",
        "origin_group",
        "dynamic_social_unit",
        "assigned_social_unit",
        "dynamic_assignment",
        "dynamic_segment_id",
        "window_start",
        "mean_latitude",
        "mean_longitude",
        "temp_group_id",
        "temp_group_size",
        "component_size",
        "n_dynamic_unit_components",
        "component_rank",
        "is_carried_night",
        "is_local_2h_supported",
        "position_support_confidence",
        "membership_confidence",
    ]
    for row in traces[trace_columns].to_dict(orient="records"):
        same_time = by_time.get(row["window_start"])
        if same_time is None:
            continue
        context = context_for_row(pd.Series(row), same_time)
        metrics = centrality_metrics_for_row(pd.Series(row), context)
        if metrics is None:
            continue
        metric_rows.append({**row, **metrics})
    if not metric_rows:
        return pd.DataFrame()
    out = pd.DataFrame(metric_rows)
    out["window_start"] = pd.to_datetime(out["window_start"])
    out["trace_start_time"] = pd.to_datetime(out["trace_start_time"])
    out["trace_end_time"] = pd.to_datetime(out["trace_end_time"])
    out["days_since_trace_start"] = (
        out["window_start"] - out["trace_start_time"]
    ).dt.total_seconds() / 86400
    out["trace_duration_days"] = out["trace_duration_hours"] / 24
    out["trace_fraction"] = np.where(
        out["trace_duration_days"].gt(0),
        out["days_since_trace_start"] / out["trace_duration_days"],
        np.nan,
    )
    return out


def summarize_traces(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    return (
        metrics.groupby(
            [
                "trace_id",
                "trace_type",
                "animal_id",
                "sex",
                "age",
                "origin_group",
                "dynamic_social_unit",
                "context_label",
            ],
            dropna=False,
            observed=True,
        )
        .agg(
            start_time=("window_start", "min"),
            end_time=("window_start", "max"),
            n_timestamps=("window_start", "nunique"),
            duration_days=("trace_duration_days", "max"),
            n_context_animals_median=("n_context_animals", "median"),
            centrality_median=("centroid_centrality_percentile", "median"),
            centrality_start=("centroid_centrality_percentile", lambda s: s.iloc[: max(1, ceil(len(s) * 0.2))].median()),
            centrality_end=("centroid_centrality_percentile", lambda s: s.iloc[-max(1, ceil(len(s) * 0.2)) :].median()),
            normalized_distance_median=("normalized_centroid_distance", "median"),
        )
        .reset_index()
    )


def fit_spline_mixedlm(
    data: pd.DataFrame,
    response: str,
    trace_type: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    use = data[
        data["trace_type"].eq(trace_type)
        & data["days_since_trace_start"].between(0, MAX_MODEL_DAYS)
    ].dropna(subset=[response, "days_since_trace_start", "trace_id"]).copy()
    use = use[np.isfinite(use[response].astype(float))]
    if use["trace_id"].nunique() < 4 or len(use) < 30:
        return pd.DataFrame(), {
            "trace_type": trace_type,
            "response": response,
            "status": "skipped_not_enough_data",
            "n_rows": int(len(use)),
            "n_traces": int(use["trace_id"].nunique()),
        }

    formula = f"bs(days_since_trace_start, df={SPLINE_DF}, degree=3, include_intercept=False)"
    x = dmatrix(formula, use, return_type="dataframe")
    y = use[response].astype(float)
    model = MixedLM(y, x, groups=use["trace_id"].astype(str))
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = model.fit(reml=True, method="lbfgs", maxiter=500, disp=False)
        status = "ok"
        warning_text = " | ".join(str(w.message) for w in caught[-5:])
    except Exception as exc:
        try:
            result = model.fit(reml=True, method="powell", maxiter=500, disp=False)
            status = "ok_powell"
            warning_text = f"lbfgs_failed: {exc}"
        except Exception as exc2:
            return pd.DataFrame(), {
                "trace_type": trace_type,
                "response": response,
                "status": f"failed: {exc2}",
                "n_rows": int(len(use)),
                "n_traces": int(use["trace_id"].nunique()),
            }

    day_grid = np.linspace(0, min(MAX_MODEL_DAYS, float(use["days_since_trace_start"].max())), 220)
    pred_base = pd.DataFrame({"days_since_trace_start": day_grid})
    pred_x = build_design_matrices([x.design_info], pred_base)[0]
    pred = np.asarray(result.predict(pred_x), dtype=float)
    predictions = pd.DataFrame(
        {
            "trace_type": trace_type,
            "response": response,
            "days_since_trace_start": day_grid,
            "prediction": pred,
            "n_model_rows": len(use),
            "n_model_traces": use["trace_id"].nunique(),
        }
    )
    summary = {
        "trace_type": trace_type,
        "response": response,
        "status": status,
        "n_rows": int(len(use)),
        "n_traces": int(use["trace_id"].nunique()),
        "n_animals": int(use["animal_id"].nunique()),
        "max_observed_day": float(use["days_since_trace_start"].max()),
        "aic": float(result.aic) if result.aic is not None else np.nan,
        "bic": float(result.bic) if result.bic is not None else np.nan,
        "scale": float(result.scale),
        "random_intercept_variance": float(result.cov_re.iloc[0, 0]) if result.cov_re.size else np.nan,
        "converged": bool(getattr(result, "converged", False)),
        "warnings": warning_text,
    }
    return predictions, summary


def fit_models(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    predictions = []
    summaries = []
    for trace_type in sorted(metrics["trace_type"].dropna().unique()):
        for response in RESPONSE_SPECS:
            pred, summary = fit_spline_mixedlm(metrics, response, trace_type)
            summaries.append(summary)
            if not pred.empty:
                predictions.append(pred)
    return (
        pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame(),
        pd.DataFrame(summaries),
    )


def binned_profile(metrics: pd.DataFrame, response: str) -> pd.DataFrame:
    use = metrics[
        metrics["days_since_trace_start"].between(0, MAX_MODEL_DAYS)
        & metrics[response].notna()
    ].copy()
    use["day_bin"] = np.floor(use["days_since_trace_start"]).astype(int)
    return (
        use.groupby(["trace_type", "day_bin"], observed=True)
        .agg(
            n_points=(response, "size"),
            n_traces=("trace_id", "nunique"),
            median=(response, "median"),
            q25=(response, lambda s: s.quantile(0.25)),
            q75=(response, lambda s: s.quantile(0.75)),
        )
        .reset_index()
    )


def plot_response(metrics: pd.DataFrame, predictions: pd.DataFrame, response: str, out_path: Path) -> None:
    trace_types = ["dispersal_join", "within_group_split_off"]
    fig, axes = plt.subplots(1, len(trace_types), figsize=(14.5, 5.2), sharey=False, constrained_layout=True)
    if len(trace_types) == 1:
        axes = [axes]
    spec = RESPONSE_SPECS[response]
    rng = np.random.default_rng(42)

    for ax, trace_type in zip(axes, trace_types):
        use = metrics[
            metrics["trace_type"].eq(trace_type)
            & metrics["days_since_trace_start"].between(0, MAX_MODEL_DAYS)
            & metrics[response].notna()
        ].copy()
        color = TRACE_COLORS.get(trace_type, "#666666")
        if use.empty:
            ax.set_title(f"{trace_type}\nno data")
            continue
        if len(use) > RAW_PLOT_SAMPLE:
            use_plot = use.iloc[rng.choice(len(use), RAW_PLOT_SAMPLE, replace=False)].copy()
        else:
            use_plot = use
        ax.scatter(
            use_plot["days_since_trace_start"],
            use_plot[response],
            s=9,
            alpha=0.08,
            color=color,
            linewidth=0,
        )
        prof = binned_profile(use, response)
        prof = prof[prof["trace_type"].eq(trace_type)]
        ax.fill_between(prof["day_bin"], prof["q25"], prof["q75"], color=color, alpha=0.16, linewidth=0)
        ax.plot(prof["day_bin"], prof["median"], color=color, linewidth=1.6, alpha=0.85, label="daily median")
        pred = predictions[
            predictions["trace_type"].eq(trace_type)
            & predictions["response"].eq(response)
        ]
        if not pred.empty:
            pred_y = pred["prediction"].astype(float).copy()
            if response == "centroid_centrality_percentile":
                pred_y = pred_y.clip(lower=0, upper=100)
            elif "distance" in response:
                pred_y = pred_y.clip(lower=0)
            ax.plot(
                pred["days_since_trace_start"],
                pred_y,
                color="#111111",
                linewidth=2.8,
                label="spline mixed model",
            )
        ax.set_title(
            f"{trace_type.replace('_', ' ')}\n"
            f"{use['trace_id'].nunique()} traces, {use['animal_id'].nunique()} animals",
            fontsize=11,
        )
        ax.set_xlabel("Days since trace start")
        ax.grid(True, color="#e8e8e8", linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if spec["ylim"] is not None:
            ax.set_ylim(*spec["ylim"])
        ax.legend(loc="best", frameon=False, fontsize=9)
    axes[0].set_ylabel(spec["label"])
    fig.suptitle(
        f"{spec['label']} after joining/splitting\n{spec['interpretation']}; random intercept per trace",
        fontsize=14,
        fontweight="bold",
    )
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metrics_path = OUT_DIR / "canonical_integration_trace_metrics.parquet"
    trace_rows_path = OUT_DIR / "canonical_integration_trace_rows.csv"
    if metrics_path.exists() and trace_rows_path.exists():
        metrics = pd.read_parquet(metrics_path)
        traces = pd.read_csv(trace_rows_path)
    else:
        data = add_split_component_flags(load_canonical())
        dispersal = build_dispersal_traces(data)
        split = build_split_traces(data)
        traces = pd.concat([dispersal, split], ignore_index=True)
        traces.to_csv(trace_rows_path, index=False)
        metrics = compute_trace_metrics(data, traces)

    metrics["log1p_normalized_centroid_distance"] = np.log1p(
        metrics["normalized_centroid_distance"].astype(float).clip(lower=0)
    )
    metrics.to_csv(OUT_DIR / "canonical_integration_trace_metrics.csv", index=False)
    metrics.to_parquet(metrics_path, index=False)

    trace_summary = summarize_traces(metrics)
    trace_summary.to_csv(OUT_DIR / "canonical_integration_trace_summary.csv", index=False)

    predictions, model_summary = fit_models(metrics)
    predictions.to_csv(OUT_DIR / "canonical_integration_gamm_predictions.csv", index=False)
    model_summary.to_csv(OUT_DIR / "canonical_integration_gamm_model_summary.csv", index=False)

    for response in RESPONSE_SPECS:
        plot_response(
            metrics,
            predictions,
            response,
            OUT_DIR / f"canonical_integration_gamm_{response}.png",
        )

    metadata = {
        "canonical": str(CANONICAL),
        "output_dir": str(OUT_DIR),
        "trace_types": {
            "dispersal_join": "dynamic_assignment == sustained_non_origin_association; context is same temp cluster and same dynamic social unit",
            "within_group_split_off": "animal is in a non-main temp_group component while its dynamic_social_unit is split across multiple components; context is all same dynamic_social_unit animals at that timestamp",
        },
        "responses": RESPONSE_SPECS,
        "model": f"MixedLM(response ~ cubic B-spline(days_since_trace_start, df={SPLINE_DF}) + random intercept by trace_id)",
        "max_gap_hours": MAX_GAP_HOURS,
        "min_split_trace_hours": MIN_SPLIT_TRACE_HOURS,
        "min_split_trace_timestamps": MIN_SPLIT_TRACE_TIMESTAMPS,
        "min_context_animals": MIN_CONTEXT_ANIMALS,
        "max_model_days": MAX_MODEL_DAYS,
        "n_trace_rows": int(len(traces)),
        "n_metric_rows": int(len(metrics)),
        "n_traces": int(metrics["trace_id"].nunique()) if not metrics.empty else 0,
    }
    (OUT_DIR / "canonical_integration_gamm_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print(f"Wrote metrics: {OUT_DIR / 'canonical_integration_trace_metrics.csv'}")
    print(f"Wrote trace summary: {OUT_DIR / 'canonical_integration_trace_summary.csv'}")
    print(f"Wrote model summary: {OUT_DIR / 'canonical_integration_gamm_model_summary.csv'}")
    print(model_summary.to_string(index=False))
    print(trace_summary.groupby("trace_type").agg(n_traces=("trace_id", "nunique"), n_animals=("animal_id", "nunique"), n_rows=("n_timestamps", "sum")).to_string())


if __name__ == "__main__":
    main()
