from __future__ import annotations

import html
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning


SOURCE_DIR = Path(
    r"C:\Users\rharel\Documents\New project\outputs\hourly_grouping_validation_shared_20260703\group_merge_mixing_dynamics"
)
HOURLY_CSV = SOURCE_DIR / "bigmerge_2min_vs_hourly_5m_no_copper_lilac_values_by_dyad_hourly_5m_metric_rows.csv"
TWOMIN_CSV = SOURCE_DIR / "bigmerge_2min_vs_hourly_5m_no_copper_lilac_2min_metric_rows.csv"
OUT_DIR = Path("outputs") / "bigmerge_5m_long_events_2min_days"

MAX_WITHIN_DAY_GAP_H = 2.5
MAX_OVERNIGHT_GAP_H = 15.0
LONG_EVENT_THRESHOLD_H = 12.0
N_SPLINE_DF = 5

METRICS = [
    ("edge_modularity_q", "Edge modularity Q", "higher = origin groups remain spatially/socially segregated"),
    ("cross_edge_fraction", "Cross-group edge fraction", "higher = more direct 5 m contacts across origin groups"),
    ("composition_entropy_norm", "Composition entropy", "higher = cluster composition is more mixed/balanced"),
    ("pair_balance", "Pair balance", "higher = the two focal groups are similarly represented"),
]


def parse_time(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, format="mixed", errors="coerce")


def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    hourly = pd.read_csv(HOURLY_CSV)
    for column in ["timestamp", "merge_start_time", "merge_end_time"]:
        if column in hourly:
            hourly[column] = parse_time(hourly[column])
    hourly = hourly.dropna(subset=["timestamp", "pair_key", "merge_episode_id"])
    hourly = hourly.drop_duplicates(["pair_key", "merge_episode_id", "timestamp"])

    two_min = pd.read_csv(TWOMIN_CSV)
    for column in ["timestamp", "bin_2min"]:
        if column in two_min:
            two_min[column] = parse_time(two_min[column])
    two_min = two_min.dropna(subset=["timestamp", "bin_2min", "pair_key", "merge_episode_id"])
    return hourly, two_min


def assign_long_events(hourly: pd.DataFrame) -> pd.DataFrame:
    chunks = []
    long_index = 0
    for pair_key, group in hourly.sort_values(["pair_key", "timestamp"]).groupby("pair_key", observed=True):
        group = group.sort_values("timestamp").copy()
        previous_time = group["timestamp"].shift()
        gap_h = (group["timestamp"] - previous_time).dt.total_seconds() / 3600.0
        same_clock_episode = gap_h.le(MAX_WITHIN_DAY_GAP_H)
        overnight_bridge = (
            previous_time.dt.hour.ge(15)
            & group["timestamp"].dt.hour.le(3)
            & gap_h.le(MAX_OVERNIGHT_GAP_H)
            & (group["timestamp"].dt.date > previous_time.dt.date)
        )
        starts_new = gap_h.isna() | ~(same_clock_episode | overnight_bridge)
        local_event = starts_new.cumsum()
        for _, event in group.groupby(local_event, sort=True):
            long_index += 1
            event = event.copy()
            start = event["timestamp"].min()
            end = event["timestamp"].max()
            event["long_merge_episode_id"] = f"LONGMERGE_{long_index:05d}"
            event["long_merge_start_time"] = start
            event["long_merge_end_time"] = end
            event["long_merge_elapsed_hours"] = (end - start).total_seconds() / 3600.0 + 1.0
            event["long_merge_elapsed_days"] = event["long_merge_elapsed_hours"] / 24.0
            event["long_merge_observed_hours"] = event["timestamp"].nunique()
            event["long_merge_is_over_12h"] = event["long_merge_elapsed_hours"].gt(LONG_EVENT_THRESHOLD_H)
            event["long_hours_since_start"] = (event["timestamp"] - start).dt.total_seconds() / 3600.0
            event["long_days_since_start"] = event["long_hours_since_start"] / 24.0
            chunks.append(event)
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()


def attach_long_events(two_min: pd.DataFrame, hourly_long: pd.DataFrame) -> pd.DataFrame:
    mapping_columns = [
        "pair_key",
        "merge_episode_id",
        "timestamp",
        "long_merge_episode_id",
        "long_merge_start_time",
        "long_merge_end_time",
        "long_merge_elapsed_hours",
        "long_merge_elapsed_days",
        "long_merge_observed_hours",
        "long_merge_is_over_12h",
    ]
    data = two_min.merge(hourly_long[mapping_columns], on=["pair_key", "merge_episode_id", "timestamp"], how="inner")
    data["long_hours_since_start"] = (
        data["bin_2min"] - data["long_merge_start_time"]
    ).dt.total_seconds() / 3600.0
    data["long_days_since_start"] = data["long_hours_since_start"] / 24.0
    data = data[data["long_hours_since_start"].ge(0)].copy()
    return data


def summarize_long_events(hourly_long: pd.DataFrame, two_min_long: pd.DataFrame) -> pd.DataFrame:
    hourly_summary = (
        hourly_long.groupby("long_merge_episode_id", observed=True)
        .agg(
            pair_key=("pair_key", "first"),
            start=("long_merge_start_time", "first"),
            end=("long_merge_end_time", "first"),
            elapsed_hours=("long_merge_elapsed_hours", "first"),
            elapsed_days=("long_merge_elapsed_days", "first"),
            observed_hours=("long_merge_observed_hours", "first"),
            old_daily_episode_count=("merge_episode_id", "nunique"),
            over_12h=("long_merge_is_over_12h", "first"),
        )
        .reset_index()
    )
    bin_summary = (
        two_min_long.groupby("long_merge_episode_id", observed=True)
        .agg(
            two_min_metric_rows=("bin_2min", "size"),
            first_2min_bin=("bin_2min", "min"),
            last_2min_bin=("bin_2min", "max"),
        )
        .reset_index()
    )
    return hourly_summary.merge(bin_summary, on="long_merge_episode_id", how="left").sort_values(
        ["elapsed_hours", "observed_hours"], ascending=False
    )


def fit_mixed_spline(data: pd.DataFrame, metric: str) -> tuple[sm.regression.mixed_linear_model.MixedLMResultsWrapper, str]:
    model_df = data[["pair_key", "long_days_since_start", metric]].dropna().copy()
    model_df = model_df.rename(columns={metric: "y"})
    formula = f"y ~ bs(long_days_since_start, df={N_SPLINE_DF}, degree=3, include_intercept=False)"
    model = sm.MixedLM.from_formula(formula, groups="pair_key", re_formula="1", data=model_df)
    fit_note = "powell"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        result = model.fit(reml=False, method="powell", maxiter=1000, disp=False)
        if any(isinstance(w.message, ConvergenceWarning) for w in caught):
            fit_note += "; convergence warning"
    return result, fit_note


def predict_fixed_effect(result: sm.regression.mixed_linear_model.MixedLMResultsWrapper, x_grid: np.ndarray) -> pd.DataFrame:
    grid = pd.DataFrame({"long_days_since_start": x_grid})
    design = patsy.build_design_matrices([result.model.data.design_info], grid)[0]
    exog = np.asarray(design)
    beta = result.fe_params.to_numpy()
    cov_beta = result.cov_params().loc[result.fe_params.index, result.fe_params.index].to_numpy()
    pred = exog @ beta
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", exog, cov_beta, exog), 0))
    return pd.DataFrame(
        {
            "long_days_since_start": x_grid,
            "predicted": pred,
            "ci_low": pred - 1.96 * se,
            "ci_high": pred + 1.96 * se,
        }
    )


def fit_all(two_min_long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = []
    predictions = []
    for metric, _, _ in METRICS:
        model_df = two_min_long[two_min_long[metric].notna()].copy()
        result, note = fit_mixed_spline(model_df, metric)
        random_var = float(result.cov_re.iloc[0, 0]) if result.cov_re.size else np.nan
        resid_var = float(result.scale)
        summaries.append(
            {
                "metric": metric,
                "n_rows": len(model_df),
                "n_dyads": model_df["pair_key"].nunique(),
                "n_long_events": model_df["long_merge_episode_id"].nunique(),
                "n_long_events_over_12h": model_df.loc[
                    model_df["long_merge_is_over_12h"], "long_merge_episode_id"
                ].nunique(),
                "x_max_days": model_df["long_days_since_start"].max(),
                "aic": result.aic,
                "bic": result.bic,
                "random_intercept_variance_dyad": random_var,
                "residual_variance": resid_var,
                "dyad_icc": random_var / (random_var + resid_var) if random_var + resid_var > 0 else np.nan,
                "fit_note": note,
            }
        )
        x_grid = np.linspace(
            model_df["long_days_since_start"].quantile(0.01),
            model_df["long_days_since_start"].quantile(0.99),
            180,
        )
        pred = predict_fixed_effect(result, x_grid)
        pred["metric"] = metric
        predictions.append(pred)
    return pd.DataFrame(summaries), pd.concat(predictions, ignore_index=True)


def plot_predictions(two_min_long: pd.DataFrame, predictions: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), sharex=True)
    axes = axes.ravel()
    for ax, (metric, label, note) in zip(axes, METRICS):
        raw = two_min_long[["long_days_since_start", metric]].dropna()
        if len(raw) > 3500:
            raw = raw.sample(3500, random_state=42)
        ax.scatter(raw["long_days_since_start"], raw[metric], s=5, color="#1f77b4", alpha=0.06, linewidth=0)
        metric_pred = predictions[predictions["metric"].eq(metric)]
        ax.plot(metric_pred["long_days_since_start"], metric_pred["predicted"], color="#1f77b4", linewidth=2.5)
        ax.fill_between(
            metric_pred["long_days_since_start"],
            metric_pred["ci_low"],
            metric_pred["ci_high"],
            color="#1f77b4",
            alpha=0.18,
        )
        ax.axvline(0.5, color="#555555", linestyle=":", linewidth=1.1)
        ax.text(0.505, 0.03, "12 h", transform=ax.get_xaxis_transform(), fontsize=9, color="#555555")
        ax.set_title(label, fontsize=13, fontweight="bold")
        ax.text(0.01, 0.98, note, transform=ax.transAxes, va="top", ha="left", fontsize=9, color="#555")
        ax.grid(True, color="#e8e8e8", linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if metric != "edge_modularity_q":
            ax.set_ylim(-0.05, 1.05)
    fig.suptitle(
        "2-minute 5 m metrics aligned to long hourly-defined merge events\n"
        "Hourly social grouping defines event continuity; 2-minute fixes define metric values",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )
    fig.supxlabel("Days since long merge event started", fontsize=12)
    fig.supylabel("Metric value", fontsize=12)
    fig.tight_layout(rect=[0.04, 0.04, 1.0, 0.92])
    out = OUT_DIR / "bigmerge_5m_long_events_2min_gamm_days.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def plot_event_durations(event_summary: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    bins = np.arange(0, np.ceil(event_summary["elapsed_hours"].max() / 6) * 6 + 6, 6)
    ax.hist(event_summary["elapsed_hours"], bins=bins, color="#4c78a8", alpha=0.85)
    ax.axvline(LONG_EVENT_THRESHOLD_H, color="#222222", linewidth=2)
    ax.text(LONG_EVENT_THRESHOLD_H + 0.5, ax.get_ylim()[1] * 0.9, "12 h", ha="left", va="top")
    ax.set_xlabel("Long hourly-defined event duration (clock hours)")
    ax.set_ylabel("Number of events")
    ax.set_title("Event durations after bridging overnight continuity")
    ax.grid(True, axis="y", color="#e8e8e8")
    out = OUT_DIR / "bigmerge_5m_long_event_duration_distribution.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def write_html(image_paths: list[Path], summary: pd.DataFrame, event_summary: pd.DataFrame) -> Path:
    out = OUT_DIR / "bigmerge_5m_long_events_2min_days_dashboard.html"
    model_rows = "\n".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in ["metric", "n_rows", "n_dyads", "n_long_events", "n_long_events_over_12h", "x_max_days", "dyad_icc", "fit_note"]
        )
        + "</tr>"
        for _, row in summary.round({"x_max_days": 3, "dyad_icc": 3}).iterrows()
    )
    event_note = (
        f"{len(event_summary):,} long events; "
        f"{int(event_summary['over_12h'].sum()):,} longer than 12 h; "
        f"max duration {event_summary['elapsed_hours'].max():.1f} h."
    )
    imgs = "\n".join(f"<img src='{html.escape(path.name)}' alt='{html.escape(path.stem)}'>" for path in image_paths)
    out.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Long event 2-min metrics</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;color:#222;background:#fff}}
.wrap{{padding:24px 30px 44px;max-width:1500px}}
.note{{line-height:1.45;color:#555;max-width:1100px}}
img{{max-width:100%;display:block;margin:18px 0 34px;border:1px solid #ddd;border-radius:6px}}
table{{border-collapse:collapse;margin:16px 0 32px;font-size:13px}}
th,td{{border:1px solid #ddd;padding:6px 8px;text-align:left}}
th{{background:#f2f2f2}}
</style></head><body><div class="wrap">
<h1>Long hourly-defined merge events with 2-minute 5 m metrics</h1>
<p class="note">Hourly social grouping defines which dyad is together and when an event remains continuous.
Events are bridged across the night only when the same dyad is present at the late-day hour and again at 03:00 the next day.
The four metric values themselves are calculated from the 2-minute fixes inside those hourly-defined clusters.</p>
<p class="note">{html.escape(event_note)}</p>
<table>
<tr><th>metric</th><th>2-min rows</th><th>dyads</th><th>long events</th><th>events &gt;12 h</th><th>max x days</th><th>dyad ICC</th><th>fit note</th></tr>
{model_rows}
</table>
{imgs}
</div></body></html>""",
        encoding="utf-8",
    )
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    hourly, two_min = read_inputs()
    hourly_long = assign_long_events(hourly)
    two_min_long = attach_long_events(two_min, hourly_long)
    event_summary = summarize_long_events(hourly_long, two_min_long)

    hourly_long.to_csv(OUT_DIR / "bigmerge_5m_hourly_long_event_assignments.csv", index=False)
    two_min_long.to_csv(OUT_DIR / "bigmerge_5m_2min_metrics_with_long_event_days.csv", index=False)
    event_summary.to_csv(OUT_DIR / "bigmerge_5m_long_event_summary.csv", index=False)
    event_summary[event_summary["over_12h"]].to_csv(OUT_DIR / "bigmerge_5m_long_events_over_12h.csv", index=False)

    summary, predictions = fit_all(two_min_long)
    summary.to_csv(OUT_DIR / "bigmerge_5m_long_events_2min_gamm_summary.csv", index=False)
    predictions.to_csv(OUT_DIR / "bigmerge_5m_long_events_2min_gamm_predictions.csv", index=False)

    image_paths = [plot_predictions(two_min_long, predictions), plot_event_durations(event_summary)]
    dashboard = write_html(image_paths, summary, event_summary)
    print(f"Wrote dashboard: {dashboard.resolve()}")
    print(f"Wrote event summary: {(OUT_DIR / 'bigmerge_5m_long_event_summary.csv').resolve()}")
    print(f"Wrote >12h events: {(OUT_DIR / 'bigmerge_5m_long_events_over_12h.csv').resolve()}")
    for path in image_paths:
        print(path.resolve())


if __name__ == "__main__":
    main()
