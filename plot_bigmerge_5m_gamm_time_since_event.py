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

PREVIOUS_OUT_DIR = Path("outputs") / "bigmerge_5m_gamm_time_since_event"
OUT_DIR = Path("outputs") / "bigmerge_5m_gamm_time_since_event_2min_only"

METRICS = [
    ("edge_modularity_q", "Edge modularity Q", "higher = origin groups remain spatially/socially segregated"),
    ("cross_edge_fraction", "Cross-group edge fraction", "higher = more direct 5 m contacts across origin groups"),
    ("composition_entropy_norm", "Composition entropy", "higher = cluster composition is more mixed/balanced"),
    ("pair_balance", "Pair balance", "higher = the two focal groups are similarly represented"),
]

SCALE_ORDER = ["2-min"]
SCALE_COLORS = {"hourly": "#222222", "2-min": "#1f77b4"}
N_SPLINE_DF = 5


def load_metric_rows() -> pd.DataFrame:
    cached_rows = PREVIOUS_OUT_DIR / "bigmerge_5m_gamm_model_input_rows.csv"
    if cached_rows.exists():
        rows = pd.read_csv(cached_rows, parse_dates=["timestamp"])
        return rows[rows["scale"].isin(SCALE_ORDER)].copy()

    hourly = pd.read_csv(HOURLY_CSV, parse_dates=["timestamp", "merge_start_time"])
    hourly["scale"] = "hourly"
    hourly["time_since_event_start_h"] = pd.to_numeric(hourly["hours_since_merge_start"], errors="coerce")

    two_min = pd.read_csv(TWOMIN_CSV, parse_dates=["timestamp", "bin_2min"])
    start_lookup = (
        hourly[["merge_episode_id", "merge_start_time"]]
        .dropna()
        .drop_duplicates("merge_episode_id")
    )
    two_min = two_min.merge(start_lookup, on="merge_episode_id", how="left")
    two_min["scale"] = "2-min"
    two_min["time_since_event_start_h"] = (
        two_min["bin_2min"] - two_min["merge_start_time"]
    ).dt.total_seconds() / 3600.0

    common_cols = [
        "scale",
        "pair_key",
        "merge_episode_id",
        "timestamp",
        "time_since_event_start_h",
        "edge_modularity_q",
        "cross_edge_fraction",
        "composition_entropy_norm",
        "pair_balance",
        "pair_n",
        "cluster_size_total",
        "total_edges",
    ]
    rows = pd.concat([hourly[common_cols], two_min[common_cols]], ignore_index=True)
    rows = rows.replace([np.inf, -np.inf], np.nan)
    rows = rows[rows["pair_key"].notna() & rows["time_since_event_start_h"].notna()]
    rows = rows[rows["time_since_event_start_h"].ge(0)]
    rows = rows[rows["scale"].isin(SCALE_ORDER)]
    return rows


def fit_mixed_spline(data: pd.DataFrame, metric: str) -> tuple[sm.regression.mixed_linear_model.MixedLMResultsWrapper, str]:
    model_df = data[["pair_key", "time_since_event_start_h", metric]].dropna().copy()
    model_df = model_df.rename(columns={metric: "y"})
    formula = f"y ~ bs(time_since_event_start_h, df={N_SPLINE_DF}, degree=3, include_intercept=False)"
    model = sm.MixedLM.from_formula(formula, groups="pair_key", re_formula="1", data=model_df)

    fit_note = "powell"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        try:
            result = model.fit(reml=False, method="powell", maxiter=1000, disp=False)
        except Exception:
            fit_note = "lbfgs fallback"
            result = model.fit(reml=False, method="lbfgs", maxiter=1000, disp=False)
        if any(isinstance(w.message, ConvergenceWarning) for w in caught):
            fit_note += "; convergence warning"
    return result, fit_note


def predict_fixed_effect(
    result: sm.regression.mixed_linear_model.MixedLMResultsWrapper,
    x_grid: np.ndarray,
) -> pd.DataFrame:
    grid = pd.DataFrame({"time_since_event_start_h": x_grid})
    design = patsy.build_design_matrices([result.model.data.design_info], grid)[0]
    exog = np.asarray(design)
    beta = result.fe_params.to_numpy()
    cov_beta = result.cov_params().loc[result.fe_params.index, result.fe_params.index].to_numpy()
    pred = exog @ beta
    se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", exog, cov_beta, exog), 0))
    return pd.DataFrame(
        {
            "time_since_event_start_h": x_grid,
            "predicted": pred,
            "ci_low": pred - 1.96 * se,
            "ci_high": pred + 1.96 * se,
        }
    )


def summarize_fit(
    data: pd.DataFrame,
    metric: str,
    scale: str,
    result: sm.regression.mixed_linear_model.MixedLMResultsWrapper,
    fit_note: str,
) -> dict[str, object]:
    model_df = data[["pair_key", "merge_episode_id", "time_since_event_start_h", metric]].dropna()
    random_var = float(result.cov_re.iloc[0, 0]) if result.cov_re.size else np.nan
    resid_var = float(result.scale)
    return {
        "metric": metric,
        "scale": scale,
        "n_rows": len(model_df),
        "n_dyads": model_df["pair_key"].nunique(),
        "n_merge_episodes": model_df["merge_episode_id"].nunique(),
        "x_min_h": model_df["time_since_event_start_h"].min(),
        "x_max_h": model_df["time_since_event_start_h"].max(),
        "aic": result.aic,
        "bic": result.bic,
        "log_likelihood": result.llf,
        "random_intercept_variance_dyad": random_var,
        "residual_variance": resid_var,
        "dyad_icc": random_var / (random_var + resid_var) if np.isfinite(random_var + resid_var) and (random_var + resid_var) > 0 else np.nan,
        "fit_note": fit_note,
    }


def fit_all(rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, object]] = []
    predictions: list[pd.DataFrame] = []
    coefficients: list[pd.DataFrame] = []

    for metric, _, _ in METRICS:
        for scale in SCALE_ORDER:
            subset = rows[rows["scale"].eq(scale)].copy()
            subset = subset[subset[metric].notna()]
            if subset["pair_key"].nunique() < 2 or len(subset) < 20:
                continue

            result, fit_note = fit_mixed_spline(subset, metric)
            summaries.append(summarize_fit(subset, metric, scale, result, fit_note))

            x_grid = np.linspace(
                subset["time_since_event_start_h"].quantile(0.01),
                subset["time_since_event_start_h"].quantile(0.99),
                160,
            )
            pred = predict_fixed_effect(result, x_grid)
            pred["metric"] = metric
            pred["scale"] = scale
            predictions.append(pred)

            coef = pd.DataFrame(
                {
                    "term": result.fe_params.index,
                    "estimate": result.fe_params.values,
                    "se": result.bse_fe.values,
                    "metric": metric,
                    "scale": scale,
                }
            )
            coefficients.append(coef)

    return (
        pd.DataFrame(summaries),
        pd.concat(predictions, ignore_index=True),
        pd.concat(coefficients, ignore_index=True),
    )


def plot_predictions(rows: pd.DataFrame, pred: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), sharex=True)
    axes = axes.ravel()

    for ax, (metric, label, note) in zip(axes, METRICS):
        raw = rows[["scale", "pair_key", "time_since_event_start_h", metric]].dropna()
        raw = raw[raw["time_since_event_start_h"].between(0, 14)]
        for scale in SCALE_ORDER:
            sub_raw = raw[raw["scale"].eq(scale)]
            if not sub_raw.empty:
                if len(sub_raw) > 2500:
                    sub_raw = sub_raw.sample(2500, random_state=42)
                ax.scatter(
                    sub_raw["time_since_event_start_h"],
                    sub_raw[metric],
                    s=5,
                    color=SCALE_COLORS[scale],
                    alpha=0.07,
                    linewidth=0,
                )

            sub_pred = pred[(pred["metric"].eq(metric)) & (pred["scale"].eq(scale))]
            if not sub_pred.empty:
                color = SCALE_COLORS[scale]
                ax.plot(
                    sub_pred["time_since_event_start_h"],
                    sub_pred["predicted"],
                    color=color,
                    linewidth=2.4,
                    label=scale,
                )
                ax.fill_between(
                    sub_pred["time_since_event_start_h"],
                    sub_pred["ci_low"],
                    sub_pred["ci_high"],
                    color=color,
                    alpha=0.16,
                    linewidth=0,
                )

        ax.set_title(label, fontsize=13, fontweight="bold")
        ax.text(0.01, 0.98, note, transform=ax.transAxes, va="top", ha="left", fontsize=9, color="#555")
        ax.grid(True, color="#e8e8e8", linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if metric != "edge_modularity_q":
            ax.set_ylim(-0.05, 1.05)

    axes[0].legend(loc="lower right", frameon=True, fontsize=10)
    fig.suptitle(
        "2-minute 5 m big-merge GAMM-like smooths over time since event start\n"
        "Cubic spline fixed effect + random intercept for group dyad; Copper-Lilac excluded",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )
    fig.supxlabel("Hours since merge episode started", fontsize=12)
    fig.supylabel("Metric value", fontsize=12)
    fig.tight_layout(rect=[0.04, 0.04, 1.0, 0.92])

    out = OUT_DIR / "bigmerge_5m_gamm_2min_time_since_event_predictions.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def plot_dyad_random_effects(rows: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.ravel()

    for ax, (metric, label, _) in zip(axes, METRICS):
        summary = (
            rows[rows["scale"].eq("2-min")]
            .groupby("pair_key", observed=True)
            .agg(mean_metric=(metric, "mean"), n=("pair_key", "size"))
            .dropna()
            .sort_values("mean_metric")
        )
        y = np.arange(len(summary))
        ax.scatter(summary["mean_metric"], y, s=np.clip(summary["n"] * 2.5, 20, 140), color="#4c78a8", alpha=0.75)
        ax.set_yticks(y)
        ax.set_yticklabels(summary.index, fontsize=8)
        ax.set_title(f"{label}\nraw 2-minute dyad means", fontsize=12, fontweight="bold")
        ax.grid(True, axis="x", color="#e8e8e8", linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Dyad-level differences behind the random intercept", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = OUT_DIR / "bigmerge_5m_gamm_2min_dyad_metric_means.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def write_html(image_paths: list[Path], summary: pd.DataFrame) -> Path:
    out = OUT_DIR / "bigmerge_5m_gamm_2min_time_since_event_dashboard.html"
    rows = "\n".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[col]))}</td>"
            for col in [
                "metric",
                "scale",
                "n_rows",
                "n_dyads",
                "n_merge_episodes",
                "x_max_h",
                "dyad_icc",
                "fit_note",
            ]
        )
        + "</tr>"
        for _, row in summary.round({"x_max_h": 2, "dyad_icc": 3}).iterrows()
    )
    imgs = "\n".join(f"<img src='{html.escape(path.name)}' alt='{html.escape(path.stem)}'>" for path in image_paths)
    out.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>5 m big-merge GAMM</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;background:#fff;color:#222}}
.wrap{{padding:24px 30px 44px;max-width:1500px}}
.note{{line-height:1.45;color:#555;max-width:1050px}}
img{{max-width:100%;display:block;margin:18px 0 34px;border:1px solid #ddd;border-radius:6px}}
table{{border-collapse:collapse;margin:16px 0 32px;font-size:13px}}
th,td{{border:1px solid #ddd;padding:6px 8px;text-align:left}}
th{{background:#f2f2f2}}
</style></head><body><div class="wrap">
<h1>2-minute 5 m big-merge GAMM-like models</h1>
<p class="note">Each model is fit for the 2-minute metric rows:
metric ~ cubic spline(hours since merge episode started) + random intercept for group dyad.
The curves are marginal fixed-effect predictions; the dyad random intercept captures stable differences among group-pairs.</p>
<table>
<tr><th>metric</th><th>scale</th><th>n rows</th><th>dyads</th><th>episodes</th><th>max x h</th><th>dyad ICC</th><th>fit note</th></tr>
{rows}
</table>
{imgs}
</div></body></html>""",
        encoding="utf-8",
    )
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_metric_rows()
    rows.to_csv(OUT_DIR / "bigmerge_5m_gamm_model_input_rows.csv", index=False)

    summary, predictions, coefficients = fit_all(rows)
    summary.to_csv(OUT_DIR / "bigmerge_5m_gamm_model_summary.csv", index=False)
    predictions.to_csv(OUT_DIR / "bigmerge_5m_gamm_predictions.csv", index=False)
    coefficients.to_csv(OUT_DIR / "bigmerge_5m_gamm_fixed_effect_coefficients.csv", index=False)

    image_paths = [plot_predictions(rows, predictions), plot_dyad_random_effects(rows)]
    dashboard = write_html(image_paths, summary)

    print(f"Wrote dashboard: {dashboard.resolve()}")
    print(f"Wrote summary: {(OUT_DIR / 'bigmerge_5m_gamm_model_summary.csv').resolve()}")
    for path in image_paths:
        print(path.resolve())


if __name__ == "__main__":
    main()
