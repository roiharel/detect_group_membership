from __future__ import annotations

import html
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


NEW_PROJECT = Path(r"C:\Users\rharel\Documents\New project")
sys.path.insert(0, str(NEW_PROJECT))

from plot_group_merge_mixing_dynamics import read_status  # noqa: E402
from plot_merge_2min_vs_hourly_5m import add_sampling_controls, load_pair_timepoints  # noqa: E402
from plot_merge_edge_radius_sensitivity import recompute_radius_metrics  # noqa: E402


OUT_DIR = (
    NEW_PROJECT
    / "outputs"
    / "hourly_grouping_validation_shared_20260703"
    / "group_merge_mixing_dynamics"
)
TWO_MIN_CSV = OUT_DIR / "bigmerge_2min_vs_hourly_5m_no_copper_lilac_2min_metric_rows.csv"
TIMEPOINTS_CSV = OUT_DIR / "hourly_merge_mixing_eps100m_minpair2_minep6h_bigmerge_medianfrac0p5_timepoints.csv"
STATUS_FILE = NEW_PROJECT / "outputs" / "rerun_20260703_full_hierarchical_sampling_filtered_eps300" / "proximity_status_full.parquet"
OUTPUT_PREFIX = "bigmerge_2min_vs_hourly_5m_no_copper_lilac_values_by_dyad"
EDGE_RADIUS_M = 5.0
EXCLUDED_PAIRS = {"Copper - Lilac"}

METRICS = [
    ("edge_modularity_q", "Edge modularity Q", "higher = more within-origin close edges / more segregated", None),
    ("cross_edge_fraction", "Cross-group edge fraction", "higher = more 5 m contacts across origin groups", (0, 1)),
    ("composition_entropy_norm", "Composition entropy", "higher = more balanced multi-group composition", (0, 1.05)),
    ("pair_balance", "Pair balance", "higher = the two focal groups are similarly represented", (0, 1.05)),
]


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    two_min = pd.read_csv(TWO_MIN_CSV, parse_dates=["timestamp", "bin_2min"])
    two_min["scale"] = "2-min"

    base = load_pair_timepoints(TIMEPOINTS_CSV, EXCLUDED_PAIRS)
    status = read_status(STATUS_FILE)
    hourly = add_sampling_controls(recompute_radius_metrics(base.copy(), status, EDGE_RADIUS_M))
    hourly["scale"] = "hourly"
    return hourly, two_min


def dyad_order(hourly: pd.DataFrame, two_min: pd.DataFrame) -> list[str]:
    counts = (
        pd.concat(
            [
                hourly.groupby("pair_key").size().rename("hourly_rows"),
                two_min.groupby("pair_key")["hour_row_id"].nunique().rename("represented_hours"),
            ],
            axis=1,
        )
        .fillna(0)
        .assign(total=lambda x: x["hourly_rows"] + x["represented_hours"])
        .sort_values(["total", "represented_hours"], ascending=False)
    )
    return counts.index.astype(str).tolist()


def summarize_2min(two_min: pd.DataFrame, metric: str) -> pd.DataFrame:
    out = two_min.copy()
    out["day_bin"] = np.floor(out["cumulative_pair_observed_days"].astype(float)).astype(int)
    return (
        out.groupby(["pair_key", "day_bin"], observed=True)
        .agg(
            x=("cumulative_pair_observed_days", "median"),
            median=(metric, "median"),
            q25=(metric, lambda s: s.quantile(0.25)),
            q75=(metric, lambda s: s.quantile(0.75)),
            n_bins=(metric, "size"),
        )
        .reset_index()
    )


def plot_metric(hourly: pd.DataFrame, two_min: pd.DataFrame, metric: str, label: str, note: str, ylim: tuple[float, float] | None) -> Path:
    pairs = dyad_order(hourly, two_min)
    ncols = 4
    nrows = int(np.ceil(len(pairs) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.25 * ncols, 2.75 * nrows),
        sharex=False,
        sharey=True if ylim is not None else False,
        squeeze=False,
    )
    for ax in axes.flat:
        ax.set_visible(False)

    two_summary = summarize_2min(two_min, metric)
    for ax, pair in zip(axes.flat, pairs):
        ax.set_visible(True)
        h = hourly[hourly["pair_key"].eq(pair)].sort_values("cumulative_pair_observed_days")
        t = two_min[two_min["pair_key"].eq(pair)]
        ts = two_summary[two_summary["pair_key"].eq(pair)].sort_values("x")

        if not t.empty:
            if len(t) > 2200:
                t_plot = t.sample(2200, random_state=42)
            else:
                t_plot = t
            ax.scatter(
                t_plot["cumulative_pair_observed_days"],
                t_plot[metric],
                s=7,
                alpha=0.08,
                color="#4c78a8",
                linewidth=0,
                label="2-min bins",
            )
        if not ts.empty:
            ax.fill_between(ts["x"], ts["q25"], ts["q75"], color="#4c78a8", alpha=0.16, linewidth=0)
            ax.plot(ts["x"], ts["median"], color="#1f5f99", linewidth=1.8, label="2-min daily median")
        if not h.empty:
            ax.plot(
                h["cumulative_pair_observed_days"],
                h[metric],
                color="#222222",
                linewidth=1.0,
                alpha=0.65,
                marker="o",
                markersize=2.3,
                label="hourly median",
            )

        ax.set_title(f"{pair}\n{len(h):,} h; {t['hour_row_id'].nunique() if not t.empty else 0:,} h w/2-min", fontsize=9)
        ax.grid(True, color="#e8e8e8", linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=8)
        if ylim is not None:
            ax.set_ylim(*ylim)

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=5, color="#4c78a8", alpha=0.5, label="2-min bins"),
        plt.Line2D([0], [0], color="#1f5f99", linewidth=2, label="2-min daily median"),
        plt.Line2D([0], [0], color="#222222", marker="o", markersize=4, linewidth=1.2, label="hourly median"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=10)
    fig.suptitle(
        f"{label} by group dyad, 5 m merge-mixing analysis\n{note}",
        fontsize=15,
        fontweight="bold",
        y=0.996,
    )
    fig.supxlabel("Cumulative observed merge days for that dyad", fontsize=11, y=0.035)
    fig.supylabel(label, fontsize=11, x=0.006)
    fig.subplots_adjust(left=0.045, right=0.995, top=0.92, bottom=0.08, hspace=0.48, wspace=0.18)

    out = OUT_DIR / f"{OUTPUT_PREFIX}_{metric}.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def write_html(paths: list[Path]) -> Path:
    out = OUT_DIR / f"{OUTPUT_PREFIX}_dashboard.html"
    rel_imgs = "\n".join(
        f"<h2>{html.escape(path.stem.replace(OUTPUT_PREFIX + '_', '').replace('_', ' '))}</h2>"
        f"<img src='{html.escape(path.name)}' alt='{html.escape(path.stem)}'>"
        for path in paths
    )
    out.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>5 m values by dyad</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;color:#222;background:#fff}}
.wrap{{padding:24px 30px 44px;max-width:1700px}}
.note{{color:#555;max-width:1100px;line-height:1.45}}
img{{max-width:100%;display:block;margin:16px 0 36px;border:1px solid #ddd;border-radius:6px}}
</style></head><body><div class="wrap">
<h1>5 m merge-mixing values by group dyad</h1>
<p class="note">Each facet is one group dyad. Blue points are 2-minute bins inside hourly merge clusters, the blue line is a daily median over cumulative observed merge time, and black points/lines are hourly-median values recomputed at the same 5 m edge radius.</p>
{rel_imgs}
</div></body></html>""",
        encoding="utf-8",
    )
    return out


def main() -> None:
    hourly, two_min = load_data()
    hourly.to_csv(OUT_DIR / f"{OUTPUT_PREFIX}_hourly_5m_metric_rows.csv", index=False)
    paths = [plot_metric(hourly, two_min, *metric_spec) for metric_spec in METRICS]
    dashboard = write_html(paths)
    print(f"Wrote dashboard: {dashboard}")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
