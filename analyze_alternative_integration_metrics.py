from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATA = Path("outputs/disperser_vs_recipient_members/disperser_recipient_member_bin_comparison.parquet")
OUT = Path("outputs/disperser_alternative_integration_metrics")
N_BOOT = 2000
SEED = 20260810


METRICS = {
    "contact_5m": ("Time within 5 m", "Proportion", True),
    "nearest_recipient_m": ("Nearest recipient", "Distance (m)", False),
    "recipient_preference": ("Recipient − origin preference", "Probability difference", True),
    "contact_day_consistency": ("Days with any 5 m contact", "Proportion of observed days", True),
    "resident_parity": ("Parity with resident members", "Focal − resident probability", True),
}


def event_week_cells(data: pd.DataFrame) -> pd.DataFrame:
    d = data[data.radius_m.eq(5)].copy()
    d["bin_2min"] = pd.to_datetime(d.bin_2min)
    d["start_time"] = pd.to_datetime(d.start_time)
    d["week"] = np.floor((d.bin_2min - d.start_time).dt.total_seconds() / (7 * 86400)).astype(int)
    d["day"] = np.floor((d.bin_2min - d.start_time).dt.total_seconds() / 86400).astype(int)
    d = d[d.week.between(0, 7)]
    daily = (d.groupby(["event_id", "week", "day"], observed=True).recipient_contact.max()
             .rename("day_contact").reset_index())
    consistency = (daily.groupby(["event_id", "week"], observed=True).day_contact.mean()
                   .rename("contact_day_consistency").reset_index())
    cells = (d.groupby(["event_id", "animal_id", "recipient_group", "week"], observed=True)
             .agg(recipient_contacts=("recipient_contact", "sum"), recipient_opportunities=("recipient_opportunity", "sum"),
                  origin_contacts=("origin_contact", "sum"), origin_opportunities=("origin_opportunity", "sum"),
                  nearest_recipient_m=("nearest_recipient_m", "median"), resident_contact_probability=("resident_contact_fraction", "mean"))
             .reset_index())
    cells["contact_5m"] = cells.recipient_contacts / cells.recipient_opportunities.replace(0, np.nan)
    cells["origin_5m"] = cells.origin_contacts / cells.origin_opportunities.replace(0, np.nan)
    cells["recipient_preference"] = cells.contact_5m - cells.origin_5m
    cells["resident_parity"] = cells.contact_5m - cells.resident_contact_probability
    return cells.merge(consistency, on=["event_id", "week"], how="left")


def summarize(cells: pd.DataFrame) -> pd.DataFrame:
    long = cells.melt(id_vars=["event_id", "week"], value_vars=list(METRICS), var_name="metric", value_name="value").dropna()
    observed = (long.groupby(["metric", "week"], observed=True)
                .agg(mean=("value", "mean"), events=("event_id", "nunique")).reset_index())
    ids = cells.event_id.unique(); rng = np.random.default_rng(SEED); draws = []
    for rep in range(N_BOOT):
        chosen = rng.choice(ids, len(ids), replace=True)
        parts = []
        for i, eid in enumerate(chosen):
            z = long[long.event_id.eq(eid)].copy(); z["sample_id"] = i; parts.append(z)
        x = pd.concat(parts, ignore_index=True).groupby(["metric", "week"], observed=True).value.mean().rename("value").reset_index()
        x["rep"] = rep; draws.append(x)
    boot = pd.concat(draws, ignore_index=True)
    ci = (boot.groupby(["metric", "week"], observed=True).value.quantile([.025, .975]).unstack()
          .reset_index().rename(columns={.025: "low", .975: "high"}))
    return observed.merge(ci, on=["metric", "week"]).query("events >= 5")


def make_plot(summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(13, 13)); axes = axes.ravel()
    for ax, (metric, (title, ylabel, higher)) in zip(axes, METRICS.items()):
        z = summary[summary.metric.eq(metric)].sort_values("week")
        ax.plot(z.week, z["mean"], marker="o", color="#3969ac", lw=2)
        ax.fill_between(z.week, z.low, z.high, color="#3969ac", alpha=.18)
        if metric in ["recipient_preference", "resident_parity"]: ax.axhline(0, color="black", ls="--", lw=1)
        ax.set(title=title, xlabel="Weeks since entry", ylabel=ylabel)
        ax.grid(alpha=.2)
        ax.text(.99, .03, "Higher = more integrated" if higher else "Lower = more integrated", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=9)
    axes[-1].axis("off")
    fig.suptitle("Five complementary measures of disperser integration\nEvent-equal means; 95% CIs resample complete dispersal events")
    fig.tight_layout(rect=[0, 0, 1, .96]); fig.savefig(OUT / "alternative_integration_metrics.png", dpi=220); plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cells = event_week_cells(pd.read_parquet(DATA)); cells.to_csv(OUT / "event_week_metric_cells.csv", index=False)
    summary = summarize(cells); summary.to_csv(OUT / "alternative_integration_metrics_summary.csv", index=False)
    make_plot(summary)
    print(summary.to_string(index=False))


if __name__ == "__main__": main()
