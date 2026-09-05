from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


BINS = Path("outputs/disperser_finescale_integration/disperser_2min_contact_rows.parquet")
EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/canonical_disperser_events.csv")
OUT = Path("outputs/disperser_integration_calendar_time")
RADII = (2.0, 5.0)
MAX_DAY = 59
MAX_WEEK = 15
MIN_EVENTS = 3
BOOTSTRAPS = 1000
SEED = 20260810


def make_cells() -> pd.DataFrame:
    bins = pd.read_parquet(BINS)
    events = pd.read_csv(EVENTS, parse_dates=["start_time"])[["event_id", "start_time"]]
    bins["bin_2min"] = pd.to_datetime(bins.bin_2min)
    bins = bins.merge(events, on="event_id", how="left", validate="many_to_one")
    elapsed_days = (bins.bin_2min - bins.start_time).dt.total_seconds() / 86400
    bins["day"] = np.floor(elapsed_days).astype(int)
    bins["week"] = np.floor(elapsed_days / 7).astype(int)
    cells = []
    for scale, maximum in [("day", MAX_DAY), ("week", MAX_WEEK)]:
        x = bins[bins[scale].between(0, maximum) & bins.radius_m.isin(RADII)]
        c = (x.groupby(["event_id", "animal_id", "recipient_group", "radius_m", scale], observed=True)
             .agg(contacts=("recipient_contact", "sum"), opportunities=("recipient_opportunity", "sum"))
             .reset_index().rename(columns={scale: "time_index"}))
        c["integration_score"] = c.contacts / c.opportunities.replace(0, np.nan)
        c["scale"] = scale
        cells.append(c.dropna(subset=["integration_score"]))
    return pd.concat(cells, ignore_index=True)


def summarize_one(cells: pd.DataFrame, group: str, rng: np.random.Generator) -> pd.DataFrame:
    c = cells if group == "All groups" else cells[cells.recipient_group.eq(group)]
    observed = (c.groupby(["scale", "radius_m", "time_index"], observed=True)
                .agg(mean=("integration_score", "mean"), events=("event_id", "nunique"), animals=("animal_id", "nunique"))
                .reset_index())
    ids = c.event_id.drop_duplicates().to_numpy()
    draws = []
    for rep in range(BOOTSTRAPS):
        selected = rng.choice(ids, len(ids), replace=True)
        # Preserve event-level weighting; repeated sampled events remain separate replicates.
        parts = []
        for sample_i, eid in enumerate(selected):
            z = c[c.event_id.eq(eid)].copy()
            z["sample_event"] = sample_i
            parts.append(z)
        sample = pd.concat(parts, ignore_index=True)
        draw = (sample.groupby(["scale", "radius_m", "time_index"], observed=True).integration_score.mean()
                .rename("value").reset_index())
        draw["rep"] = rep
        draws.append(draw)
    boot = pd.concat(draws, ignore_index=True)
    ci = (boot.groupby(["scale", "radius_m", "time_index"], observed=True).value
          .quantile([.025, .975]).unstack().reset_index().rename(columns={.025: "low", .975: "high"}))
    out = observed.merge(ci, on=["scale", "radius_m", "time_index"], how="left")
    out = out[out.events.ge(MIN_EVENTS)].copy()
    out["group"] = group
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cells = make_cells()
    cells.to_parquet(OUT / "disperser_event_calendar_cells.parquet", index=False)
    rng = np.random.default_rng(SEED)
    groups = ["All groups"] + sorted(cells.recipient_group.unique().tolist())
    summaries = []
    for group in groups:
        n_events = cells.event_id.nunique() if group == "All groups" else cells.loc[cells.recipient_group.eq(group), "event_id"].nunique()
        if n_events >= MIN_EVENTS:
            summaries.append(summarize_one(cells, group, rng))
    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(OUT / "disperser_integration_by_day_week.csv", index=False)
    coverage = (cells.groupby(["scale", "time_index"], observed=True)
                .agg(events=("event_id", "nunique"), animals=("animal_id", "nunique"), groups=("recipient_group", "nunique"))
                .reset_index())
    coverage.to_csv(OUT / "calendar_time_coverage.csv", index=False)
    print(f"events={cells.event_id.nunique()} animals={cells.animal_id.nunique()} groups={cells.recipient_group.nunique()}")
    print(summary[summary.group.eq('All groups')].to_string(index=False))


if __name__ == "__main__":
    main()
