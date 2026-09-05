from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/canonical_disperser_events.csv")
MEMBERSHIP = Path(r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_local_2h_support\canonical_hourly_membership.parquet")
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
OUT = Path("outputs/disperser_finescale_integration")
RADII = (2.0, 5.0, 10.0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Measure fine-scale origin and recipient integration of dispersing animals.")
    p.add_argument("--events-file", type=Path, default=EVENTS)
    p.add_argument("--membership-file", type=Path, default=MEMBERSHIP)
    p.add_argument("--gps-file", type=Path, default=GPS)
    p.add_argument("--output-dir", type=Path, default=OUT)
    p.add_argument("--bootstrap-replicates", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20260808)
    return p.parse_args()


def distance_m(lon: float, lat: float, other_lon: np.ndarray, other_lat: np.ndarray) -> np.ndarray:
    radius = 6_371_000.0
    lon1, lat1 = np.radians(lon), np.radians(lat)
    lon2, lat2 = np.radians(other_lon), np.radians(other_lat)
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return radius * 2 * np.arctan2(np.sqrt(a), np.sqrt(np.maximum(0, 1 - a)))


def bootstrap(events: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    ids = events["event_id"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(seed)
    draws = []
    for replicate in range(replicates):
        chosen = rng.choice(ids, len(ids), replace=True)
        sample = pd.concat([events[events["event_id"].eq(x)] for x in chosen], ignore_index=True)
        for radius, frame in sample.groupby("radius_m", observed=True):
            origin = frame["origin_contacts"].sum() / frame["origin_opportunity_bins"].sum()
            recipient = frame["recipient_contacts"].sum() / frame["recipient_opportunity_bins"].sum()
            draws.append({"replicate": replicate, "radius_m": radius, "origin_contact_probability": origin,
                          "recipient_contact_probability": recipient, "recipient_minus_origin": recipient - origin})
    return pd.DataFrame(draws)


def plot(summary: pd.DataFrame, animal: pd.DataFrame, draws: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    x = np.arange(len(RADII)); width = 0.36
    for offset, metric, label, color in [
        (-width / 2, "origin_contact_probability", "Origin group", "#3b7a57"),
        (width / 2, "recipient_contact_probability", "Recipient group", "#7b3294"),
    ]:
        values = summary.set_index("radius_m").reindex(RADII)[metric]
        low = draws.groupby("radius_m")[metric].quantile(0.025).reindex(RADII)
        high = draws.groupby("radius_m")[metric].quantile(0.975).reindex(RADII)
        errors = np.vstack([
            np.maximum(0, values.to_numpy() - low.to_numpy()),
            np.maximum(0, high.to_numpy() - values.to_numpy()),
        ])
        axes[0].bar(x + offset, values, width, label=label, color=color,
                    yerr=errors, capsize=4)
    axes[0].set_xticks(x); axes[0].set_xticklabels([f"{int(r)} m" for r in RADII])
    axes[0].set_ylabel("Contact probability given simultaneous coverage")
    axes[0].set_title("All dispersal events")
    axes[0].legend(frameon=False)

    five = animal[animal["radius_m"].eq(5)].sort_values("recipient_minus_origin")
    labels = [f"{r.animal_id}: {r.origin_group} → {r.recipient_group}" for r in five.itertuples()]
    colors = np.where(five["recipient_minus_origin"].ge(0), "#7b3294", "#3b7a57")
    axes[1].barh(labels, five["recipient_minus_origin"], color=colors)
    axes[1].axvline(0, color="#444444", lw=1)
    axes[1].set_xlabel("Recipient minus origin contact probability")
    axes[1].set_title("Individual affiliation at 5 m")
    fig.suptitle("Fine-scale integration of dispersing animals\nEvent-bootstrap 95% CIs on pooled estimates")
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    events = pd.read_csv(args.events_file, parse_dates=["start_time", "end_time"])
    membership = pd.read_parquet(args.membership_file, columns=["animal_id", "window_start", "dynamic_social_unit"])
    membership["window_start"] = pd.to_datetime(membership["window_start"])
    event_hours = []
    relevant_members = []
    for event in events.itertuples(index=False):
        hours = pd.date_range(event.start_time, event.end_time, freq="h")
        event_hours.append(pd.DataFrame({"event_id": event.event_id, "focal_animal": event.animal_id,
                                         "origin_group": event.origin_group, "recipient_group": event.dynamic_social_unit, "hour": hours}))
        subset = membership[
            membership["window_start"].isin(hours)
            & (membership["dynamic_social_unit"].isin([event.origin_group, event.dynamic_social_unit])
               | membership["animal_id"].eq(event.animal_id))
        ].copy()
        relevant_members.append(subset)
    event_hours = pd.concat(event_hours, ignore_index=True)
    relevant_members = pd.concat(relevant_members, ignore_index=True).drop_duplicates(["animal_id", "window_start"])
    animals = sorted(relevant_members["animal_id"].astype(str).unique())
    start = event_hours["hour"].min().tz_localize("UTC")
    end = (event_hours["hour"].max() + pd.Timedelta(hours=1)).tz_localize("UTC")
    gps = pd.read_parquet(args.gps_file, columns=["animal_id", "timestamp", "location.long", "location.lat"],
                          filters=[("animal_id", "in", animals), ("timestamp", ">=", start), ("timestamp", "<", end)])
    gps["timestamp"] = pd.to_datetime(gps["timestamp"], utc=True).dt.tz_localize(None)
    gps["hour"] = gps["timestamp"].dt.floor("h")
    needed_hours = set(event_hours["hour"])
    gps = gps[gps["hour"].isin(needed_hours)].rename(columns={"location.long": "longitude", "location.lat": "latitude"})
    gps["bin_2min"] = gps["timestamp"].dt.floor("2min")
    fixes = gps.groupby(["hour", "bin_2min", "animal_id"], observed=True, as_index=False).agg(
        longitude=("longitude", "median"), latitude=("latitude", "median"))
    fixes = fixes.merge(relevant_members.rename(columns={"window_start": "hour"}), on=["animal_id", "hour"], how="inner")

    bin_rows = []
    for event in events.itertuples(index=False):
        hours = event_hours.loc[event_hours["event_id"].eq(event.event_id), "hour"]
        frame = fixes[fixes["hour"].isin(hours)]
        for (hour, bin_time), bin_frame in frame.groupby(["hour", "bin_2min"], observed=True):
            focal = bin_frame[bin_frame["animal_id"].eq(event.animal_id)]
            if focal.empty: continue
            focal = focal.iloc[0]
            others = bin_frame[~bin_frame["animal_id"].eq(event.animal_id)].copy()
            if others.empty: continue
            others["distance_m"] = distance_m(focal.longitude, focal.latitude, others.longitude.to_numpy(), others.latitude.to_numpy())
            for radius in RADII:
                origin = others[others["dynamic_social_unit"].eq(event.origin_group)]
                recipient = others[others["dynamic_social_unit"].eq(event.dynamic_social_unit)]
                bin_rows.append({
                    "event_id": event.event_id, "animal_id": event.animal_id, "origin_group": event.origin_group,
                    "recipient_group": event.dynamic_social_unit, "hour": hour, "bin_2min": bin_time, "radius_m": radius,
                    "origin_opportunity": not origin.empty, "recipient_opportunity": not recipient.empty,
                    "origin_contact": bool((origin["distance_m"] <= radius).any()),
                    "recipient_contact": bool((recipient["distance_m"] <= radius).any()),
                    "nearest_origin_m": origin["distance_m"].min() if not origin.empty else np.nan,
                    "nearest_recipient_m": recipient["distance_m"].min() if not recipient.empty else np.nan,
                    "n_origin_observed": len(origin), "n_recipient_observed": len(recipient),
                })
    bins = pd.DataFrame(bin_rows)
    event_summary = bins.groupby(["event_id", "animal_id", "origin_group", "recipient_group", "radius_m"], observed=True).agg(
        bins=("bin_2min", "size"), origin_opportunity_bins=("origin_opportunity", "sum"),
        recipient_opportunity_bins=("recipient_opportunity", "sum"), origin_contacts=("origin_contact", "sum"),
        recipient_contacts=("recipient_contact", "sum"), median_nearest_origin_m=("nearest_origin_m", "median"),
        median_nearest_recipient_m=("nearest_recipient_m", "median"),
    ).reset_index()
    event_summary["origin_contact_probability"] = event_summary["origin_contacts"] / event_summary["origin_opportunity_bins"].replace(0, np.nan)
    event_summary["recipient_contact_probability"] = event_summary["recipient_contacts"] / event_summary["recipient_opportunity_bins"].replace(0, np.nan)
    event_summary["recipient_minus_origin"] = event_summary["recipient_contact_probability"] - event_summary["origin_contact_probability"]
    overall = event_summary.groupby("radius_m", observed=True).agg(
        events=("event_id", "nunique"), animals=("animal_id", "nunique"), origin_contacts=("origin_contacts", "sum"),
        origin_opportunity_bins=("origin_opportunity_bins", "sum"), recipient_contacts=("recipient_contacts", "sum"),
        recipient_opportunity_bins=("recipient_opportunity_bins", "sum"),
    ).reset_index()
    overall["origin_contact_probability"] = overall["origin_contacts"] / overall["origin_opportunity_bins"]
    overall["recipient_contact_probability"] = overall["recipient_contacts"] / overall["recipient_opportunity_bins"]
    overall["recipient_minus_origin"] = overall["recipient_contact_probability"] - overall["origin_contact_probability"]
    animal = event_summary.groupby(["animal_id", "origin_group", "recipient_group", "radius_m"], observed=True).agg(
        events=("event_id", "nunique"), origin_contacts=("origin_contacts", "sum"), origin_opportunity_bins=("origin_opportunity_bins", "sum"),
        recipient_contacts=("recipient_contacts", "sum"), recipient_opportunity_bins=("recipient_opportunity_bins", "sum"),
    ).reset_index()
    animal["origin_contact_probability"] = animal["origin_contacts"] / animal["origin_opportunity_bins"].replace(0, np.nan)
    animal["recipient_contact_probability"] = animal["recipient_contacts"] / animal["recipient_opportunity_bins"].replace(0, np.nan)
    animal["recipient_minus_origin"] = animal["recipient_contact_probability"] - animal["origin_contact_probability"]
    draws = bootstrap(event_summary, args.bootstrap_replicates, args.seed)
    intervals = draws.groupby("radius_m").quantile([0.025, 0.5, 0.975]).reset_index()
    bins.to_parquet(args.output_dir / "disperser_2min_contact_rows.parquet", index=False)
    event_summary.to_csv(args.output_dir / "disperser_event_integration_summary.csv", index=False)
    animal.to_csv(args.output_dir / "disperser_animal_integration_summary.csv", index=False)
    overall.to_csv(args.output_dir / "disperser_overall_integration_summary.csv", index=False)
    intervals.to_csv(args.output_dir / "disperser_event_bootstrap_intervals.csv", index=False)
    plot(overall, animal, draws, args.output_dir / "disperser_finescale_integration.png")
    metadata = {"events_file": str(args.events_file.resolve()), "events_with_simultaneous_coverage": int(event_summary.event_id.nunique()),
                "animals": int(event_summary.animal_id.nunique()), "two_minute_bins": int(bins[["event_id", "bin_2min"]].drop_duplicates().shape[0]),
                "radii_m": RADII, "opportunity_definition": "focal and at least one member of the comparison social unit observed in the same 2-minute bin",
                "bootstrap_cluster_unit": "dispersal event", "caution": "contact is spatial proximity, not a directed social interaction"}
    (args.output_dir / "disperser_finescale_integration_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__": main()
