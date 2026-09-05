from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path(
    r"C:\Users\rharel\Documents\New project\outputs"
    r"\canonical_robust_hourly_membership_local_2h_support"
    r"\canonical_hourly_membership.parquet"
)
DEFAULT_OUTPUT = Path("outputs/chartreuse_social_reorganization")

STATE_ORDER = [
    "with_chartreuse_core",
    "within_group_splinter",
    "isolated",
    "dispersed_to_other_group",
    "uncertain",
]
STATE_COLORS = {
    "with_chartreuse_core": "#3b7a57",
    "within_group_splinter": "#e69f00",
    "isolated": "#999999",
    "dispersed_to_other_group": "#7b3294",
    "uncertain": "#d9d9d9",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an identity-resolved origin-centred panel for social reorganization."
    )
    parser.add_argument("--membership-file", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--origin-group", default="Chartreuse")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-gap-hours", type=float, default=2.5)
    parser.add_argument("--phase-window-days", type=int, default=14)
    return parser.parse_args()


def haversine_m(lon1: pd.Series, lat1: pd.Series, lon2: pd.Series, lat2: pd.Series) -> np.ndarray:
    radius = 6_371_000.0
    lon1r, lat1r, lon2r, lat2r = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2r - lon1r
    dlat = lat2r - lat1r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2) ** 2
    return radius * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def add_event_ids(active: pd.DataFrame, max_gap_hours: float) -> pd.DataFrame:
    active = active.sort_values(["animal_id", "window_start"]).copy()
    gap = active.groupby("animal_id", observed=True)["window_start"].diff().dt.total_seconds().div(3600)
    new_event = gap.isna() | gap.gt(max_gap_hours)
    active["event_number"] = new_event.groupby(active["animal_id"]).cumsum().astype(int)
    active["event_id"] = (
        "chartreuse_separation__"
        + active["animal_id"].astype(str)
        + "__"
        + active["event_number"].astype(str)
    )
    return active.drop(columns="event_number")


def build_panel(membership: pd.DataFrame, origin_group: str) -> pd.DataFrame:
    membership = membership.copy()
    membership["window_start"] = pd.to_datetime(membership["window_start"])
    focal = membership[membership["origin_group"].eq(origin_group)].copy()
    if focal.empty:
        raise ValueError(f"No rows found for origin_group={origin_group!r}")

    core_candidates = focal[
        focal["dynamic_social_unit"].eq(origin_group)
        & ~focal["dynamic_assignment"].eq("sustained_non_origin_association")
    ].copy()
    focal_cluster = (
        focal.groupby(["window_start", "temp_group_id"], observed=True)
        .agg(
            focal_cluster_size=("animal_id", "nunique"),
            cluster_lon=("longitude", "median"),
            cluster_lat=("latitude", "median"),
        )
        .reset_index()
    )
    core_counts = (
        core_candidates.groupby(["window_start", "temp_group_id"], observed=True)["animal_id"]
        .nunique().rename("core_candidate_size").reset_index()
    )
    core_rank = core_counts.sort_values(
        ["window_start", "core_candidate_size", "temp_group_id"],
        ascending=[True, False, True],
    )
    core_rank["focal_cluster_rank"] = (
        core_rank.groupby("window_start", observed=True).cumcount() + 1
    )
    core = core_rank[core_rank["focal_cluster_rank"].eq(1)].merge(
        focal_cluster,
        on=["window_start", "temp_group_id"],
        how="left",
    )[["window_start", "temp_group_id", "cluster_lon", "cluster_lat", "core_candidate_size"]].rename(
        columns={
            "temp_group_id": "core_temp_group_id",
            "cluster_lon": "core_longitude",
            "cluster_lat": "core_latitude",
            "core_candidate_size": "core_focal_size",
        }
    )
    focal = focal.merge(focal_cluster, on=["window_start", "temp_group_id"], how="left")
    focal = focal.merge(core, on="window_start", how="left")
    focal["is_core_component"] = focal["temp_group_id"].eq(focal["core_temp_group_id"])
    focal["distance_to_core_m"] = haversine_m(
        focal["longitude"], focal["latitude"], focal["core_longitude"], focal["core_latitude"]
    )

    dispersed = focal["dynamic_assignment"].eq("sustained_non_origin_association")
    isolated = focal["social_context"].eq("isolated") | focal["temp_group_size"].eq(1)
    supported = focal["temp_group_id"].notna()
    focal["social_state"] = np.select(
        [dispersed, isolated, supported & ~focal["is_core_component"], supported & focal["is_core_component"]],
        ["dispersed_to_other_group", "isolated", "within_group_splinter", "with_chartreuse_core"],
        default="uncertain",
    )
    focal["recipient_group"] = focal["dynamic_target_group"].where(dispersed)
    focal["is_separated"] = focal["social_state"].isin(
        ["within_group_splinter", "isolated", "dispersed_to_other_group"]
    )
    focal["outsider_count_in_cluster"] = (
        focal["temp_group_size"] - focal["focal_cluster_size"]
    ).clip(lower=0)
    focal["outsider_fraction_in_cluster"] = (
        focal["outsider_count_in_cluster"] / focal["temp_group_size"].replace(0, np.nan)
    )
    focal["interaction_context"] = np.select(
        [
            focal["is_mixed_temp_group"].fillna(False) & dispersed,
            focal["is_mixed_temp_group"].fillna(False),
            focal["social_state"].eq("isolated"),
        ],
        ["mixed_with_recipient", "mixed_with_other_group", "alone_or_low_support"],
        default="same_origin_only",
    )
    return focal


def build_events(panel: pd.DataFrame, max_gap_hours: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    active = add_event_ids(panel[panel["is_separated"]].copy(), max_gap_hours)
    events = (
        active.groupby("event_id", observed=True)
        .agg(
            animal_id=("animal_id", "first"),
            sex=("sex", "first"),
            age=("age", "first"),
            start_time=("window_start", "min"),
            end_time=("window_start", "max"),
            observed_hours=("window_start", "nunique"),
            dominant_state=("social_state", lambda x: x.value_counts().index[0]),
            n_splinter_hours=("social_state", lambda x: int((x == "within_group_splinter").sum())),
            n_isolated_hours=("social_state", lambda x: int((x == "isolated").sum())),
            n_dispersed_hours=("social_state", lambda x: int((x == "dispersed_to_other_group").sum())),
            recipient_group=("recipient_group", lambda x: x.dropna().mode().iloc[0] if not x.dropna().empty else None),
            max_splinter_size=("focal_cluster_size", "max"),
            mean_splinter_size=("focal_cluster_size", "mean"),
            max_distance_to_core_m=("distance_to_core_m", "max"),
            median_distance_to_core_m=("distance_to_core_m", "median"),
            mean_outsider_fraction=("outsider_fraction_in_cluster", "mean"),
            n_mixed_hours=("is_mixed_temp_group", "sum"),
            n_local_2h_supported=("is_local_2h_supported", "sum"),
            n_carried_night=("is_carried_night", "sum"),
        )
        .reset_index()
    )
    events["duration_hours"] = (
        (events["end_time"] - events["start_time"]).dt.total_seconds().div(3600) + 1
    )
    events["duration_days"] = events["duration_hours"] / 24

    lookup = panel[["animal_id", "window_start", "social_state", "dynamic_social_unit"]].copy()
    for label, offset in [("pre", pd.Timedelta(hours=-1)), ("post", pd.Timedelta(hours=1))]:
        key = events[["event_id", "animal_id", "start_time", "end_time"]].copy()
        key["lookup_time"] = key["start_time"] + offset if label == "pre" else key["end_time"] + offset
        key = key.merge(
            lookup,
            left_on=["animal_id", "lookup_time"],
            right_on=["animal_id", "window_start"],
            how="left",
        )
        events = events.merge(
            key[["event_id", "social_state", "dynamic_social_unit"]].rename(
                columns={"social_state": f"{label}_state", "dynamic_social_unit": f"{label}_social_unit"}
            ),
            on="event_id",
            how="left",
        )
    events["immediate_outcome"] = np.select(
        [
            events["post_state"].eq("with_chartreuse_core"),
            events["post_state"].eq("dispersed_to_other_group"),
            events["post_state"].isin(["within_group_splinter", "isolated"]),
            events["post_state"].isna(),
        ],
        ["returned_to_core", "remained_or_moved_to_recipient", "continued_separation", "unobserved_after_event"],
        default="other_transition",
    )
    return active, events


def build_phase_panel(panel: pd.DataFrame, events: pd.DataFrame, window_days: int) -> pd.DataFrame:
    pieces = []
    delta = pd.Timedelta(days=window_days)
    for event in events.itertuples(index=False):
        rows = panel[
            panel["animal_id"].eq(event.animal_id)
            & panel["window_start"].between(event.start_time - delta, event.end_time + delta)
        ].copy()
        if rows.empty:
            continue
        rows["event_id"] = event.event_id
        rows["event_start"] = event.start_time
        rows["event_end"] = event.end_time
        rows["event_phase"] = np.select(
            [rows["window_start"].lt(event.start_time), rows["window_start"].gt(event.end_time)],
            ["before", "after"],
            default="during",
        )
        rows["relative_hour"] = np.where(
            rows["event_phase"].eq("after"),
            (rows["window_start"] - event.end_time).dt.total_seconds() / 3600,
            (rows["window_start"] - event.start_time).dt.total_seconds() / 3600,
        )
        pieces.append(rows)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def plot_timeline(panel: pd.DataFrame, output: Path) -> None:
    daily = panel.assign(date=panel["window_start"].dt.floor("D"))
    daily = (
        daily.groupby(["animal_id", "date"], observed=True)["social_state"]
        .agg(lambda x: x.value_counts().index[0])
        .unstack("date")
    )
    activity = panel.groupby("animal_id", observed=True)["is_separated"].mean()
    daily = daily.loc[activity.sort_values(ascending=False).index]
    code = {state: i for i, state in enumerate(STATE_ORDER)}
    image = daily.apply(lambda column: column.map(code)).astype(float).to_numpy()
    cmap = plt.matplotlib.colors.ListedColormap([STATE_COLORS[s] for s in STATE_ORDER])
    fig, ax = plt.subplots(figsize=(16, max(7, 0.27 * len(daily))))
    ax.imshow(image, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.5, vmax=len(code) - 0.5)
    ax.set_yticks(range(len(daily)))
    ax.set_yticklabels(daily.index, fontsize=7)
    ticks = np.linspace(0, len(daily.columns) - 1, min(10, len(daily.columns))).astype(int)
    ax.set_xticks(ticks)
    ax.set_xticklabels([daily.columns[i].strftime("%Y-%m-%d") for i in ticks], rotation=35, ha="right")
    handles = [plt.Line2D([0], [0], color=STATE_COLORS[s], lw=8, label=s.replace("_", " ")) for s in STATE_ORDER]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=3, frameon=False)
    ax.set_title("Chartreuse-origin animals: dominant daily social state")
    ax.set_xlabel("Date")
    ax.set_ylabel("Animal ID (ordered by fraction of separated hours)")
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    columns = [
        "animal_id", "window_start", "tag_id", "sex", "age", "origin_group", "longitude", "latitude",
        "is_observed", "is_carried_night", "is_local_2h_supported", "position_support_type",
        "temp_group_id", "temp_group_size", "temp_group_origin_counts", "is_mixed_temp_group",
        "social_context", "dynamic_social_unit", "dynamic_assignment", "dynamic_target_group",
    ]
    membership = pd.read_parquet(args.membership_file, columns=columns)
    panel = build_panel(membership, args.origin_group)
    active, events = build_events(panel, args.max_gap_hours)
    phase = build_phase_panel(panel, events, args.phase_window_days)

    panel.to_parquet(args.output_dir / "chartreuse_hourly_identity_panel.parquet", index=False)
    events.to_csv(args.output_dir / "chartreuse_separation_events.csv", index=False)
    phase.to_parquet(args.output_dir / "chartreuse_event_phase_panel.parquet", index=False)
    state_summary = (
        panel.groupby(["animal_id", "sex", "age", "social_state"], dropna=False, observed=True)
        .size().rename("hours").reset_index()
    )
    state_summary.to_csv(args.output_dir / "chartreuse_animal_state_summary.csv", index=False)
    plot_timeline(panel, args.output_dir / "chartreuse_daily_social_state_timeline.png")

    metadata = {
        "membership_file": str(args.membership_file.resolve()),
        "origin_group": args.origin_group,
        "definition_core": "largest temporary cluster among animals whose dynamic social unit remains the origin group; ties use temp_group_id order",
        "definition_separated": list(STATE_ORDER[1:4]),
        "max_gap_hours_for_event_stitching": args.max_gap_hours,
        "phase_window_days": args.phase_window_days,
        "panel_rows": len(panel),
        "animals": int(panel["animal_id"].nunique()),
        "events": len(events),
        "state_counts": panel["social_state"].value_counts(dropna=False).to_dict(),
        "event_dominant_state_counts": events["dominant_state"].value_counts(dropna=False).to_dict(),
        "recipient_counts": events["recipient_group"].value_counts(dropna=False).to_dict(),
        "important_limitation": "temporary-cluster membership is structural association, not 5 m proximity integration",
    }
    (args.output_dir / "chartreuse_social_reorganization_metadata.json").write_text(
        json.dumps(metadata, indent=2, default=str), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, default=str))


if __name__ == "__main__":
    main()
