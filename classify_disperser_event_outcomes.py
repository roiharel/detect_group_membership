from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


CANONICAL = Path(r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_shared_full_20260722\canonical_hourly_membership_with_association_events.parquet")
EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/canonical_disperser_events.csv")
CELLS = Path("outputs/disperser_integration_calendar_time/disperser_event_calendar_cells.parquet")
OUT = Path("outputs/disperser_integration_calendar_time")
POST_DAYS = 14


def plot_event_lines(cells: pd.DataFrame, output: Path) -> None:
    cells = cells[cells.radius_m.isin([2.0, 5.0])].copy()
    colors = {"Stayed in recipient group": "#16857b", "Left — returned to origin": "#d95f02",
              "Left — isolated": "#7b3294", "Left — moved to another group": "#3f62a8",
              "Censored / no later data": "#777777"}
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharey="row")
    for row, radius in enumerate([2.0, 5.0]):
        for col, (scale, maximum) in enumerate([("day", 13), ("week", 7)]):
            ax = axes[row, col]
            z = cells[(cells.radius_m.eq(radius)) & cells.scale.eq(scale) & cells.time_index.le(maximum)]
            for event_id, g in z.groupby("event_id", observed=True):
                g = g.sort_values("time_index"); fate = g.fate.iloc[0]
                ax.plot(g.time_index, g.integration_score, color=colors.get(fate, "#777777"),
                        alpha=.72, lw=1.25, marker="o", ms=2.5)
            ax.set_xlabel("Days since entry" if scale == "day" else "Weeks since entry")
            ax.set_ylabel(f"Proportion within {int(radius)} m" if col == 0 else "")
            ax.set_title(f"{int(radius)} m integration — {scale}ly")
            ax.set_ylim(-.015, max(.5, z.integration_score.max() * 1.05))
            ax.grid(alpha=.2)
    handles = [plt.Line2D([0], [0], color=color, lw=2, marker="o", ms=4, label=label)
               for label, color in colors.items() if label in set(cells.fate)]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(.5, .955))
    fig.suptitle("Disperser integration trajectories by event and subsequent dynamic-group outcome\n"
                 "Outcome uses the focal animal’s modal dynamic social unit during the 14 days after event end",
                 y=.995)
    fig.tight_layout(rect=[0, 0, 1, .91])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    events = pd.read_csv(EVENTS, parse_dates=["start_time", "end_time"])
    c = pd.read_parquet(CANONICAL, columns=["window_start", "animal_id", "dynamic_social_unit", "observed_cluster_id"])
    c["window_start"] = pd.to_datetime(c.window_start)
    cluster_sizes = (c.groupby(["window_start", "observed_cluster_id"], observed=True).animal_id
                     .transform("nunique"))
    c["observed_cluster_size"] = cluster_sizes
    rows = []
    for e in events.itertuples(index=False):
        post = c[(c.animal_id.astype(str).eq(str(e.animal_id))) &
                 (c.window_start.gt(e.end_time)) &
                 (c.window_start.le(e.end_time + pd.Timedelta(days=POST_DAYS)))].copy()
        if post.empty:
            fate, destination = "Censored / no later data", ""
            n_post, isolated_fraction = 0, float("nan")
        else:
            n_post = post.window_start.nunique()
            isolated_fraction = post.observed_cluster_size.eq(1).mean()
            modal_unit = str(post.dynamic_social_unit.mode().iloc[0])
            if isolated_fraction >= 0.5:
                fate, destination = "Left — isolated", modal_unit
            elif modal_unit == str(e.dynamic_social_unit):
                fate, destination = "Stayed in recipient group", modal_unit
            elif modal_unit == str(e.origin_group):
                fate, destination = "Left — returned to origin", modal_unit
            else:
                fate, destination = "Left — moved to another group", modal_unit
        rows.append({"event_id": e.event_id, "animal_id": e.animal_id, "origin_group": e.origin_group,
                     "recipient_group": e.dynamic_social_unit, "end_time": e.end_time, "fate": fate,
                     "destination_group": destination, "post_observed_hours": n_post,
                     "post_isolated_fraction": isolated_fraction})
    outcomes = pd.DataFrame(rows)
    outcomes.to_csv(OUT / "disperser_event_outcomes_14d.csv", index=False)
    cells = pd.read_parquet(CELLS).merge(outcomes, on=["event_id", "animal_id", "recipient_group"], how="left")
    cells.to_csv(OUT / "disperser_event_lines_with_outcome.csv", index=False)
    plot_event_lines(cells, OUT / "disperser_event_lines_by_outcome.png")
    print(outcomes[["event_id", "animal_id", "recipient_group", "fate", "destination_group", "post_observed_hours"]].to_string(index=False))
    print("\n", outcomes.fate.value_counts().to_string())


if __name__ == "__main__":
    main()
