from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CELLS = Path("outputs/disperser_network_centrality/event_week_centrality_cells.csv")
OUTCOMES = Path("outputs/disperser_integration_calendar_time/disperser_event_outcomes_14d.csv")
OUT = Path("outputs/disperser_network_centrality")
METRICS = ["degree", "strength", "eigenvector", "betweenness", "harmonic_closeness"]
SEED = 20260810


def main() -> None:
    cells = pd.read_csv(CELLS)
    outcomes = pd.read_csv(OUTCOMES)[["event_id", "fate"]]
    outcomes["outcome"] = np.select(
        [outcomes.fate.eq("Stayed in recipient group"), outcomes.fate.str.startswith("Left", na=False)],
        ["Stayed", "Left"], default="Censored")
    cells = cells.merge(outcomes, on="event_id", how="left")
    confirmed = cells[cells.outcome.isin(["Stayed", "Left"])].copy()
    observed = (confirmed.groupby(["metric", "week", "outcome"], observed=True)
                .agg(events=("event_id", "nunique"), mean_percentile=("focal_percentile", "mean"),
                     minimum=("focal_percentile", "min"), maximum=("focal_percentile", "max"))
                .reset_index())
    # Event-bootstrap CI for stayers only; with two leavers their observed range is more honest.
    stay = confirmed[confirmed.outcome.eq("Stayed")]; ids = stay.event_id.unique(); rng = np.random.default_rng(SEED); draws = []
    for rep in range(2000):
        chosen = rng.choice(ids, len(ids), replace=True)
        s = pd.concat([stay[stay.event_id.eq(e)] for e in chosen], ignore_index=True)
        q = s.groupby(["metric", "week"], observed=True).focal_percentile.mean().rename("value").reset_index(); q["rep"] = rep; draws.append(q)
    ci = (pd.concat(draws).groupby(["metric", "week"], observed=True).value.quantile([.025, .975]).unstack()
          .reset_index().rename(columns={.025: "low", .975: "high"}))
    observed = observed.merge(ci, on=["metric", "week"], how="left")
    observed.to_csv(OUT / "centrality_stayed_vs_left_summary.csv", index=False)
    confirmed.to_csv(OUT / "centrality_stayed_vs_left_event_cells.csv", index=False)
    labels = {"degree": "Degree", "strength": "Weighted strength", "eigenvector": "Eigenvector",
              "betweenness": "Betweenness", "harmonic_closeness": "Harmonic closeness"}
    fig, axes = plt.subplots(3, 2, figsize=(13, 13)); axes = axes.ravel()
    for ax, metric in zip(axes, METRICS):
        z = observed[observed.metric.eq(metric)]
        s = z[z.outcome.eq("Stayed")].sort_values("week"); l = z[z.outcome.eq("Left")].sort_values("week")
        ax.plot(s.week, s.mean_percentile, marker="o", color="#16857b", label="Stayed")
        ax.fill_between(s.week, s.low, s.high, color="#16857b", alpha=.16)
        ax.plot(l.week, l.mean_percentile, marker="s", color="#d95f02", label="Left")
        ax.fill_between(l.week, l.minimum, l.maximum, color="#d95f02", alpha=.14)
        ax.axhline(.5, color="black", ls="--", lw=1); ax.set(title=labels[metric], xlabel="Weeks since entry", ylabel="Focal percentile", ylim=(-.03, 1.03)); ax.grid(alpha=.2)
    axes[0].legend(frameon=False); axes[-1].axis("off")
    fig.suptitle("Network integration of dispersers that stayed versus left\nStayers: event-bootstrap 95% CI; leavers: observed range of two events")
    fig.tight_layout(rect=[0, 0, 1, .96]); fig.savefig(OUT / "centrality_stayed_vs_left.png", dpi=220); plt.close(fig)
    print(observed.to_string(index=False))


if __name__ == "__main__": main()
