"""Does the betweenness decline survive correction for collar count, and how does
the integration signal vary with spatial scale?

Normalised betweenness falls as a network grows, because more nodes means more
alternative routes. Collar count roughly doubles at the 2025-07-31 deployment, so
an uncorrected decline is confounded with sample size. `rar_*` columns hold
betweenness recomputed on repeated subsamples of exactly 5 collars per origin
group, which removes that confound. Comparing raw against rarefied says which
findings were real.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

PROJECT = Path(__file__).resolve().parent
DEPLOYMENT = pd.Timestamp("2025-08-01")
RAW, RAR = "#9ec5f4", "#0d366b"
COPPER, LILAC, BROKER = "#2a78d6", "#eb6834", "#4a3aa7"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.labelcolor": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
})


def trend(s: pd.DataFrame, col: str):
    d = s.dropna(subset=[col])
    if len(d) < 12:
        return None
    t = (d.week - d.week.min()).dt.days / 365.25
    fit = sm.OLS(d[col].astype(float).to_numpy(),
                 sm.add_constant(pd.DataFrame({"t": t.to_numpy()}))).fit(
        cov_type="HAC", cov_kwds={"maxlags": 4})
    ci = fit.conf_int()
    return fit.params["t"], ci.loc["t", 0], ci.loc["t", 1], fit.pvalues["t"]


def dress(ax, t, s, yl):
    ax.set_title(t, loc="left", fontweight="bold", color=INK, pad=26, fontsize=10.5)
    ax.text(0, 1.035, s, transform=ax.transAxes, fontsize=8.2, color=MUTED)
    ax.set_ylabel(yl, fontsize=8.5)
    ax.grid(axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, length=3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=PROJECT / "outputs" / "copper_lilac_weekly_network_metrics_2026-09-03")
    a = ap.parse_args()
    d = pd.read_csv(a.dir / "weekly_network_metrics.csv", parse_dates=["week"])
    radii = sorted(d.radius_m.unique())

    fig, axes = plt.subplots(1, 3, figsize=(14.6, 5.0))
    fig.subplots_adjust(top=0.72, bottom=0.15, left=0.06, right=0.985, wspace=0.30)

    # A: raw vs rarefied brokerage at 5 m
    ax = axes[0]
    s = d[d.radius_m == 5.0].sort_values("week")
    for col, c, lab in [("between_cross_mean", RAW, "raw (all collars)"),
                        ("rar_cross", RAR, "rarefied to 5+5 collars")]:
        ax.plot(s.week, s[col].rolling(5, center=True, min_periods=2).mean(),
                color=c, lw=2.3, zorder=3, label=lab)
    ax.axvline(DEPLOYMENT, color=LILAC, lw=1.3, ls=(0, (4, 3)), zorder=2)
    ax.legend(loc="upper right", fontsize=8.2, labelcolor=INK2)
    ax.xaxis.set_major_locator(mpl.dates.MonthLocator(interval=5))
    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter("%Y-%m"))
    dress(ax, "A  Brokerage collapse is real", "Copper<->Lilac betweenness at 5 m, 5-week mean",
          "normalised betweenness")

    # B: raw vs rarefied trend, per metric, at 5 m
    ax = axes[1]
    labels = [("between_within_copper_mean", "rar_within_copper", "within Copper"),
              ("between_within_lilac_mean", "rar_within_lilac", "within Lilac"),
              ("between_cross_mean", "rar_cross", "Cop<->Lil brokerage")]
    y = np.arange(len(labels))
    ax.axvline(0, color=AXIS, lw=1.2, zorder=1)
    for i, (rawc, rarc, lab) in enumerate(labels):
        for col, c, off in [(rawc, RAW, 0.17), (rarc, RAR, -0.17)]:
            r = trend(s, col)
            if r:
                est, lo, hi, p = r
                ax.errorbar(est, y[i] + off, xerr=[[est - lo], [hi - est]], fmt="o",
                            color=c, ms=6.5, lw=2, capsize=3, zorder=3)
                ax.text(hi + 0.004, y[i] + off, f"p={p:.3f}", fontsize=7.2,
                        color=INK if p < 0.05 else MUTED, va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, _, l in labels], fontsize=8.5)
    ax.set_ylim(-0.6, len(labels) - 0.25)
    ax.invert_yaxis()
    ax.set_xlabel("trend per year (5 m)", fontsize=8.5)
    ax.grid(axis="y", lw=0)
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=7, color=RAW, label="raw"),
                       plt.Line2D([], [], marker="o", ls="", ms=7, color=RAR,
                                  label="rarefied to 5+5")],
              loc="lower left", fontsize=8.2, labelcolor=INK2, ncol=2)
    dress(ax, "B  Only brokerage survives the correction",
          "trend with 95% CI, HAC(4); raw vs collar-controlled", "")

    # C: scale gradient of the integration signal
    ax = axes[2]
    xs = np.arange(len(radii))
    for col, c, mk, lab in [("Q_origin", "#0d366b", "o", "Q_origin"),
                            ("ARI_vs_origin", "#008300", "s", "ARI vs origin"),
                            ("rar_cross", BROKER, "^", "brokerage (rarefied)")]:
        est, lo, hi = [], [], []
        for r in radii:
            t = trend(d[d.radius_m == r].sort_values("week"), col)
            est.append(t[0] if t else np.nan)
            lo.append(t[1] if t else np.nan)
            hi.append(t[2] if t else np.nan)
        est, lo, hi = np.array(est), np.array(lo), np.array(hi)
        ax.errorbar(xs, est, yerr=[est - lo, hi - est], fmt=mk + "-", color=c, ms=6,
                    lw=1.8, capsize=3, zorder=3, label=lab)
    ax.axhline(0, color=AXIS, lw=1.2, zorder=1)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{r:g} m" for r in radii], fontsize=8.5)
    ax.set_xlabel("proximity radius", fontsize=8.5)
    ax.legend(loc="lower right", fontsize=8.2, labelcolor=INK2)
    dress(ax, "C  The signal fades with scale",
          "trend per year by radius, 95% CI", "trend per year")

    fig.suptitle("Correcting betweenness for collar count, and how integration varies with scale",
                 x=0.06, y=0.965, ha="left", va="top", fontsize=13.5,
                 fontweight="bold", color=INK)
    fig.text(0.06, 0.885,
             "Normalised betweenness falls as a network grows, and collar count roughly doubles at the deployment - so an uncorrected decline is confounded with sample size.\n"
             "`Rarefied` recomputes it on repeated random subsamples of exactly 5 collars per group (60 draws/week), holding node count fixed across the whole series.",
             ha="left", va="top", fontsize=8.5, color=INK2, linespacing=1.55)
    out = a.dir / "betweenness_rarefaction_by_scale.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print("wrote", out.name)

    print(f"\n{'metric':<26}" + "".join(f"{r:>12g} m" for r in radii))
    for col in ["Q_origin", "ARI_vs_origin", "rar_cross", "rar_within_copper",
                "rar_within_lilac"]:
        line = f"{col:<26}"
        for r in radii:
            t = trend(d[d.radius_m == r].sort_values("week"), col)
            line += f"{t[0]:>+11.4f}{'*' if t and t[3] < 0.05 else ' '}" if t else f"{'--':>12}"
        print(line)


if __name__ == "__main__":
    main()
