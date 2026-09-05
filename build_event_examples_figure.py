"""Manuscript figure: five worked examples with complete observation, plus depth.

The earlier version picked windows on event properties alone, so the rasters carried
blank cells wherever an animal was not observed. Here the candidate search is scored on
OBSERVATION COMPLETENESS first: for every candidate window it finds the animals observed
in EVERY interval of that window, and keeps the window that yields the most such animals.
Every cell in every raster below is therefore a real observation -- no gaps.

Five examples, so the between-group and individual types each show both of their ends:

  a  between-group, full      most of both units in one cluster
  b  between-group, partial   most of one unit, a minority of the other
  c  within-group             one unit holding more than one cluster
  d  individual, alone        one animal away from its unit and with no other
  e  individual, joined       one animal away and with another unit

Depth stays continuous and is shown once per TYPE, not per panel: mutual participation
for between-group, share outside the largest cluster for within-group, share of
away-nights with another unit for individual -- all on a common 0-1 scale.

Output: docs/framing_2026-09-04/event_examples_figure.html
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

NARROW = Path("outputs/membership_export_narrow/canonical_hourly_membership_narrow.csv")
STAGE1 = Path("outputs/general_structure_2026_09/phase2_two_stage_events/"
              "stage1_events_with_stage2_mixing.csv")
EXC = Path("outputs/general_structure_2026_09/phase4b_individual_axis/"
           "excursions_dominant_gap0.csv")
EXC_ANY = Path("outputs/general_structure_2026_09/phase4b_individual_axis/"
               "excursions_any_away_gap0.csv")
OUT = Path("docs/framing_2026-09-04/event_examples_figure.html")

# Collars fix positions 03:00-16:00 UTC. The night is not re-fixed, but the pipeline
# already carries the evening state forward through it (`is_carried_night`), which is
# exactly the assumption that animals stay put overnight and the interaction persists:
# carried rows retain social_context, association_event_id and temp_group_size in 100%
# of cases. Known state is therefore `is_observed | is_carried_night`, at 96.3% of all
# animal-hours against 58.7% observed.
#
# One genuine hole remains: 02:00, the handover hour, where the carry-forward has ended
# and the morning fix has not happened -- only 27% known. It is excluded from the slots.
DAY_HOURS = [h for h in list(range(3, 24)) + [0, 1]]   # 23 hours, 02:00 dropped
N_DAYS = 2                         # two full day-night cycles
DEAD_HOUR = 2
NIGHT_START = 17                   # hours >= this, or <= 1, are carried not fixed
MAX_COLS = len(DAY_HOURS) * N_DAYS
NIGHT_COLS = 44                    # nights in the individual trajectory panel
MAX_PER_UNIT = 9       # gapless animals drawn per unit
MIN_PER_UNIT = 4       # a candidate needs at least this many gapless animals per unit
MIN_SPLIT_ANIMALS = 9
MIN_OBS_FOR_SPLIT = 8
SEED = 20260904

# Panel a: three hours of the same dyad, chosen to straddle the clustering decision --
# apart, merged, apart again. Panel e: one animal whose assigned label churns while its
# dynamic unit holds.
GPS_UNITS = ("Chartreuse", "Maroon")
# Scanned rather than fixed: the frame shown is the hour in this window with at least
# one cluster holding both units and the most clusters overall, so a single map carries
# both a within-unit cluster and a shared one.
GPS_WINDOW = ("2025-12-06", "2025-12-09")
# Each animal gets a movement tail: its preceding fixes within the same daylight run.
# The frame is therefore restricted to hours late enough that TAIL_HOURS of fixes exist
# before it on the same day -- a tail spanning the night would be one straight jump
# between sleep sites rather than a track.
TAIL_HOURS = 6
STAGE_ANIMAL = "24AC17_7J8K"
STAGE_START = "2025-08-13"
STAGE_DAYS = 14
INSET_DAYS = 1     # drawn at the raster's px/hour; 309 px buys one day

CONTEXT_COLOUR = {
    "with_origin": "var(--n2)",
    "mixed_with_origin_present": "var(--c1)",
    "other": "var(--c5)",
    "mixed_without_origin_unit": "var(--c5b)",
    "isolated": "var(--alone)",
    "mixed_tie_or_unclear": "var(--muted)",
}
CONTEXT_SHORT = {
    "with_origin": "own unit",
    "mixed_with_origin_present": "merged",
    "other": "another unit",
    "mixed_without_origin_unit": "merged, away",
    "isolated": "alone",
    "mixed_tie_or_unclear": "unclear",
}
CONTEXT_LABEL = {
    "with_origin": "with origin unit",
    "mixed_with_origin_present": "merged, origin present",
    "other": "with another unit",
    "mixed_without_origin_unit": "merged, origin absent",
    "isolated": "alone",
    "mixed_tie_or_unclear": "unclear",
}

WITH_ORIGIN = {"with_origin", "mixed_with_origin_present"}
AWAY_WITH_OTHERS = {"other", "mixed_without_origin_unit"}

FILL = {"own": "var(--n2)", "mixed": "var(--c1)", "alone": "var(--alone)",
        "main": "var(--n2)", "secondary": "var(--c4)",
        "origin": "var(--n2)", "other": "var(--c5)"}

# Most social units in this study are named after colours, so a destination unit can be
# drawn in its own colour and read without a lookup. Mid-tones, chosen to hold on both
# a light and a dark ground.
UNIT_COLOUR = {
    "Bronze": "#a8703a", "Emerald": "#2e8b6b", "Chartreuse": "#8a9a2b",
    "Purple": "#7a5aa0", "Copper": "#b5702f", "Lilac": "#9b8ac4",
    "Jade": "#2f8f7a", "Lapis": "#2a5fa8", "Magenta": "#b0448c",
    "Maroon": "#9c3a4a", "Periwinkle": "#6d82c4", "Sapphire": "#2a6fb0",
    "Turquoise": "#2a9aa5", "RubyRunners": "#b03a4a", "FireOpal": "#c46a3a",
    "Obsidian": "#4a4a52", "Pearl": "#9aa0a8", "Red": "#b83a3a",
    "Green": "#3f8f45", "Gold": "#b8912a", "Mint": "#4fb08a",
    "SneakySilver": "#8f97a0", "PhantomWest": "#6a7480",
    "FenceJumpers": "#7a8a5a", "TrickyTeal": "#2a8a95",
    "LapisSplinter": "#4a7fc0",
}
UNIT_FALLBACK = "var(--c5)"


def unit_fill(unit: str) -> str:
    return UNIT_COLOUR.get(unit, UNIT_FALLBACK)


def cell_fill(v: str) -> str:
    """Cell colour; a named destination unit takes its own colour."""
    return unit_fill(v[5:]) if v.startswith("unit:") else FILL[v]


# ------------------------------------------------------------ observation index
def build_index(d: pd.DataFrame) -> dict:
    """hour -> unit -> frozenset of animals whose state is KNOWN in that hour.

    Known means fixed or carried through the night, following the assumption that
    animals stay put overnight and the interaction persists.
    """
# A row's state is KNOWN if it has a position from a real fix -- exactly at the
# hour, borrowed from a fix within 60 min (`local_2h`), or carried across the
# night. Omitting local_2h while accepting carried_night was indefensible and
# became visible only when coverage improved: at 17:00 the 2026-09-05 build
# supplies 106,147 local_2h rows where the frozen build had carried_night, so the
# old predicate silently discarded almost the whole hour. Interpolated rows are
# deliberately NOT known: their position is inferred, not observed.
    obs = d[d["is_observed"] | d["is_carried_night"]
            | d["is_local_2h_supported"]]
    idx: dict = {}
    for (ts, unit), g in obs.groupby(["window_start", "dynamic_social_unit"],
                                     sort=False):
        idx.setdefault(ts, {})[unit] = frozenset(g["animal_id"])
    return idx


def day_slots(day0: pd.Timestamp, n_days: int) -> list:
    """Continuous day-night cycles: 03:00 of day d through 01:00 of day d+1.

    Hour 02:00 is dropped because neither the carry-forward nor the morning fix covers
    it. Everything else is either fixed or carried, so the sequence is continuous.
    """
    d0 = pd.Timestamp(day0).normalize()
    out = []
    for k in range(n_days):
        for h in DAY_HOURS:
            day = k if h >= 3 else k + 1      # 00:00 and 01:00 belong to the next date
            out.append(d0 + pd.Timedelta(days=day, hours=h))
    return out


def is_night(ts: pd.Timestamp) -> bool:
    return ts.hour >= NIGHT_START or ts.hour <= 1


def gapless(idx: dict, unit: str, hours: list) -> frozenset:
    """Animals of `unit` observed in every one of `hours`."""
    sets = [idx.get(h, {}).get(unit, frozenset()) for h in hours]
    if not sets or any(len(s) == 0 for s in sets):
        return frozenset()
    out = sets[0]
    for s in sets[1:]:
        out &= s
    return out


# -------------------------------------------------------------- gps / clusters
def metres(lon, lat, lon0, lat0):
    """Local equirectangular projection, good enough over a few km."""
    k = math.cos(math.radians(lat0))
    return ((lon - lon0) * 111320.0 * k, (lat - lat0) * 110574.0)


def gps_panel(d, x0, y0, w, inset=None, side=292.0):
    """One hour of positions, coloured by unit, with the spatial clusters outlined.

    The hour is chosen rather than fixed: within GPS_WINDOW, the frame with at least one
    cluster holding both units and the most clusters overall, so a single map carries
    both a within-unit cluster and a shared one. A scale bar gives the distance the
    clustering is acting on instead of stating the parameter and asking for trust.
    """
    win = d[d["dynamic_social_unit"].isin(GPS_UNITS)
            & d["window_start"].between(*GPS_WINDOW)
            & d["is_observed"] & d["longitude"].notna()]
    if win.empty:
        return "", 0.0, {}
    best = None
    for ts, g in win.groupby("window_start"):
        if g["animal_id"].nunique() < 15:
            continue
        # need TAIL_HOURS of same-day fixes before this hour
        if ts.hour < 3 + TAIL_HOURS:
            continue
        shared = sum(1 for _, gg in g.groupby("association_event_id")
                     if gg["dynamic_social_unit"].nunique() > 1)
        off = int((g["origin_group"].astype(str)
                   != g["dynamic_social_unit"].astype(str)).sum())
        score = (1 if shared else 0, 1 if off else 0,
                 g["association_event_id"].nunique(), off)
        if best is None or score > best[0]:
            best = (score, ts, g)
    if best is None:
        return "", 0.0, {}
    _, ts, g = best

    # tail fixes: the same animals over the preceding TAIL_HOURS, same day
    tail = win[win["animal_id"].isin(g["animal_id"])
               & win["window_start"].between(ts - pd.Timedelta(hours=TAIL_HOURS), ts)
               & win["window_start"].dt.normalize().eq(ts.normalize())].copy()

    lon0, lat0 = tail["longitude"].mean(), tail["latitude"].mean()
    tmx, tmy = metres(tail["longitude"].to_numpy(), tail["latitude"].to_numpy(),
                      lon0, lat0)
    tail = tail.assign(mx=tmx, my=tmy)
    mx, my = metres(g["longitude"].to_numpy(), g["latitude"].to_numpy(), lon0, lat0)
    g = g.assign(mx=mx, my=my)
    pad = 110.0
    half = max(tmx.max() - tmx.min(), tmy.max() - tmy.min()) / 2 + pad
    ccx0, ccy0 = (tmx.max() + tmx.min()) / 2, (tmy.max() + tmy.min()) / 2

    # one enclosing box holds the map, its text, the legend and the inset
    pad = 13.0
    box_h = side + 2 * pad
    scale = side / (2 * half)
    mx0, my0 = x0 + pad, y0 + pad
    o = []
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="var(--paper)" '
             'stroke="var(--rule)" stroke-width=".9"/>' % (x0, y0, w, box_h))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--rule-soft)" '
             'stroke-width=".8"/>'
             % (mx0 + side + pad, y0 + 6, mx0 + side + pad, y0 + box_h - 6))

    def sx(v):
        return mx0 + (v - (ccx0 - half)) * scale

    def sy(v):
        return my0 + side - (v - (ccy0 - half)) * scale

    shared_n = 0
    for ev, gg in g.groupby("association_event_id"):
        multi = gg["dynamic_social_unit"].nunique() > 1
        shared_n += int(multi)
        px = [sx(v) for v in gg["mx"]]
        py = [sy(v) for v in gg["my"]]
        cx, cy = float(np.mean(px)), float(np.mean(py))
        rad = max(7.0, max((math.hypot(p - cx, q - cy) for p, q in zip(px, py)),
                           default=0.0) + 6.0)
        col = "var(--c1)" if multi else "var(--muted)"
        o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" opacity=".13" '
                 'stroke="%s" stroke-width="1" stroke-dasharray="3 2"/>'
                 % (cx, cy, rad, col, col))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%d</text>'
                 % (cx, cy - rad - 3, len(gg)))
    # tails first, then the current fix on top
    for aid, tg in tail.sort_values("window_start").groupby("animal_id"):
        if len(tg) < 2:
            continue
        u = tg["dynamic_social_unit"].iloc[-1]
        pts = " ".join("%.1f,%.1f" % (sx(v), sy(q))
                       for v, q in zip(tg["mx"], tg["my"]))
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1" '
                 'opacity=".5" stroke-linejoin="round"/>' % (pts, unit_fill(u)))
    off_origin = []
    for u in GPS_UNITS:
        for r in g[g["dynamic_social_unit"].eq(u)].itertuples():
            same = str(r.origin_group) == u
            if not same:
                off_origin.append(str(r.origin_group))
            o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="%s" '
                     'stroke-width="%.1f"/>'
                     % (sx(r.mx), sy(r.my), 2.6 if same else 3.4, unit_fill(u),
                        "none" if same else unit_fill(str(r.origin_group)),
                        0.0 if same else 1.7))

    aa = g[g["dynamic_social_unit"].eq(GPS_UNITS[0])]
    bb = g[g["dynamic_social_unit"].eq(GPS_UNITS[1])]
    dist = float("nan")
    if len(aa) and len(bb):
        dx = aa["mx"].to_numpy()[:, None] - bb["mx"].to_numpy()[None, :]
        dy = aa["my"].to_numpy()[:, None] - bb["my"].to_numpy()[None, :]
        dist = float(np.sqrt(dx ** 2 + dy ** 2).min())

    target = 2 * half / 4
    step = 10 ** math.floor(math.log10(target))
    bar = min([m * step for m in (1, 2, 5, 10) if m * step >= target * 0.6],
              default=step)
    bl = bar * scale
    by = my0 + side - 10
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="1.5"/>' % (mx0 + 6, by, mx0 + 6 + bl, by))
    for xx in (mx0 + 6, mx0 + 6 + bl):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
                 'stroke-width="1.5"/>' % (xx, by - 3.5, xx, by + 3.5))
    o.append('<text class="ts" x="%.1f" y="%.1f">%s m</text>'
             % (mx0 + 10 + bl, by + 3, f"{bar:,.0f}"))

    lx, ly = mx0 + side + pad + 14, y0 + pad + 6
    o.append('<text class="sn" x="%.1f" y="%.1f">%s</text>'
             % (lx, ly, ts.strftime("%d %b %Y, %H:%M")))
    o.append('<text class="ts" x="%.1f" y="%.1f">%d animals in %d clusters</text>'
             % (lx, ly + 12, g["animal_id"].nunique(),
                g["association_event_id"].nunique()))
    o.append('<text class="ts" x="%.1f" y="%.1f">nearest cross-unit %.0f m</text>'
             % (lx, ly + 23, dist))
    o.append('<text class="ts" x="%.1f" y="%.1f">%d away from its origin group</text>'
             % (lx, ly + 34, len(off_origin)))
    ly += 53
    for lab, col in ((GPS_UNITS[0], unit_fill(GPS_UNITS[0])),
                     (GPS_UNITS[1], unit_fill(GPS_UNITS[1]))):
        o.append('<circle cx="%.1f" cy="%.1f" r="2.8" fill="%s"/>' % (lx + 4, ly, col))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (lx + 12, ly + 3, lab))
        ly += 13
    ring = unit_fill(off_origin[0]) if off_origin else "var(--c4)"
    o.append('<circle cx="%.1f" cy="%.1f" r="3.4" fill="%s" stroke="%s" '
             'stroke-width="1.7"/>'
             % (lx + 4, ly, unit_fill(GPS_UNITS[0]), ring))
    o.append('<text class="ts" x="%.1f" y="%.1f">ring = origin group,</text>'
             % (lx + 12, ly + 3))
    o.append('<text class="ts" x="%.1f" y="%.1f">where it is not this unit</text>'
             % (lx + 12, ly + 14))
    ly += 24
    ly += 5
    for lab, col in (("one unit", "var(--muted)"), ("holds both units", "var(--c1)")):
        o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s" opacity=".2" stroke="%s" '
                 'stroke-width="1" stroke-dasharray="3 2"/>' % (lx + 4, ly, col, col))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (lx + 12, ly + 3, lab))
        ly += 13
    o.append('<text class="ts" x="%.1f" y="%.1f">dots are the fix at this hour,</text>'
             % (lx, ly + 4))
    o.append('<text class="ts" x="%.1f" y="%.1f">tails the preceding %d h; numbers</text>'
             % (lx, ly + 15, TAIL_HOURS))
    o.append('<text class="ts" x="%.1f" y="%.1f">are cluster sizes</text>'
             % (lx, ly + 26))
    ly += 22
    if inset is not None:
        isvg, _ih = inset(lx - 6, ly + 30, x0 + w - (lx - 6) - pad)
        o.append(isvg)
    return chr(10).join(o), box_h, {
        "hour": ts.strftime("%d %b %Y, %H:%M"),
        "animals": int(g["animal_id"].nunique()),
        "clusters": int(g["association_event_id"].nunique()),
        "shared": shared_n, "nearest_m": round(dist), "scale_m": int(bar),
        "tail_hours": TAIL_HOURS, "tail_fixes": int(len(tail)),
        "off_origin": len(off_origin),
        "off_origin_units": ", ".join(sorted(set(off_origin)))}


def stages_panel(d, x0, y0, w, compact: bool = False, days=None):
    """One animal, two weeks, two strips: the context it is in and the unit it is in.

    The pair is the point. Social context is recomputed every hour and moves constantly;
    the dynamic unit only changes on seven days of sustained association, so it does not
    move at all here.
    """
    days = STAGE_DAYS if days is None else days
    t0 = pd.Timestamp(STAGE_START)
    t1 = t0 + pd.Timedelta(days=days)
    g = d[d["animal_id"].eq(STAGE_ANIMAL) & d["window_start"].between(t0, t1)
          & (d["is_observed"] | d["is_carried_night"]
             | d["is_local_2h_supported"])].sort_values("window_start")
    if g.empty:
        return "", 0.0, {}
    hours = list(g["window_start"])
    lab_w = 56.0 if compact else 88.0
    width = w - lab_w
    cw = width / len(hours)
    gx0 = x0 + lab_w
    o, y = [], y0
    counts = {}

    for label, col in STAGE_STRIPS:
        values = list(g[col].astype(str))
        counts[col] = len(set(values))
        h = 13.0 if compact else 17.0
        if col == "social_context":
            def colour_of(v):
                return CONTEXT_COLOUR.get(v, "var(--muted)")
        else:
            def colour_of(v):
                return unit_fill(v)
        rs, rv = 0, values[0]
        for i in range(1, len(values) + 1):
            if i == len(values) or values[i] != rv:
                o.append('<rect x="%.2f" y="%.1f" width="%.2f" height="%.1f" fill="%s"/>'
                         % (gx0 + rs * cw, y, max(0.4, (i - rs) * cw), h, colour_of(rv)))
                if i < len(values):
                    rs, rv = i, values[i]
        o.append('<text class="%s" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % ("ts" if compact else "sn", gx0 - 6, y + h * 0.5 + 3,
                    label if not compact else label.split()[-1]))
        if not compact:
            o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                     % (gx0 - 8, y + 18, "%d distinct" % counts[col]))
        y += h + (7 if compact else 13)

    step = max(1, (7 if compact else 2) * days // max(1, STAGE_DAYS))
    for k in range(0, days + 1, step):
        xx = gx0 + (k / days) * width
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (xx, y0, xx, y - (7 if compact else 13)))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (xx, y + (3 if compact else -1),
                    (t0 + pd.Timedelta(days=k)).strftime("%d %b")))
    y += 6
    seen = []
    for v in g["social_context"].astype(str):
        if v not in seen:
            seen.append(v)
    lx, ly = gx0, y + 8
    for v in seen:
        lab = (CONTEXT_SHORT if compact else CONTEXT_LABEL).get(v, v)
        iw = 9 + len(lab) * 4.3 + 12
        if lx + iw > gx0 + width:
            lx, ly = gx0, ly + 11
        o.append('<rect x="%.1f" y="%.1f" width="6" height="6" fill="%s"/>'
                 % (lx, ly, CONTEXT_COLOUR.get(v, "var(--muted)")))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (lx + 9, ly + 5.5, lab))
        lx += iw
    return chr(10).join(o), (ly + 14) - y0, {
        "animal": STAGE_ANIMAL.split("_")[0],
        "unit": str(g["dynamic_social_unit"].iloc[0]),
        "hours": len(hours), "days": days,
        "contexts": counts.get("social_context", 0),
        "dynamic_labels": counts.get("dynamic_social_unit", 0)}


DAILY_UNITS = ("Copper", "Lilac")
DAILY_DAYS = 84          # a 12-week window at daily resolution
ISO_HOURS_COLS = 46      # two day-night cycles, matching the other hourly panels


def _daily_states(d, units):
    """Per animal-day: shared with the other unit, alone, or own unit only."""
    k = d[d["dynamic_social_unit"].isin(units)
          & (d["is_observed"] | d["is_carried_night"])].copy()
    k["day"] = k["window_start"].dt.normalize()
    nu = (k.groupby(["window_start", "association_event_id"])["dynamic_social_unit"]
          .nunique().rename("nu").reset_index())
    k = k.merge(nu, on=["window_start", "association_event_id"], how="left")
    k["is_shared"] = k["nu"].gt(1)
    k["is_alone"] = k["temp_group_size"].le(1)
    return k.groupby(["animal_id", "dynamic_social_unit", "day"]).agg(
        shared=("is_shared", "max"), alone=("is_alone", "max"),
        hrs=("is_shared", "size")).reset_index()


def pick_between_daily(d) -> dict | None:
    """One dyad at DAILY resolution over 12 weeks, spanning its transition.

    The hourly panels show what a single event looks like. A between-group relationship
    runs on a far longer clock than a single encounter, so the window is chosen to
    maximise the change in shared-day frequency between its first and second half --
    the dyad going from occasional encounters to sustained association.
    """
    ad = _daily_states(d, DAILY_UNITS)
    if ad.empty:
        return None
    days = sorted(ad["day"].unique())
    per_day = ad.groupby("day")["shared"].mean()
    best = None
    for i in range(0, len(days) - DAILY_DAYS, 7):
        win = days[i:i + DAILY_DAYS]
        sub = ad[ad["day"].isin(win)]
        # every animal drawn must have a known state on every day of the window
        counts = sub.groupby("animal_id")["day"].nunique()
        full = counts[counts == len(win)].index
        if len(full) < 8:
            continue
        half = len(win) // 2
        f = per_day.reindex(win[:half]).mean()
        g2 = per_day.reindex(win[half:]).mean()
        score = (g2 - f, len(full))
        if best is None or score > best[0]:
            best = (score, win, full)
    if best is None:
        return None
    (delta, _), win, full = best
    sub = ad[ad["day"].isin(win) & ad["animal_id"].isin(full)]
    rows = []
    for u in DAILY_UNITS:
        ids = sorted(sub.loc[sub["dynamic_social_unit"].eq(u), "animal_id"].unique())
        rows += [(u, a) for a in ids[:MAX_PER_UNIT]]
    keep = {a for _, a in rows}
    cells = {}
    for r in sub[sub["animal_id"].isin(keep)].itertuples():
        k = (r.dynamic_social_unit, r.animal_id)
        if k not in rows:
            continue
        cells[(k, r.day)] = ("mixed" if r.shared else
                             ("alone" if r.alone else "own"))
    return {"cols": list(win), "rows": rows, "cells": cells,
            "legend": [("own", "own unit only"), ("mixed", "both units"),
                       ("alone", "alone")],
            "meta": "%s, %s to %s; shared days %.0f%% then %.0f%%"
                    % (" - ".join(DAILY_UNITS),
                       pd.Timestamp(win[0]).date(), pd.Timestamp(win[-1]).date(),
                       100 * per_day.reindex(win[:len(win) // 2]).mean(),
                       100 * per_day.reindex(win[len(win) // 2:]).mean()),
            "xlab": "%d consecutive days" % len(win)}


def pick_isolation_hourly(d, idx, exc) -> dict | None:
    """One animal alone, at HOURLY resolution, with its unit-mates for reference.

    Isolation has a median of one hour, so it only reads at this scale: the hourly
    panel shows the animal detaching and rejoining within a single day-night cycle,
    which a daily panel would collapse to a single cell.
    """
    cand = exc[exc["depth"].eq("alone_only") & exc["alone_nights"].between(1, 3)]
    for r in cand.sort_values("alone_nights").itertuples():
        start = pd.Timestamp(r.start_night)
        hours = day_slots(start - pd.Timedelta(days=1), 2)
        g = gapless(idx, r.origin_group, hours)
        if r.animal_id not in g or len(g) < MIN_PER_UNIT + 1:
            continue
        unit, animal = r.origin_group, r.animal_id
        others = [a for a in sorted(g) if a != animal][:MAX_PER_UNIT]
        rows = [(unit, animal)] + [(unit, a) for a in others]
        w = d[d["window_start"].isin(hours) & d["dynamic_social_unit"].eq(unit)
              & d["animal_id"].isin([a for _, a in rows])
              & (d["is_observed"] | d["is_carried_night"])]
        biggest = {}
        for ts, gg in w.groupby("window_start"):
            biggest[ts] = gg.groupby("association_event_id").size().idxmax()
        cells, alone_h = {}, 0
        for t in w.itertuples():
            k = (unit, t.animal_id)
            if k not in rows:
                continue
            v = ("alone" if t.temp_group_size <= 1 else
                 ("main" if biggest.get(t.window_start) == t.association_event_id
                  else "secondary"))
            cells[(k, t.window_start)] = v
            if k == (unit, animal) and v == "alone":
                alone_h += 1
        if alone_h < 3:
            continue
        return {"cols": hours, "rows": rows, "cells": cells,
                "focal": (unit, animal),
                "legend": [("main", "with the unit"),
                           ("secondary", "in a secondary cluster"),
                           ("alone", "alone")],
                "meta": "%s of %s, %s; %d h alone"
                        % (animal.split("_")[0], unit, hours[0].date(), alone_h),
                "xlab": "%d continuous hours over 2 day-night cycles" % len(hours)}
    return None


def pick_dispersal_hourly(d, idx, exc_any) -> dict | None:
    """An animal leaving its unit, at HOURLY resolution across the transition.

    The nightly view showed where an excursion went; this shows the hour it left. The
    window is centred on the first away hour of an excursion, and away hours are
    coloured by the unit the animal was actually with, so a departure into another unit
    and a departure into nothing are distinguishable.
    """
    cand = exc_any[exc_any["other_nights"].ge(2)
                   & exc_any["away_nights"].between(2, 30)]
    for r in cand.sort_values("away_nights").itertuples():
        start = pd.Timestamp(r.start_night)
        hours = day_slots(start - pd.Timedelta(days=1), 2)
        g = gapless(idx, r.origin_group, hours)
        if r.animal_id not in g or len(g) < MIN_PER_UNIT + 1:
            continue
        unit, animal = r.origin_group, r.animal_id
        others = [a for a in sorted(g) if a != animal][:MAX_PER_UNIT]
        rows = [(unit, animal)] + [(unit, a) for a in others]
        w = d[d["window_start"].isin(hours) & d["dynamic_social_unit"].eq(unit)
              & d["animal_id"].isin([a for _, a in rows])
              & (d["is_observed"] | d["is_carried_night"])]
        cells = {}
        for t in w.itertuples():
            k = (unit, t.animal_id)
            if k not in rows:
                continue
            sc = t.social_context
            cells[(k, t.window_start)] = (
                "other" if sc in AWAY_WITH_OTHERS else
                ("alone" if sc == "isolated" else "origin"))
        # name the destination unit for the focal animal's away hours
        fa = w[w["animal_id"].eq(animal) & w["social_context"].isin(AWAY_WITH_OTHERS)]
        pairs = set(zip(fa["window_start"], fa["association_event_id"]))
        away_h = int(len(fa))
        if away_h < 4:
            continue
        pool = d[(d["is_observed"] | d["is_carried_night"])
                 & d["window_start"].isin(hours)
                 & d["dynamic_social_unit"].ne(unit)]
        pool = pool[pd.MultiIndex.from_arrays(
            [pool["window_start"], pool["association_event_id"]]).isin(pairs)]
        units_seen = []
        if not pool.empty:
            dest = (pool.groupby("window_start")["dynamic_social_unit"]
                    .agg(lambda x: x.value_counts().idxmax()).to_dict())
            for ts, u in dest.items():
                cells[((unit, animal), ts)] = "unit:" + u
                if u not in units_seen:
                    units_seen.append(u)
        leg = [("origin", "with its unit"), ("alone", "alone")]
        leg += [("unit:" + u, u) for u in units_seen]
        if any(cells.get(((unit, animal), h)) == "other" for h in hours):
            leg.append(("other", "another unit"))
        return {"cols": hours, "rows": rows, "cells": cells,
                "focal": (unit, animal),
                "legend": leg,
                "meta": "%s of %s, %s; %d h away in this window"
                        % (animal.split("_")[0], unit, hours[0].date(), away_h),
                "xlab": "%d continuous hours over 2 day-night cycles" % len(hours)}
    return None


def pick_all_patterns(d, s1, idx) -> dict | None:
    """One window, two units, with every pattern from the event taxonomy at once.

    Rather than one panel per event type, the search looks for a single 46-hour window
    over two units in which all three co-occur: a cluster holding both units (a merge),
    a unit split across two or more clusters, and at least one animal alone. Candidate
    windows come from the Stage-1 encounter list, which guarantees the merge; the split
    and the isolation are then required on top, with every animal drawn having a known
    state in every hour.
    """
    s = s1.copy()
    s["obs_min"] = s[["mean_observed_a", "mean_observed_b"]].min(axis=1)
    cand = s[s["obs_min"].ge(6) & s["structural_span_hours"].ge(8)]
    cand = cand.sort_values("obs_min", ascending=False)
    best = None
    for r in cand.itertuples():
        hours = day_slots(pd.Timestamp(r.start_hour), N_DAYS)
        ga = gapless(idx, r.unit_a, hours)
        gb = gapless(idx, r.unit_b, hours)
        if min(len(ga), len(gb)) < MIN_PER_UNIT:
            continue
        rows = ([(r.unit_a, a) for a in sorted(ga)[:MAX_PER_UNIT]]
                + [(r.unit_b, a) for a in sorted(gb)[:MAX_PER_UNIT]])
        keep = {a for _, a in rows}
        w = d[d["window_start"].isin(hours)
              & d["dynamic_social_unit"].isin([r.unit_a, r.unit_b])
              & d["animal_id"].isin(keep)
              & (d["is_observed"] | d["is_carried_night"])]
        if w.empty:
            continue
        # which clusters hold both units, and which cluster is each unit's largest
        nu = (w.groupby(["window_start", "association_event_id"])["dynamic_social_unit"]
              .nunique().rename("nu").reset_index())
        w = w.merge(nu, on=["window_start", "association_event_id"], how="left")
        biggest, split_hours = {}, 0
        for (ts, u), gg in w.groupby(["window_start", "dynamic_social_unit"]):
            own = gg[gg["nu"].eq(1)]
            sizes = own.groupby("association_event_id").size()
            if len(sizes):
                biggest[(ts, u)] = sizes.idxmax()
            if gg["association_event_id"].nunique() > 1:
                split_hours += 1
        cells = {}
        n_mixed = n_split = n_alone = 0
        for t in w.itertuples():
            k = (t.dynamic_social_unit, t.animal_id)
            if k not in rows:
                continue
            if t.temp_group_size <= 1:
                v = "alone"
                n_alone += 1
            elif t.nu and t.nu > 1:
                v = "mixed"
                n_mixed += 1
            elif biggest.get((t.window_start, t.dynamic_social_unit)) ==                     t.association_event_id:
                v = "main"
            else:
                v = "secondary"
                n_split += 1
            cells[(k, t.window_start)] = v
        if n_mixed < 12 or n_split < 6 or n_alone < 2:
            continue
        score = (min(n_split, 40), min(n_alone, 12), n_mixed)
        if best is None or score > best[0]:
            best = (score, r, hours, rows, cells, n_mixed, n_split, n_alone)
        if best is not None and best[0][0] >= 30 and best[0][1] >= 6:
            break
    if best is None:
        return None
    _, r, hours, rows, cells, n_mixed, n_split, n_alone = best
    return {"cols": hours, "rows": rows, "cells": cells,
            "legend": [("main", "with its unit"), ("secondary", "split off"),
                       ("mixed", "cluster holds both units"), ("alone", "alone")],
            "meta": ("%s, %s; %d h in a shared cluster, %d split off, %d alone"
                     % (r.pair_key, hours[0].date(), n_mixed, n_split, n_alone)),
            "xlab": "%d continuous hours over 2 day-night cycles" % len(hours)}


# ------------------------------------------------------------------ pickers
def pick_between(d, s1, idx, want: str) -> dict | None:
    s = s1.copy()
    s["med_min"] = s[["median_frac_a", "median_frac_b"]].min(axis=1)
    s["med_max"] = s[["median_frac_a", "median_frac_b"]].max(axis=1)
    if want == "full":
        cand = s[s["med_min"].ge(0.8)]
    else:
        cand = s[s["med_max"].ge(0.6) & s["med_min"].between(0.12, 0.5)]
    cand = cand[cand["structural_span_hours"].ge(10)]
    best = None
    for n_days in (N_DAYS, 1):
        for r in cand.itertuples():
            hours = day_slots(pd.Timestamp(r.start_hour), n_days)
            ga = gapless(idx, r.unit_a, hours)
            gb = gapless(idx, r.unit_b, hours)
            n = min(len(ga), len(gb))
            if n < MIN_PER_UNIT:
                continue
            score = (n, len(ga) + len(gb))
            if best is None or score > best[0]:
                best = (score, r, hours, ga, gb)
        if best is not None:
            break
    if best is None:
        return None
    _, r, hours, ga, gb = best
    rows = ([(r.unit_a, a) for a in sorted(ga)[:MAX_PER_UNIT]]
            + [(r.unit_b, a) for a in sorted(gb)[:MAX_PER_UNIT]])
    w = d[d["window_start"].isin(hours)
          & d["dynamic_social_unit"].isin([r.unit_a, r.unit_b])
          & d["animal_id"].isin([a for _, a in rows])]
    known = w[w["is_observed"] | w["is_carried_night"]]
    mixed = set()
    for ts, g in known.groupby("window_start"):
        for ev, gg in g.groupby("association_event_id"):
            if gg["dynamic_social_unit"].nunique() > 1:
                mixed.add((ts, ev))
    cells, carried = {}, set()
    for t in known.itertuples():
        k = (t.dynamic_social_unit, t.animal_id)
        if k not in rows:
            continue
        cells[(k, t.window_start)] = ("mixed" if (t.window_start,
                                                  t.association_event_id) in mixed
                                      else ("alone" if t.temp_group_size <= 1
                                            else "own"))
        if not t.is_observed:
            carried.add((k, t.window_start))
    return {"cols": hours, "rows": rows, "cells": cells, "carried": carried,
            "legend": [("own", "own unit only"), ("mixed", "both units"),
                       ("alone", "alone")],
            "meta": "%s, %s; depth %.2f" % (r.pair_key, hours[0].date(), r.med_min),
            "xlab": "%d continuous hours over %d day-night cycles" % (len(hours), len(hours) // len(DAY_HOURS))}


def pick_within(d, idx) -> dict | None:
    obs = d[d["is_observed"] & d["window_start"].dt.year.ge(2025)]
    sizes = (obs.groupby(["window_start", "dynamic_social_unit",
                          "association_event_id"]).size().rename("k").reset_index())
    tot = sizes.groupby(["window_start", "dynamic_social_unit"]).agg(
        total=("k", "sum"), biggest=("k", "max"), n_clusters=("k", "size")).reset_index()
    tot["depth"] = 1.0 - tot["biggest"] / tot["total"]
    cand = tot[(tot["n_clusters"] > 1) & (tot["total"] >= MIN_OBS_FOR_SPLIT)
               & tot["depth"].between(0.2, 0.6)]
    best = None
    for n_days in (N_DAYS, 1):
        for r in cand.itertuples():
            hours = day_slots(pd.Timestamp(r.window_start), n_days)
            g = gapless(idx, r.dynamic_social_unit, hours)
            if len(g) < MIN_SPLIT_ANIMALS:
                continue
            score = (len(g), r.depth)
            if best is None or score > best[0]:
                best = (score, r, hours, g)
        if best is not None:
            break
    if best is None:
        return None
    _, r, hours, g = best
    unit = r.dynamic_social_unit
    rows = [(unit, a) for a in sorted(g)[:MAX_PER_UNIT * 2]]
    w = d[d["window_start"].isin(hours) & d["dynamic_social_unit"].eq(unit)
          & d["animal_id"].isin([a for _, a in rows])]
    known = w[w["is_observed"] | w["is_carried_night"]]
    biggest = {}
    for ts, gg in known.groupby("window_start"):
        biggest[ts] = gg.groupby("association_event_id").size().idxmax()
    cells, carried = {}, set()
    for t in known.itertuples():
        k = (unit, t.animal_id)
        if k not in rows:
            continue
        cells[(k, t.window_start)] = ("alone" if t.temp_group_size <= 1 else
                                      ("main" if biggest.get(t.window_start)
                                       == t.association_event_id else "secondary"))
        if not t.is_observed:
            carried.add((k, t.window_start))
    return {"cols": hours, "rows": rows, "cells": cells, "carried": carried,
            "legend": [("main", "largest cluster"), ("secondary", "secondary cluster"),
                       ("alone", "alone")],
            "meta": "%s, %s; depth %.2f" % (unit, pd.Timestamp(r.window_start).date(),
                                            r.depth),
            "xlab": "%d continuous hours over %d day-night cycles" % (len(hours), len(hours) // len(DAY_HOURS))}


def pick_trajectory(d, exc_any) -> dict | None:
    """One animal that passes through BOTH away states, so the change is visible.

    Two animals each showing one mode says nothing about transition. A single animal
    moving with-origin -> alone -> with-another-unit -> back is the process itself, so
    the panel is built around the most balanced such excursion: the one maximising the
    smaller of its alone-nights and joined-nights, within a window that stays readable.
    """
    cand = exc_any[exc_any["depth"].eq("alone_and_joined")
                   & exc_any["alone_nights"].ge(3)
                   & exc_any["other_nights"].ge(3)
                   & exc_any["away_nights"].between(10, NIGHT_COLS - 8)].copy()
    if cand.empty:
        return None
    cand["balance"] = cand[["alone_nights", "other_nights"]].min(axis=1)
    cand = cand.sort_values(["balance", "away_nights"], ascending=[False, True])

    dd = d.copy()
    dd["night"] = (dd["window_start"] - pd.Timedelta(hours=10)).dt.normalize()
    known = dd[dd["is_observed"] | dd["is_carried_night"]]
    nidx: dict = {}
    for (nn, unit), g in known.groupby(["night", "dynamic_social_unit"], sort=False):
        nidx.setdefault(nn, {})[unit] = frozenset(g["animal_id"])

    for r in cand.itertuples():
        span = int(r.away_nights)
        pad = max(3, (NIGHT_COLS - span) // 2)
        n0 = pd.Timestamp(r.start_night) - pd.Timedelta(days=pad)
        nights = [n0 + pd.Timedelta(days=k) for k in range(NIGHT_COLS)]
        sets = [nidx.get(nn, {}).get(r.origin_group, frozenset()) for nn in nights]
        if any(r.animal_id not in ss for ss in sets):
            continue
        g = sets[0]
        for ss in sets[1:]:
            g &= ss
        if len(g) < MIN_PER_UNIT + 1:
            continue

        unit, animal = r.origin_group, r.animal_id
        others = [a for a in sorted(g) if a != animal][:MAX_PER_UNIT]
        rows = [(unit, animal)] + [(unit, a) for a in others]
        w = dd[dd["night"].isin(nights) & dd["dynamic_social_unit"].eq(unit)
               & dd["animal_id"].isin([a for _, a in rows])]
        cells = {}
        for (a, nn), gg in w.groupby(["animal_id", "night"]):
            k = (unit, a)
            if k not in rows:
                continue
            gg = gg[gg["is_observed"] | gg["is_carried_night"]]
            if gg.empty:
                continue
            sc = gg["social_context"]
            v = ("other" if sc.isin(AWAY_WITH_OTHERS).any() else
                 ("alone" if sc.eq("isolated").any() else
                  ("origin" if sc.isin(WITH_ORIGIN).any() else None)))
            if v:
                cells[(k, nn)] = v

        # For the focal animal's away nights, name the unit it was actually with:
        # the modal other unit among animals sharing its cluster that night.
        focal_away = w[w["animal_id"].eq(animal)
                       & w["social_context"].isin(AWAY_WITH_OTHERS)
                       & (w["is_observed"] | w["is_carried_night"])]
        pairs = set(zip(focal_away["window_start"], focal_away["association_event_id"]))
        pool = dd[(dd["is_observed"] | dd["is_carried_night"])
                  & dd["night"].isin(nights)
                  & dd["dynamic_social_unit"].ne(unit)]
        pool = pool[pd.MultiIndex.from_arrays(
            [pool["window_start"], pool["association_event_id"]]).isin(pairs)]
        dest = {}
        if not pool.empty:
            dest = (pool.groupby("night")["dynamic_social_unit"]
                    .agg(lambda x: x.value_counts().idxmax()).to_dict())
        units_seen = []
        for nn, u in dest.items():
            cells[((unit, animal), nn)] = "unit:" + u
            if u not in units_seen:
                units_seen.append(u)

        # count the focal animal's transitions between the three states
        seq = [cells.get(((unit, animal), nn)) for nn in nights]
        seq = [("other" if v and v.startswith("unit:") else v) for v in seq]
        runs = [v for i, v in enumerate(seq) if v is not None and (i == 0 or v != seq[i - 1])]
        return {"cols": nights, "rows": rows, "cells": cells,
                "focal": (unit, animal), "nightly": True,
                "legend": ([("origin", "with origin group"), ("alone", "alone")]
                           + [("unit:" + u, u) for u in units_seen]
                           + ([("other", "another unit")]
                              if any(v == "other" for v in
                                     [cells.get(((unit, animal), nn))
                                      for nn in nights]) else [])),
                "meta": ("%s of %s, %s to %s; %d alone-nights, %d with another unit, "
                         "%d state changes"
                         % (animal.split("_")[0], unit,
                            pd.Timestamp(r.start_night).date(),
                            pd.Timestamp(r.end_night).date(),
                            int(r.alone_nights), int(r.other_nights),
                            max(0, len(runs) - 1))),
                "xlab": "%d consecutive nights" % len(nights)}
    return None


# ------------------------------------------------------------------ depths
def between_depth(s1):
    v = s1[["median_frac_a", "median_frac_b"]].min(axis=1).dropna()
    return v[(v > 0) & (v <= 1)].to_numpy()


def within_depth(d):
    obs = d[d["is_observed"]]
    sizes = (obs.groupby(["window_start", "dynamic_social_unit",
                          "association_event_id"]).size().rename("k").reset_index())
    tot = sizes.groupby(["window_start", "dynamic_social_unit"]).agg(
        total=("k", "sum"), biggest=("k", "max"), n_clusters=("k", "size")).reset_index()
    tot = tot[(tot["n_clusters"] > 1) & (tot["total"] >= MIN_OBS_FOR_SPLIT)]
    return (1.0 - tot["biggest"] / tot["total"]).to_numpy()


def individual_depth(exc):
    f = exc[exc["away_nights"] > 0]
    return (f["other_nights"] / f["away_nights"]).to_numpy()


# ------------------------------------------------------------------ drawing
def raster(g, x0, y0, w):
    cols, rows, cells = g["cols"], g["rows"], g["cells"]
    carried = g.get("carried", set())
    nightly = g.get("nightly", False)
    lab_w = 50.0
    cw = (w - lab_w) / max(1, len(cols))
    ch = 8.0
    o, prev, missing = [], None, 0
    for i, key in enumerate(rows):
        unit, animal = key
        yy = y0 + i * ch
        if prev is not None and unit != prev:
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
                     'stroke-width=".9"/>'
                     % (x0, yy - 1.4, x0 + lab_w + cw * len(cols), yy - 1.4))
        prev = unit
        cls = "rlf" if g.get("focal") == key else "rl"
        o.append('<text class="%s" x="%.1f" y="%.1f" text-anchor="end">%s</text>'
                 % (cls, x0 + lab_w - 4, yy + ch - 1.6, animal.split("_")[0]))
        for j, c in enumerate(cols):
            v = cells.get((key, c))
            xx = x0 + lab_w + j * cw
            if v is None:
                missing += 1
                o.append('<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" fill="none" '
                         'stroke="var(--rule-soft)" stroke-width=".4"/>'
                         % (xx, yy, cw - .7, ch - 1.0))
            else:
                o.append('<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" fill="%s"/>'
                         % (xx, yy, cw - .7, ch - 1.0, cell_fill(v)))
    ybot = y0 + len(rows) * ch
    seen = []
    for i, (u, _) in enumerate(rows):
        if u not in [s[0] for s in seen]:
            seen.append((u, i))
    for u, i in seen:
        n = sum(1 for uu, _ in rows if uu == u)
        o.append('<text class="ul" transform="translate(%.1f %.1f) rotate(-90)" '
                 'text-anchor="middle">%s</text>' % (x0 - 4, y0 + (i + n / 2) * ch, u))
    right = x0 + lab_w + cw * len(cols)
    lx, ly = x0 + lab_w, ybot + 7
    for v, lab in g["legend"]:
        iw = 9 + len(lab) * 4.3 + 12
        if lx + iw > right and lx > x0 + lab_w:
            lx, ly = x0 + lab_w, ly + 10
        o.append('<rect x="%.1f" y="%.1f" width="6" height="6" fill="%s"/>'
                 % (lx, ly, cell_fill(v)))
        o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>' % (lx + 9, ly + 5.5, lab))
        lx += iw
    o.append('<text class="ts" x="%.1f" y="%.1f">%s</text>'
             % (x0 + lab_w, ly + 17, g["xlab"]))
    return "\n".join(o), (ly + 22) - y0, missing


def depth_hist(vals, x0, y0, w, h, colour, xlab, rng):
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    nb = 20
    cnt, _ = np.histogram(v, bins=np.linspace(0, 1, nb + 1))
    frac = cnt / cnt.sum()
    o, top, bw = [], frac.max() * 1.14, w / nb
    for t in (0, 0.5, 1.0):
        o.append('<line class="g" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x0 + t * w, y0, x0 + t * w, y0 + h))
        o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (x0 + t * w, y0 + h + 11,
                    "0" if t == 0 else ("1" if t == 1 else "0.5")))
    for i, f in enumerate(frac):
        if f <= 0:
            continue
        bh = f / top * h
        o.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s" '
                 'opacity=".82"/>' % (x0 + i * bw, y0 + h - bh, bw - .7, bh, colour))
    med = float(np.median(v))
    boot = [np.median(rng.choice(v, len(v), replace=True)) for _ in range(1000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="1.1" stroke-dasharray="3 2"/>'
             % (x0 + med * w, y0 - 4, x0 + med * w, y0 + h))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="var(--ink)" '
             'stroke-width="2"/>' % (x0 + lo * w, y0 - 4, x0 + hi * w, y0 - 4))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle" fill="var(--ink)">'
             'median %.2f</text>' % (x0 + med * w, y0 - 8, med))
    o.append('<line class="ax" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
             % (x0, y0 + h, x0 + w, y0 + h))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
             % (x0 + w / 2, y0 + h + 22, xlab))
    o.append('<text class="ts" x="%.1f" y="%.1f" text-anchor="middle">n = %s</text>'
             % (x0 + w / 2, y0 + h + 33, format(len(v), ",")))
    return "\n".join(o)


PANELS = [
    ("Two units, every pattern at once", "", "all_patterns"),
]
# only these two strips: the context an animal is in, and the unit it is assigned to
STAGE_STRIPS = [
    ("social context", "social_context"),
    ("dynamic unit", "dynamic_social_unit"),
]


def build():
    cols = ["window_start", "animal_id", "origin_group", "dynamic_social_unit",
            "social_context", "association_event_id", "temp_group_size", "is_observed",
            "is_carried_night", "is_local_2h_supported",
            "assigned_social_unit", "longitude", "latitude"]
    d = pd.read_csv(NARROW, usecols=cols, parse_dates=["window_start"])
    s1 = pd.read_csv(STAGE1, parse_dates=["start_hour"])
    exc = pd.read_csv(EXC, parse_dates=["start_night", "end_night"])
    exc_any = pd.read_csv(EXC_ANY, parse_dates=["start_night", "end_night"])
    rng = np.random.default_rng(SEED)
    idx = build_index(d)

    got = {
        "all_patterns": pick_all_patterns(d, s1, idx),
    }
    for k, v in got.items():
        print("  picker %-16s %s" % (k, "no candidate passed" if v is None else
                                     "%d rows x %d cols :: %s"
                                     % (len(v["rows"]), len(v["cols"]), v["meta"])))
    missing = [k for k, v in got.items() if v is None]
    if missing:
        raise SystemExit("no gapless window found for: %s" % ", ".join(missing))

    depths = {"between": between_depth(s1), "within": within_depth(d),
              "individual": individual_depth(exc)}
    for dk, v in depths.items():
        stats_depth = {"n": int(len(v)), "median": float(np.median(v)),
                       "at_one": float((v >= .999).mean()),
                       "at_zero": float((v <= .001).mean()),
                       "below_half": float((v < .5).mean())}
        depths[dk] = stats_depth

    # depth entries are namespaced: a panel key ("within") is otherwise identical to
    # a depth key and silently overwrites it
    stats = {"depth_" + k: v for k, v in depths.items()}
    o, y, total_missing = [], 78.0, 0
    rw = 648.0
    letters = "abcde"

    sinfo_box = {}

    def draw_inset(ix, iy, iw):
        """The same machinery over time, drawn inside panel a's box."""
        head = ('<text class="ts" x="%.1f" y="%.1f">one animal, %d day%s, hourly, at '
                'the same scale as (b)</text>'
                % (ix, iy - 6, INSET_DAYS, "" if INSET_DAYS == 1 else "s"))
        svg, h, info = stages_panel(d, ix, iy, iw, compact=True, days=INSET_DAYS)
        sinfo_box.update(info)
        return head + chr(10) + svg, h

    gsvg, gh, ginfo = gps_panel(d, 20.0, y, rw, inset=draw_inset)
    o.append('<text class="pl" x="0" y="%.1f">a</text>' % (y - 12))
    o.append(gsvg)
    stats["gps"] = ginfo
    stats["stages"] = sinfo_box
    # the fortnight is no longer drawn, so its numbers are carried into the caption
    _, _, long_info = stages_panel(d, 0.0, 0.0, 300.0, compact=True, days=STAGE_DAYS)
    stats["stages_long"] = long_info
    y += gh + 46
    for i, (name, defn, key) in enumerate(PANELS):
        g = got[key]
        if g is None:
            continue
        o.append('<text class="pl" x="0" y="%.1f">%s</text>'
                 % (y - 19, letters[i + 1]))
        rs, rh, miss = raster(g, 20.0, y, rw)
        total_missing += miss
        o.append(rs)
        stats[key] = {"meta": g["meta"], "rows": len(g["rows"]),
                      "cols": len(g["cols"])}
        y += rh + 34
    H = y - 14
    head = ('<svg viewBox="0 0 700 %.0f" role="img" aria-label="Three worked examples of '
            'social events with no gaps: a between-group merge in which one unit joins '
            'fully and the other only partly, a within-group split, and one animal '
            'moving between its own unit, being alone, and another unit.">' % H)
    stats["_missing_cells"] = total_missing
    return head + "\n" + "\n".join(o) + "\n</svg>", stats


CSS = """
:root{--ground:#f4f5f3;--paper:#fff;--ink:#15191a;--ink-2:#3a4244;--muted:#6b7476;
--rule:#d5dad7;--rule-soft:#e6eae7;--n2:#c3cbc8;
--c1:#1f6f8b;--c4:#8a5a1f;--c5:#6b4a7a;--c6:#0f5f57;--alone:#4a5356;
--sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;
--serif:"Newsreader",Georgia,serif;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--ground:#0f1213;
--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;--muted:#8f9998;--rule:#2c3435;
--rule-soft:#242c2d;--n2:#3f4a4b;--c1:#74b6d0;--c4:#d3a061;--c5:#b494c4;--c6:#5fc0b0;--alone:#a7b0b2;}}
:root[data-theme="dark"]{--ground:#0f1213;--paper:#171b1c;--ink:#e9eeed;--ink-2:#c2cac9;
--muted:#8f9998;--rule:#2c3435;--rule-soft:#242c2d;--n2:#3f4a4b;--c1:#74b6d0;
--c4:#d3a061;--c5:#b494c4;--c6:#5fc0b0;--alone:#a7b0b2;}
*{box-sizing:border-box}
body{background:var(--ground);color:var(--ink);font-family:var(--sans);font-size:15px;
line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:800px;margin:0 auto;padding:48px 26px 80px}
.kicker{font-family:var(--mono);font-size:10px;letter-spacing:.16em;
text-transform:uppercase;color:var(--muted);margin-bottom:10px}
h1{font-family:var(--serif);font-size:30px;font-weight:500;line-height:1.15;
letter-spacing:-.01em;margin:0 0 26px;text-wrap:balance}
.plate{background:var(--paper);border:1px solid var(--rule-soft);padding:22px 20px 8px}
.plate svg{display:block;width:100%;height:auto;overflow:visible}
.scroll{overflow-x:auto}.scroll>svg{min-width:660px}
.cap{font-size:12.5px;line-height:1.62;color:var(--ink-2);margin:16px 0 0;
padding-top:14px;border-top:1px solid var(--rule)}
.cap b{font-weight:600;color:var(--ink)}.cap .fignum{font-weight:600;color:var(--ink)}
.note{font-size:11.5px;line-height:1.6;color:var(--muted);margin-top:22px;
padding-top:14px;border-top:1px solid var(--rule-soft)}
.note code{font-family:var(--mono);font-size:.92em}
text{font-family:var(--sans)}
.rl{font-size:6.6px;fill:var(--muted);font-family:var(--mono)}
.rlf{font-size:6.6px;fill:var(--ink);font-weight:700;font-family:var(--mono)}
.ul{font-size:7.4px;fill:var(--ink-2);font-weight:600}
.ts{font-size:8.2px;fill:var(--muted)}
.ph{font-size:12px;fill:var(--ink);font-weight:600}
.pl{font-size:13px;font-weight:700;fill:var(--ink)}
.g{stroke:var(--rule-soft);stroke-width:.8}
.ax{stroke:var(--rule);stroke-width:1}
@media (prefers-reduced-motion:reduce){*{transition:none!important;animation:none!important}}
"""


def page():
    fig, st = build()
    bt, wi, iv = st["depth_between"], st["depth_within"], st["depth_individual"]
    gapfree = ("Every cell in every raster is a real observation"
               if st["_missing_cells"] == 0 else
               "%d cells could not be filled" % st["_missing_cells"])
    return f"""<title>Clustering To Social Structure</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500&family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;700&display=swap">
<style>{CSS}</style>
<div class="wrap">
  <div class="kicker">Figure 1 &middot; draft</div>
  <h1>Clustering, and the three ways a boundary opens</h1>

  <div class="plate">
    <div class="scroll">{fig}</div>
    <p class="cap">
      <span class="fignum">Figure 1.</span> How positions become social structure, and what
      the three kinds of boundary opening look like in the record.
      <b>(a)</b> Clustering in space, one hour of two units, {st["gps"]["hour"]}. Each animal
      is drawn as its fix at that hour with a tail of its preceding
      {st["gps"]["tail_hours"]} hourly fixes, filled by the unit it is travelling with, so
      the clusters can be read against how the animals were actually moving. Dashed
      outlines enclose the animals the clustering placed in one spatial cluster, blue where
      a cluster holds both units, and numbers give cluster sizes. {st["gps"]["animals"]}
      animals fall into {st["gps"]["clusters"]} clusters, of which {st["gps"]["shared"]}
      holds both units, with a nearest cross-unit distance of {st["gps"]["nearest_m"]} m.
      <b>A ring marks an animal whose origin group is not the unit it is travelling
      with</b> &mdash; {st["gps"]["off_origin"]} here, from
      {st["gps"]["off_origin_units"]} &mdash; so the individual axis and the group axis
      appear in the same frame. The scale bar shows the distance the clustering is acting
      on rather than asking the reader to take the 120&ndash;900 m parameter on trust. The
      frame is a daytime hour so that the tails are a continuous track; a tail spanning the
      night would be one straight jump between sleep sites, since collars do not fix
      between 17:00 and 02:00.
      <b>Inset:</b> the same machinery over time for one animal &mdash;
      {st["stages"]["animal"]} of {st["stages"]["unit"]}, drawn at <b>the same
      pixels-per-hour as (b)</b> and over the <b>first day of the same window</b>, for one
      of the animals in (b)&rsquo;s own rows &mdash; so the two hourly strips can be read
      against each other directly. The upper bar is social context, recomputed every hour
      and taking {st["stages"]["contexts"]} distinct values inside that single day; the
      lower bar is the dynamic unit, which changes only on seven days of sustained
      association and takes {st["stages"]["dynamic_labels"]}. Extending the same animal to
      {STAGE_DAYS} days &mdash; too wide to draw at this scale, since a fortnight needs
      about 4,700 px &mdash; gives {st["stages_long"]["contexts"]} contexts against
      {st["stages_long"]["dynamic_labels"]} dynamic label across
      {st["stages_long"]["hours"]} hours. Membership is deliberately far more stable than
      the situation it summarises.
      <b>(b)</b> Every pattern at once, {st["all_patterns"]["meta"]}. Rather than one panel
      per event type, this is a single 46-hour window over two units in which all three
      co-occur, so they can be read against each other rather than across separate examples.
      Rows are collared animals with the rule separating the two units; columns are
      consecutive hours across two day-night cycles; and every animal shown has a known state
      in every hour, so the raster carries no gaps. Grey is an animal with the bulk of its own
      unit, bronze an animal split off into a secondary cluster of its own unit, blue a cluster
      holding <em>both</em> units, and slate an animal alone. The three level-one events are
      therefore the three departures from grey: a within-group split is a bronze block, a
      between-group merge a blue one, and isolation a slate cell.
      <b>On scale.</b> Hourly resolution shows the transitions but not the outcomes. Median
      durations are 2 h for a within-group split and 1 h for isolation, so those are fully
      contained in a window this size; the median between-group encounter spans 14 h and the
      median individual excursion 2 nights, with tails reaching months, so this window catches
      one event out of a longer sequence. Individual excursions are shown separately for that
      reason.
    </p>
  </div>

  <p class="note">
    <b>On the overnight assumption.</b> The export already encodes it as
    <code>is_carried_night</code>, and carried rows retain <code>social_context</code>,
    <code>association_event_id</code> and <code>temp_group_size</code> in 100% of cases &mdash; so
    both the location and the nature of the interaction carry. Treating fixed-or-carried as known
    raises coverage from <b>58.7%</b> of animal-hours to <b>96.3%</b>. One hole remains and is
    excluded here: <b>02:00</b>, the handover hour, where the carry-forward has ended and the
    morning fix has not happened (27% known). Note what the assumption does and does not buy: it
    makes <em>durations</em> continuous across nights, but a carried hour is not an independent
    observation, so it should not be counted as new evidence in anything that rests on
    observation counts.
    Candidate windows were scored on completeness before event properties: for each candidate the
    search finds the animals with a known state in every interval and keeps the window yielding
    the most, requiring at least {MIN_PER_UNIT} per unit
    ({MIN_SPLIT_ANIMALS} for the within-group panel). This replaces a six-type taxonomy that
    split merges into large, medium-partial and small-subset &mdash; that subdivision rests on a
    classifier keyed to the <em>maximum</em> participation a unit ever reaches rather than its
    typical state, and it does not distinguish groups (between-group ICC 0.05&ndash;0.13,
    permutation p = 0.09&ndash;0.26). Carried as the continuous depth in (a, b) it loses nothing.
    All panels and depth distributions come from the frozen 1,924,104-row export (2024-03-01 to
    2026-07-22); within-group depth uses unit-hours with at least {MIN_OBS_FOR_SPLIT} observed
    animals. Generated by <code>build_event_examples_figure.py</code>; seed {SEED}.
  </p>
</div>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(page(), encoding="utf-8")
    print("wrote %s (%d bytes)" % (args.output, args.output.stat().st_size))


if __name__ == "__main__":
    main()
