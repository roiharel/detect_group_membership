# Extending panel 2d from hours to 2 minutes

## What it would take, and what the data will actually support

Panel 2d currently spans **1 hour to ~13,000 hours** because every state in it is assigned
hourly. The proposal is to reach **2 minutes**, giving a single axis from 2 minutes to a
year. That is possible, but not as one continuous measurement, and the reason is in the
data rather than in the code.

---

## 1. The hard constraint: eleven hours a day have no fixes at all

Measured on the canonical file, two weeks of 2026-02 (922,218 fixes, 241 animals):

| quantity | value |
|---|---:|
| 2-minute bins in the span | 10,021 |
| bins with at least one fix | **5,544 (55.3%)** |
| median animals per bin, inside the fixing window | **153** of 241 |
| 90th percentile animals per bin | 184 |
| bins with ≥12 animals | 54.6% |

Median animals per 2-minute bin, by hour UTC:

```
00:0 01:0 02:0  03:147 04:161 05:165 06:166 07:165 08:166 09:163
10:165 11:166 12:166 13:166 14:166 15:164  16:0 17:0 18:0 ... 23:0
```

**Inside 03:00–15:00 UTC the density is excellent** — 153 of 241 animals present in a
typical 2-minute bin, which is more than enough to cluster. **Outside it there is
nothing.** Not sparse: zero. The hourly pipeline covers those eleven hours with
`carried_night`, which is a defensible assumption at hour scale and an indefensible one
at 2-minute scale — it would assert 2-minute state transitions from a position that was
last measured before dusk.

**Therefore the floor is not 2 minutes across the board. It is 2 minutes for states that
begin and end inside one daylight window (at most ~13 h), and one hour beyond that.**

## 2. The existing 2-minute products cannot be reused

`canonical_5m_shared_history/canonical_5m_2min_merge_metric_rows.csv` looks like the right
input but is not:

| | |
|---|---:|
| rows | 159,979 |
| distinct 2-minute bins | 116,396 |
| **dyads** | **12** |
| span | 2024-07-15 → 2026-06-12 |

It is built *from Stage-1 encounter windows* with `--max-pairs 20`, so it exists only for
dyads that already had a detected encounter, and only during those encounters. It is a
contact-depth product, not a membership product, and it cannot yield durations for
arbitrary units, dyads or animals.

## 3. What has to be built

**(a) A 2-minute membership pass, daylight only.** The same adaptive rule
(k=2, factor 1.65, clipped 120–900 m) applied per 2-minute bin instead of per hour, over
03:00–15:00 UTC. Scale: ~366,000 populated bins over the record against 18,101 hours, so
roughly **20x the clustering work** of the hourly build — hours, not minutes, and it
should be run as an overnight job. A full 2-minute membership table would be on the order
of **56 million rows** (~2 GB as parquet, ~10 GB as CSV), so it should be written as
parquet and, better, **not persisted at all** — the pipeline should emit bout intervals
directly and discard the per-bin table.

**(b) A dwell rule, because GPS error will manufacture short bouts.** Position error is
roughly 5–15 m against a linkage threshold near 198 m, so two animals near the threshold
will flicker in and out of the same cluster from noise alone. A state must persist for
*k* consecutive bins before it counts as a transition. `k` should be chosen from the
noise floor in (c), not picked by hand.

**(c) An empirical noise floor.** Re-cluster each frame with positions jittered within
their error, repeatedly, and measure the bout-length distribution that arises from error
alone. Report durations only above the length at which the observed distribution
separates from the jittered one. Without this the short end of the axis is unreadable:
a 2-minute "bout" is not obviously a biological event.

**(d) Right-censoring, and this is the part that decides whether the figure is honest.**
Every bout that is still running when the fixes stop at ~15:00 is **truncated, not
finished**. Treating those as completed bouts would pile up spurious mass at ~12 h and
bias every short-end statistic downwards. The short arm must be summarised with a
survival estimator (Kaplan–Meier, censoring at the daylight boundary), not with a raw
empirical CDF like the current panel uses.

**(e) A cohort statement.** The 22 collars that fix every two hours have a fix in roughly
one bin in sixty. They cannot take part in a 2-minute pass, and interpolating them to
2-minute resolution would fabricate ~60 positions per gap. **The short arm therefore
describes a different cohort from the long arm** — the ~311 dense collars — and the long
arm should be recomputed on that same cohort for the comparison in (f) to mean anything.

**(f) An overlap band, measured both ways.** Between about **1 hour and 13 hours** both
resolutions can measure the same bouts. Build the figure only if the two agree there.
That overlap is what licenses splicing them onto one axis; without it, the joint axis is
two different measurements drawn as though they were one.

## 4. What the figure would then show

A duration axis from 2 minutes to 1 year, drawn as **two arms with a marked seam**:

* **2 min – 13 h** — 2-minute resolution, daylight only, dense collars, Kaplan–Meier with
  censoring at dusk, shaded below the empirical noise floor.
* **1 h – 1 year** — the current hourly product, unchanged.
* **1 h – 13 h** — both, plotted together as the agreement check.

The seam must be drawn, not hidden. A reader who cannot see where the resolution changes
will read the short end as though it had the same standing as the long end, and it does
not: the short end is daylight-only, dense-collar-only, and censored.

## 5. Recommended order of work

1. Build the 2-minute pass for **one month** first, not the whole record. Enough to
   produce the noise floor, choose the dwell rule, and test the overlap band. Cheap, and
   it answers whether the rest is worth running.
2. If the overlap band agrees, run the full record overnight, emitting bout intervals
   only.
3. Rebuild panel 2d as a two-arm figure with the seam and the censoring shown.

If the overlap band does **not** agree, that is the result, and the honest figure keeps
the hourly axis and states why 2-minute durations were not adopted.
