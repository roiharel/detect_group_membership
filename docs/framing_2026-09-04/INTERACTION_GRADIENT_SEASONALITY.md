# Is the *type* of interaction seasonal? The full radius gradient

Prepared 4 September 2026. Corrects the framing in [SEASONALITY_AND_TRENDS.md](SEASONALITY_AND_TRENDS.md) §3.

Script: [analyze_seasonality_interaction_gradient.py](../../analyze_seasonality_interaction_gradient.py).
Outputs: `outputs/general_structure_2026_09/phase4c_seasonality/interaction_gradient/`.

---

## What was wrong with the earlier answer

I wrote that the encounter season is "geography, not sociality," and treated that as settling the question. It does not. It settles one question — *do two groups end up in the same place* — and then I collapsed everything downstream of a merge into a single number, the edge-weighted cross-group contact fraction at 5 m.

That is a thin reading. Once two groups occupy the same spatial cluster, at a scale of hundreds of metres, what happens inside spans a gradient, and each level is a different social question:

| Scale | Question | Events | Dyads |
|---|---|---:|---:|
| ~10² m (cluster) | how much of **each** group actually joins? | 1,480 | 66 |
| 5 m | cross-group contact rate; does the merged mass hold two origin communities or one? | 716 | 12 |
| 2 m | the same, where only real proximity survives | 654 | 12 |
| 2 m / 5 m | how much of the 5 m structure survives tightening? | 323 | 12 |
| betweenness | who brokers, and how concentrated is it? | 1 dyad, 77 weeks, 4 radii | 1 |

A merge can be large but structurally segregated, or small but thoroughly mixed. A single 5 m cross-fraction cannot distinguish those. **Power also decays down the gradient**, which is itself a result: a null at the cluster scale rests on 66 dyads, a null at 2 m on 12, and betweenness on one.

---

## What the gradient shows

Annual harmonic, trend term, duration control, GEE clustered on dyad, well-observed window from 2024-12-01.

| Level | Response | n | Dyads | Amplitude | Peak | Joint *p* |
|---|---|---:|---:|---:|---|---:|
| **cluster** | **mutual participation** — logit min(frac A, frac B) | 1,480 | 66 | 1.481 | Apr | **0.0033** |
| **cluster** | **large merge rather than partial** | 1,480 | 66 | 0.714 | Mar | **0.0029** |
| cluster | log largest joint cluster size | 1,480 | 66 | 0.078 | Jan | 0.184 |
| 5 m | mixture deficit vs composition expectation | 716 | 12 | 0.401 | Sep | 0.320 |
| 5 m | *same, duration held out* | 716 | 12 | 0.399 | Aug | 0.084 |
| 5 m | **within-merge modularity** | 716 | 12 | 0.032 | Oct | **0.066** |
| 5 m | within-merge entropy | 716 | 12 | 0.032 | Sep | 0.178 |
| 5 m | within-merge balance | 716 | 12 | 0.048 | Sep | 0.300 |
| 5 m | fraction of bins with any cross-contact | 716 | 12 | 0.017 | Jul | 0.327 |
| 2 m | mixture deficit | 654 | 12 | 0.070 | Jul | 0.911 |
| 2 m | **within-merge modularity** | 654 | 12 | 0.031 | Oct | **0.061** |
| 2 m | entropy / balance / contact fraction | 654 | 12 | 0.028–0.055 | Sep–Jun | 0.123–0.180 |
| ratio | log(2 m contact / 5 m contact) | 323 | 12 | 0.093 | Oct | 0.196 |

So the picture is **not** a flat null. Two things at the coarse end look seasonal at *p* < 0.005, and within-merge modularity is borderline at *p* ≈ 0.06 **at both radii with the same October phase** — consistency across radii being the kind of thing that argues against noise.

Descriptively, mutual participation runs from **0.95 in May–June to 0.62–0.69 in August and November** (dyad-balanced: 0.949, 0.924 → 0.689, 0.743). That is a large swing, and it survives a balanced panel of the 11 dyads spanning the window (*p* = 0.0036).

---

## And then it fails, for a reason that is about the record, not the animals

The well-observed window runs 2024-12 to 2026-07, holding **1.63 annual cycles**. That is not evenly distributed across the calendar:

| Sampled in two or more years | Sampled in one year only |
|---|---|
| Jan, Feb, Mar, Apr, May, Jun, Jul, Dec | **Aug, Sep, Oct, Nov** |

August to November exist **only in 2025** — and those are precisely the low-participation months. So the fitted harmonic's trough is a single autumn.

Two checks, both fatal:

| Check | Result |
|---|---|
| Refit on months with 2+ years of data (Jan–Jul, Dec) | *n* = 1,193, 60 dyads, **joint *p* = 0.511** |
| Season demeaned within dyad | participation *p* = 0.169; large-merge *p* = 0.116 |
| Same calendar month, between-year spread in participation | up to **0.236** — as large as the apparent seasonal range |

December makes the point on its own: mean participation 0.643 in 2024 against 0.880 in 2025, the same month one year apart.

**The participation signal is a single-year late-2025 feature.** With Aug–Nov sampled once, no harmonic fitted to this record can separate an annual cycle from one unusual autumn. The borderline modularity signals at 5 m and 2 m share the same October phase and therefore the same limitation.

---

## The corrected bottom line

Six statements, ordered down the gradient:

1. **Whether two groups end up in the same place varies through the year — and that is geography.** The signal decays monotonically under distance restriction (amplitude 0.285 → 0.087; *p* 0.010 → 0.629) and dies at 5 km.
2. **Whether they meet, given they are already close, has no seasonal pattern.** This is the strongest null: 14,775 dyad-days, 115 dyads.
3. **How much of each group actually joins a merge varies a great deal — from about 62% to 98% — but the variation is not established as seasonal.** The low months come from one autumn, and the same month differs between years by as much as the apparent cycle. Report it as documented variation in social context, not as a season.
4. **How closely they mix once merged, at 5 m or 2 m, shows no seasonal pattern**, and neither does how much of the 5 m structure survives tightening to 2 m.
5. **Whether the merged mass stays structurally two groups is borderline** (*p* ≈ 0.06 at both radii, October phase) and rests on the same single-year window. An open signal, not a result.
6. **Betweenness within a merge is measurable in one dyad only** (Copper–Lilac, 38 animals, 4 radii). Its weekly series looks strongly seasonal at every radius, with the amplitude falling as the radius widens (1.74 at 1 m to 0.53 at 20 m) — but it is one dyad passing through a fusion, clustered on 7 quarter-blocks, and the source metadata records that its fusion hours saturate after 2025-08-01. A secular transition would alias with season here. **Not reportable as seasonality.**

### So the earlier one-liner needs replacing

Not: *opportunity is seasonal, the crossing decision is not.*

Instead:

> Where groups are follows the year. Whether they meet once close does not. **How much of each group joins a merge varies widely, and this record cannot say whether that variation is seasonal** — the months where participation is lowest are sampled in a single year. How closely they mix once merged shows no seasonal pattern at any radius tested.

The middle clause is the honest change, and it is the user's point: the type of interaction is a real and separate question from the fact of meeting, it does vary, and collapsing it into one 5 m number hid that.

---

## What would settle it

1. **A third and fourth annual cycle.** August–November need at least two years each. Nothing else in this list matters as much.
2. **A second harmonic**, for bimodal rainfall — untestable at 1.6 cycles.
3. **More fine-scale dyads.** 12 of 68 have 2-minute products, so every 5 m and 2 m row above is clustered on 12 units. Gap 7.5.
4. **Betweenness beyond Copper–Lilac.** The within-merge brokerage question is currently a one-dyad case study. Extending it needs the same fine-scale regeneration as (3).
5. **A within-group cohesion control on participation.** Mutual participation correlates with joint cluster size at *r* = 0.364, and low participation could reflect a group being internally split rather than declining to join. Axis B found no seasonality in splitting, which is reassuring, but the two have not been modelled together.
