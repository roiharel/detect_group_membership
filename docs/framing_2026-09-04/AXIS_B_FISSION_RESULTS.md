# Axis B: fission, modularity, and what predicts it

Prepared 4 September 2026. Amends [STABLE_GROUPS_FLUID_BOUNDARIES.md](STABLE_GROUPS_FLUID_BOUNDARIES.md) and completes the third axis alongside [INDIVIDUAL_AXIS_MERGE.md](INDIVIDUAL_AXIS_MERGE.md).

Script: [analyze_axis_b_fission_structure.py](../../analyze_axis_b_fission_structure.py).
Outputs: `outputs/general_structure_2026_09/phase4b_axis_b_fission/`.

Three questions: is within-group modularity real once coverage is controlled and is it a group property or a state; does general co-sitting or immediately-preceding proximity predict who ends up on which side of a split; and do environmental conditions predict when a group splits and for how long.

Answers, in order: **real and a state; general co-sitting, not pre-split proximity; no.**

> **Source caveat throughout.** The weekly network metrics derive from the legacy hourly source filtered to 2025-01-01 onward, while axes A and C use the frozen 1,924,104-row export. Cross-axis pooling is still blocked on gap 7.3. NDVI (Sentinel-2, per group per day) runs 2024-02-29 to 2026-04-29, so Part 3 stops three months short of the membership record. `n_animals_observed` is collar coverage, never group size.

---

## Part 1 — Modularity is real above about twelve collars, and it is a state

### Detectability scales with collar count, so nothing can be compared unmatched

| Collars observed | Group-weeks | Units | Weeks modular | Weeks with a split | p90 modularity | Median communities |
|---|---:|---:|---:|---:|---:|---:|
| 1–4 | 491 | 24 | **1.2%** | 28.7% | 0.000 | 1 |
| 5–7 | 306 | 19 | 1.3% | 61.4% | 0.000 | 1 |
| 8–9 | 184 | 12 | 3.3% | 65.8% | 0.000 | 1 |
| 10–11 | 123 | 15 | 9.8% | 82.1% | 0.000 | 1 |
| 12–14 | 175 | 12 | 12.6% | 76.6% | 0.007 | 1 |
| 15–19 | 107 | 6 | 31.8% | 81.3% | 0.026 | 1 |
| 20+ | 32 | 2 | **46.9%** | 78.1% | 0.064 | **2** |

Spearman correlation with collar count: **0.225** for modularity, **0.489** for split-timestamp fraction. Half the variation in *split detection* is collar count. That is the axis-B analogue of axis A's fabricated negatives, and it means the 4,570-event inventory has a coverage-driven detection rate that no unadjusted rate can survive.

### At matched coverage, groups genuinely differ

Restricting to group-weeks with **12–16** collars observed, so every group is measured with a comparable number of nodes:

| Unit | Group-weeks | Median collars | Nominal size | Weeks modular | Max modularity | Median largest-community fraction |
|---|---:|---:|---:|---:|---:|---:|
| Lilac | 39 | 13 | 60 | **53.8%** | 0.458 | 0.750 |
| Chartreuse | 20 | 15 | 74 | **45.0%** | 0.065 | 0.967 |
| FireOpal | 7 | 12 | 62 | 14.3% | 0.006 | 1.000 |
| Copper | 8 | 13.5 | 50 | 12.5% | 0.003 | 1.000 |
| PhantomWest | 9 | 13 | 78 | 11.1% | 0.080 | 0.917 |
| Maroon | 9 | 12 | 42 | 11.1% | 0.018 | 1.000 |
| Lapis | 45 | 13 | 74 | 2.2% | 0.007 | 1.000 |
| Magenta | 7 | 13 | 56 | **0.0%** | ~0 | 1.000 |
| Purple | 34 | 14 | 62 | **0.0%** | ~0 | 1.000 |
| RubyRunners | 27 | 13 | 64 | **0.0%** | ~0 | 1.000 |

Lilac is modular in more than half its well-covered weeks and reaches a modularity of 0.458. Purple, RubyRunners and Magenta are modular in *none* of theirs, at the same collar count. **That is a real between-group difference and it is reportable.**

This revises the earlier framing, which said modularity should not be reported as biology. That judgement was correct for the full 20-group set — where median modularity is zero in 19 of 20 groups purely because most are tracked with 5–8 collars — and wrong for the well-covered subset. The correct statement is narrower and stronger: **above about twelve collars, groups differ in internal cohesion; below it, the measure returns a single cluster regardless.**

### But modularity is mostly a state, not a group property

Random-intercept model on well-covered group-weeks (314 group-weeks, 12 units, ≥12 collars):

| | Value |
|---|---:|
| Between-unit variance | 0.000319 |
| Residual (within-unit, over time) variance | 0.001835 |
| **ICC, between-unit** | **0.148** |
| Collar-count coefficient | 0.00064 (p = 0.43) |

Only **15% of the variance in modularity is between groups; 85% is within group across time.** And once above the well-covered threshold, collar count no longer predicts modularity at all (p = 0.43) — which is the check that licenses this subset.

Week-to-week persistence, on units with at least 30 consecutive week-pairs:

| Unit | Consecutive week-pairs | Lag-1 autocorrelation |
|---|---:|---:|
| Lilac | 37 | **0.413** |
| Chartreuse | 73 | **0.331** |
| Purple | 33 | 0.132 |
| RubyRunners | 38 | −0.027 |
| Lapis | 75 | −0.018 |

So modularity is **episodic but sticky**: a group enters a modular phase and stays in it for some weeks, rather than fluctuating week to week. It is a state with duration, which is exactly the shape axes A and C found.

**This is Result 5's third leg.** Extreme permeability is a state, not an entity — now demonstrated on all three axes, with three different measurements:

| Axis | The state | Evidence |
|---|---|---|
| A | sustained association | 16 of 1,705 encounters, 5 of 68 dyads; contrast +1.376 [+0.338, +1.543] |
| B | modular phase | ICC 0.148 — 85% of variance within group over time; lag-1 *r* 0.33–0.41 in the two groups that enter it |
| C | terminal excursion | 61 of 338 excursions, 16 of 91 animals reach settlement |

---

## Part 2 — The split-composition signal is a bond, not a pre-split cue

### What was tested

The saved analysis (mean event AUC 0.568) uses only the cohesive hours **immediately preceding** each split onset. That measures the animals' current positioning, not their relationship. The obvious alternative was never tested: does a dyad's *general* co-sitting rate predict which side of a split it ends up on?

General bond, per event-dyad: the dyad's 5 m co-sitting rate over **all other cohesive hours in its group, excluding a ±7-day window around the focal split onset**. A temporal holdout, so the bond measure cannot contain the event it predicts.

### Both predict, equally, and they are not redundant

Identical retained sample throughout — 461 events, 15 units, 27,812 event-dyads, minimum 10 co-observed 2-minute bins per dyad:

| Predictor | Mean event AUC | 95% CI |
|---|---:|---|
| Immediately-preceding proximity | 0.5702 | 0.5601 – 0.5816 |
| **General bond** (±7-day holdout) | **0.5742** | 0.5617 – 0.5873 |
| Paired difference (general − preceding) | +0.0039 | −0.0093 – +0.0179 |

The preceding-proximity figure reproduces the saved 0.568. The general bond does **at least as well**, and the paired difference straddles zero.

Both in one model, each z-scored then demeaned within event so every comparison is between dyads of the same split — Binomial GEE clustered on event, 26,058 event-dyads in 465 events:

| Term | OR per SD, together | 95% CI | OR per SD, alone |
|---|---:|---|---:|
| General bond | **1.286** | 1.215 – 1.360 | 1.347 |
| Preceding proximity | **1.194** | 1.137 – 1.254 | 1.269 |

Each retains roughly 90% of its solo effect with the other in the model, so they carry partly independent information — but the general bond is the larger of the two, and it is measured on weeks that have nothing to do with the split.

### What this changes

The earlier statement was "animals that co-sat within 5 m before a split are more likely to end up on the same side of it." That is true but misattributed. The effect is **not specifically a pre-split signal** — a dyad's ordinary, long-run co-sitting rate does the same work. Splits follow standing social structure; they are not anticipated by immediate positioning.

The effect also remains weak, and that has not changed: AUC ≈ 0.57, and the saved mean event ARI of **0.100** still says split composition is mostly not predicted. Report the reframing and the ARI together.

---

## Part 3 — Environment does not predict when groups split, or for how long

### The occurrence response, tested eight ways

Group-weeks are axis B's opportunity table: a group-week with observation support and no split is an observed negative. 836 group-weeks in 19 units retained (requires a successful NDVI week for that group and ≥5 collars).

Two NDVI terms are kept separate, following the axis-A design rule: `ndvi_unit_mean` is the unit's average range greenness — a **between-group** term — and `ndvi_deviation` is that week's departure from the unit's own norm, the **within-group** term. Only the second can answer "does *this group* split more when conditions are greener."

| Subset | Response | n | Base rate | Between-group NDVI | **Within-group NDVI** | Collar coverage |
|---|---|---:|---:|---|---|---|
| all | any split timestamp | 836 | 0.713 | 1.596 [1.051, 2.426] | 1.029 [0.902, 1.174] | 1.543 [0.936, 2.545] |
| all | split in >50% of timestamps | 836 | 0.246 | 1.340 [0.887, 2.025] | 1.088 [0.874, 1.356] | 2.083 [1.327, 3.267] |
| all | any multi-animal split | 836 | 0.392 | 2.194 [1.482, 3.247] | 1.005 [0.833, 1.214] | 1.530 [1.046, 2.237] |
| all | multi-animal split >25% | 836 | 0.106 | 2.906 [1.891, 4.466] | **0.768 [0.592, 0.997]** | 2.080 [1.380, 3.135] |
| ≥12 collars | any split timestamp | 284 | 0.799 | 4.476 [1.217, 16.465] | 0.926 [0.759, 1.130] | 1.190 [0.745, 1.903] |
| ≥12 collars | split in >50% of timestamps | 284 | 0.433 | 2.429 [1.919, 3.074] | 0.975 [0.716, 1.328] | 1.192 [0.729, 1.949] |
| ≥12 collars | any multi-animal split | 284 | 0.525 | 4.504 [2.080, 9.753] | 1.117 [0.829, 1.505] | 1.131 [0.882, 1.449] |
| ≥12 collars | multi-animal split >25% | 284 | 0.239 | 5.490 [2.067, 14.579] | **0.703 [0.545, 0.907]** | 0.828 [0.339, 2.024] |

Odds ratios per standard deviation, Binomial GEE clustered on unit. Three readings:

1. **The within-group NDVI term is null in six of eight specifications**, and in the two strictest it turns *negative* with the interval barely excluding 1 — greener-than-usual weeks have marginally **fewer** substantial multi-animal splits. There is no support for an environmental trigger, and what little signal exists points toward scarcity fragmenting, not abundance permitting.
2. **The between-group NDVI term is always large** (OR 1.34 to 5.49) and is not evidence of anything ecological. It is range greenness across 19 units, inseparable from group identity, size and range. This is the same failure mode as axis A's seasonal term, which fell from OR 0.62 [0.44, 0.86] to 0.86 [0.69, 1.07] once restricted to within-opportunity comparisons.
3. **Collar coverage predicts split detection in the full sample** (OR 1.53–2.08) and stops predicting it above twelve collars (OR 0.83–1.19). That validates the well-covered restriction from Part 1 and confirms the detection confound is real below it.

The base rate of the originally reported response is 0.713 — a group-week with thirteen collars almost always contains *some* timestamp with animals in more than one cluster. The stricter responses are the meaningful ones, and they are where the within-group term goes negative.

### Extent, given a split week

Gaussian GEE on split-timestamp fraction among the 596 split weeks: both NDVI terms null (between-group 1.030 [0.966, 1.099]; within-group 1.009 [0.969, 1.050]), and the only term with an interval excluding zero is **collar coverage** (1.132 [1.063, 1.206]). Split *extent*, as currently defined, measures observation rather than behaviour.

### Duration

655 events with NDVI on the onset day, 14 units. Median persistent split duration **2 h**, p90 **13 h**, max **240 h**. Log-duration GEE clustered on unit:

| Term | Coefficient | 95% CI |
|---|---:|---|
| Within-group NDVI deviation | +0.0505 | +0.0021, +0.0989 |
| Between-group NDVI mean | −0.0456 | −0.2179, +0.1267 |
| Animals observed at onset | **+0.2367** | +0.0536, +0.4197 |

The within-group NDVI term is positive and marginal — its lower bound is 0.002 — and it points the **opposite way** from the occurrence result: greener weeks give slightly *longer* splits but slightly *fewer* substantial ones. Two marginal effects with opposite signs are not a finding. The largest term is again collar coverage.

### Verdict

**No robust environmental effect on fission — not on whether a group splits, not on how much of the week it spends split, and not on how long a split lasts.** Every apparently strong environmental coefficient is a between-group range term that cannot be separated from group identity, and the two within-group signals are marginal and mutually contradictory.

This retires the symmetric prediction in the old Phase 5 plan — *scarcity fragments, abundance connects* — on the fission side, as written. Report it as a tested and unsupported prediction, which is more useful than leaving it as an unstated ambition.

---

## What this settles, and what it leaves

**Settled.**

| Question | Answer |
|---|---|
| Is within-group modularity measurable? | Yes, above ~12 collars. Below that the estimator returns one cluster regardless — the 19-of-20-groups-at-zero figure is a measurement artefact |
| Do groups differ in internal cohesion? | Yes, at matched coverage. Lilac 53.8% and Chartreuse 45.0% of weeks modular against Purple, RubyRunners and Magenta at 0.0%, all at 13–15 collars |
| Is modularity a group property or a state? | **A state.** ICC 0.148 — 85% of variance is within group over time — but sticky, with lag-1 *r* of 0.33–0.41 in the two groups that enter it |
| Does pre-split proximity predict split composition? | Yes, weakly — AUC 0.570 |
| Does *general* co-sitting predict it? | Yes, at least as well — AUC 0.574, paired difference +0.004 [−0.009, +0.018]. **The effect is a bond, not a pre-split cue** |
| Are the two redundant? | No. Together: bond OR 1.286 [1.215, 1.360], preceding OR 1.194 [1.137, 1.254] per SD, within event |
| Does environment predict split timing? | **No.** Within-group NDVI null in 6 of 8 specifications, weakly negative in the strictest two |
| Does environment predict split duration? | **No.** Marginal positive, opposite in sign to occurrence, smaller than the coverage term |

**Left open.**

1. **Cohort re-derivation (gap 7.3) is now the last structural blocker on this axis.** Parts 1 and 3 run on the legacy source from 2025-01-01; axes A and C are on the frozen export from 2024-03-01. Re-deriving the weekly network metrics on the frozen source is required before any cross-axis table.
2. **Split extent needs a coverage-independent definition.** As it stands the response is predicted by collar count and nothing else. A rate over group-weeks with matched coverage, or a modularity-based extent, would fix it.
3. **The modular phases themselves are unexamined.** Lilac and Chartreuse enter modular states that persist for weeks. Whether those phases precede a permanent fission — LapisSplinter exists, so at least one group in this population did divide — is the obvious next question, and it connects axis B to axis C directly.
4. **Environmental predictors beyond NDVI are untested here.** The Sentinel-2 product carries NDWI, LSWI and EVI. Water availability is a different mechanism from vegetation greenness and may behave differently; NDVI alone is not a complete test of the ecological hypothesis, only of the version that was stated.
