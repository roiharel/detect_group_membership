# Seasonality and trends, on all three axes

Prepared 4 September 2026. Completes the framing in [STABLE_GROUPS_FLUID_BOUNDARIES.md](STABLE_GROUPS_FLUID_BOUNDARIES.md).

Script: [analyze_seasonality_and_trends.py](../../analyze_seasonality_and_trends.py).
Outputs: `outputs/general_structure_2026_09/phase4c_seasonality/`.

**Short answer.** The environment in this system is strongly seasonal. Boundary crossing is not. The only seasonal signal anywhere is in **range overlap** — where groups are — and it disappears once you condition on the groups already being close. No secular trend survives its effort control on any axis.

---

## 1. Observation effort rises 121-fold, so read no trend before this

| | Month | Collared animals | Units | Nightly resolvable |
|---|---|---:|---:|---:|
| First | 2024-02 | 2 | 1 | 0.00 |
| Peak | 2026-03 | **242** | 23 | **0.95** |

Resolvability crosses 0.87 in November 2024 and stays above 0.87 thereafter, so the **well-observed window is 2024-12-01 onward** — and that window contains only **1.53–1.63 seasonal cycles**. That is the binding limit on everything below, and no interval reported here fixes it.

The good news is that the confound is *separable* rather than fatal. In the retained window:

| Diagnostic | Axis A | Axis B | Axis C |
|---|---:|---:|---:|
| Seasonal cycles in window | 1.63 | 1.53 | 1.63 |
| corr(season, trend) | 0.34 / −0.30 | 0.21 / −0.09 | 0.26 / −0.31 |
| corr(season, effort) | 0.03 / 0.07 | 0.01 / 0.12 | 0.09 / −0.01 |
| Max VIF across terms | 1.25 | 1.10 | 1.21 |

Season and effort are essentially uncorrelated, and variance inflation is negligible. So a seasonal coefficient here is not an effort artefact — it is simply estimated from 1.6 cycles, which is thin. Trends are a different matter, and §3 shows why.

---

## 2. The environmental season is real and strong

NDVI deviation from each unit's own mean, 7,423 group-days across 23 units, same window and same harmonic:

| | Value |
|---|---|
| Seasonal amplitude | **0.046** NDVI units |
| Peak | **day 199 ≈ 18 July** |
| Joint Wald *p* | **< 0.0001** |

This matters because it removes the easy explanation for a null. The absence of seasonality on axes B and C below is not because there is no season for the animals to respond to.

---

## 3. Axis A — the season is in geometry, not in interaction

### The signal looks real at first

Encounter occurrence on the opportunity table, 67,501 dyad-days, 317 dyads, Binomial GEE clustered on dyad:

| Term | OR per SD | 95% CI |
|---|---:|---|
| season_sin | 0.967 | 0.838 – 1.116 |
| season_cos | **0.820** | 0.710 – 0.946 |
| trend_months | **0.748** | 0.615 – 0.909 |
| shared observed hours | 1.039 | 0.955 – 1.129 |

Seasonal amplitude **0.285** logit, peak **day 193 ≈ 12 July**, joint Wald **p = 0.0103**. Dropping the effort term changes nothing (0.822 vs 0.820), and the peak lands within **six days** of NDVI's own peak. It is a tempting result.

### It fails the restriction that killed the Phase 4 seasonal term

A seasonal encounter signal can mean two different things: groups' ranges overlap seasonally (**geometry**), or groups already near each other meet seasonally (**interaction**). Restricting to dyad-days where the units were already close removes the first.

| Restriction | Dyad-days | Dyads | Encounter rate | Seasonal amplitude | Peak | Joint *p* |
|---|---:|---:|---:|---:|---|---:|
| all supported dyad-days | 67,501 | 317 | 0.037 | 0.285 | Jul | **0.0103** |
| min centroid distance ≤ 10 km | 30,794 | 165 | 0.080 | 0.250 | Jul | **0.0074** |
| min centroid distance ≤ 5 km | 14,775 | 115 | 0.160 | 0.178 | Jul | 0.231 |
| min centroid distance ≤ 2 km | 5,102 | 82 | 0.433 | 0.087 | Oct | 0.629 |

The amplitude decays monotonically — 0.285 → 0.250 → 0.178 → 0.087 — and stops being distinguishable from zero at 5 km. At 2 km the phase itself flips to October, which is what a fitted harmonic does when there is no harmonic to fit. This is the same restriction, and the same failure mode, that took the Phase 4 Stage-1 seasonal term from OR 0.62 [0.44, 0.86] to 0.86 [0.69, 1.07].

With the balanced dyad panel *and* the distance restriction applied together:

| Restriction | Dyad-days | Dyads | Seasonal amplitude | Peak | Joint *p* |
|---|---:|---:|---:|---|---:|
| balanced panel, all distances | 26,258 | 45 | 0.348 | Aug | 0.045 |
| balanced panel and ≤ 5 km | 9,149 | 36 | 0.248 | Aug | 0.158 |
| balanced panel and ≤ 2 km | 3,628 | 31 | 0.107 | Oct | 0.494 |

Same decay, same crossing point.

**Verdict.** There is an annual cycle in *where groups are* — ranges shift with the vegetation cycle, which is why the unrestricted peak tracks NDVI to within a week. There is no annual cycle in *whether groups interact once co-located*. The NDVI phase match was a property of the geometry, not of the behaviour.

This is the Stage-1/Stage-2 dissociation once more: geometry has a season, the meeting decision does not. It belongs with Result 4, not as a separate ecological finding.

### The declining trend is denominator composition

The unrestricted trend of OR 0.748 [0.615, 0.909] looks like encounter probability falling across the record. On a **balanced panel of the 45 dyads observed throughout**, it is OR **1.100 [0.884, 1.368]** — null, and the sign flips.

The mechanism is visible in the denominator:

| | Value |
|---|---:|
| Dyads present throughout the window | 45 of 317 |
| Dyads in the denominator, first month → last month | 45 → **210** |
| Median closest approach, first half of window | 8,670 m |
| Median closest approach, second half | **11,437 m** |

As units are collared, the denominator fills with pairs that are further apart and less likely ever to meet. The apparent decline is that composition change, not a behavioural trend. **No secular trend in encounter survives.**

---

## 4. Axis B — no season, and the effort confound demonstrated cleanly

Multi-animal split in more than 25% of a group-week's timestamps, 1,148 group-weeks, 26 units:

| Term | With effort | Without effort |
|---|---|---|
| season_sin | 0.832 [0.541, 1.280] | 0.862 [0.543, 1.367] |
| season_cos | 0.973 [0.725, 1.305] | 1.066 [0.862, 1.318] |
| trend_months | 1.152 [0.873, 1.520] | **1.389 [1.022, 1.889]** |
| collars observed | **2.667 [1.055, 6.745]** | — |
| **Seasonal joint *p*** | **0.644** | 0.778 |

Two things to take from this. **No season** — joint *p* = 0.64, and the fitted peak wanders between September and October depending on specification, which is the signature of noise. And a textbook effort confound: the trend term is significant when effort is omitted (OR 1.389, CI excluding 1) and null once collars are in the model (OR 1.152). Collar coverage itself carries OR 2.667. That is the same lesson as §3, on a different axis.

Modularity in well-covered group-weeks (308 weeks, 12 units, ≥12 collars): **no season** (joint *p* = 0.615). The trend term is OR 1.440 [1.040, 1.993] and survives the effort control — but it **cannot be verified**. Only **2 of 12** well-covered units (Chartreuse, Lapis) span the window with 20 or more weeks, so a balanced-panel fit would be clustered on two units, where cluster-robust standard errors are badly anti-conservative. Report the modularity trend as *unverifiable*, not as confirmed. (A two-cluster fit does return OR 2.913 with a seasonal *p* of 0.0000; that is an artefact of two-cluster inference and is not reported as a result.)

---

## 5. Axis C — no season, no trend, on the cleanest denominator

Resolvable animal-nights, 71,827 rows, 22 units:

| Response | season_sin | season_cos | trend_months | Seasonal joint *p* |
|---|---|---|---|---:|
| Away from origin group tonight | 1.016 [0.750, 1.376] | 0.914 [0.729, 1.145] | 1.159 [0.873, 1.537] | **0.730** |
| **Excursion starts tonight** | 0.900 [0.683, 1.186] | 0.829 [0.584, 1.176] | 0.820 [0.578, 1.164] | **0.573** |
| Excursion starts, effort dropped | 0.899 [0.686, 1.177] | 0.827 [0.595, 1.149] | 0.814 [0.614, 1.079] | 0.523 |

Nothing, on either response, with or without the effort term. This is the most informative null of the three, because axis C has the cleanest denominator in the project — 82,205 animal-nights in four explicit states, with the sparse-collar artefact already removed. **Individual departure from the origin group has no annual rhythm** in this record.

---

## 6. What to report

| Claim | Evidence | Verdict |
|---|---|---|
| Observation effort rises 121-fold across the record | 2 → 242 collared animals; resolvability 0.00 → 0.95 | **report** — it licenses and forbids every trend |
| The environment is strongly seasonal | NDVI deviation amplitude 0.046, peak 18 July, *p* < 0.0001 | **report** as the reference cycle |
| Range overlap between groups is seasonal | encounter amplitude 0.285, peak 12 July, *p* = 0.010 unrestricted | **report as geometry** |
| Groups meet more in some seasons, once co-located | amplitude falls 0.285 → 0.087 and *p* 0.010 → 0.629 under distance restriction | **do not report** — fails the within-5 km test, as in Phase 4 |
| Encounter probability declines over the record | OR 0.748 [0.615, 0.909] unrestricted; **1.100 [0.884, 1.368]** on a balanced panel | **do not report** — denominator composition |
| Fission is seasonal | joint *p* = 0.644 with effort, 0.778 without | **report as tested and null** |
| Fission is increasing over the record | OR 1.389 [1.022, 1.889] without effort; **1.152 [0.873, 1.520]** with it | **do not report** — effort confound |
| Modularity is increasing | OR 1.440 [1.040, 1.993], survives effort | **unverifiable** — only 2 of 12 well-covered units span the window |
| Individual excursions are seasonal | away *p* = 0.730; onset *p* = 0.573 | **report as tested and null** |

### The general statement this supports

> The environment cycles strongly and the animals' ranges cycle with it, but boundary crossing itself has no season. What is seasonal is *where groups are*; what is not seasonal is *whether they cross a boundary once the opportunity exists* — at any of the three scales: between groups, within groups, or for individuals.

That sits directly alongside Result 4. Frequency and depth already had different predictors; now the same dissociation appears in time. Opportunity is seasonal; the crossing decision is not.

### Limits carried forward

1. **1.53–1.63 seasonal cycles.** Every seasonal estimate here rests on less than two annual cycles. The nulls are the more robust half of these results, because a null from a short window is a weak claim of absence, whereas a *positive* from a short window is easily an artefact — and the one positive found (axis A, unrestricted) turned out to be geometry.
2. **A single annual harmonic.** No second harmonic, no bimodal rainfall structure. Systems with two wet seasons need a second harmonic before an annual null is conclusive.
3. **Axis B still runs on the legacy source** filtered to 2025-01-01, so its window is shorter still and gap 7.3 applies to everything in §4.
4. **NDVI stops on 2026-04-29** while membership runs to 2026-07-22, so the environmental reference cycle and the behavioural windows are not exactly co-terminous.
5. **Effort is entered as a covariate, not as a design.** A properly effort-matched subsample across the whole record would be stronger than an effort term, and is feasible now that all three axes have denominators.
