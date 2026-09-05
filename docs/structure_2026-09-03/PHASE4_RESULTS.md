# Phase 4: the two-stage models

Prepared 4 September 2026. Continues [PHASE0_3_RESULTS.md](PHASE0_3_RESULTS.md).

Script: [phase4_two_stage_models.py](../../phase4_two_stage_models.py). Outputs: `outputs/general_structure_2026_09/phase4_two_stage_models/`.

---

## What was fitted, and what it is called

| | Method | Why |
|---|---|---|
| Stage 1 primary | Binomial GEE, independence working correlation, robust SEs clustered on dyad | 317 dyads support population-averaged inference |
| Stage 1 random effects | Mean-field **variational Bayes** mixed GLM (`BinomialBayesMixedGLM.fit_vb`) | recovers variance components; it is *not* MCMC and no MCMC diagnostics apply |
| Stage 1 validation | 5-fold cross-validation with **whole dyads held out** | tests whether the model transfers to an unseen dyad |
| Stage 2 primary | Binomial GLM on cross-edge counts with `logit(composition expectation)` as **offset**, plus dyad fixed effects | the intercept is then the mean logit deficit, and 12 dyads cannot support between-dyad inference |
| Stage 2 uncertainty | Cluster bootstrap over encounters within dyads, 2,000 draws, seed 20260904 | overdispersion is 12–27×, so model SEs are not usable |

`pymc` is installed but has no C++ compiler in this environment, so MCMC would run in Python fallback mode. No MCMC was attempted, and nothing here is described as a hierarchical Bayesian fit.

### Design rules enforced in code

- **Covariates strictly precede the outcome.** Same-day distance is endogenous during an encounter, so the within-dyad distance term is a **lag-1** deviation and the event-level distance term is measured over the 7 days before the event starts.
- **Between-dyad and within-dyad proximity are separate terms.** "These two ranges overlap" and "they were unusually close yesterday" are different claims.
- **One effect per group, applying at either end of the dyad** — a 26-column design where each column is 1 whichever side the group occupies, replacing alphabetical `group_a` / `group_b` terms.
- **Collar counts are named collar coverage** and treated as observation covariates, never as group size.
- **One retained sample per stage**, so the nested comparisons are valid.

---

## Stage 1 — whether two groups meet

Sample: 72,724 dyad-days, 317 dyads, 2,852 encounters (3.92%). Dropped from the 73,353 Phase 1 rows: 174 with insufficient support, 317 with no previous day, 138 whose previous observation was more than 3 days earlier.

Odds ratios per standard deviation, full model:

| Term | OR | 95% CI |
|---|---:|---|
| Between-dyad typical distance | 0.143 | 0.115 – 0.177 |
| Within-dyad distance yesterday | 0.346 | 0.310 – 0.386 |
| Prior encounters (strictly preceding) | **1.762** | **1.486 – 2.088** |
| Collar coverage, smaller unit | 0.887 | 0.666 – 1.181 |
| Collar coverage, total | 1.326 | 0.974 – 1.805 |
| Shared observed hours | 0.934 | 0.891 – 0.980 |
| Season (sin) | 0.616 | 0.439 – 0.863 |
| Season (cos) | 0.824 | 0.655 – 1.037 |

### Three results

**1. History predicts encounters beyond geometry.** Prior encounters carry OR 1.76 after both distance components are in the model. The strong test is the restricted sample: among dyad-days where the groups were already within 5 km the previous day (12,156 rows, 103 dyads, 20.1% encounter rate), the effect holds at **OR 1.57 [1.27, 1.93]**. Repeated contact between particular group pairs is therefore not a restatement of spatial overlap — which is the first real support for Result 5 in the manuscript plan.

**2. Encounter propensity is a dyad property, not a group property.** Variance components from the mixed model, across three prior scales:

| Prior scale | Dyad effect SD | Shared group effect SD |
|---:|---:|---:|
| 0.5 | 1.516 | 0.120 |
| 1.0 | 1.578 | 0.044 |
| 2.0 | 1.589 | 0.022 |

The dyad component is large and stable across priors. The group component is small and shrinks as the prior tightens, which means it is not identified — there is no evidence that some groups are simply more sociable. Whatever drives repeated encounters lives in the pair.

**3. Collar coverage does not manufacture encounters.** Both coverage terms have intervals spanning 1, in the full sample and in the restricted sample. This is the single most important robustness check in the project, and it passes.

### Transfer to unseen dyads

| | Value |
|---|---|
| Out-of-sample AUC, whole dyads held out | **0.880** |
| Calibration slope | 0.540 |
| Calibration intercept | −0.963 |
| Mean predicted / observed rate | 0.0365 / 0.0392 |

Discrimination transfers well; calibration does not. A slope of 0.54 means predictions are too extreme for a dyad the model has never seen, which is what should happen when a large dyad-level variance component cannot be estimated for a new pair. Report it as a limit rather than smoothing it over.

Seasonality is the one effect that does **not** survive restriction: season (sin) is 0.62 [0.44, 0.86] in the full sample but 0.86 [0.69, 1.07] within 5 km. The apparent seasonal signal is mostly between-dyad range structure, not a within-opportunity effect.

---

## Stage 2 — how much they mix, given an encounter

Sample: 738 discrete encounters across 12 dyads (746 including sustained associations). Observed cross-edge fraction 0.0589 against an edge-weighted composition expectation of 0.4504.

| Term | Estimate | 95% bootstrap CI |
|---|---:|---|
| Intercept = mean logit deficit | **−3.073** | **−3.433 – −2.765** |
| log structural span | **+0.624** | **+0.302 – +0.828** |
| Pre-event distance (7 days before) | +0.005 | −0.080 – +0.083 |
| log supported exposure | −0.095 | −0.329 – +0.188 |
| Collar coverage | −0.015 | −0.182 – +0.144 |
| Prior encounters | −0.256 | −0.355 – −0.174 |

An intercept of −3.07 means the odds of a cross-group link during a discrete encounter are about **4.6% of the composition expectation**. That is the Phase 3 gradient restated as a model, and it is the paper's central number.

### One robust predictor, and two coefficients that did not survive

**Duration is the only robust predictor of mixing.** log span is +0.62 [0.30, 0.83] and holds in every specification: without dyad effects (+0.61), with calendar time (+0.54). Longer encounters mix more, within dyads.

**Prior encounters does not survive calendar time.** Within a dyad, prior encounters largely indexes time (r = 0.285 with elapsed days). Adding an explicit calendar-time term moves it from −0.256 [−0.355, −0.174] to **−0.082 [−0.289, +0.060]**, and the time term takes −0.246 [−0.423, +0.008]. So there is no support for mixing declining with encounter history; what little signal exists is a weak temporal trend. Do not report the history effect at Stage 2.

**Distance and coverage effects at Stage 2 are between-dyad confounds.** Both are null within dyad and both look real when dyad effects are dropped:

| Term | With dyad effects | Without dyad effects |
|---|---:|---:|
| Pre-event distance | +0.005 [−0.080, +0.083] | −0.123 [−0.216, −0.018] |
| Collar coverage | −0.015 [−0.182, +0.144] | −0.389 [−0.487, −0.295] |

Overdispersion also halves when dyad effects are included (Pearson χ²/df 12.8 versus 27.1). This is the clearest available demonstration that a mixing model without dyad structure reports dyad identity as if it were a covariate.

### The sustained-association contrast, confirmed with a caveat

Within dyad, and unconditioned:

| Model | `is_sustained` | 95% bootstrap CI |
|---|---:|---|
| A2 — stratum contrast alone | **+1.376** | **+0.338 – +1.543** |
| A1 — with exposure and other terms | −0.109 | −0.789 – +0.462 |

A2 is the honest estimate. A sustained association mixes about 1.4 logit units closer to composition expectation than the *same dyads'* discrete encounters. A1's null is a second instance of the collinearity trap: sustained events are by construction the highest-exposure ones, so `log supported exposure` absorbs the stratum, exactly as `log span` did before it. Both cases are recorded in the script's collinearity audit.

Two caveats stand. The contrast rests on **16 sustained events in 3 dyads**, and the bootstrap interval is strongly asymmetric — its lower tail comes from resamples in which Copper–Lilac's sustained events are under-represented. And the model cannot separate the stratum from duration or exposure, because the stratum is defined by duration. The Phase 3 per-dyad description remains the primary evidence; this is a confirmation, not an independent test.

---

## What Phase 4 establishes for the manuscript

| Result | Status after Phase 4 |
|---|---|
| 1 — identity persists while membership varies | supported, descriptive |
| 2 — reconfiguration is routine and multi-scale | supported, pending cohort re-derivation |
| 3 — shared space rarely means mixing | supported and now modelled: mixing is 4.6% of expectation |
| 4 — encounter and mixing have different predictors | **supported.** Geometry and history drive whether groups meet; neither predicts how much they mix. Duration is the only robust mixing predictor. |
| 5 — permeability structures the population | **first real support.** History predicts encounters beyond geometry, and the effect is a dyad property with no identifiable group component. |

Result 4 is now a positive finding rather than an absence: the two stages genuinely have different predictors, which is what justified separating them.

## Limits carried forward

- Stage 2 rests on 12 dyads with fine-scale support, out of 68 with structural encounters. Between-dyad statements are unavailable; every Stage 2 result above is within-dyad.
- Calibration does not transfer to unseen dyads (slope 0.54).
- No ecological covariate is in either model. NDVI, rainfall and shared concentrated resources belong to Phase 5, and the seasonal terms here are placeholders that already fail the within-5 km test.
- Demographic group size is still absent. Nothing above may be read as a group-size effect.
- The mixed model is variational Bayes. Its variance components are informative and prior-sensitivity is reported, but they are a mean-field approximation.

## Phase 5 entry conditions

Both stages now have a frozen sample, a specification that survives its own sensitivity checks, and a documented set of dead coefficients. Phase 5 can add the environmental covariates to exactly these samples and to a matched fission model, with the requirement that any ecological effect must survive the within-5 km restriction that already killed the seasonal term.
