# Stable groups, fluid boundaries

## General framing, organised by axis of fluidity and by what population coverage supports

Prepared 4 September 2026. Successor to [PAPER_STRUCTURE_GENERAL.md](../structure_2026-09-03/PAPER_STRUCTURE_GENERAL.md), which is retained for its phase plan.

The earlier framing was correct about the central problem — a case study had to become a distribution — but it made **one** axis of fluidity carry the paper: between-group encounter and mixing. Three other axes were listed and then demoted (within-group fission, isolation, dispersal). This document restates the argument with four co-equal axes, and makes the coverage argument explicit: for each axis, what a 350-animal, 26-group, 29-month cohort does and does not license.

Every number is read from saved outputs in this project. No model was rerun.

> **Amended 4 September 2026 — three axes, not four.** See [INDIVIDUAL_AXIS_MERGE.md](INDIVIDUAL_AXIS_MERGE.md). Axes C (isolation) and D (dispersal) are one process: an animal away from its origin group is either alone or with another group, and those are graded depths of a single individual-level trip. Held apart they gave a 285-animal isolation inventory that was ~76% observation artefact plus a 17-animal dispersal case series with no denominator. Merged, and filtered by a **located-reference-cluster** rule that removes 100% of the 1- and 2-collar contexts, axis C becomes **338 excursions, 91 animals, 18 groups over 82,205 animal-nights in four explicit states** — the second axis able to carry a two-stage model. Read §4 (axes C and D), §5's table row, and §7.1–7.4 as superseded by that document.
>
> **Amended again — axis B is now analysed.** See [AXIS_B_FISSION_RESULTS.md](AXIS_B_FISSION_RESULTS.md). Three corrections to §4 (axis B): (1) modularity *is* reportable **above about twelve collars** — at matched coverage Lilac is modular in 53.8% of weeks and Chartreuse 45.0%, against Purple, RubyRunners and Magenta at 0.0%; below that threshold the estimator returns one cluster regardless, so the "19 of 20 groups at zero" figure is a measurement artefact, not biology. (2) Modularity is a **state**, not a group property — between-unit ICC **0.148**, so 85% of the variance is within group over time, and it is sticky (lag-1 *r* 0.33–0.41 in the two groups that enter it). This gives Result 5 its third leg. (3) The split-composition signal is a **bond, not a pre-split cue**: a dyad's general co-sitting rate, measured with a ±7-day holdout around the split, predicts sides at least as well as the immediately preceding proximity (AUC 0.574 vs 0.570, paired difference +0.004 [−0.009, +0.018]). And **NDVI does not predict fission** — within-group deviation null in six of eight specifications — which retires the old Phase 5 "scarcity fragments" prediction on the fission side.

>
> **Amended again — seasonality and trends.** See [SEASONALITY_AND_TRENDS.md](SEASONALITY_AND_TRENDS.md). Observation effort rises **121-fold** across the record (2 → 242 collared animals, resolvability 0.00 → 0.95), so no raw trend is interpretable; the well-observed window starts 2024-12-01 and holds only **1.6 seasonal cycles**. NDVI itself is strongly seasonal (amplitude 0.046, peak 18 July, *p* < 0.0001), so nulls cannot be blamed on an absent environmental cycle. Only axis A shows a season — amplitude 0.285 logit peaking 12 July, within six days of NDVI — and it **decays to nothing under distance restriction** (0.285 → 0.087; *p* 0.010 → 0.629), the same test that killed the Phase 4 seasonal term. The season is in *where groups are*, not in whether they interact once co-located. Axis A's apparent secular decline (OR 0.748) is null on a balanced dyad panel (OR 1.100) — denominator composition. Axes B and C show no season at all. **Opportunity is seasonal; the crossing decision is not.**

>
> **Amended again — cohort closure, two corrections, and rare patterns.** See [COHORT_CLOSURE_AND_RARE_PATTERNS.md](COHORT_CLOSURE_AND_RARE_PATTERNS.md). Gap 7.3 is closed for axis B's weekly metrics, rebuilt on the frozen export — and the check caught two errors. (1) **The split-detection confound was overstated**: Spearman with collar count is 0.163 on the frozen cohort, not 0.489. (2) **Chartreuse's modular rate does not replicate** (45.0% of 20 weeks → 14.7% of 34); Lilac's does (55.3% vs 53.8%), as does the ICC (0.169 vs 0.148) and the stickiness. (3) **LapisSplinter is not a documented fission**: six of its seven animals carry `origin_group = LapisSplinter` and were first observed on the day they were collared, so the split predates tracking — the fission-precedes-permeability question is *untestable* here, not negative. (4) **Seasonality was tested on occurrence only**; depth and duration are now tested too (seven responses, joint *p* 0.234–0.826) and are also null, with one open question — axis A's mixing depth reaches *p* = 0.084 once the duration term is dropped. (5) **Editorial rule changed**: rare patterns are reported *as variation* with their frequency stated, not suppressed. Measurement failures are still withheld. The two are different, and 26 units over 29 months is a subset — the tail's existence is itself the result.
>
> **Amended 5 September 2026 - the encounter set contains dispersal.** A `detected_encounter` fires whenever two units share an association cluster, and it does not ask how many animals each unit contributed. **639 of the 2,867 detected-encounter dyad-days (22.3%) have exactly one animal on the smaller side**, spread over 43 of the 68 dyads, and 9 dyads are *only ever* single-animal. Those are not two groups meeting; they are one animal visiting - the axis C process, counted on axis A. The geometry says so independently: on those days the two group centroids sit a median **2,460 m** apart, against **44 m** on a genuine group encounter. Requiring **two or more animals on both sides** leaves **2,228 dyad-days across 59 dyads** (77.7% of the days, 87% of the dyads). Wherever axis A and axis C are compared - and the level-2 figure overlays them on one node layout - the restricted set is the one to use, or the comparison is partly circular. Under it, **44 dyads meet without ever exchanging an animal, 1 exchanges without a recorded group encounter, and 15 do both**. The unrestricted 68-dyad figure stands wherever the question is *contact of any kind*, including S5.3's 16-of-1,705 sustained-association result, which is measured on encounter events rather than dyad-days; the two counts are not interchangeable and should be labelled wherever they appear.
>
> **Amended again 5 September 2026 - hourly resolution, the crossing ledger, and how far the two networks are independent.** See `derive_level2_inputs.py`.
> 
> *(i) An encounter is an hourly fact, and the dyad-day table was hiding that.* Rebuilt from the narrow export as contiguous hours in which both units had two or more animals in one association cluster (gaps of up to 2 h bridged, because 02:00 is a known coverage hole), the population has **1,445 group-encounter bouts across 58 dyads**. Median **14 h**; **1,375 of 1,445 (95%) end within a day** and 565 (39%) within six hours. The dyad-day table floors every encounter at 24 h, so it was overstating the typical encounter by roughly a factor of two and hiding the entire short mode. The tail is real but thin: the longest bout is **8,372 h** (Copper-Lilac, a pair that effectively stopped separating).
> 
> *(ii) Repeat single-animal crossings on one dyad are separate dispersals, not one merge and not one edge.* A dyad can carry group encounters and single-animal crossings at once, and can carry several crossings months apart. The per-dyad ledger (`dyad_ledger.csv`) now records both sides. On the **audited** axis-C set - dominant nightly rule, named destination - there are **182 crossing episodes by 38 animals over 19 dyads, and 14 of those 19 dyads carry more than one episode**; 7 of them involve more than one animal. Maroon-Sapphire alone carries 44 episodes by 5 animals. On the **broad** set - every resolved away-night run, before the excursion rule - it is 882 episodes by 146 animals over 51 dyads, with 38 dyads repeating. The broad set is not filtered by the located-reference-cluster rule and should not be substituted for the audited one; it is reported here only to show that the audited set is a lower bound.
> 
> *(iii) The two networks are mostly, but not entirely, independent evidence.* Of the 58 encounter dyads and 19 audited crossing dyads, **16 carry both**, 42 carry encounters only, and 3 carry crossings only - among them **Emerald-PhantomWest, which has zero group-encounter bout-hours and 4 audited crossings**; its 12 unrestricted "encounter" dyad-days were all single-animal, all one disperser, and all removed by the two-per-side rule. Asking how much of the group-encounter record sits inside a crossing's stay: **158 of 848 bouts (18.6%)** on the shared dyads, and **18.2% of encounter hours once Copper-Lilac is set aside**. Unfiltered by hours the figure is 54.6%, but that is one 8,372-hour bout dominating the sum, not a general confounding. **Three of the 16 shared dyads are not clean** - Copper-Lilac (84% of encounter hours inside a stay) and Maroon-Sapphire (91%) most of all, and both are heavy edges in the level-2 figure. Those two should not be cited as evidence that encounters and transfers are separate channels; the other 13 can be.
>
> **Amended again 5 September 2026 - all three states at hourly resolution, and a resolution floor.** See `derive_level2_inputs.py`. The three axes were each measured at whatever resolution their pipeline happened to use - dyad-days for A, weeks for B, nights for C - which made them incomparable and floored two of them above the scale on which they actually happen. All three are now also derived hourly from the same export.
> 
> **Axis B has an hourly face and it is short.** A unit with 12+ animals present occupying two or more clusters of two or more animals gives **914 split bouts across 12 units, median 1 h, 96% inside a day**. This does not contradict the weekly modularity: it explains it. *A group can be modular for seven consecutive weeks without ever being apart for a whole day.* Modularity is a property of the week's association network, not a period of separation, and the two should never be described in the same words.
> 
> **Axis C's hourly curve is not a dispersal duration.** Contiguous hours of non-origin social context gives 3,583 bouts for the audited 91-animal cohort, median 2 h - but a months-long absence is cut every time the hourly context returns to the origin, so the longest hourly run is 5,478 h where the nightly excursion reaches 500 nights. **The nightly dominant rule stays the measure of how long leaving lasts;** the hourly view only adds the bottom of the distribution. Unrestricted, the hourly excursion state fires for **319 animals**, the sparse-collar population the located-reference-cluster rule exists to exclude - a reminder that the rule is doing real work.
> 
> **A resolution floor applies to all of it.** **50% of split bouts and 48% of excursion bouts are exactly one hour, against 16% of encounter bouts.** One animal's cluster membership flickers where a whole group's co-occurrence does not, so anything read below about two hours on the individual and within-group axes is assignment noise. The level-2 figure shades that strip rather than filtering it away.
> 
> **The three axes are now one figure.** `build_figure3.py` composes them as `docs/framing_2026-09-04/figure3.html`, a two-by-two: (a) one group's weekly modularity as a line through time with two of its own weeks inset as association networks, (b) group meetings and animal transfers **overlaid on one node layout**, (c) all four concentration curves on shared axes, (d) all three states' durations at hourly resolution. The overlay in (b) is the comparison worth making, with one caveat stated in the caption: 16 dyads carry both chord types, and on Copper-Lilac and Maroon-Sapphire a doubled chord is one finding rather than two. `build_level2_figure.py` and `build_modularity_figure.py` remain as the sources of the individual panels.
> 
> **Modularity now has its own figure.** `build_modularity_figure.py` - one group, Lilac, week by week, and four of its own weeks drawn as association networks at Q = 0.45, 0.28, 0.14 and 0.00. The pictures make two points a score cannot: **the middle of the scale is a heavily bridged community structure, not a group in two halves**, and at Q = 0.28 the partition is into *three* communities, not two. One unit produces the whole range inside nine months, which is the case for reading modularity as a state rather than a trait - but Lilac is not a typical unit. Across the 7 well-covered units, 5 reach modularity in some week and 2 never leave zero, and **Lilac and Chartreuse alone hold 37 of the 42 modular weeks (88%)**. This range is what a unit *can* do, not what most units do.
>
> **Amended again 5 September 2026 - the clustering rule, the GPS source, and a second pipeline.** See [CLUSTERING_AND_SOURCE_DECISIONS.md](CLUSTERING_AND_SOURCE_DECISIONS.md). Six method questions were opened and closed. **(1)** The adaptive clip is not what makes the rule permissive - the median edge threshold actually used is **198 m**, only 0.3% of dyad-thresholds exceed 600 m, and tightening 900 m to 600 m changes **4 dyad-hours in ~250,000** (ARI 1.00). **(2)** HDBSCAN is ruled out: it finds 25 clusters where the adaptive rule finds 10, with silhouette 0.36 against 0.94, stability 0.56 and persistence 0.43 - it reads the density gradient inside a foraging group as a boundary. **(3)** The clustering *rule* barely matters: this project's adaptive rule and the parallel pipeline's DBSCAN at 300 m give **identical partitions in the median hour (ARI 1.00)**. The *input* does matter - collars fix every 2 minutes and there is a median of **30 shared 2-minute bins per dyad-hour** that an hourly median discards; feeding the same rule the fine-scale matrix raises silhouette 0.941 to 0.955 and persistence 0.978 to 1.000. Banked for the next rebuild, not applied, because it invalidates every `association_event_id`. **(4) The real blocker is the GPS vintage:** three sources are in circulation with no matching key spaces, and **the frozen export runs to 2026-07-22, later than any raw file available locally**, so it was built from a file we do not have; 25AA07_4S5T's last record differs three ways (2026-03-16, -03-18, -03-23). **(5)** Crosswalking the two pipelines over **1,075,208 shared animal-hours** gives 84.7% concordance - and **independently confirms the single-animal correction**: their `DISPERSAL_WITH_OTHER_GROUP` state has the animal alone on its own side in **100%** of its 13,417 rows, 95.6% of them this project's `mixed_without_origin_unit`. The residual 15% is *not* a labelling artefact (origin labels agree on 99.9% of rows) and, since the clustering agrees at ARI 1.00, must come from the support / carry-forward / sparse-veto layer around it. **(6)** The three-way partial merge - pure A, pure B and a mixed A+B in one hour - occurs in **534 hours, 166 events, 19 pairs**, but **the mixed side is a single animal in 56.1% of them** and three dyads carry 86% of the total, two involving LapisSplinter. Reportable as a restricted observation, not as a population-level phenomenon.

---

## 1. The thesis

A group is not a container. It is a **baseline that animals return to**, and the interesting quantity is not whether the baseline exists but how, how often, and how deeply animals leave it.

On 82.7% of 81,695 animal-nights an animal sleeps with its origin group. That is the stability. On the other 17.3% it is somewhere else, and on 33.2% of nights the *configuration* it sleeps in is not its whole group — it is merged with another group, split from part of its own, or alone. Stability and fluidity are properties of the same population at the same time, measured on the same rows.

The paper's job is to turn that pair of facts into a general, quantified statement about **which routes across a group boundary are used, how far each one goes, and what a population-scale dataset can actually establish about each.**

---

## 2. Four routes across the same boundary

A group boundary can become fluid in exactly four ways, and the frozen cohort measures all four on the same rows.

| Axis | The route | What crosses | Nightly share |
|---|---|---|---:|
| **A. Between-group interaction and mixture** | Two groups occupy the same space; individuals may or may not mix within it | whole groups, then individuals inside the merged mass | 25.1% between-group merge |
| **B. Within-group fission and modularity** | A group divides into subunits that may or may not have persistent membership | internal structure | 4.6% within-group split |
| **C. Isolation** | An animal is with no group at all | one animal, temporarily out of the system | 3.5% isolated / low support |
| **D. Dispersal** | An animal leaves one group and joins another | one animal, permanently | (subset of the 14.0% with another group) |

These are not four separate studies. They are four measurements of one quantity — **boundary permeability** — taken at different grain, and they interact:

- LapisSplinter is a unit **created by** axis B (a Lapis fission product). It then supplies the *most* permeable dyad on axis A: Lapis–LapisSplinter has the smallest deficit of the twelve discrete-encounter dyads (−0.295), LapisSplinter–Magenta the fourth smallest (−0.356), and LapisSplinter appears in two of the five dyads that ever entered sustained association. It also receives an axis-D transfer (animal `24AC02_3C4D`, origin Lapis, 655 observed hours in LapisSplinter).
- Sapphire is the only group with substantial median within-group modularity (0.278, axis B) and it is one of three dyads with a sustained association on axis A (Maroon–Sapphire).
- Axis C is the residual of the other three: 5,828 single-animal separation events (301 of 350 animals) are the raw material from which both a brief absence and a permanent dispersal are drawn. The classifier that separates them is the weakest link in the project.

**The framing claim.** Fission products and internally divided groups are the most permeable units in the population. That is one sentence, it is cross-axis, and it is currently an *observation from saved tables* — not a tested claim. Section 7 says what would test it.

---

## 3. What population-level coverage buys

This is the general argument, and it is worth stating as a method contribution rather than burying it.

Population coverage does three specific things that no case study can do, and each has already changed a conclusion in this project.

### 3.1 It converts a case into a distribution with a documented tail

Before: Copper–Lilac fused, therefore group boundaries are permeable. After: twelve measurable dyads all mix far below composition expectation (dyad-weighted −0.405 [−0.442, −0.367], 12 of 12 negative), and Copper–Lilac's *encounters* rank tenth of twelve. The general result is stronger than the case ever was, and the case became the labelled endpoint.

### 3.2 It separates an entity from a state

This is the most valuable thing coverage does, and it generalises across all four axes.

With one dyad you cannot tell whether an extreme value is a property of that pair or of a condition the pair happened to be in. With twelve you can. Copper–Lilac's discrete encounters give −0.339 (inside the pack); its sustained association gives −0.033 (near expectation). **The outlier is a state, not a dyad** — and five dyads entered that state.

The same test is available, and not yet run, on the other three axes:

| Axis | The entity-versus-state question | Available? |
|---|---|---|
| A | Is high mixing a dyad property or a state? | **Answered: state.** Sustained association, 5 dyads |
| B | Is high modularity a group property or a state some groups enter? | Sample exists (4,570 events, 24 groups); test not run |
| C | Is isolation an animal property or a state most animals pass through? | Strongly suggestive: 285 of 350 animals have ≥1 isolated event, but the per-animal count is median 4, p90 37, **max 364**. Not yet decomposed |
| D | Is "disperser" a type or a phase? | Suggestive: 71 segments in 17 animals, of which 46 permanent and 25 returned — animals have *multiple* segments with *different* outcomes. Not yet modelled |

If the answer is "state" on all four axes, the paper has a single general finding rather than four descriptive chapters.

### 3.3 It makes negatives real

Coverage is what allows an *observed absence* to be distinguished from a *missing row*. The Phase 1 opportunity table has 70,312 observed-no-encounter dyad-days with a median 24 shared observed hours and a median closest approach of 10.8 km, versus 2,867 detected-encounter days with a median closest approach of 66 m. The saved pipeline had zero-filled 2,406 real encounters across 53 dyads as non-events, and had no exposure measure behind any zero.

The same problem is unfixed on axes B, C and D: there is no opportunity table for a group failing to split, an animal failing to be alone, or an animal failing to leave. Every rate on those axes is currently an event count without a denominator.

---

## 4. The coverage ledger

The four axes are not equally supported. Scaling each claim to its coverage is the single most important editorial decision in the paper.

| | **A. Between-group** | **B. Within-group** | **C. Isolation** | **D. Dispersal** |
|---|---|---|---|---|
| Units | 26 social units, 317 co-observed dyads | 24 groups with split events; 20 with weekly networks | 350 animals | 350 animals |
| Units with the phenomenon | 68 dyads with encounters | 24 groups | 285 animals (81.4%) | 17 animals (4.9%) |
| Events | 1,705 Stage-1 encounters | 4,570 splits | 4,293 isolated; 5,828 single-animal separations | 71 segments / 28 events (unreconciled) |
| Denominator | **exists** — 73,353 dyad-days, four explicit states | absent | absent | absent |
| Depth measurement | 5 m / 2 m cross-group contact, 12 dyads, 751 events | weekly modularity + split fraction, 20 groups | duration only | centrality trajectory, 60 segments |
| **Coverage-limited claim level** | population rate + within-dyad model | population inventory; **modularity is coverage-bound** | population inventory; per-animal decomposition available | **case series** |

### Axis A — best supported, and the model is fitted

Two-stage separation works. Stage 1 (72,724 dyad-days, 2,852 encounters, 3.92%): whether groups meet is driven by geometry and by *history* — prior encounters carry OR 1.762 [1.486, 2.088] after both between-dyad and lag-1 within-dyad distance are in the model, and OR 1.57 [1.27, 1.93] in the restricted sample where the groups were already within 5 km yesterday. Encounter propensity is a **dyad** property: dyad-effect SD 1.58 versus shared-group-effect SD 0.044, stable across three prior scales. Out-of-sample AUC with whole dyads held out is 0.880; the calibration slope is 0.540, which is a stated limit.

Stage 2 (738 discrete encounters, 12 dyads): the odds of a cross-group link during a discrete encounter are **4.6% of the composition expectation** (intercept −3.073 [−3.433, −2.765]). Duration is the only robust predictor (log span +0.624 [+0.302, +0.828]). Nothing that predicts meeting predicts mixing.

**What coverage forbids:** any between-dyad statement at Stage 2. Twelve of 68 structural dyads have fine-scale support. Every Stage 2 result is within-dyad, and the paper must say so in the results text, not only in limitations.

### Axis B — broad event coverage, and modularity is coverage-bound

The event side is genuinely general: 4,570 split events across 24 groups, median duration 2 h, p90 14 h, median maximum split fraction 0.25, p90 0.50. Splits are overwhelmingly binary (2 clusters in 3,829 events, 83.8%; 3 in 628; ≥4 in 113) and shallow — the typical split detaches a quarter of the observed animals for two hours.

The modularity side is where coverage bites, and this must be stated as a measurement limit rather than a biological result. Median tracked animals per group per week is 8, against nominal group sizes of 23–108. Coverage runs from **7.4%** (Jade: 8 of 108) to 25.8% (Emerald: 8 of 31), median ≈ 17%.

Median weekly within-group modularity is **exactly 0 in 19 of 20 groups**. That is not evidence that baboon groups have no internal community structure. It is what a modularity estimator returns when it is handed 5–18 nodes drawn from a group of 50: a single cohesive cluster. The one exception (Sapphire, 0.278) rests on 9 weeks and 4 tracked animals and should not be reported as the group that is different.

The honest positive result on this axis is the *anticipation* result, because it uses dyads rather than communities: across 466 events in 15 groups, animals that co-sat within 5 m before a split are more likely to end up on the same side of it — mean event AUC 0.568 [0.558, 0.578], permutation p = 0.002, mean same-minus-different-side 5 m contact +0.0092 [0.0074, 0.0111]. It is a small, real, replicated effect, and the mean event ARI of 0.100 says correctly that split composition is *mostly* not predicted. The matched Chartreuse–Purple modularity comparison is a clean null (difference 0.0015 [−0.025, 0.023]) and should be reported as one.

**What coverage forbids:** any claim about within-group community structure, subgroup persistence, or modularity differences between groups. Present the event inventory and the pairwise anticipation result; present modularity only as a coverage-limited descriptive with the collar fraction on the same axis.

### Axis C — the most general axis in the dataset, and the least analysed

4,293 isolated events across **285 of 350 animals and all 26 origin groups**. This is broader coverage than any other axis, including axis A. And the shape is informative: median duration 1 h, p90 13 h, max 814 h. Being briefly alone is near-universal and near-instantaneous; being alone for a month happens.

The per-animal distribution is the interesting one: median 4 events, p90 37, **max 364**. That range spans two orders of magnitude and is exactly the entity-versus-state question from §3.2. It has three candidate explanations and the data can separate them:

1. a genuine behavioural type (some animals are peripheral);
2. a state (some animals pass through a peripheral phase — pre-dispersal, illness, injury);
3. an observation artefact (an animal whose group-mates are poorly collared appears isolated).

Explanation 3 must be excluded first, and the cohort supports it: isolation counts can be regressed on the same collar-coverage variables that Stage 1 already showed do **not** manufacture encounters. That check is the entry condition for anything else on this axis.

**What coverage allows that is not yet used:** isolation is the one axis where 350 animals is a large sample. A per-animal decomposition — how much of the variance is animal identity, how much is a time-varying state, how much is coverage — is feasible now on the frozen cohort and would make this a result rather than a footnote.

### Axis D — a case series, and it should be labelled one

Two irreconcilable inventories are in circulation: 28 disperser events in 9 animals (with fates: 16 stayed, 10 censored, 2 returned) and 71 segments in 17 animals (46 permanent, 25 returned; 60 with centrality). The discrepancy is not a bug to be split — the two products use different definitions, different sources, and different censoring rules, and one of them must be chosen and documented.

Even after reconciliation, 17 animals with repeated segments and heavy censoring is a case series. The population-scale contribution of this axis is not the dispersal rate; it is the demonstration that **outcome is not a fixed property of the animal**: the same animals contribute both permanent-dispersal and returned-to-origin segments.

**What coverage forbids:** dispersal rates, sex or age effects, and any statement about which animals disperse. State the reconciliation, report the trajectories as a labelled case series, and let axis C carry the individual-level generality.

---

## 5. The general pattern the four axes share

Read across, the four axes tell one story with the same shape at every scale.

> **Boundary crossing is common, brief, partial, and rarely deep — at every scale it is measured.**

| Axis | Common | Brief | Partial | Rarely deep |
|---|---|---|---|---|
| A | 68 of 317 dyads; 25.1% of animal-nights in a merged configuration | median structural span 14 h | median active contact 0.17 h within a 14 h span | mixing at 4.6% of composition expectation; 16 of 1,705 events sustained |
| B | 4,570 events, 24 groups | median 2 h | median split fraction 0.25; 83.8% binary | split composition barely anticipated (AUC 0.568, ARI 0.100) |
| C | 285 of 350 animals | median 1 h | one animal at a time | 3.5% of animal-nights, and the long tail is rare |
| D | 5,828 separation events feed it | most reverse | one animal at a time | 46 permanent segments in 17 animals |

Three implications follow, and they are the paper's general contributions.

**5.1 The unit of analysis is the state, not the entity.** Axis A already demonstrates it: a dyad's extreme mixing value belongs to a condition (sustained association) that five dyads entered, not to the dyad. The same decomposition is what axes B, C and D need, and running it on all four is what turns a descriptive survey into a single finding.

**5.2 Depth and frequency are separately regulated.** On axis A this is a positive, fitted result: geometry and history predict *whether* groups meet; neither predicts *how much* they mix, and duration is the only thing that does. That dissociation is the reason the two-stage structure exists, and the corresponding dissociation on axis B (splits happen often; who ends up with whom is weakly determined, ARI 0.100) is the same pattern at a smaller scale. Two axes, two grains, same dissociation — that is a general claim about fission–fusion regulation and it is the strongest broad-audience angle in the project.

**5.3 Depth requires duration, and duration is rare.** The only route to near-complete mixing observed anywhere in the population is sustained association: +1.376 logit units [+0.338, +1.543] closer to expectation than the same dyads' discrete encounters. It occurred in 16 of 1,705 encounters, in 5 of 68 dyads, and only one dyad carried it to near-completion. Permeability is therefore not a property that groups have or lack; it is a threshold that almost nothing crosses.

---

## 6. Result architecture

Six results. Each carries its coverage in its own statement.

| # | Result | Coverage that carries it | Status |
|---|---|---|---|
| 1 | Group identity persists while configuration varies constantly | 350 animals, 26 groups, 81,695 animal-nights; two orthogonal state distributions | frozen; needs coverage/ambiguity weighting |
| 2 | Boundary crossing happens by four routes, and all four are common, brief and partial | 1,705 encounters, 4,570 splits, 4,293 isolations, 71 dispersal segments, one cohort | assembled here; needs cohort re-derivation for B/C/D |
| 3 | Groups share space far more often than they mix, and mixing runs at 4.6% of composition expectation | 12 dyads, 738 encounters, 12 of 12 below expectation | supported and modelled |
| 4 | Frequency and depth have different predictors, at two scales | Stage 1 (317 dyads) vs Stage 2 (12 dyads); split occurrence vs split composition | supported on axis A; axis B parallel needs stating |
| 5 | Extreme permeability is a state, not an entity | 5 dyads in sustained association; contrast +1.376 | supported on axis A; **extend to B, C, D** |
| 6 | Repeated crossing structures the population, and the structure is dyadic | prior-encounter OR 1.57 within 5 km; dyad SD 1.58 vs group SD 0.044 | first real support; needs time-respecting connectivity |

Result 5 is the new spine. It is what the four-axis framing buys, and it is the result that no single-axis paper could state.

Supporting, not main: Copper–Lilac as the labelled endpoint of Result 3; Chartreuse–Purple as a reported null; individual dispersal centrality as a case series; activity controls (VeDBA on 101/108 events, GPS movement on 97/108; r ≈ 0.28 and 0.21 against saved 5 m integration).

### Figure plan

1. **Cohort and coverage.** Animals × time ribbon, and — critically — collar fraction per group against nominal size. This figure is what licenses and forbids the rest, so it goes first and is not a supplement.
2. **Four routes, one cohort.** Four panels sharing a duration axis: encounter span, split duration, isolation duration, dispersal segment length, each with its event count and unit count. The "common, brief, partial" pattern should be visible in one glance.
3. **The mixing gradient.** Per-dyad observed-minus-expected with support counts and intervals; radius as a second dimension; Copper–Lilac labelled at the tail.
4. **Frequency versus depth.** Stage 1 predicted encounter probability and Stage 2 mixing conditional on encounter, side by side, with the axis-B analogue (split occurrence versus split-composition AUC) as a third panel.
5. **State, not entity.** The entity-versus-state decomposition on all four axes: within-entity variance against between-entity variance. This is Result 5 and it is the figure that does not exist yet.
6. **Population network.** Dyadic encounter network against an opportunity-adjusted null.

---

## 7. What this framing requires that does not exist

The four-axis structure exposes gaps the single-axis structure hid. In priority order.

**7.1 Opportunity tables for axes B, C and D.** Axis A has 73,353 dyad-days with four explicit states, and it is the reason axis A has a model. Axes B, C and D have event counts with no denominator. Each needs its own: group-weeks with observation support for splits; animal-days with a resolvable group context for isolation; animal-windows at risk for dispersal. Until these exist, results 2, 5 and 6 are inventories, not rates.

**7.2 The entity-versus-state decomposition on all four axes** (Result 5, Figure 5). Axis A's version is done. Axis C is next and cheapest: 350 animals is a real sample, and the collar-coverage confound is testable with the check that already passed at Stage 1. Axis B follows once 7.1 gives it group-weeks. Axis D may only support a description.

**7.3 Cohort re-derivation for axes B, C and D.** Axes B and C run on the legacy 1,703,133-row source filtered to 2025-01-01 onward; axis A runs on the frozen 1,924,104-row export from 2024-03-01. Result 2 puts all four axes in one table, so it currently mixes denominators. This is the single largest correctness risk in this document and it must be fixed before Result 2 is written.

**7.4 The dispersal reconciliation.** 28 events / 9 animals versus 71 segments / 17 animals. Pick one definition, document the censoring, and label the product a case series.

**7.5 Fine-scale extension on axis A.** Twelve of 68 structural dyads have 2-minute products. Regenerating that product for more dyads in `New project` is the only route to a between-dyad Stage 2, and it was not attempted in Phases 0–4.

**7.6 Demographic group size.** Absent on every axis. Collar counts are correctly labelled as coverage throughout Phases 2–4 and must not be promoted to group size in the cross-axis tables. Axis B needs this most: split fraction is a fraction of *observed* animals.

**7.7 The cross-axis prediction in §2.** Fission products and internally divided groups are the most permeable units. LapisSplinter and Sapphire both point that way, from four independent saved tables. It is a directional prediction that can fail, and testing it needs 7.1 and 7.3 first. Do not assert it before then; do state it as the prediction the framing generates.

---

## 8. Claim, coverage, verdict

| Claim | Coverage behind it | Verdict |
|---|---|---|
| Identity persists while configuration varies | 81,695 animal-nights; 82.7% origin, 33.2% not whole-group | **report** |
| All four routes are common, brief and partial | 4 inventories, 1 cohort (denominators unreconciled) | **report after 7.3** |
| Shared space rarely means mixing; 4.6% of expectation | 12 dyads, 738 encounters, 12/12 negative | **report** |
| Frequency and depth are separately regulated | Stage 1 vs Stage 2; split occurrence vs composition | **report** (axis A fitted; axis B descriptive) |
| Extreme permeability is a state, not an entity | 5 dyads sustained; contrast +1.376 [+0.338, +1.543] on 16 events / 3 dyads | **report on axis A; extend before generalising** |
| Repeated crossing is dyadic and beyond geometry | OR 1.57 [1.27, 1.93] within 5 km; dyad SD 1.58 vs group SD 0.044 | **report** |
| Split composition is partly anticipated | 466 events, 15 groups, AUC 0.568 [0.558, 0.578], ARI 0.100 | **report with the ARI** |
| Groups differ in internal modularity | median modularity 0 in 19/20 groups at 17% collar coverage | **do not report as biology** — coverage-bound |
| Isolation is an animal type | 285/350 animals; per-animal count median 4, max 364 | **not tested** — needs 7.2 |
| Dispersal rates / who disperses | 17 animals, repeated segments, heavy censoring | **do not report** — case series only |
| Fission products are the most permeable units | LapisSplinter: smallest discrete deficit, 2 of 5 sustained dyads, receives a transfer | **state as prediction**, not result |
| Ecology regulates permeability | one model with known input faults; seasonal term already fails the within-5 km test | **Phase 5** |

---

## 9. Where this leaves the phase plan

Phases 0–4 of the earlier plan are complete and unaffected — they built axis A, and axis A is the model the other three axes should be built to match. The revision is to the ordering of what remains:

| Old | New | Why |
|---|---|---|
| Phase 5 — symmetric ecological test | **deferred** | ecology needs denominators on axes B and C, which do not exist |
| Phase 6 — connectivity | keep, unchanged | Result 6 |
| Phase 7 — case studies and hygiene | keep, unchanged | includes 7.4 |
| — | **new Phase 4b — axis B/C/D opportunity tables and cohort re-derivation** (7.1, 7.3) | prerequisite for Result 2 |
| — | **new Phase 4c — entity-versus-state decomposition on all four axes** (7.2) | Result 5, Figure 5 |

Phases 0–4 plus 4b and 4c produce the paper this framing describes: one cohort, four routes, four coverage-scaled claims, and one general finding — that permeability at every scale is a state a few units enter briefly, and almost none sustain.
