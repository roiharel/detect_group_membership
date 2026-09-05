# General structure and continuation plan

Prepared 3 September 2026. Companion to `docs/handoff_2026-09-03/CLAUDE_HANDOFF.md`.

Purpose: restate the project as a **population-level** argument in which no single dyad carries a claim, and give an ordered plan to finish it. The Copper–Lilac trajectory is retained, but as one endpoint on a general gradient rather than as the paper's spine.

All counts below were read from saved outputs in this project. No model was rerun.

---

## 1. Why the structure needed to change

Two facts, both from saved tables, argue against a case-led paper.

**Copper–Lilac dominates the fine-scale evidence.** In `all_supported_2m_merge_dyad_summary.csv`, Copper–Lilac supplies 27,891 of 39,975 2-minute rows — **69.8%**. Any statistic pooled across those rows is largely a Copper–Lilac statistic. Cross-dyad claims must therefore be dyad-weighted, not row-weighted.

**Copper–Lilac is the outlier, not the exemplar.** All 12 supported dyads have a negative mean observed-minus-shuffled cross-group edge fraction, but the spread is wide:

| Dyad | mean obs − shuffle | 2-min rows |
|---|---:|---:|
| Copper – Lilac | −0.048 | 27,891 |
| Maroon – Sapphire | −0.236 | 624 |
| Lapis – LapisSplinter | −0.304 | 3,187 |
| LapisSplinter – Magenta | −0.349 | 930 |
| Lapis – Periwinkle | −0.356 | 1,121 |
| Bronze – Lilac | −0.384 | 374 |
| Bronze – Magenta | −0.400 | 677 |
| Emerald – Lilac | −0.455 | 568 |
| Copper – Magenta | −0.477 | 1,069 |
| Chartreuse – Purple | −0.482 | 777 |
| Lilac – Magenta | −0.485 | 1,008 |
| Jade – Purple | −0.486 | 1,749 |

Read generally, this is the paper's strongest simple result: **when two groups share space, close-range mixing is almost always far below what composition alone would produce.** Copper–Lilac is the rare dyad that approaches its own composition expectation. That makes it the demonstration of what full fusion looks like — the top of a gradient — and it is more persuasive in that role than as the lead result.

Structurally, this converts a case study into a **distribution with a documented tail**.

> **Revised by Phase 3** — see [PHASE0_3_RESULTS.md](PHASE0_3_RESULTS.md). Separating discrete encounters from sustained association sharpens this considerably. Copper–Lilac's *encounters* give a 5 m deficit of −0.339, tenth of twelve — inside the pack. Only its *fused* period reaches −0.033. The outlier is therefore a **state, not a dyad**: encounters barely mix everywhere (dyad-weighted −0.405 [−0.442, −0.367], 12 of 12 below expectation), and near-expectation mixing belongs to sustained association, a state five dyads entered.

---

## 2. The general claim architecture

Three nested levels. Each level's claim is supported by a population-scale dataset, and each level constrains the next.

### Level 1 — Membership is a state, not a label

An origin label does not describe an animal's current social setting. Measured over 350 animals, 26 origin groups and 29 months (1,924,104 animal-hours; 81,695 animal-nights):

| Nightly membership | Animal-nights | Share |
|---|---:|---:|
| With origin group | 67,561 | 82.7% |
| With another group | 11,418 | 14.0% |
| Isolated | 2,499 | 3.1% |
| Mixed / unclear | 217 | 0.3% |

| Nightly association structure | Animal-nights | Share |
|---|---:|---:|
| Whole group | 54,572 | 66.8% |
| Between-group merge | 20,468 | 25.1% |
| Within-group split | 3,789 | 4.6% |
| Isolated / low support | 2,866 | 3.5% |

**Claim.** Group identity persists — animals are with their origin group on roughly five nights in six — yet a quarter of animal-nights are spent in a configuration that includes another group. Persistence and permeability are simultaneous properties of the same population, not competing descriptions.

This level is the paper's foundation and is already supported by a frozen, validated data product. It needs uncertainty and coverage treatment, not new analysis.

### Level 2 — Encounter and mixing are different processes measured at different scales

The event inventory in `outputs/canonical_group_merge_scale_log_scatter/` contains 13,615 events across five types:

| Event type | Events |
|---|---:|
| Single-animal separation | 5,828 |
| Within-group split | 4,570 |
| Large merge | 1,420 |
| Small subset merge | 1,132 |
| Medium partial merge | 665 |

The 3,217 merge events span **74 dyads**. Copper–Lilac contributes 615 of them (19%) — substantial, not dominant.

Against that, the fine-scale meeting model in `daily_interaction_hurdle_daily_event_rows.csv` has **56,128 candidate dyad-days across 207 dyads but only 117 positive days in 21 dyads — a 0.21% positive rate** — and 56,000 of those days carry zero eligible 2-minute bins.

The two numbers are not in conflict; they measure different things. But it means the encounter response currently used for the ecological models is a **rare-event indicator of detected close contact**, not an indicator of groups meeting. That definitional gap is the largest structural problem in the project, and it is general — it does not depend on any dyad.

**Claim.** Groups meet at broad scale far more often than they mix at close range, and the two must be modeled as separate stages: whether a structural encounter occurs, and — given one — how much animals actually mix.

### Level 3 — A permeable population can retain internal structure

Candidate routes connecting groups: repeated dyadic encounters (74 dyads with merges, concentrated in a much smaller set), individual departures (5,828 separation events; 71 dispersal segments from 17 animals), and within-group fission that redistributes companions (4,570 split events; 466 analysed events from 15 groups, preceding-proximity AUC 0.568, bootstrap 0.558–0.578, mean ARI 0.100).

**Claim, currently unproven.** These routes generate a structured, non-random network among groups. This is the level most likely to justify a broad-audience venue and the level with the least completed work. Treat it as the paper's stated ambition, and assert it only once time-respecting connectivity is quantified.

---

## 3. The one decision that unblocks the rest

**Define the encounter unit at the structural scale, and measure mixing inside it.**

Two-stage definition:

1. **Structural encounter.** Two social units occupy a shared spatial cluster (the existing `association_event_id` / merge-scale machinery), detected independently of close-range contact. Stage-1 sample: on the order of 3,217 events across 74 dyads.
2. **Mixing within the encounter.** Within each detected encounter, retain **all** eligible fine-scale bins with a valid denominator, including bins with zero cross-group contact. Report a rate with observation, distinct from an undefined rate with no denominator.

The saved pipeline currently does the opposite: `build_positive_episodes()` filters to `has_cross_edge == True` before stitching, so the reported integration fraction describes only the contact-positive part of an episode. Every ecological coefficient inherits that conditioning.

Adopting the two-stage definition:

- raises the stage-1 sample from 21 dyads to about 74, making a general model feasible;
- removes the rare-event pathology from the occurrence model;
- makes "shared space without mixing" a measured outcome rather than a missing row;
- places the Copper–Lilac trajectory on the same axis as every other dyad.

Everything in the plan below is ordered around this.

---

## 4. Definitions to freeze before any refit

| Quantity | Current problem | Required form |
|---|---|---|
| Encounter opportunity | Candidate days built from both groups having any position that day; unmatched rows zero-filled | Dyad-window table with shared observation duration and four explicit states: observed no encounter, detected encounter, insufficient support, excluded |
| Encounter occurrence | Positive 5 m cross-group contact | Structural encounter (stage 1) |
| Mixing | Cross-edge fraction over contact-positive bins only | Cross-edge rate over all supported bins within the encounter, against an origin- or composition-preserving reference |
| Duration | Elapsed span bridging positive bins up to 14 h apart | Report separately: structural span, supported exposure, active-contact time |
| Group size | `group_size_total = group_a_animals + group_b_animals`, i.e. observed collars, but labeled and plotted as group size | Demographic size as a covariate with provenance; collar coverage retained as a separate observation variable |
| Membership field | Scripts mix `origin_group`, `assigned_social_unit`, `dynamic_social_unit` | One documented choice per question; origin labels preserved wherever mixing is being tested |
| Cohort | At least three membership products in circulation (1,924,104 / 1,703,133 / 1,070,061 rows; 26 vs 21 groups; different start dates) | One frozen cohort with an event-ID crosswalk; older results re-derived or explicitly labeled legacy |

The event inventory in §2 comes from the 1,703,133-row legacy source, filtered to 2025-01-01 onward, while the frozen membership export covers 1,924,104 rows from 2024-03-01. Reconciling those two is a prerequisite for reporting population rates, because the denominator differs.

---

## 5. General manuscript structure

Working title: **Persistent groups with permeable boundaries: how often animal groups meet, and how little they mix.**

| # | Result | Population-scale evidence | Status |
|---|---|---|---|
| 1 | Group identity persists while membership varies | 350 animals, 26 groups, 81,695 animal-nights; state and association distributions | Data frozen; needs coverage and ambiguity treatment |
| 2 | Social reconfiguration is common and structured by type and scale | 13,615 events in five types; size and duration distributions across 74 dyads | Exists on the older source; needs cohort reconciliation |
| 3 | Shared space rarely means close mixing | 12 supported dyads, all below composition expectation, spread −0.05 to −0.49 | Strongest general result; needs dyad-weighting, support-aware uncertainty, radius sensitivity |
| 4 | Encounter and mixing have different predictors | Two-stage model over about 74 dyads | Blocked on §3 and §4 |
| 5 | Permeability structures the population network | Repeated dyads, transfers, fission-mediated redistribution | Ambition; least complete |

Supporting analyses, not main results: the Copper–Lilac trajectory (endpoint case for Result 3), the Chartreuse–Purple matched comparison (a null: modularity difference 0.0015, interval −0.025 to 0.023), individual dispersal centrality trajectories, and activity controls (VeDBA on 101/108 events, GPS movement on 97/108; r ≈ 0.28 and 0.21 against saved 5 m integration).

### Figure plan

1. **Cohort and states.** Animals × time coverage, state ribbon, nightly state distribution. Establishes scale and honesty about observation.
2. **Event inventory.** Size versus duration by event type, five panels, log axes. Shows that reconfiguration is routine and multi-scale.
3. **The mixing gradient.** Per-dyad observed minus expected cross-group contact, ordered, with support counts and intervals; radius as a second dimension. Copper–Lilac appears as the labeled right-hand tail.
4. **Two-stage encounter model.** Structural-encounter probability, and mixing conditional on encounter — predicted responses, not per-group coefficient forests.
5. **Population network.** Dyadic encounter network against an opportunity-adjusted null.

Copper–Lilac's own radius-by-time trajectory becomes a supplementary figure supporting panel 3.

---

## 6. Continuation plan

Each phase has a completion criterion. Do not start a phase before its predecessor's criterion is met. Put new work in `outputs/general_structure_2026_09/`; do not overwrite existing outputs.

### Phase 0 — Cohort ledger and definition sheet (prerequisite)

Build a source-to-analysis ledger for all membership products (rows, animals, groups, date range, builder parameters, and which saved analyses consumed it) plus an event-ID crosswalk between the 1,924,104-row canonical export and the legacy 1,703,133-row inventory source. Write one definitions sheet covering every row of §4.

*Complete when:* every saved result in §5 maps to exactly one cohort, and every quantity in §4 has one written definition with a worked example.

### Phase 1 — Opportunity table

Replace zero-filled candidate days with a dyad-window opportunity table carrying shared observation duration and the four explicit states. Validate against known examples, including the Copper–Lilac days currently coded as 747 zero-encounter days.

*Complete when:* the number of dyad-windows in each state is reported, no state is inferred from a missing row, and supported negatives are distinguishable from unobserved windows.

### Phase 2 — Two-stage encounter unit

Implement stage 1 (structural encounter) and stage 2 (mixing over all supported bins within it). Produce an event table with structural span, supported exposure and active-contact time as separate columns, plus one worked event showing which bins are retained, which are dropped, and why.

*Complete when:* the event table reproduces a hand-checked example, and stage-1 counts are reported per dyad with observation support.

### Phase 3 — The general mixing gradient (first publishable general result)

Apply one observation-aware metric across all supported dyads with origin labels preserved. Weight by dyad, not by bin. Report cross-origin and within-origin rates alongside their ratio, so a falling reference is visible. Sensitivity: radius, same-bin timing tolerance, GPS error, partner availability, stable-collar subsets.

*Complete when:* Figure 3 exists with per-dyad support and intervals, and the pooled statement is demonstrably not a restatement of one dyad.

### Phase 4 — Two-stage model with sound inference

Fit stage 1 and stage 2 on the repaired samples. Use a likelihood matched to the response, uncertainty matched to repeated dyads, groups and animals, and a shared group effect for unordered pairs rather than alphabetical `group_a` / `group_b` terms. Rename or replace the current normal-approximation and conjugate-sampling fits; the columns named `p_value` in saved Bayesian tables are posterior sign proportions, not frequentist p-values. Use genuinely preceding covariates, and separate between-dyad proximity from within-dyad temporal change.

*Complete when:* the same retained sample is used across model comparisons, prior sensitivity and calibration are reported, and no result is attributed to an added predictor without a matched-sample comparison.

### Phase 5 — Symmetric ecological test

Fission and intergroup encounter share one environmental framework, with directional predictions that can fail: scarcity fragments, abundance connects. Add the missing environmental model of fission, control for season and shared concentrated resources, and separate demographic size from collar coverage.

*Complete when:* both processes are modeled on the frozen cohort with the same covariates, and the two-sided prediction is either supported or rejected in print.

### Phase 6 — Connectivity

Quantify time-respecting connections among groups from encounters and transfers, against an opportunity-adjusted null.

*Complete when:* the network claim in Result 5 is either quantified or explicitly withdrawn from the main text.

### Phase 7 — Case studies and supplement hygiene

Recast Copper–Lilac as the labeled endpoint of Phase 3, with phase-contrast uncertainty and alternative phase boundaries. Resolve the dispersal-outcome reconciliation (71 segments from 17 animals versus 28 events from 9 animals; repeated segments per animal; censoring). Fix the two confirmed export faults: leaver bootstrap intervals mislabeled in `compare_centrality_stayed_left.py`, and outcome-based selection (more than 100 positive 5 m bins) in `plot_permanent_5m_colocation_cadence_filtered.py`.

*Complete when:* every supplementary claim states its selection rule, and no figure selects on the measured outcome without saying so.

### Ordering note

Phases 0–3 are the critical path and produce a defensible paper on their own: persistence with permeability, a multi-scale event inventory, and a general mixing gradient. Phase 5 decides whether the headline is *structural* (how much groups mix) or *ecological* (what makes boundaries open). Recommendation: build 0–3 first and keep ecology as a test appended in Phase 5, rather than promising an ecological mechanism the current denominators cannot support.

---

## 7. What to stop doing

The project has roughly 20 single-purpose scripts targeting individual dyads and animals (Chartreuse variants, disperser cadence, permanent-5 m variants). They were useful for discovery. Continuing to add them raises the reconciliation burden without moving any of the five results forward. Freeze them, and add new analysis only through the shared Phase 2 event table, so that every dyad is measured the same way.

---

## 8. Claim–evidence–gap summary

| Claim | Evidence now | Gap |
|---|---|---|
| Identity persists while membership varies | 82.7% of 81,695 animal-nights with origin group; 14.0% with another group | Coverage and ambiguity weighting; reconcile 26 vs 21 groups |
| Reconfiguration is routine and multi-scale | 13,615 events, five types, 74 dyads | Older cohort; needs re-derivation on the frozen source |
| Shared space rarely means mixing | 12 of 12 dyads below composition expectation, −0.05 to −0.49 | Dyad-weighted estimate, support-aware intervals, radius sensitivity |
| Full fusion is possible | Copper–Lilac 5 m ratio 0.162 → 0.911 across phases | Phase-contrast uncertainty; position on the general gradient |
| Encounter and mixing differ in predictors | Saved hurdle models | Blocked: denominator, exclusions, positive-bin conditioning, size covariate |
| Fission composition is partly anticipated | 466 events, 15 groups, AUC 0.568 (0.558–0.578) | Sensitivity to overlapping events and window length |
| Permeability structures the population | Repeated dyads, 71 dispersal segments | Not yet quantified; do not assert |
| Ecology regulates permeability | Positive NDVI association in a model with known input faults | Repair inputs; add a fission model; test symmetrically |
