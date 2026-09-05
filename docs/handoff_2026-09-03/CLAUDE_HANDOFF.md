# Baboon group membership: argument, evidence and next work

Prepared 3 September 2026 for continuing this project with Claude.

Project root: `C:\Users\rharel\Documents\group_mebership`. Source paths below are relative to this root unless an absolute path is given. This document is self-contained for an initial scientific discussion. Claude needs access to the project files to inspect or change the analyses.

This handoff combines earlier project framing with a fresh inspection of saved local metadata, result tables, selected scripts and two figures. The accompanying `evidence_snapshot.json` records the numerical checks and SHA-256 hashes of 31 inspected sources. No scientific models were rerun, source GPS data revalidated, remote files accessed, or literature claims checked. The source analyses remain unchanged.

## The main argument

**Baboon groups can retain recognizable social identities while their boundaries remain permeable. Temporary fragmentation, intergroup association and individual movement create connections between groups, but occupying the same space and mixing at close range are distinct processes.**

The strongest current manuscript would explain the *degree and scale of social permeability*. Its contribution is the distinction between structural association, close proximity and sustained membership. Ecological explanations can strengthen that argument once their models are repaired and validated.

The argument develops in four steps:

1. **Membership varies through time.** An origin label alone cannot represent an animal's current social setting. Animals can stay with an origin group, form a splinter, associate with another group, return, or enter a longer association. The membership pipeline makes these states measurable; a final population-wide description still needs a frozen cohort and uncertainty audit.
2. **Association is graded.** A shared spatial cluster does not imply equal close-range mixing across original groups. Compare cross-origin contact with an appropriate within-origin or composition-preserving reference, conditional on tracking opportunity.
3. **Integration can develop differently across scales.** Copper–Lilac is the clearest example: early broad-scale association coexists with a large close-range deficit, and that deficit subsequently narrows. This is consistent with broad spatial fusion preceding close mixing; a lag or change-point estimate has not yet established the timing formally.
4. **A connected population can retain internal structure.** Repeated intergroup proximity and individual transfers are candidate routes connecting groups. Quantifying how much they change population connectivity is a remaining analysis, rather than a demonstrated population-level consequence.

Keep three biological questions separate: **whether groups meet; how long a meeting lasts; how much animals mix during it.** The existing modeling architecture recognizes these questions, but its operational definitions currently need work. In particular, “meeting” in the saved hurdle pipeline means detected positive 5 m cross-group contact, not every broad spatial encounter.

Suggested working title: **Permeable social boundaries: group membership and spatial integration in baboons.** This is a framing suggestion, not a literature-verified novelty claim.

## Evidence that can carry the argument

### Copper–Lilac: the strongest longitudinal example

The effort-corrected analysis defines each pair's contact rate as the fraction of simultaneously observed 2-minute bins within a distance threshold, during canonical fusion hours. Pair rates are averaged within individuals, individuals equally within original groups, and the two original groups equally. The integration ratio is the balanced cross-origin rate divided by the mean within-origin reference. **A ratio of 1 means equal cross- and within-origin contact density; it does not mean every animal is in contact.**

Saved phase means of monthly integration ratios are:

| Radius | Early mixing | Transition | High merge |
|---|---:|---:|---:|
| 2 m | 0.155 | 0.873 | 0.884 |
| 5 m | 0.162 | 0.884 | 0.911 |
| 20 m | 0.277 | 0.930 | 0.962 |
| 100 m | 0.567 | 0.964 | 0.989 |
| 200 m | 0.799 | 0.983 | 0.997 |
| 400 m | 0.940 | 0.996 | 0.999 |

Early mixing is May 2024–July 2025; transition is August 2025–February 2026; high merge is March–June 2026. The phase algorithm also identifies March–April 2024 as low contact, absent from the fusion-conditioned ratio summary. Phases were chosen from monthly large/whole merge probabilities, not specified independently before seeing those data.

At 5 m, the mean monthly cross-origin rate increases from about **1.06% to 3.58%**, while the within-origin reference decreases from **6.70% to 3.93%**. Thus ratio convergence reflects both increasing cross-origin proximity and changing within-origin proximity. These are means of monthly quantities; the mean monthly ratio is not exactly the ratio of the phase-mean rates.

The saved analysis covers 38 animals and 5,785 fusion hours. It requires at least 60 simultaneous bins per pair-month and two partners per individual-month. Monthly intervals resample individuals separately within original group and month. The phase summary itself has no phase-contrast interval. Claims about “complete integration,” statistical equivalence to 1, or interactions involving uncollared animals would exceed these results.

Sources: `outputs/copper_lilac_effort_corrected_integration/copper_lilac_phase_integration_by_radius.csv`, its metadata and monthly summary; `outputs/copper_lilac_phase_alluvial_from_202403/copper_lilac_data_driven_phases.csv`; `analyze_copper_lilac_effort_corrected_integration.py`.

### Across group pairs: shared space does not erase origin structure

The saved 2 m analysis includes 12 supported dyads, 1,863 event units, 5,034 unique hours and 39,975 2-minute rows. All 12 dyad point estimates for mean observed-minus-shuffled cross-edge fraction are negative. Copper–Lilac is much closer to its shuffle reference (−0.0476) than, for example, Chartreuse–Purple (−0.4819) or Jade–Purple (−0.4864).

This supports heterogeneous close-range mixing among observed dyads. The shuffle comparison and the Copper–Lilac cross/within rate ratio have different denominators and weighting, so their numbers are not interchangeable. GPS accuracy and simultaneous sampling are especially consequential at 2 m. The 1,863 event units also follow a different pipeline from the 108 episodes in the hurdle analysis.

Sources: `outputs/all_supported_2m_group_merges/all_supported_2m_merge_dyad_summary.csv`, `all_supported_2m_merge_event_bootstrap_intervals.csv` and metadata.

### Fission: a modest link to preceding association structure

Across 466 eligible split events from 15 groups, preceding 5 m proximity discriminates pairs that enter the same versus different split subgroups with mean event AUC **0.568** (saved event-bootstrap interval **0.558–0.578**). The corresponding within-event label-permutation result is approximately 0.002; subgroup sizes are preserved. Mean adjusted Rand index for recovered communities is only **0.100**.

The useful claim is that preceding proximity contains **modest information** about subsequent subgroup composition. This is not yet an externally validated forecast, a strong predictor, or evidence of deliberate partner choice. Overlapping events, repeated groups, spatial persistence and the short preceding window need sensitivity checks. The plots call the metric “co-sitting,” but the code uses GPS proximity; sitting behavior has not been established by that label.

Sources: `outputs/presplit_cositting_prediction/prediction_summary.csv`, `event_prediction_metrics.csv`; `analyze_presplit_cositting_predicts_split.py`.

### Individual integration: useful examples, unresolved outcome inference

The earlier fine-scale dataset has 28 events from nine animals. Its 14-day modal-membership classification yields 16 “stayed,” two “returned,” and ten censored events. The later established-outcome table has **71 segments from 17 animals**, labeled 46 “Permanent dispersal” and 25 “Returned to origin.” Only **60 segments** contribute centrality estimates: 35 and 25 respectively.

These are different event definitions, not successive sample sizes from one analysis. Several segments belong to the same animal or departure episode. Departure from an origin group also does not establish permanent residence in one particular recipient group. Use observed follow-up and censoring explicitly. Calendar-time plots are preferable for statements about days or weeks since entry; normalized event progress answers a different question.

Sources: `outputs/disperser_integration_calendar_time/disperser_event_outcomes_14d.csv`; `outputs/established_dispersal_centrality/established_events.csv`, `event_week_centrality_cells.csv`; associated scripts. These analyses are best treated as supporting examples until event identity and outcome definitions are reconciled.

## What the current evidence does not establish

“Resource limitation fragments groups whereas favorable conditions connect groups” is a useful hypothesis. It is not the current demonstrated conclusion. The saved meeting model shows a positive NDVI association, but the input issues below prevent using it as a reliable ecological result. A matched model of environmental conditions and fission is also needed to support the two-sided hypothesis.

Proximity does not distinguish affiliative contact, tolerance, aggression or coincident use of a resource. Information exchange, disease transmission, fitness benefits and resource competition are possible implications to discuss with evidence from elsewhere, not measured outcomes here.

VeDBA exists for 101/108 episodes and GPS movement for 97/108. Their Pearson correlations with the saved 5 m integration fraction are approximately 0.278 and 0.209. These are descriptive associations. GPS-controlled models exist, but modest correlations do not establish that activity has been ruled out as an explanation; comparable-sample adjusted effects and the role of activity in the causal question need inspection.

The matched Chartreuse–Purple comparison also resists an overly tidy story: its modularity difference is about 0.0015 with a saved bootstrap interval spanning zero (−0.0254 to 0.0232). It does not establish that another group's presence causes fragmentation. Source: `outputs/chartreuse_modularity_purple_matched/chartreuse_modularity_purple_matched_metadata.json`.

## Missing pieces, in priority order

### 1. Repair the meeting-model denominator and exclusions before interpreting coefficients

**Confirmed issue.** The raw encounter source explicitly excludes Copper–Lilac, while candidate days retain that dyad. The saved daily table has **747 Copper–Lilac days, all with zero encounters and zero saved eligible bins**. Of these, **703 enter the NDVI-complete daily model**. The weekly equivalents are **109 candidate weeks and 103 fitted weeks**, also all zero.

This is an inconsistent exclusion, not evidence that Copper–Lilac never meets. Either include its positive encounter data under the same rule as other pairs or exclude that dyad consistently from the candidate population and fitted sample.

A broader issue needs tracing: **56,000 of 56,128 candidate days** have zero eligible bins in this saved metric table. Candidates are formed from two groups having positions somewhere in a day; unmatched encounter rows are filled with zeros. Some absent rows may represent defensible non-encounters, while others may reflect missing overlap or upstream filtering. Zero saved bins alone does not tell us which. Simply dropping every zero-bin day would also risk selecting only encounters.

**Next deliverable:** a dyad-day opportunity table with shared observation duration and explicit states for observed no encounter, detected encounter, insufficient support and excluded. Validate known examples, including Copper–Lilac, and retain supported negative observations.

Code: `fit_daily_interaction_hurdle_model.py`, `load_daytime_candidate_dyad_days()` and `build_model_rows()` (currently lines 147 and 434). Evidence: daily/weekly CSVs in `outputs/dynamic_social_unit_merge_gamm/daily_interaction_hurdle/` and the accompanying snapshot.

### 2. Separate demographic group size from collar coverage

**Confirmed issue.** `group_size_total` equals `group_a_animals + group_b_animals` in every saved candidate day and week. Those animal counts are observed collars. Event size predictors are also constructed from averages of observed collar counts, although plots label them “Total group size.”

The nominal demographic sizes used by other plotting scripts have not thereby become demographic covariates in this model. Existing size effects cannot be interpreted as effects of whole-group population size.

**Next deliverable:** join independently sourced group-size estimates, with dates and provenance where available; retain collar coverage as a separate observation variable. If only static demographic estimates exist, state that limitation. Refit only after the observation denominator is sound.

Code: `fit_daily_interaction_hurdle_model.py`, `build_model_rows()`, `aggregate_to_weekly_rows()` and `build_event_duration_rows()`.

### 3. Align “encounter” and “integration” with the scientific question

**Confirmed definition.** `build_positive_episodes()` first filters to `has_cross_edge == True` (currently line 341), stitches those bins, and sums their edges. Thus the saved event integration fraction omits supported zero-cross-contact bins, including such bins between positive contacts. It describes contact-positive portions of an episode.

To ask how integrated animals are *during a broad encounter*, detect that encounter independently at the structural scale, then retain all eligible fine-scale bins, including zero cross-contact bins. Distinguish a zero rate with observation from an undefined rate with no denominator. If positive-contact episodes remain the target, name and interpret them accordingly.

Also report **elapsed episode span** separately from **observed active-contact time**. The saved rule bridges positive bins up to 14 hours apart. This is not evidence of continuous overnight interaction. Test reasonable shorter/longer gap rules and account for bouts already underway or still ongoing at coverage boundaries.

**Next deliverable:** one operational definition sheet, a worked event showing retained and omitted bins, and an event table with structural span, supported exposure and active-contact time.

### 4. Freeze the cohort and validate membership decisions

Multiple products are in use:

| Product | Saved scope | Role |
|---|---|---|
| Full narrow hourly export | 1,924,104 rows; 350 animals; 26 origin groups; March 2024–July 2026 | Membership data product |
| Nightly export from that source | 81,695 animal-nights; 2,162 ambiguous nights | Evening-labeled biological nights |
| Older dynamic proximity-status product | 1,070,061 rows; 354 animals; 21 analysis groups | Input to existing encounter models |
| Older weekly spline/GEE | 96 dyad-weeks/21 dyads; NDVI subset 75/18 | Exploratory link-fraction models |
| Current saved meeting models | 51,582 days or 7,605 weeks; 190 dyads after covariate filtering | Occurrence models needing repair |
| Current saved event models | Duration: 88 events/20 dyads with NDVI; integration: 108/21 without NDVI | Different fitted samples |

These counts are not contradictory once their definitions are stated. A newer export does not mean earlier results were rerun against it. The older dynamic-status metadata also records **96,349 unmatched status rows**, whose membership fallback requires review.

Keep `origin_group` (historical identity), `assigned_social_unit` (immediate assignment) and `dynamic_social_unit` (conservative longer-term assignment) distinct. Select the field for the question and document it; origin comparisons intentionally preserve origin labels. Many current scripts use `dynamic_social_unit`, so verify the choice rather than relabeling their results as immediate membership.

Do not equate `temp_group_id`, an hourly spatial cluster, with `association_event_id`, a continuing association event. Nightly inference must retain coverage, ambiguity, exact-hour, local-2h and carried-night support separately.

The Emerald/Bronze case in `docs/emerald_bronze_sparse_case_24AA11.md` is a concrete validation target: sparse fixes from a same-origin companion weaken a confident lone-animal transfer interpretation. This handoff verifies the note exists; it does not establish that its recommended change was applied to every downstream product.

**Next deliverable:** a source-to-analysis cohort ledger, event-ID crosswalk and small manually reviewed set of clear and ambiguous membership transitions. Prioritize exact-hour versus local-support sensitivity and collar-loss cases.

### 5. Establish statistical reliability after the input repairs

The weekly “GAMM” files use Python spline-basis binomial GEE with group fixed effects. The displayed equivalent `mgcv` formula is not an executed R fit. Other trajectory scripts use MixedLM and should be described separately.

The hurdle script's Bayesian logistic models use a normal approximation around a penalized optimum. Its Gaussian event models use conjugate posterior sampling with fixed prior scales for group and dyad indicators. They do not estimate a conventional set of group/dyad variance components with hyperpriors. Columns named `p_value` in these Bayesian tables are computed from posterior sign proportions, not classical frequentist p-values.

Reassess the model rather than just renaming it. Use an uncertainty model appropriate to repeated dyads, groups and animals, and a likelihood appropriate to the chosen rate or duration. Fine-scale edges and bins are not independent replicates. For unordered group pairs, consider a shared group effect contributing at either end of the dyad; separate alphabetical `group_a`/`group_b` effects need justification.

Inspect prior sensitivity, predictive calibration, residual dependence and influence of individual dyads. If using MCMC, save its convergence diagnostics; those diagnostics would not validate the current normal-approximation fit. Compare effects on the same retained samples before attributing changes to added predictors.

For temporal explanations, use genuinely preceding covariates. Existing occurrence distance is a pair's mean across its available days; event predictors use a window extending one week before through one week after the event. A before/after distance comparison already exists, but its size and NDVI covariates still use the wider window. Separate between-pair proximity from within-pair temporal changes.

### 6. Test how general the spatial-scale result is

Copper–Lilac is one important trajectory. Apply a common, observation-aware metric across other supported dyads, with original labels preserved when testing origin mixing. Examine sensitivity to radius, same-bin timing differences, GPS error, partner availability and stable collar subsets.

For Copper–Lilac, report uncertainty for phase contrasts and continuous calendar-time change, and test alternative phase boundaries. Individual resampling within each month does not by itself address shared dyads and serial dependence across months. Show the cross- and within-origin rates alongside their ratio so a falling reference is visible. Large-radius rates approach a ceiling, limiting discrimination.

**Next deliverable:** a radius-by-time figure with support and uncertainty, plus a restrained assessment of which other dyads show similar or different trajectories. Extend population connectivity claims only if time-respecting connections are actually quantified.

### 7. Resolve dispersal outcomes and selective follow-up

Reconcile departure episodes, recipient segments, returns and observation endpoints before comparing outcomes. Cluster repeated segments at the animal/departure level as appropriate. Inspect whether apparent change over weeks is instead caused by the changing set of animals still observed.

Two concrete problems should be retained in Claude's audit:

- In `compare_centrality_stayed_left.py`, bootstrap intervals are calculated for stayers but merged onto both outcome categories in the CSV. The figure uses the leavers' observed range, but the CSV's leaver `low`/`high` values are not leaver intervals. Correct the exported table before reuse.
- `plot_permanent_5m_colocation_cadence_filtered.py` retains animals with **more than 100 positive 5 m contact bins**. This is selection on the measured outcome, not just tracking opportunity. These selected examples cannot estimate typical integration among all dispersers. Cadence filtering and positive-contact selection must be described separately.

The earlier 14-day classifier treats any nonempty follow-up as classifiable; adequate support needs an explicit rule. “No observed return by the endpoint” is the defensible observation when permanent fate is unknown.

### 8. Decide how much ecology and mechanism the paper needs

The core argument can be developed through dynamic membership, heterogeneity in mixing and the Copper–Lilac scale/time result. Fission organization and individual entry trajectories provide supporting strands.

To make ecological opportunity the main explanation, add the missing environmental analysis of fission, repair the meeting pipeline, account for seasonal/shared-resource alternatives, and make directional predictions that could fail. A literature review is still needed to place the operational measures and novelty claims. Weather effects, affiliative behavior, demographic benefits and population transmission are not established by the artifacts inspected here.

## A practical manuscript and figure sequence

1. **Define the social states and observation coverage.** Show representative membership timelines and a cohort diagram. An event-size/duration overview exists under `outputs/canonical_group_merge_scale_log_scatter/`, but uses an older source and needs reconciliation before publication.
2. **Show that broad association permits different degrees of close mixing.** Start from `outputs/all_supported_2m_group_merges/all_supported_2m_group_merge_integration.png`; retain support counts, uncertainty and a distance-sensitivity companion.
3. **Show the Copper–Lilac trajectory across scales.** Start from `outputs/copper_lilac_effort_corrected_integration/copper_lilac_integration_scale_summary.png` and `copper_lilac_effort_corrected_integration_5m.png`. The scale plot was visually inspected for this handoff. Add temporal support and phase-contrast uncertainty before treating it as a final inferential figure.
4. **Choose the final supporting result after scope is settled.** The fission figure `outputs/presplit_cositting_prediction/presplit_cositting_prediction.png` was visually inspected and supports a modest association. A repaired occurrence/duration/integration figure could instead anchor an ecology-focused paper. Keep sparsely supported dispersal trajectories supplementary unless their inferential problems are resolved.

The discussion should explain what these results reveal about permeable boundaries, distinguish the mechanisms still open, and state the limits imposed by observing a changing subset of animals.

## How Claude should continue

Start with a short claim–evidence–gap table and an agreed core paper scope. Then inspect the three highest-priority issues: the meeting denominator/exclusion, collar counts labeled as group size, and positive-bin conditioning of integration. Produce concrete corrected definitions and small worked examples before a broad rerun. Preserve current source data and place candidate results in a new dated output directory.

The first successful continuation should deliver a reproducible cohort/observation audit, a proposed final event definition, and an ordered analysis plan tied to manuscript claims. A polished abstract is premature if it promotes the current ecological coefficients to settled results.

Useful entry points are `fit_daily_interaction_hurdle_model.py`, `analyze_copper_lilac_effort_corrected_integration.py`, `analyze_presplit_cositting_predicts_split.py`, and the tables named above. Upstream membership builders and larger GPS caches are mainly under `C:\Users\rharel\Documents\New project`; this project folder is not a standalone reproduction bundle.

To refresh this handoff's saved-output checks, run `python docs\handoff_2026-09-03\audit_saved_outputs.py` from the project root. It rewrites only `evidence_snapshot.json`, does not refit models, and records the inspected source hashes. The prose is a dated assessment and should be updated deliberately if later outputs change.
