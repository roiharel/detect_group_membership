# Phases 0–3: results

Prepared 3 September 2026. Executes Phases 0–3 of [PAPER_STRUCTURE_GENERAL.md](PAPER_STRUCTURE_GENERAL.md).

Scripts: [phase0_cohort_ledger.py](../../phase0_cohort_ledger.py), [phase1_opportunity_table.py](../../phase1_opportunity_table.py), [phase2_two_stage_events.py](../../phase2_two_stage_events.py), [phase3_mixing_gradient.py](../../phase3_mixing_gradient.py).
Outputs: `outputs/general_structure_2026_09/`. No existing output was modified.

Frozen source for all four phases:
`New project/outputs/canonical_robust_hourly_membership_shared_full_20260722/canonical_hourly_membership_with_association_events.parquet`
(1,924,104 animal-hours, 350 animals, 26 origin groups, 2024-03-01 to 2026-07-22).

---

## Part 1 — What the two definitions actually measure

### The meeting model

`fit_daily_interaction_hurdle_model.py` builds its response in two unconnected halves.

**Denominator** (`load_daytime_candidate_dyad_days`): for each calendar day, form every unordered pair of dynamic social units that each had at least one daytime position. Nothing about proximity enters. The result is 56,128 dyad-days across 207 dyads, and the centroid separation on those days is:

| | 5% | 25% | 50% | 75% | 95% | max |
|---|---:|---:|---:|---:|---:|---:|
| Centroid distance (m) | 2,102 | 5,939 | 10,676 | 15,701 | 23,857 | 37,906 |

The median candidate "opportunity to meet" is two groups **10.7 km apart**.

**Numerator** (`load_positive_interactions`): a dyad-day counts as a meeting if it appears in one fine-scale 2-minute file with at least one bin containing a 5 m cross-group edge. On positive days the centroid separation is median **440 m** (max 2,930 m).

The two halves are joined with a left merge and `fillna(0)`. That is where the artefact enters:

- only **128** of 56,128 candidate dyad-days appear in the fine-scale file at all;
- **117** of those 128 are positive;
- the remaining **56,000** days are zero-filled, and their `eligible_2min_bins` is 0 — there is no exposure measure behind the zero;
- the fine-scale file is `bigmerge_..._no_copper_lilac_2min_metric_rows.csv`. The exclusion is in the filename, but the excluded dyad's candidate days stay in the denominator as zeros.

So the fitted question is: *given that two collared groups were each tracked somewhere in a 38 km landscape today, was a 5 m cross-group contact recorded in a file that was itself built from already-detected large merges and that drops one dyad by construction?* The 0.21% positive rate is a property of that construction, not of baboons.

### The structural merge

`plot_canonical_group_merge_scale_log_scatter.py` asks a different question at a different scale. A merge exists in an hour when animals from two social units are members of **the same spatial cluster** — adaptive clustering with a 120–900 m edge range, so the scale is hundreds of metres, not 5 m. Consecutive hours are stitched into events, and events are classified by what fraction of each group's observed animals is present.

There is no denominator, because it is an inventory rather than a model: 3,217 merge events across 74 dyads.

### The three differences

| | Meeting model | Structural merge |
|---|---|---|
| Spatial scale | 5 m inter-individual contact | cluster co-membership, order 10² m |
| Denominator | every co-tracked dyad-day, median 10.7 km apart | none — an inventory |
| Provenance of the positive class | a filtered downstream product that excludes one dyad | the membership table itself |
| Absence means | zero-filled, exposure unknown | not applicable |

Both are legitimate measurements. Using one as the response and the other as the concept is what makes the ecological coefficients uninterpretable.

---

## Part 2 — Phase 0: cohort ledger

20 products inventoried across four tiers, in [cohort_ledger.csv](../../outputs/general_structure_2026_09/phase0_cohort_ledger/cohort_ledger.csv), with an identifier crosswalk in [event_id_crosswalk.csv](../../outputs/general_structure_2026_09/phase0_cohort_ledger/event_id_crosswalk.csv).

Three corrections to figures previously in circulation:

1. The legacy inventory source (`canonical_robust_hourly_membership_local_2h_support`) holds **1,703,133 rows and 347 animals**. The 1,516,740 figure quoted from its metadata is `rows_used` after filtering to 2025 onward, not the source size.
2. The saved candidate table spans **207 dyads over 21 groups** because it derives from `proximity_status_dynamic_social_unit.parquet`, whose unit set includes `SneakySilver` but omits Gold, Green, Mint, Red and TrickyTeal. The frozen source has 26 units and 317 co-observed dyads.
3. The 2 m fine-scale product's file is *named* `canonical_5m_shuffle_expectation_2min_rows.csv`. Radius is set by the directory, not the filename. Anything that resolves this file by name alone will silently analyse the wrong radius.

The crosswalk records which identifier answers which question: `temp_group_id` (one hour, one cluster), `association_cluster_id` (the join key for fine-scale bins), `association_event_id` (Stage-1 encounter identity), and the legacy `merge_episode_id` / positive `episode_id`, both of which should now be used only to reproduce saved results.

---

## Part 3 — Phase 1: the opportunity table

[opportunity_dyad_day.csv](../../outputs/general_structure_2026_09/phase1_opportunity/opportunity_dyad_day.csv) — 73,353 dyad-days, 317 dyads, 26 social units, 873 days. Every row carries shared observed hours, hours with exact-hour support for both units, minimum and median **hourly** centroid distance (not a pair-level mean), and one explicit state.

| State | Dyad-days | Dyads | Median shared hours | Median of min hourly distance |
|---|---:|---:|---:|---:|
| observed_no_encounter | 70,312 | 316 | 24 | 10,783 m |
| detected_encounter | 2,867 | 68 | 24 | 66 m |
| insufficient_support | 174 | 119 | 7 | 11,123 m |

The state assignment separates cleanly on a variable it never used: encounter days have a median closest approach of **66 m**, non-encounter days **10.8 km**.

### The fabricated negatives are general

Cross-tabulating against the saved candidate table:

| Phase 1 state | Saved call | Dyad-days |
|---|---|---:|
| detected_encounter | saved as zero meetings | **2,406** |
| detected_encounter | absent from saved table | 344 |
| detected_encounter | saved as positive | 117 |
| observed_no_encounter | saved as zero | 52,848 |
| observed_no_encounter | absent from saved table | 17,464 |
| insufficient_support | saved as zero | 8 |

The 2,406 fabricated negatives span **53 dyads**. Copper–Lilac contributes 592 of them — the largest single share, but under a quarter. Copper–Magenta (200), Lilac–Magenta (155), Lapis–LapisSplinter (121), Jade–Purple (114) and 48 further dyads carry the rest. Repairing the Copper–Lilac exclusion alone would not fix the response.

All 749 saved rows absent from Phase 1 are `SneakySilver` dyads, a unit that does not exist in the frozen source's dynamic units.

---

## Part 4 — Phase 2: the two-stage encounter unit

[stage1_events_with_stage2_mixing.csv](../../outputs/general_structure_2026_09/phase2_two_stage_events/stage1_events_with_stage2_mixing.csv) — **1,705 Stage-1 events across 68 dyads** (3-hour stitching gap), each with three separately reported durations:

- `structural_span_hours` — last encounter hour minus first, plus one;
- `supported_exposure_hours` — eligible fine-scale bins × 2 min;
- `active_contact_hours` — contact-positive bins × 2 min.

Median structural span 14 h; median supported exposure 2.3 h; median active contact 0.17 h where fine-scale data exist. Those three numbers are not interchangeable, and the saved pipeline reported one figure for all three roles.

### Discrete encounters versus sustained associations

Sixteen Stage-1 runs last ≥ 168 h (7 days) with observed-hour fraction ≈ 1. These are not encounters; they are sustained associations, and the 168 h threshold matches the membership builder's own 7-day dynamic-reassignment rule.

They occur in **five dyads**, not one: Bronze–Magenta, Copper–LapisSplinter, Copper–Lilac, LapisSplinter–Lilac, Maroon–Sapphire. Copper–LapisSplinter and LapisSplinter–Lilac share identical windows in August–September 2025, which is a three-way fusion, not two dyadic events.

The remaining 1,689 events are discrete encounters with a median span of 14 h.

### Fine-scale coverage and what the saved pipeline dropped

751 events (12 dyads) have 5 m bins; 668 have 2 m bins. Across the retained 5 m exposure, **42.5% of eligible bins carry zero cross-group contact** — precisely the bins `build_positive_episodes()` removes before measuring integration.

Worked example, [worked_example_bins.csv](../../outputs/general_structure_2026_09/phase2_two_stage_events/worked_example_bins.csv): Lapis–LapisSplinter event `E0080`, 19–21 May 2026, structural span 42 h, **517 eligible bins**, of which the saved rule keeps **53**. It discards **89.8%** of the observed exposure and then reports the remainder as the encounter.

Gap-rule sensitivity — the dyad count never moves:

| Gap rule | Events | Dyads | Discrete | Sustained | Median span |
|---:|---:|---:|---:|---:|---:|
| 2 h | 1,812 | 68 | 1,796 | 16 | 14 h |
| 3 h | 1,705 | 68 | 1,689 | 16 | 14 h |
| 6 h | 1,547 | 68 | 1,531 | 16 | 14 h |
| 14 h | 1,109 | 68 | 1,093 | 16 | 6 h |

---

## Part 5 — Phase 3: the general mixing gradient

Estimator: edge-weighted deficit = `sum(cross_edges)/sum(total_edges)` minus the edge-weighted composition expectation carried in the fine-scale product. Uncertainty: percentile bootstrap over Stage-1 encounters within a dyad, 2,000 draws, seed 20260903. Pooling by dyad, with the bin-weighted value reported alongside.

Figure: [mixing_gradient_5m.png](../../outputs/general_structure_2026_09/phase3_mixing_gradient/mixing_gradient_5m.png).

### Discrete encounters — all twelve dyads, none near expectation

| Dyad | Observed | Expected | Deficit | 95% interval | Bins | Events |
|---|---:|---:|---:|---|---:|---:|
| Jade – Purple | 0.0225 | 0.5083 | −0.4858 | −0.4926, −0.4793 | 7,960 | 119 |
| Chartreuse – Purple | 0.0047 | 0.4875 | −0.4828 | −0.4959, −0.4700 | 3,033 | 24 |
| Lilac – Magenta | 0.0363 | 0.5168 | −0.4805 | −0.4929, −0.4639 | 5,533 | 57 |
| Copper – Magenta | 0.0354 | 0.5045 | −0.4692 | −0.4805, −0.4554 | 7,934 | 128 |
| Emerald – Lilac | 0.0438 | 0.4961 | −0.4523 | −0.4729, −0.4247 | 3,849 | 46 |
| Bronze – Magenta | 0.0616 | 0.5082 | −0.4467 | −0.4617, −0.4270 | 3,032 | 48 |
| Bronze – Lilac | 0.0535 | 0.4353 | −0.3818 | −0.4186, −0.3451 | 1,691 | 22 |
| LapisSplinter – Magenta | 0.1386 | 0.4946 | −0.3560 | −0.3869, −0.3283 | 6,559 | 86 |
| Lapis – Periwinkle | 0.0142 | 0.3565 | −0.3423 | −0.3653, −0.3232 | 4,089 | 65 |
| Copper – Lilac | 0.1972 | 0.5364 | −0.3392 | −0.3642, −0.3141 | 6,205 | 47 |
| Maroon – Sapphire | 0.0213 | 0.3483 | −0.3270 | −0.3396, −0.3153 | 1,209 | 16 |
| Lapis – LapisSplinter | 0.0564 | 0.3518 | −0.2954 | −0.3270, −0.2713 | 10,803 | 82 |

Dyad-weighted mean **−0.405 [−0.442, −0.367]**; 12 of 12 below expectation. The largest dyad now contributes 22.8% of the bins, so the pooled value is no longer one dyad's statistic.

### Sustained associations — the state, not the dyad

| Dyad | Observed | Expected | Deficit | 95% interval | Bins | Events |
|---|---:|---:|---:|---|---:|---:|
| Copper – Lilac | 0.4719 | 0.5048 | −0.0329 | −0.3306, −0.0242 | 93,763 | 8 |
| Maroon – Sapphire | 0.2085 | 0.3544 | −0.1458 | — | 1,614 | 1 |
| Bronze – Magenta | 0.0576 | 0.4129 | −0.3552 | — | 2,077 | 2 |

### The reframing this forces

**Copper–Lilac is not an outlier dyad. It is an ordinary dyad that entered an extraordinary state.**

Its discrete encounters give −0.339, which ranks tenth of twelve — inside the pack. Its sustained-association value is −0.033, near composition expectation. The earlier pooled figure of −0.048 was simply its fused period dominating its own record: Copper–Lilac holds 93,763 of its 99,968 bins in the fused stratum, and 85.9% of every dyad's bins pooled together.

The general claim is therefore stronger and no longer rests on one case:

> Groups that encounter each other barely mix, and they barely mix everywhere — the deficit is large and negative in all twelve measurable dyads. Mixing to near composition expectation is a property of sustained association, a state five dyads entered and only one carried to near-completion.

### Weighting, demonstrated

| Stratum | Dyad-weighted | Bin-weighted | Largest dyad's share of bins |
|---|---:|---:|---:|
| Discrete encounters | −0.405 | −0.421 | 22.8% (Lapis–LapisSplinter) |
| Sustained associations | −0.178 | −0.057 | 99.4% (Copper–Lilac) |
| All encounters pooled | −0.369 | −0.198 | 85.9% (Copper–Lilac) |

Pooling by bin instead of by dyad halves the estimate. The choice is not cosmetic.

### Radius behaves consistently across dyads

Using the 5 m eligible bin set as a common denominator, so the 2 m rate is not conditioned on a 2 m contact existing:

| Stratum | 2 m / 5 m contact-rate ratio |
|---|---|
| Discrete encounters, 12 dyads | 0.086 – 0.194 (median ≈ 0.14) |
| Copper–Lilac, sustained | 0.262 |
| Maroon–Sapphire, sustained | 0.141 |

Across every measurable dyad, tightening the radius from 5 m to 2 m removes roughly the same proportion of contact. Dyads differ in *how much* they mix, not in how their mixing is distributed across distance. Only sustained fusion shifts that ratio.

---

## What Phases 0–3 did not settle

- **Fine-scale coverage is 12 of 68 structural dyads.** The mixing gradient rests on the dyads for which an upstream 2-minute product exists. Extending it means regenerating that product for more dyads in `New project`, which was not attempted here.
- **Eligible bins are bins with at least one within-radius edge** in a shared cluster. A fully unconditioned denominator needs the 2-minute GPS join. The radius-nested comparison bounds how much this matters; it does not remove the conditioning.
- **The composition-preserving shuffle is the reference used throughout.** The within-origin reference used in the Copper–Lilac analysis is a different estimator and the two are not interchangeable.
- **Demographic group size is still absent.** Phase 2 and 3 use observed collar counts, correctly labeled as such. Phase 4 must not reuse them as group size.
- **Nothing here is a model.** These are repaired descriptive quantities and their uncertainty. Phase 4 fits the two stages.

## Phase 4 entry conditions, now met

The two-stage samples exist and are auditable: 73,353 dyad-days with explicit states for Stage 1, and 1,705 events with three separated durations plus dyad-level mixing for Stage 2. Phase 4 can fit occurrence on the opportunity table and mixing on the event table, with a shared group effect for unordered pairs, encounters as the clustering unit, and discrete encounters and sustained associations kept in separate strata.
