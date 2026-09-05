# Gap 7.3 closed, two corrections, and how to treat rare patterns

Prepared 4 September 2026. Amends [AXIS_B_FISSION_RESULTS.md](AXIS_B_FISSION_RESULTS.md) and [SEASONALITY_AND_TRENDS.md](SEASONALITY_AND_TRENDS.md).

Scripts: [rederive_axis_b_on_frozen_cohort.py](../../rederive_axis_b_on_frozen_cohort.py), [analyze_seasonality_depth_duration.py](../../analyze_seasonality_depth_duration.py).
Outputs: `outputs/general_structure_2026_09/phase4d_axis_b_frozen/`, `.../phase4c_seasonality/depth_duration/`.

Four things: the cohort gap is closed and it caught two errors; the seasonality claim was overstated and is now properly tested; the LapisSplinter story is not what I said it was; and the editorial rule for rare patterns needs changing.

---

## 1. Gap 7.3 — closed for axis B

### What the gap was

Axis B ran on the legacy hourly source (1,703,133 rows, 21 units, filtered to 2025-01-01 onward). Axes A and C run on the frozen export (1,924,104 rows, 26 units, from 2024-03-01). Result 2 puts all three axes in one table, so it was mixing denominators, and no axis-B number could be pooled with A or C.

### How it was closed

The weekly within-group network metrics were rebuilt from the frozen narrow export, matching the legacy `combined_1100_1600` definition exactly: two sampling hours per day (11:00 and 16:00 UTC, so 14 timestamps in a full week — verified against the legacy file's own `n_timestamps`), observed positions only, cluster identity from `association_event_id`, and modularity of a greedy-modularity partition of the weighted co-clustering graph. 1,268 group-weeks in 23 units, against the legacy 1,418 in 26.

### What survives, and what does not

| | Frozen | Legacy | |
|---|---:|---:|---|
| Group-weeks / units | 1,268 / 23 | 1,418 / 26 | |
| Spearman(modularity, collars) | **0.222** | 0.225 | holds |
| **Spearman(split fraction, collars)** | **0.163** | **0.489** | **does not hold** |
| ICC between-unit | **0.169** | 0.148 | holds |
| ICC group-weeks / units | 368 / 12 | 314 / 12 | |
| Collar term *p* in the ICC model | 0.641 | 0.428 | holds (null both ways) |
| Lilac, % well-covered weeks modular | **55.3%** | 53.8% | holds |
| Lilac, max modularity | **0.449** | 0.458 | holds |
| **Chartreuse, % weeks modular** | **14.7%** (34 wks) | **45.0%** (20 wks) | **does not hold** |
| Lilac lag-1 autocorrelation | 0.417 | 0.413 | holds |
| Chartreuse lag-1 autocorrelation | 0.399 | 0.331 | holds |
| Purple / RubyRunners / Turquoise | 0.0% | 0.0% | holds |

**Correction 1 — I overstated the split-detection confound.** I wrote that "half the variation in split detection is collar count," from a Spearman of 0.489. On the frozen source it is **0.163**. The detectability ladder is correspondingly much flatter: split-weeks run 38.3% → 58.6% across the collar bands on the frozen source, against 28.7% → 78.1% on the legacy one. Split detection *is* coverage-sensitive, but not nearly to the degree I claimed, and the claim that no unadjusted split rate can survive is too strong. The modularity/collar relationship, by contrast, is unchanged (0.222 vs 0.225), so the ≥12-collar threshold for modularity stands.

**Correction 2 — Chartreuse's modularity is cohort-sensitive; Lilac's is not.** On the frozen source, Chartreuse is modular in 14.7% of its 34 well-covered weeks rather than 45.0% of 20. Its *stickiness* survives (lag-1 0.399, higher than legacy) and its maximum modularity is similar (0.071 vs 0.065), so it still enters structured phases — just far less often than the legacy source implied. **Lilac replicates almost exactly** (55.3% vs 53.8%, max 0.449 vs 0.458, lag-1 0.417 vs 0.413). So the robust statement is:

> Lilac is the one unit in this dataset with a well-replicated, recurrent modular phase. Chartreuse shows sticky structure at a much lower rate than first reported. Purple, RubyRunners and Turquoise are measured and cohesive.

**What is unaffected.** The ICC — 0.169 frozen against 0.148 legacy, both on 12 units — so *modularity is a state, not a group property* holds on the frozen cohort, which was the finding that mattered most. And the collar term stays null above the threshold on both sources, which is what licenses the well-covered subset.

Everything else in [AXIS_B_FISSION_RESULTS.md](AXIS_B_FISSION_RESULTS.md) Part 3 (the NDVI models) still runs on the legacy weekly file and has **not** been re-derived. That is now the remaining piece of 7.3, and it is smaller: the environmental conclusion was a null in six of eight specifications, and nulls are less likely to be cohort-artefacts than positives.

---

## 2. LapisSplinter is not a documented fission — correction

I twice used LapisSplinter as evidence that "at least one group in this population actually divided," and built the cross-axis prediction on it. That was wrong, and the frozen export shows it plainly:

| Animal | `origin_group` | First appears in LapisSplinter | First observed anywhere |
|---|---|---|---|
| 24AE10_8X9Y | **LapisSplinter** | 2024-11-11 | 2024-11-11 |
| 24AE09_6V7W | **LapisSplinter** | 2024-11-11 | 2024-11-11 |
| 24AE11_021A | **LapisSplinter** | 2024-11-11 | 2024-11-11 |
| 24AE12_3B4C | **LapisSplinter** | 2024-11-11 | 2024-11-11 |
| 24AE14_7F8G | **LapisSplinter** | 2024-11-11 | 2024-11-11 |
| 24AE17_3L4M | **LapisSplinter** | 2024-11-11 | 2024-11-11 |
| 24AC02_3C4D | Lapis | 2025-04-28 | 2024-07-13 |

Six of the seven animals carry **`origin_group = LapisSplinter`** and were all first observed on the same day — the day they were collared. Only one animal actually transferred from Lapis, 289 days into its own tracking record. Lapis's observed collars jump from 6 to 14 on that same date.

**So the unit name records field knowledge of a split that happened before, or outside, the tracking record.** The dataset documents that a splinter unit *exists* and that it is unusually permeable with its parent; it does not document the fission. Lapis's 18 pre-2024-11-11 group-weeks carry 4–6 collars — below the 12-collar detectability threshold — so the question "did Lapis show modular phases before it divided?" is **not testable here**. Not negative: untestable.

What still stands, and is worth reporting: **Lapis–LapisSplinter is the most permeable of the twelve measurable dyads** (deficit −0.295), LapisSplinter appears in two of the five dyads that ever entered sustained association, and it receives the one genuine inter-unit transfer in that lineage. A parent and its splinter mix more than any other pair. That is an observation about a splinter unit's permeability, not about fission dynamics, and the earlier cross-axis prediction should be restated accordingly.

---

## 3. Seasonality — I overstated it, and here is the properly scoped version

You asked whether "opportunity is seasonal, the crossing decision is not" holds for all social aspects. It did not, as stated: every seasonal test I had run used an **occurrence** response — does an encounter happen, does a split happen, does an animal leave. **Depth and duration had not been tested on any axis.** That was a real gap, and there was a specific reason to expect a difference: on two axes duration is the thing that carries depth, so season could reach depth through duration even with occurrence flat.

Now tested, same annual harmonic, same clustering:

| Axis | Response | n | Clusters | Amplitude | Joint *p* |
|---|---|---:|---:|---:|---:|
| A | encounter duration (log structural span) | 1,480 | 66 | 0.084 | 0.589 |
| A | mixing depth given an encounter | 716 | 12 | 0.401 | 0.320 |
| A | mixing depth, duration term dropped | 716 | 12 | 0.399 | **0.084** |
| B | split duration (log persistent hours) | 857 | 16 | 0.113 | 0.708 |
| C | excursion duration (log away-nights) | 322 | 18 | 0.104 | 0.778 |
| C | depth: reached another unit | 322 | 18 | 0.514 | 0.234 |
| C | depth: reached settlement | 322 | 18 | 0.229 | 0.826 |

**No season in depth or duration either.** So the statement survives — but it is now a claim about seven responses across three axes rather than three, and the depth and duration samples are one to two orders of magnitude smaller than the occurrence tables, so these are weaker nulls.

One thing is worth naming rather than burying. Axis A's mixing depth reaches *p* = 0.084 when the duration term is dropped and 0.320 when it is included — and duration itself carries a large effect on depth (+1.286 [1.127, 1.445] per SD). That is exactly the mechanism to be suspicious of: if season acts on depth anywhere, it acts through duration. On 716 events in 12 dyads there is not enough to resolve it. Report it as the open question it is, not as a null.

One more from these fits: axis C excursion duration carries a **trend** of +0.235 [0.033, 0.437] — excursions getting longer across the record. Given that later excursions are better observed, and that observation effort rises 121-fold, treat that as effort-suspect until it survives a balanced-panel check.

---

## 4. Rare patterns are the finding, not a problem — revised editorial rule

You are right, and this changes several verdicts. I had been applying one filter — *is this well enough covered to claim?* — and using it to suppress patterns that are simply rare. Those are two different things, and conflating them throws away the most informative part of a population-scale dataset.

### The distinction

**Measurement failure — do not report.** The number does not measure what it claims. It would not become true with more data of the same kind; it is wrong now.

- Modularity below 12 collars. The estimator returns one cluster whatever the animals are doing, so the "19 of 20 groups at zero" figure is a property of the algorithm.
- The 285-animal isolation inventory. About 76% is artefact; 1- and 2-collar contexts contribute 15,908 animal-hours of which **zero** survive a located-reference test.
- Any per-animal isolation ranking. The top animal by raw count has 2,880 isolated hours and zero supported ones, because it is the only collar in its group.
- A trend fitted on two clusters.
- Split extent as currently defined, which is predicted by collar count and nothing else.

**Rare but measured — report as documented variation.** The observation is sound; the state is uncommon. Suppressing it hides the range of the system.

| Pattern | Frequency in this sample | Why it matters |
|---|---|---|
| Sustained association | 16 of 1,705 encounters, 5 of 68 dyads | the only route to near-complete mixing observed anywhere |
| Full fusion | 1 dyad (Copper–Lilac) reaching −0.033 against expectation | the ceiling of between-group permeability |
| Recurrent modular phase | 1 of 11 well-covered units (Lilac, 55% of its weeks) | shows a group can hold internal community structure for weeks |
| Settlement excursion | 61 of 338 excursions, 16 of 91 animals | dispersal as an outcome of an ordinary trip |
| Alone-then-joined excursion | 11 of 338 (3%) | the classical float-then-transfer pathway does occur, rarely |
| Splinter-unit permeability | 1 lineage (Lapis–LapisSplinter, deficit −0.295) | the most permeable dyad measured |

### The argument to make in print

> This is 26 social units over 29 months in one landscape — a subset of a population and a subset of the environmental conditions that population experiences. Within that subset, every axis of boundary permeability shows states that only a few units enter: one dyad fuses almost completely, one group holds a recurrent modular phase, sixteen animals settle elsewhere. The frequencies are low, and they are reported as such. But the presence of a tail in a 26-unit sample is itself the result: **the population contains variation in social context, and the sample is too small to have exhausted it.** A larger set of groups, across a wider range of environmental conditions, would sample that tail better — it would not make it disappear.

That reframes the rare states from underpowered claims into the paper's evidence for variation, which is a stronger and more honest position than the coverage-forbids framing I had been using. It also puts the nulls in their place: the seasonal and environmental nulls say that *within the range this sample covers*, no annual or vegetation signal organises boundary crossing. They do not say no such signal exists across a wider environmental range.

### What changes in the verdict table

Three-way, not two-way:

- **report** — supported at population scale, with its coverage stated
- **report as variation** — rare but soundly measured; give the frequency, do not suppress
- **do not report** — measurement failure

Rows moving to *report as variation*: full fusion; Lilac's modular phase; settlement excursions; alone-then-joined excursions; splinter-unit permeability. Rows staying at *do not report*: modularity below 12 collars; the 285-animal isolation figure; per-animal isolation type; the two-cluster modularity trend; split extent as defined.

---

## 5. What remains

1. **Re-derive the axis-B NDVI models on the frozen weekly file.** The last piece of 7.3. The conclusion there was a null in six of eight specifications, so it is less likely to move than the two positives that did.
2. **Balanced-panel check on axis C's excursion-duration trend** (+0.235 [0.033, 0.437]), which is effort-suspect.
3. **Axis A mixing depth and season, with duration held out** — the *p* = 0.084 hint. Needs more fine-scale dyads (12 of 68 currently), which is gap 7.5.
4. **The fission-precedes-permeability question is untestable here** and should be stated as such rather than framed as a prediction the data can adjudicate. It needs a group tracked with 12+ collars *through* a division.
5. **Second seasonal harmonic** for bimodal rainfall, and the untested indices (NDWI, LSWI, EVI). The well-observed window holds 1.6 annual cycles, which limits what any harmonic can say.
