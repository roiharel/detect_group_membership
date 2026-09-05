# Merging isolation and dispersal into one individual axis

Prepared 4 September 2026. Amends [STABLE_GROUPS_FLUID_BOUNDARIES.md](STABLE_GROUPS_FLUID_BOUNDARIES.md), which treated isolation (axis C) and dispersal (axis D) as separate.

Script: [analyze_individual_axis_identifiability.py](../../analyze_individual_axis_identifiability.py).
Outputs: `outputs/general_structure_2026_09/phase4b_individual_axis/`.
Source: the frozen narrow export — 1,924,104 animal-hours, 350 animals, 26 origin groups, 2024-03-01 to 2026-07-22.

Two problems, one fix each. They are independent and both are settled below.

---

## Part 1 — The unification

### The state space

Isolation is not a phenomenon parallel to dispersal; it is one of the two things a dispersing animal can be doing. An animal's position relative to its origin group at any moment is one of four states, and they are exhaustive:

| State | Meaning | Depth |
|---|---|---:|
| `with_origin` | in a cluster containing its origin group — including a between-group merge, where the origin group is still present | 0 |
| `alone` | away from the origin group, with no group | 1 |
| `with_other` | away from the origin group, in a cluster with another unit | 2 |
| `unresolvable` | the data cannot distinguish — see Part 2 | — |

`unresolvable` is a first-class state, not a gap. That is the axis-A lesson (four explicit opportunity states, no state inferred from a missing row) applied to individuals.

### The unit: an excursion, not an event

The current products cut two incompatible inventories out of this space — 4,293 "isolated events" and 28-or-71 "disperser events/segments" — which is why they cannot be joined. The unit that joins them is the **excursion**: a run of consecutive away-from-origin nights, ending when the animal is back with its origin group.

An excursion carries a **depth** (the deepest state it reached) and an **outcome** (returned to origin, still away and censored, or settled). Dispersal is then not a separate category but a *terminal excursion*: one whose `with_other` phase passes the builder's own 7-night sustained threshold. "Disperser" stops being a type of animal and becomes a property of a trip.

### This is axis A's two-stage structure, at the individual scale

That parallel is what makes the merge worth doing rather than merely tidy.

| | Axis A — groups | Axis C — individuals |
|---|---|---|
| Stage 1 | does the dyad **meet**? | does the animal **leave**? |
| Stage 1 sample | 73,353 dyad-days, four states | 82,205 animal-nights, four states |
| Stage 2 | given a meeting, **how much do they mix**? | given a departure, **how deep does it go** — alone, joined, settled? |
| Depth reference | composition expectation | the origin-group baseline |
| Result so far | mixing at 4.6% of expectation; duration is the only robust predictor | see below |

### What the merged axis yields

Excursions built from resolvable nights, dominant-state nightly rule, no gap bridging:

| | Value |
|---|---:|
| Excursions | **338** |
| Animals | **91** |
| Origin groups | **18** |
| Median away-nights | 2 |
| Max away-nights | 500 |

Depth composition, and the finding that falls straight out of it:

| Depth | Excursions | Median away-nights |
|---|---:|---:|
| `alone_only` — left, joined nothing | 163 | **1** |
| `joined_only` — left, went straight to another unit | 164 | **4** |
| `alone_and_joined` — both, in one trip | 11 | 3 |

**Depth requires duration, on this axis too.** Excursions that only ever reach `alone` last a single night at the median; excursions that reach another group last four times as long; settlement requires seven. That is the same dissociation axis A found — where duration was the only robust predictor of mixing, and sustained association was the only route to near-complete mixing.

Settlement candidates (≥7 nights with another unit): **61 excursions in 16 animals** — which recovers the old dispersal sample almost exactly (17 animals), now nested inside a denominator instead of standing alone as a case series.

### One result that qualifies the premise

The premise behind merging was that a dispersing individual either moves alone or joins another group — two phases of one process. The data says they are mostly **not sequential**: only 11 of 338 excursions (3%) contain both states, rising to just 15 of 196 under generous gap bridging, and 361 of 3,534 (10%) under the permissive nightly rule. Most excursions are one thing or the other.

This does not undermine the merge; it vindicates it. On separate axes the question could not be asked at all. Merged, the first answer is that **floating is rarely the route to transfer in this population** — animals that end up in another group mostly appear there directly. That is a real, reportable finding, and it is only sayable once the two are on one axis.

---

## Part 2 — The sparse-collar artefact

### The problem, as raised

With few collars, an animal can be scored `isolated` by accident: when it is the only tracked animal in its group, or when a group with two collars splits and both sides look alone. Neither is isolation, and neither should feed a dispersal analysis.

### The rule

**A social state is assignable only when the data could have shown a different state.** For `alone`, that means the origin group has to be *locatable somewhere else* at the same moment.

An animal-hour is scored `alone` (or `with_other`) only if:

1. the focal animal's own position is **observed**, not carried; and
2. at least **2 collared animals of its origin group are observed in the same hour, co-clustered with each other**, in a cluster the focal animal is not in.

Condition 2 is the whole fix. It is a *located reference group*, not a partner count — "I can see your group over there, and you are not in it" rather than "I cannot see anyone with you."

### It removes exactly the cases it was built for

| Observation context | Animal-hours | Animals | Groups | Passing the rule |
|---|---:|---:|---:|---:|
| Only the focal animal observed | 15,613 | 52 | 21 | **0** (0.0%) |
| Exactly 2 animals of the group observed | 295 | 29 | 10 | **0** (0.0%) |

Both cases go to zero, by construction rather than by tuning: with one collar there is never a reference cluster, and with two collars apart neither side reaches size 2. A group-mate count would not do this — a plain "≥1 observed group-mate" test passes the two-collar split, which is precisely the context you flagged.

### What survives

| Rule (each adds to the one above) | Animal-hours | Share of raw | Animals | Groups |
|---|---:|---:|---:|---:|
| Raw isolated animal-hours | 48,775 | 100% | 307 | 26 |
| Focal position observed | 27,608 | 56.6% | 280 | 26 |
| + ≥1 observed group-mate anywhere | 11,995 | 24.6% | 274 | 23 |
| + ≥3 observed group-mates anywhere | 11,210 | 23.0% | 265 | 21 |
| **+ located reference cluster ≥2** | **11,668** | **23.9%** | **267** | **21** |
| + located reference cluster ≥3 | 11,061 | 22.7% | 263 | 21 |
| + located reference cluster ≥4 | 9,779 | 20.1% | 258 | 21 |

Three things to take from this table:

1. **About 76% of the raw isolation record is an observation artefact.** Nearly half of it (21,167 hours) is not even an observed position — the focal animal's location was carried forward.
2. **The rule is insensitive to its own threshold** — 23.9% / 22.7% / 20.1% for a reference cluster of 2 / 3 / 4. The result does not depend on where the cut is placed, which is what a support rule should look like.
3. **Five origin groups drop to zero supported isolated hours** — TrickyTeal, Green, Gold, Red, Mint. Their entire isolation record was artefact. Twenty-six groups can produce an isolation count; only 21 can support an isolation claim.

Survival is strongly structured by coverage, which is the point:

| Group | Observed isolated hours | Supported | Surviving |
|---|---:|---:|---:|
| Maroon | 1,033 | 1,033 | 100% |
| Lapis | 1,022 | 1,017 | 99.5% |
| Magenta | 370 | 364 | 98.4% |
| LapisSplinter | 972 | 913 | 93.9% |
| Chartreuse | 4,447 | 3,120 | 70.2% |
| Bronze | 3,612 | 763 | 21.1% |
| Copper | 2,157 | 389 | 18.0% |
| Emerald | 2,980 | 394 | 13.2% |
| SneakySilver | 1,005 | 16 | 1.6% |
| TrickyTeal | 1,683 | 0 | **0%** |
| Green | 1,058 | 0 | **0%** |

And the per-animal ranking reshuffles completely. `25AA22_2W3X` (TrickyTeal) has 2,880 raw isolated hours and **zero** supported ones — it is the only collar in its group, so it was scored isolated by construction for its entire record. `24AA10_4R7W` (Emerald) goes from 3,035 raw to 1. Any per-animal isolation statistic built on the raw field is measuring collar deployment.

### The resulting ledger

Four states, all 1,924,104 animal-hours, nothing dropped:

| State | Animal-hours | Share | Animals |
|---|---:|---:|---:|
| `with_origin` | 1,029,679 | 53.5% | 342 |
| `unresolvable` | 850,875 | 44.2% | 350 |
| `with_other` | 31,882 | 1.7% | 269 |
| `alone` | 11,668 | 0.6% | 267 |

The 44.2% unresolvable is honest, not new: 795,467 of those hours are carried positions — the existing `is_observed` field — and only 53,008 are observed hours without a located reference cluster. At nightly resolution, where a night aggregates roughly fifteen hours, resolvability rises to **92.0%** (dominant rule) or **95.3%** (permissive rule).

---

## Part 3 — Two definitional choices that must be pinned

Neither is a defect; both are places where a default would silently decide a result.

**The nightly-state rule moves the excursion count by an order of magnitude.**

| Nightly rule | Animal-nights resolvable | Excursions | Animals | Groups | Median away-nights |
|---|---:|---:|---:|---:|---:|
| `dominant` — the night takes its most frequent hourly state | 92.0% | 338 | 91 | 18 | 2 |
| `any_away` — the night is away if any hour was resolvably away | 95.3% | 3,541 | 303 | 21 | 1 |

These answer different questions. `dominant` asks where the animal spent the night; `any_away` asks whether it was ever away during it. 71% of `any_away` excursions are single nights, so that rule is largely counting brief daytime departures. **Recommendation: `dominant` as primary**, with `any_away` reported as the sensitivity — and, following axis A's three separated durations, report span, away-nights and deepest-state nights as separate columns rather than one "duration".

**Gap tolerance matters under the dominant rule and not under the permissive one.** Bridging unresolvable nights:

| Gap bridged (nights) | 0 | 1 | 2 | 3 | 7 | 14 |
|---|---:|---:|---:|---:|---:|---:|
| Excursions, `dominant` | 338 | 234 | 210 | 199 | 196 | 196 |
| Excursions, `any_away` | 3,541 | 3,538 | 3,534 | 3,534 | 3,534 | 3,534 |
| Animals, either rule | 91 / 303 | 91 / 303 | 91 / 303 | 91 / 303 | 91 / 303 | 91 / 303 |

Under `dominant`, bridging merges what were separate trips: 338 excursions collapse to 196, stabilising from gap 3 onward. The animal count never moves, on either rule or any gap — so the *sample* is stable and only the *event parsing* is sensitive. Report the ladder, take gap 3 as primary by analogy with axis A's 3-hour rule, and note that settlement candidates fall from 61 to 22 as trips merge.

---

## Part 4 — What this changes in the framing

**Three axes, not four.**

| Axis | What crosses | Units with the phenomenon | Denominator | Claim level |
|---|---|---|---|---|
| **A** Between-group interaction and mixture | whole groups, then individuals inside them | 68 of 317 dyads | 73,353 dyad-days, 4 states | population rate + within-dyad model |
| **B** Within-group fission and modularity | internal structure | 24 groups | absent | inventory; modularity coverage-bound |
| **C** Individual excursion — isolation *and* dispersal | one animal, to nothing or to another group | 91 animals, 18 groups, 338 excursions | 82,205 animal-nights, 4 states | **population rate + depth model** |

The merge is a straight upgrade at both ends. Axis D was a 17-animal case series that could not carry a rate; axis C's raw inventory was 285 animals of which roughly three-quarters was artefact. Merged and filtered, the axis has **91 animals, 18 groups and 338 excursions over an explicit 82,205-animal-night denominator** — a real population sample with a graded depth outcome, and the only axis besides A that can support a two-stage model.

**It also closes gap 7.1 for this axis.** The earlier document listed "opportunity tables for axes B, C and D" as the top blocker. The nightly four-state ledger *is* that table for the individual axis. Axis B still has none.

**The general claim gets a second independent leg.** "Depth requires duration, and duration is rare" was an axis-A statement resting on 16 sustained encounters in 3 dyads. It now holds at the individual scale too, on a different sample and a different measurement: alone-only excursions median 1 night, joined-only 4, settlement ≥7. Two scales, two mechanisms, one pattern.

**And one new reportable finding:** floating alone is rarely the route into another group — 3% of excursions contain both states, 10% under the permissive rule. Worth stating plainly, because the intuitive model of dispersal (leave, float, join) is not what these trips look like.

### Remaining work on axis C

1. **Reconcile the night definition.** This script uses a `window_start − 10h` calendar cut, giving 82,205 animal-nights; the canonical product uses a 16:00–06:00 biological night and reports 81,695. Before any rate is published these must agree.
2. **Pin the two rules** in Part 3 and report the ladder, not a single number.
3. **Terminal outcomes and censoring.** Excursion outcome (returned / settled / censored at the end of tracking) is not yet assigned, and censoring is heavy — the 500-night Chartreuse excursion is still open at the end of the record.
4. **Retire the old products.** `canonical_isolated_events.csv` (4,293 events) and both dispersal inventories should be labelled legacy. Do not report the 285-animal isolation figure again.
5. **Stage-1 model.** Whether an animal leaves, on the animal-night table, with the same design rules as axis A Stage 1: strictly preceding covariates, animal and group effects, and collar coverage entered as an observation covariate — where it should now be *expected* to matter, unlike at axis A Stage 1, and must be shown not to drive the result.
