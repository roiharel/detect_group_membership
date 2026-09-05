# Clustering rule, GPS source, and the two pipelines

## A decision record, 5 September 2026

Prompted by a parallel `New project` thread that built an independent hourly-grouping
table, six questions were opened at once: whether to tighten the adaptive clip, whether
HDBSCAN should replace the current rule, whether to cluster on all fixes instead of hourly
medians, which GPS file is canonical, whether the two pipelines agree, and whether the
"three-way partial merge" case is worth reporting. All six are settled or bounded below.

Everything in sections 1–3 is measured on the canonical GPS file
(`data/gps_v1_canonical_20260905.parquet`, see section 4) and **replicated in six windows
spanning the whole record**, because collar coverage rises roughly 120-fold across it and
a kNN-scaled rule is exactly the kind that can misbehave when an animal's 2nd-nearest
neighbour is kilometres away:

| window | animals in window | median animals/hour scored | adaptive silhouette | fine-scale silhouette |
|---|---:|---:|---:|---:|
| 2024-03 → 06 | **11** | 4–5 | **0.509** | 0.489 |
| 2024-09 → 12 | 97 | 9 | 0.888 | 0.926 |
| 2025-03 → 06 | 83 | 10 | 0.941 | 0.955 |
| 2025-09 → 12 | 216 | 19 | 0.874 | **0.957** |
| 2026-01 → 05 | 297 | 22 | 0.931 | 0.947 |
| 2026-06 → 09 | 307 | 24 | 0.938 | 0.950 |

**The classifier is sound from about 2024-09 onward and is a different animal before
that.** In the first quarter of the record only 11 animals were collared, and there the
adaptive rule separates poorly (silhouette 0.509 against 0.87–0.94 everywhere else) and
the fine-scale input is the only place in the whole record where it does not help. That
period sits below the coverage floor the paper already excludes — the seasonality work
starts its well-observed window at 2024-12-01 — but it should be excluded *because it was
measured*, not by assumption.

The first version of this comparison ran on a stale local copy in a single window; both
faults are corrected, and none of the conclusions moved.

Scripts: `analyze_clustering_options.py`, `build_clustering_disagreement_figure.py`,
`audit_gps_sources.py`, `crosswalk_state_taxonomies.py`,
`analyze_three_way_partial_merge.py`.

---

## 1. The upper clip: leave it at 900 m — 600 m is a no-op

The adaptive rule sets each animal's scale from its 2nd-nearest neighbour and clips to
[120, 900] m, then links at `1.65 x min(scale_i, scale_j)`. Measured on the canonical file
over 120 sampled hours in each of **two independent windows** — the second deliberately
chosen for much better collar coverage:

| quantity | win 1 (2025-03→06, 83 animals) | win 2 (2026-01→05, **297 animals**) |
|---|---:|---:|
| median edge threshold actually used | **198 m** | **198 m** |
| dyad-thresholds above 600 m | 0.30% | 3.02% |
| dyad-thresholds pinned at the 900 m ceiling | 0.23% | 2.99% |
| dyad-hours where 600 m would have cut an edge 900 m kept | 4 | 24 |
| ARI, 900 m vs 600 m | **1.00** | **1.00** |
| group-pair call rate, 900 m / 600 m | 0.0526 / 0.0526 | 0.0412 / 0.0412 |

The ceiling binds ten times more often at high coverage (3.0% of dyad-thresholds versus
0.3%) — more collars means more animals whose 2nd-nearest neighbour is far away — but it
still **changes no partition**: ARI is 1.00 in both windows and the call rates are
identical to four decimals.

**With one exception, at the sparse end.** In 2024-03→06, with 11 collared animals, the
clip is *not* a no-op: the group-pair call rate is **0.0166 at 900 m against 0.0083 at
600 m — a factor of two** — and the median cluster count moves 4.5 to 5.0. When almost
nobody is collared, the 900 m ceiling is what links animals that a tighter rule would
leave apart, and it doubles the apparent encounter rate. This is the sparse-collar
artefact in its purest form. It does not change the recommendation for the analysed
period, but it is a reason the early record must stay excluded rather than merely
down-weighted. **The ceiling is not what makes this rule permissive** — the
`factor x kNN scale` term is, and its median sits at 198 m in both windows. Changing 900
to 600 would force a full re-derivation and change nothing. If the rule is ever to be
tightened, `ADAPT_FACTOR` is the knob.

## 2. HDBSCAN: ruled out, and not because it was handicapped

Both windows, `k` = median clusters per hour:

| method | input | k₁ | k₂ | silh₁ | silh₂ | stab₁ | stab₂ | pers₁ | pers₂ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| adaptive_900 *(canonical)* | hourly median | 10 | 22 | 0.941 | 0.931 | 1.00 | 0.996 | 0.978 | 0.983 |
| adaptive_600_2m | fine-scale 2 min | 10 | 21 | 0.955 | 0.947 | 1.00 | 1.00 | 1.000 | 0.989 |
| dbscan_300 *(New project)* | hourly median | 10 | 21 | 0.949 | 0.944 | 1.00 | 1.00 | 0.990 | 0.991 |
| dbscan_300_2m | fine-scale 2 min | 10 | 21 | 0.956 | 0.950 | 1.00 | 1.00 | 0.998 | 0.994 |
| dbscan_500 | hourly median | 9 | 20 | 0.960 | 0.951 | 1.00 | 1.00 | 1.000 | 0.995 |
| dbscan_100 | hourly median | 12 | 25 | 0.859 | 0.869 | 0.978 | 0.973 | 0.938 | 0.942 |
| **hdbscan** | hourly median | **25** | **75** | **0.361** | **0.274** | **0.563** | **0.465** | **0.431** | **0.342** |
| hdbscan_2m | fine-scale 2 min | 15 | 44 | 0.782 | 0.713 | 0.825 | 0.762 | 0.739 | 0.662 |

**The ranking replicates exactly, and HDBSCAN gets worse as coverage improves.** At 297
collared animals it returns **75 clusters an hour** where every other method returns
20–25, and its cluster count tracks the collar count at **r = 0.938** — it is, almost
literally, measuring the collars rather than the baboons. Its ARI against the adaptive
rule falls from 0.418 to 0.343. Across all six windows its cluster count is 4, 11, 25, 55,
75, 71 as coverage goes 11 → 307 animals, and its k~n correlation is 0.89–0.97 in every
window with 97 animals or more.

**HDBSCAN is only competitive in the one window where nothing is:** at 11 collared animals
it has the best silhouette of any method (0.631 against the adaptive rule's 0.509). With
four animals in an hour there is no meaningful partition to find, and that is a statement
about the coverage, not a case for the method.

**One caution about the coverage criterion itself.** `k~n` rises for *every* method in the
better-covered window (0.72–0.94 against 0.37–0.60), so it is a within-window comparator
only; cluster counts must never be compared across periods with different collar counts.

Criteria, none of which need a ground truth: *silhouette* of the hourly partition scored
against the fine-scale (2-minute) distance matrix; *stability* as ARI against the
partition from a random half of the hour's fixes; *persistence* as ARI against the next
hour on animals present in both; *k~n bias* as the correlation between cluster count and
collar count.

`allow_single_cluster` makes **no difference at all** (ARI 1.00 between settings), so the
earlier agreement run's choice was not the cause. HDBSCAN genuinely shatters groups: in
the worst hour, 68 animals, the adaptive rule finds 14 clusters and HDBSCAN finds **35**.
`docs/framing_2026-09-04/clustering_disagreement.html` draws the three worst
disagreements with identical points in every column and only the outlines changing — the
pattern is consistent, HDBSCAN reads the density gradient *inside* a foraging group as a
boundary. Right for finding substructure, wrong for deciding who is together.

**A bug worth remembering.** `sklearn.cluster.HDBSCAN` builds its mutual-reachability
graph **in place** on a precomputed distance matrix, silently corrupting the caller's
distances including the diagonal. The first run of `analyze_clustering_options.py` was
wrong because of it; both label functions now copy. `analyze_clustering_method_agreement.py`
fits on coordinates rather than a precomputed matrix, so its published numbers are
unaffected.

## 3. The rule barely matters; the input does

Median ARI between **adaptive_900 and dbscan_300 is 1.00** in window 1 and **0.993** in
window 2 — the canonical pipeline and the `New project` pipeline produce essentially
identical partitions despite looking entirely different on paper. adaptive vs dbscan_500
is 0.969.

What does move the criteria is the **input**. Collars fix every ~2 minutes, up to 30 times
an hour, and the fixes are not synchronised across animals; there is a median of **30
shared 2-minute bins per dyad per hour**, all of which an hourly median discards. Feeding
the same rule the fine-scale dyadic distance matrix instead improves every criterion:

| criterion, adaptive rule | median input | fine-scale input |
|---|---:|---:|
| silhouette, win 1 / win 2 | 0.941 / 0.931 | **0.955 / 0.947** |
| stability, win 1 / win 2 | 1.00 / 0.996 | **1.00 / 1.00** |
| persistence, win 1 / win 2 | 0.978 / 0.983 | **1.000 / 0.989** |
| collar-count bias, win 1 / win 2 | 0.54 / 0.81 | **0.45 / 0.73** |

The improvement is in the same direction on every criterion in both windows, and the same
holds for DBSCAN 300 m. Note the median shared 2-minute bins per dyad falls from 30 to
14.5 between the windows — more collars means more dyads that overlap only partly — so the
fine-scale advantage is not an artefact of a period with unusually complete overlap.

**Decision: banked, not applied.** Switching the input invalidates every
`association_event_id` and every product derived from them, including all three current
figures. It should happen at the next full rebuild, together with item 4 — not as a
standalone change.

## 4. The GPS source: settled. One canonical file, updated in place

**Resolved.** `New project/outputs/canonical_robust_hourly_membership_shared_full_20260722/canonical_hourly_membership_metadata.json`
records the frozen export's own source:

```
"gps_file": "\\10.126.19.90\EAS_shared\baboon\working\data\processed\2025\gps\v1_cleaned\gps_v1.parquet"
```

That is the same path the parallel hourly-grouping pipeline uses. **There was never a
second canonical source.** The file is updated in place, so every apparently competing
source was a different vintage of one file, and the confusion was entirely about dates.

The share is reachable from PowerShell (Git Bash cannot resolve the UNC here, and neither
can a Python process launched from it — hence `census_shared_gps.py`, which must be run
via `powershell python census_shared_gps.py`).

| | rows | animals | groups | span | modified |
|---|---:|---:|---:|---|---|
| **share, current** | **30,143,804** | **392** | 27 | 2024-03-01 → **2026-09-05** | 2026-09-05 09:32 |
| local `network_v1_cleaned_gps_v1.parquet` | 24,153,071 | 372 | 27 | 2024-03-01 → 2026-06-18 | 2026-06-18 09:03 |

The share has **5,990,733 more fixes, 20 more animals, and runs 79 days later**.

**The frozen export is right and the local raw copy is the outlier.** Probe animals
against the current canonical file:

| animal | canonical share | frozen export | local raw copy |
|---|---|---|---|
| 25AA07_4S5T | 2026-03-16 **06:58** | 2026-03-16 **07:00** ✓ | 2026-03-23 13:26 ✗ |
| 24AA01_5O8B | 2026-07-10 **05:56** | 2026-07-10 **06:00** ✓ | 2026-06-17 00:00 ✗ |
| 24AC17_7J8K | 2025-08-29 16:00 | 2025-08-30 06:00 | 2025-08-29 16:00 |

Our frozen export matches the canonical file to the hour. The local copy holds **2,915
fixes for 25AA07_4S5T that the canonical cleaning has since retracted** — which is why the
other thread saw this animal move in and out of the disperser set. The discrepancy was
never a pipeline disagreement; it was one stale local copy.

**Two consequences.**

*(a) The clustering comparison in sections 1–3 was re-run on the canonical copy, and
every number is identical.* It had first been measured on the stale local file. Repointed
at `data/gps_v1_canonical_20260905.parquet` (30,143,804 fixes, 392 animals, verified
against the share) the March–June 2025 window yields **2,036,598 fixes and 83 animals as
before, and all four criteria agree to three decimals for all eleven methods**. The
retraction was confined to 2026, so sections 1–3 stand as measured. The local copy is
also byte-equivalent in content to the share: same rows, animals, groups, span, and
probe-animal endpoints to the second.

*(b) The frozen export is six weeks stale.* It was built 2026-07-22; the share now carries
data to 2026-09-05 and **20 animals that the frozen cohort has never seen**. That is a
decision, not a bug — the cohort was frozen deliberately — but it should be a stated
decision with a date, not an accident of when someone last ran the builder.

**Naming rule adopted:** any local copy of a mutable upstream file is stored with the
upstream mtime in its name (`gps_v1_canonical_20260905.parquet`). An undated copy of a
file that changes in place is how this went wrong.

**One incidental clarification.** The builder's metadata records `"fixed_eps_m": 500.0`
alongside `"method": "adaptive"` — the 500 m DBSCAN radius is an available-but-unused
option that the run metadata logs regardless. That is the source of the "the manuscript
uses DBSCAN eps = 500 m" note in `analyze_clustering_method_agreement.py`. Nothing in this
project's canonical chain has ever used it.

## 5. The two pipelines agree on 84.7% of animal-hours — and the disagreement is not labelling

`crosswalk_state_taxonomies.py` joins the two tables on (animal, hour): **1,075,208
animal-hours, 342 animals**, of which 1,030,314 are observed here. Concordance under the
expected mapping is **84.7%**.

**The single-animal boundary is independently confirmed.** Their
`DISPERSAL_WITH_OTHER_GROUP` — the state created because *"5 Lilac + 1 Jade is not a
between-group merge but a dispersal event"* — covers 13,417 observed animal-hours, and the
animal is **alone on its own side in 100% of them**. 95.6% of those rows are this
project's `mixed_without_origin_unit`. Two pipelines, two GPS vintages, two clustering
rules, same boundary. That is as strong as this kind of confirmation gets, and it retires
any doubt about stripping the 639 single-animal dyad-days out of the encounter network.

**One definitional gap remains.** 16,404 observed rows have the animal as the only member
of its group in a mixed cluster; they label 13,417 of them dispersal but **2,498
BETWEEN_GROUP_MERGE**. The reason is that their state is a property of the *cluster*
(2+ animals from 2+ groups anywhere in it) while the restriction applied here is a
property of the *dyad* (2+ on both sides of the pair). A cluster of Lilac:1, Jade:5,
Copper:3 is a merge by their rule, and the lone Lilac animal inherits that label. Both
rules are defensible; they are not the same rule, and 2,498 rows sit in the gap.

**The residual 15% is not a labelling artefact.** Origin labels agree on 99.9% of joined
rows, and the disagreeing rows are no more likely to have `dynamic_social_unit != origin_group`
than the agreeing ones (1.6% vs 1.5% for `other → FULL_ORIGIN_GROUP`; 0.0% vs 0.6% for
`mixed_with_origin_present → FULL_ORIGIN_GROUP`). Since the clustering step itself agrees
at ARI 1.00, the disagreement must come from **everything around the clustering step** —
this project's 60-minute support radius, `is_local_2h_supported`, the overnight
carry-forward, and the 500 m sparse-origin veto, none of which the other pipeline has.
That is a specific, bounded place to look, and it is where the next investigation belongs.

## 6. Three-way partial merges exist, but most are one animal

`analyze_three_way_partial_merge.py`, on the frozen export: same hour, both units with
more than 5 animals observed, a cluster of 2+ A and no B, a cluster of 2+ B and no A, and
a cluster with at least one of each.

| | this project | New project (reported) |
|---|---:|---:|
| hour-cases | **534** | 207 |
| events (consecutive hours collapsed) | **166** | 137 |
| unit pairs | **19** | 12 |
| median event length | 2 h | – |
| longest event | 14 h | – |

**The qualification matters more than the count. The mixed cluster holds a single animal
from one side in 56.1% of cases, and the median is 1.** So the majority of these are not
group-level partial merges at all — they are the dispersal channel again, seen from a
third angle. And the signal is concentrated: LapisSplinter + Lilac accounts for 272 of the
534 hours, Copper + LapisSplinter 108, Copper + Lilac 81 — three dyads carry 86% of it,
two of them involving **LapisSplinter, a unit already established not to be a documented
fission** (its animals were collared after any split, see the cohort-closure amendment).

**Decision: do not promote this to a headline.** It is worth a sentence in the axis-A
discussion, restricted to the 44% of cases where both sides contribute two or more
animals, and it should not be presented as a population-level phenomenon on the strength
of three dyads.

---

## What to do next, in order

1. ~~Resolve the GPS vintage.~~ **Done.** One canonical file,
   `\10.126.19.90\EAS_shared\...1_cleaned\gps_v1.parquet`, updated in place;
   copied locally as `data/gps_v1_canonical_20260905.parquet`. The frozen export is
   faithful to it; the local `network_v1_cleaned_gps_v1.parquet` is a stale divergent copy
   and **should not be used again**. `analyze_clustering_options.py` has been repointed,
   and `build_clustering_disagreement_figure.py` inherits the path from it.
2. ~~Re-run sections 1–3 on the canonical copy.~~ **Done, unchanged** —
   `analyze_clustering_options.py` now reads the dated local copy and every criterion
   reproduces to three decimals.
3. **Decide the cohort freeze explicitly.** The share now carries 392 animals and data to
   2026-09-05; the frozen cohort is 350 animals to 2026-07-22. Either re-freeze at a
   stated date or record 2026-07-22 as the deliberate cut. Right now it is neither.
4. ~~Validate sections 1–3 in a second window.~~ **Done** — 2026-01→05, 297 animals,
   7.7M fixes. Ranking identical; HDBSCAN worse at higher coverage; the clip still a
   no-op. Run as `python analyze_clustering_options.py --start ... --end ... --tag ...`.
5. **Investigate the residual 15%.** The clustering agrees at ARI 1.00 and the labels agree
   at 99.9%, so the disagreement is in the support / carry-forward / sparse-veto layer that
   only this pipeline has. Crosswalk one well-covered week with `is_local_2h_supported`,
   `is_carried_night` and the veto flags included as predictors of disagreement.
6. **At the next rebuild, switch the clustering input to the fine-scale 2-minute matrix**
   and keep the adaptive rule and its clip as they are. Now supported in two windows.
7. Adopt the `DISPERSAL_WITH_OTHER_GROUP` distinction explicitly in this project's
   vocabulary, and decide cluster-level versus dyad-level for the 2,498 ambiguous rows.
8. Leave HDBSCAN out. Keep `clustering_disagreement.html` as the record of why.

## Operating rules adopted from this

* **Date every local copy of a mutable upstream file.** The whole confusion came from an
  undated copy of a file that changes in place.
* **The UNC share is invisible to Git Bash and to Python launched from it.** Anything that
  must read it runs via `powershell python <script>`; anything else reads the dated local
  copy.
* **Copy from the share with `robocopy /Z /J`, not `Copy-Item`.** `Copy-Item` fails part
  way through a 1.34 GB transfer with an SMB `Invalid Signature` IOException. Also note
  robocopy's exit code 1 means *files were copied* — a wrapper that treats non-zero as
  failure will report a successful copy as failed.
* **Never pass a distance matrix to `sklearn.cluster.HDBSCAN` without copying it.** It
  mutates the input in place.
