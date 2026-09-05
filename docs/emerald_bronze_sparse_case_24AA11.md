# Emerald/Bronze Sparse-Data Case: 24AA11_9A7D

## Question

Was the June-July 2024 event a true move of `24AA11_9A7D` from Emerald to Bronze, or a sparse-data artifact / small Emerald fragment moving near Bronze?

## Current canonical behavior

The current canonical local-2h membership output assigns `24AA11_9A7D` to dynamic social unit `Bronze` from:

- `2024-06-05 12:00` through `2024-07-10 13:00`

This is encoded as `dynamic_assignment = sustained_non_origin_association`, alternating between:

- `mixed_without_origin_unit`, usually `mixed_origin_group|Bronze:2-3;Emerald:1`
- `isolated`, for hours when the focal is not clustered with a multi-animal group

## Sparse proximity evidence

The sparse proximity diagnostics show this is not a clean Bronze-only move.

For the suspected period, the near-context summary was:

| window | category | hours | fraction |
|---|---:|---:|---:|
| strict +/-15 min | Bronze only <=500m | 78 | 0.331 |
| strict +/-15 min | Emerald only <=500m | 27 | 0.114 |
| strict +/-15 min | Both <=500m | 4 | 0.017 |
| strict +/-15 min | Neither <=500m | 127 | 0.538 |
| broad +/-60 min | Bronze only <=500m | 54 | 0.229 |
| broad +/-60 min | Emerald only <=500m | 55 | 0.233 |
| broad +/-60 min | Both <=500m | 28 | 0.119 |
| broad +/-60 min | Neither <=500m | 99 | 0.419 |

The other Emerald collar, `24AA10_4R7W`, is highly informative in the sparse-neighbor table:

- 132 sparse neighbor observations with `24AA11_9A7D` from `2024-06-06` to `2024-07-20`
- median distance to `24AA11_9A7D`: 37.9 m
- 75th percentile distance: 102.5 m
- all were within 900 m; most were within 250 m

However, `24AA10_4R7W` has no canonical membership rows during the key June 2024 interval; its canonical rows resume on `2024-07-26`. This means the canonical builder cannot fully represent the Emerald social context during the suspected Bronze episode.

## Interpretation

This does not look like strong evidence that `24AA11_9A7D` alone moved to Bronze.

The better interpretation is:

1. Around `2024-06-06` to `2024-06-14`, `24AA11_9A7D` was often close to Bronze collars, while sparse `24AA10_4R7W` fixes also show Emerald proximity. This looks like a small Emerald fragment / dyad in a Bronze-contact context, not necessarily a single-animal transfer.
2. From about `2024-06-15` onward, the sparse evidence favors `24AA11_9A7D` remaining close to `24AA10_4R7W` while Bronze is usually far.
3. The full `2024-06-05` to `2024-07-10` canonical Bronze assignment is therefore too strong.

## Recommendation

Do not treat this as a confident dynamic social-unit switch of `24AA11_9A7D` to Bronze.

Keep `dynamic_social_unit = Emerald` for this event, or flag the period as:

- `sparse_supported_origin_pair`
- `mixed_bronze_contact`
- `low_confidence_dynamic_switch`

The canonical table should preserve the Bronze-contact information in a separate context field, for example:

- `contact_context = Bronze`
- `contact_context_confidence = low/moderate`
- `sparse_origin_support = true`
- `sparse_origin_support_animals = 24AA10_4R7W`

This case argues for using sparse data as evidence against reassignment when a same-origin companion is repeatedly close, even if that companion does not satisfy the stricter canonical hourly membership rules.

