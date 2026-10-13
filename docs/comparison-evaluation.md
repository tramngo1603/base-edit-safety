# Mock comparison evaluator

This command checks evaluation arithmetic on mock candidates and outcomes. It accepts `data_kind: artificial`; other input kinds are unsupported.

Create an environment and install the package first; `PYTHONPATH` alone does not provide distribution metadata. From a full source checkout with Python 3.10 or newer:

```sh
python3 -m venv .venv-evaluation
.venv-evaluation/bin/python -m pip install .
```

Then run one command into a new directory:

```sh
.venv-evaluation/bin/python -I -m base_edit_safety.population.comparison_evaluation evaluate \
  --candidates tests/fixtures/comparison/evaluation/candidates.json \
  --labels tests/fixtures/comparison/evaluation/labels.json \
  --methods tests/fixtures/comparison/evaluation/methods.json \
  --protocol tests/fixtures/comparison/evaluation/protocol.json \
  --splits tests/fixtures/comparison/evaluation/splits.json \
  --output /tmp/bes-artificial-evaluation
```

Choose a different output path for a repeat; outputs are never overwritten. The output contains `evaluation.json` (output format version 3) and `receipt.json`. Input and split formats remain version 1. The result is deterministic for the same inputs and code; the receipt separately records actual elapsed evaluator time, timestamp, Python/package/source identity, command arguments, input file checksums and result checksum. `declared_method_runtime_seconds_total` is caller-supplied metadata for the entire method run, repeated in each stratum for context. The mock fixture supplies zero for this declared value.

## Input boundary

Five ordinary JSON files keep label-free grouping separate from outcomes:

| File | Required contents |
| --- | --- |
| candidates | Unique candidate IDs; canonical study, guide and family IDs; editor and assay; explicit eligibility boolean |
| labels | Exactly one of the shared four adjudicated states per candidate, an evidence reference and adjudication reason; untested candidates remain explicit |
| methods | Unique method ID, method-specific enumerated candidate IDs, scored candidate/score pairs, and declared total runtime seconds |
| protocol | Format version 1, `artificial` data kind, protocol ID, integer budget per editor/assay, sensitivity targets, bootstrap seed and replicate count, baseline method ID |
| splits | Saved assignments from label-free candidate metadata, including its digest, split seed and fractions |

The mock labels provide fixed expected outcomes for arithmetic tests. Use the evidence workflow to import and review measured observations.

Candidate IDs identify distinct assay units within the declared editor/assay context. Before real use, repeated measurements of the same biological candidate must be reconciled under a frozen rule, not assigned new IDs to inflate counts. The fixture has one adjudicated label per unit. Study, guide and family IDs are stable identities across studies, not reused local labels such as `guide_1`. Unknown linkage is not assumed to mean independence: missing and common unknown sentinels are rejected. Real pseudonymous linkage and context definitions remain pending.

Every scored candidate must belong to that method's enumeration and the full candidate universe. Enumerated candidates may lack a score. Missing scores are recorded as omissions, never filled with labels or an assumed low score. Scores, sensitivity targets and split fractions must use JSON integers or finite decimal strings. Fractional or exponent JSON numbers (including `0.5`, `1.0` and `5e-1`) are rejected locally because the shared JSON reader converts them to binary floats before exact evaluation. Write `"0.5"`, not `0.5`; integers such as `0` and `1` remain accepted. This rule also applies to direct Python calls. Runtime seconds are descriptive metadata and may still use finite fractional JSON numbers. Scores are higher-is-better, and must be comparable across guides within an editor/assay stratum. Exact decimal ordering is retained for the accepted strings and integers. For example, target `"0.500000000000000000000000000001"` with four positives requires three positives; its numeric JSON form fails before output rather than rounding into a two-positive target. Duplicate candidates, duplicate method scores, unknown IDs, incomplete label coverage, nonfinite scores and invalid outcomes fail before output.

## Grouped splits and leakage checks

All candidate rows enter grouping before eligibility or editor/assay filtering. A shared study **or** guide **or** family connects rows; connections are transitive. Thus an ineligible bridge or a related family in another study cannot allow relatives into different partitions. The connected component is both the split and resampling unit.

The split command does not read labels, scores or outcomes:

```sh
.venv-evaluation/bin/python -I -m base_edit_safety.population.comparison_evaluation split \
  --candidates candidate-metadata.json --seed 17 \
  --training-fraction 0.6 --development-fraction 0.2 \
  --output /tmp/bes-label-free-split
```

Components are ordered by a hash of the seed and component identity. Training and development receive the floor of their fraction times the component count; evaluation receives the remainder. Fractions allocate groups, not rows; imbalanced sizes are not corrected using outcomes. The evaluator reconstructs and verifies the entire manifest and metadata digest before exclusions. Arbitrary reassignment, partial membership and altered metadata fail. The mock fixture assigns all three groups to evaluation.

## Budget, denominators and ties

Each method gets the same **global budget within each editor/assay stratum**, not a separate budget for every guide. One nominated candidate assay costs one unit. The model uses equal unit costs. `assays_used` and matched costs count simulated testing slots, including currently untested candidates. `completed_assays` is null with reason `not_recorded_in_evaluation_inputs` because the inputs do not record completed assays. Every ranked selection consumes a unit, including candidates with inconclusive or untested labels. If too few scores exist, budget underfill is reported; no fallback nomination is invented.

Sort by descending score, then lexical candidate ID. A tie at the budget boundary selects only the remaining number of units; it does not admit an entire tie block for free. The identical order defines every precision–recall prefix and matched-sensitivity cost. Candidate IDs and tie rules must be fixed before outcomes are seen.

The eligible positive denominator is all `observed_positive` candidates in the evaluation partition and current editor/assay stratum, including positives absent from enumeration or scoring. Training/development and explicitly ineligible candidates do not enter metrics. An ineligible candidate still affects grouping.

- Enumeration sensitivity = enumerated eligible positives / all eligible positives.
- Recall at budget = selected eligible positives / all eligible positives.
- **Precision among selected conclusive records** = selected positives / (selected positives + selected tested-no-signal). The denominator and selected inconclusive/untested counts accompany it. This is not precision over all selected candidates or all performed assays; selective outcome availability can bias it. Unresolved labels are never inferred negatives.
- Conclusive outcome coverage = selected conclusive records / simulated testing slots. It appears next to conclusive precision at the budget and every prefix. Zero selections gives null coverage; selections with only unresolved outcomes give zero coverage and null precision. Read precision and coverage together.
- The precision–recall output gives every charged prefix, the same positive denominator, its conclusive count and unresolved counts. It includes the zero-selection prefix with undefined precision.

With no eligible positives, sensitivity and recall are null. With no selected conclusive records, precision is null. A null is not zero accuracy or a favorable result.

## Cluster uncertainty and matched sensitivity

For each bootstrap replicate, draw as many evaluation dependency groups as were present, with replacement. One set of group multiplicities is shared across every method and editor/assay stratum, preserving cross-stratum dependencies. Only groups with eligible evaluation candidates enter the draw. Each occurrence repeats the entire group's eligible candidates; copies are not deduplicated. Rerank that resampled population and apply the same fixed budget. Repeated candidates consume repeated assay units, including when the budget cuts through copies. This estimates variability of the selection procedure on resampled groups, not uncertainty for one fixed shortlist.

The output gives nearest-rank 2.5th and 97.5th percentile endpoints from the seeded replicates. It withholds an interval if the original metric is undefined/unreachable, fewer than three relevant independent groups exist, or any replicate has an undefined denominator/unreachable target. Counts of defined and undefined replicates remain visible; undefined replicates are not silently removed. For recall, enumeration and matched costs, relevant groups contain eligible positives. For selected-conclusive precision, they contain an originally selected conclusive candidate. Three groups is only an operational minimum; it does not establish reliable coverage or adequate scientific sample size. Highly connected studies can leave only one effective group, which is reported rather than treated as independent rows.

For each target sensitivity, require `ceil(target * all eligible positives)` positives and walk the complete ranked list until that count is reached, charging every preceding selection. Decimal target arithmetic is exact. This curve is separate from the fixed-budget result. A positive target with no positives is undefined; target zero requires zero assays. Missed/unscored positives can make a target unreachable. The comparison reports method cost minus baseline cost only if both reach the original target. A negative difference would mean fewer simulated units at the same validated-site sensitivity. Paired uncertainty uses the same group draws; an unreachable original comparison can never gain an interval from favorable resamples.

## Paired recall improvement at the fixed budget

Evaluator output `format_version: "3"` (input and split versions remain 1) adds
`recall_at_budget_comparisons` to each editor/assay stratum, keyed by each
non-baseline method. Each entry records `baseline_method`, `status`
(`comparable` or `not_comparable`), `budget`, `eligible_positives`,
`method_recall`, `baseline_recall`,
`recall_difference_method_minus_baseline` and an `uncertainty` interval using
the existing interval fields and reasons. Positive values favor the method.
The difference subtracts the selected positive counts first and then divides
by the shared per-stratum eligible-positive denominator; the two marginal
interval bounds are never subtracted. Counts use all eligible validated
positives, including missed and unscored ones. A zero denominator gives a null
difference, not a zero improvement. The interval reuses exactly the existing
paired transitive study/guide/family draws, repeated-copy budget charges and
nearest-rank 95% rules: zero stays inside an interval unless the actual
arithmetic excludes it. Success criteria are not part of this calculation. The tiny `tests/fixtures/comparison/recall-improvement/`
cases exercise a hand-calculated signed difference and the zero-denominator
rule through the same CLI using mock inputs.

## Hand-calculated fixture

There are eight eligible candidates and four positives across three independent mock studies. The proposed method enumerates three positives, leaves one of those unscored, and misses the fourth. Its order is inconclusive, positive, untested, positive, no-signal, no-signal. The baseline order is no-signal, positive, no-signal, positive, inconclusive, positive, untested, positive.

| Metric | Proposed | Baseline |
| --- | --- | --- |
| Enumeration sensitivity | 3/4 | 4/4 |
| Recall at four unit assays | 2/4 | 2/4 |
| Selected-conclusive precision at four | 2/2 | 2/4 |
| Conclusive outcome coverage at four | 2/4 | 4/4 |
| Completed assays | not recorded | not recorded |
| Inconclusive / untested selections at four | 1 / 1 | 0 / 0 |
| Assays for 50% sensitivity | 4 | 4 |
| Assays for 75% sensitivity | unreachable | 6 |
| Assays for 100% sensitivity | unreachable | 8 |

Focused regressions also cover transitive ineligible bridges, linked families across studies, exclusions, deterministic split reconstruction, multiplicity, cross-stratum paired draws, empty denominators, score ties, long decimal boundaries, original unreachable targets with favorable resamples, malformed inputs and output receipts.
