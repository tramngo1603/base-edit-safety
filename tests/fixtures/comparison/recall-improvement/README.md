# Mock paired recall checks

These mock cases check the difference in recall between two methods at the same
budget. They use the `comparison_evaluation` CLI and output `format_version: "3"`.
The expected results are calculated by hand below.

## `paired-difference/`

One editor/assay stratum, three dependency groups, one observed positive per
group (three eligible positives), budget 1.

- `proposed` ranks each group’s positive above its non-positive, so its top
  pick is always a positive: method recall = 1/3.
- `baseline` ranks `n1` (tested no signal) first: baseline recall = 0/3 = 0.
- Hand-calculated signed difference: `(1 - 0) / 3 = 1/3`
  (`recall_difference_method_minus_baseline`, positive favors the method).

In the paired resamples, baseline picks a positive instead of `n1` exactly when
the draw gives group 1 zero copies: `p2` if group 2 is present, otherwise `p3`.
Resampled differences are only `0` or `1/3`; the fixed seed's 100 draws give
interval `[0, 1/3]`. `splits.json` was generated with
the unchanged `split` command (`--training-fraction 0 --development-fraction 0`,
seed 17) and assigns every group to the evaluation partition.

## `zero-denominator/`

One stratum with no eligible positives (tested no signal, inconclusive and
untested). Both methods have undefined recall; the comparison must report
`status: not_comparable` with a **null** difference, never a zero improvement.
