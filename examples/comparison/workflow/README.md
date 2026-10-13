# Complete offline comparison

This mock example compares two editing designs in two diploid samples. It searches
two loci, imports a measurement, creates a carrier/noncarrier question and plans
the treatments, readouts and controls. It runs offline with the installed Python
package.

From the repository root after installing the package:

```sh
BES_COMPARE_CONFIG=examples/comparison/workflow/study.json
BES_COMPARE_RUN=/tmp/bes-comparison-offline
bes-compare validate --config "$BES_COMPARE_CONFIG"
bes-compare search --config "$BES_COMPARE_CONFIG" --output "$BES_COMPARE_RUN/search"
bes-compare observations --config "$BES_COMPARE_CONFIG" --output "$BES_COMPARE_RUN/observations"
bes-compare evidence --config "$BES_COMPARE_CONFIG" \
  --search "$BES_COMPARE_RUN/search" --observations "$BES_COMPARE_RUN/observations" \
  --output "$BES_COMPARE_RUN/evidence"
bes-compare compare --config "$BES_COMPARE_CONFIG" \
  --search "$BES_COMPARE_RUN/search" --evidence "$BES_COMPARE_RUN/evidence" \
  --output "$BES_COMPARE_RUN/comparison"
bes-compare questions --config "$BES_COMPARE_CONFIG" \
  --comparison "$BES_COMPARE_RUN/comparison" --output "$BES_COMPARE_RUN/questions"
bes-compare plan --config "$BES_COMPARE_CONFIG" \
  --questions "$BES_COMPARE_RUN/questions" --evidence "$BES_COMPARE_RUN/evidence" \
  --output "$BES_COMPARE_RUN/plan"
bes-compare report --config "$BES_COMPARE_CONFIG" \
  --comparison "$BES_COMPARE_RUN/comparison" --questions "$BES_COMPARE_RUN/questions" \
  --plan "$BES_COMPARE_RUN/plan" --output "$BES_COMPARE_RUN/report"
```

Read `comparison-report.md` in the report directory. Its accompanying TSV and
JSON files contain the source records, actual target sequences, requirements,
operations and unique resource charges. Output directories use mode `0700` and
files use `0600`. Detailed exports may contain identifiers; the default readable
report omits direct donor identifiers and local paths.

[The example report](example-report.md) shows the expected narrative.
Its export links refer to files generated together in your report directory.
`expected.json` records the expected counts and cost. Each generated `run-manifest.json`
records that run's code, input and output hashes; it stays in the output directory.
Historical example reports and receipts are retained outside the public repository.

The plan completes one question using ten operations: four treatments, four
treated panels and two donor mocks. Each panel measures the intended locus T
and candidate locus L, with its declared assay controls and quality checks.
Twenty requirement slots are filled. Total cost is 18 simulation units: setup
2, donor activation 2, treatments 4, treated panels 8 and mocks 2.

At L, one sample has an A>G substitution on one copy; the other has two reference
copies. Both have resolved intended sequences at T. This sequence difference
provides the reason for the planned comparison. The imported
historical measurement has a different cell context, unresolved donor DNA and
unknown controls. Its exact value, `9.8000000000000007` percent, remains visible;
it supplies no usable donor evidence for the planned contrast. Its extraction map
is a legacy map with no explicit outcome adjudication, so the import records the
measurement as inconclusive with the documented legacy reason; the report shows
that adjudication, the verbatim source outcome and the stated detection limit.

## Change the budget

The second configuration changes only the budget. Its search, observations,
evidence, comparison and question stages can reuse the first run's output files.

```sh
BES_COMPARE_SHORT=examples/comparison/workflow/study-budget17.json
bes-compare search --config "$BES_COMPARE_SHORT" --output "$BES_COMPARE_RUN/search" --reuse
bes-compare plan --config "$BES_COMPARE_SHORT" \
  --questions "$BES_COMPARE_RUN/questions" --evidence "$BES_COMPARE_RUN/evidence" \
  --output "$BES_COMPARE_RUN/plan17"
bes-compare report --config "$BES_COMPARE_SHORT" \
  --comparison "$BES_COMPARE_RUN/comparison" --questions "$BES_COMPARE_RUN/questions" \
  --plan "$BES_COMPARE_RUN/plan17" --output "$BES_COMPARE_RUN/report17"
```

The 17-unit plan cannot complete the required contrast. It reports
`infeasible_proven` for this fully enumerated, tiny declared catalog and
`meets_request: false`. Both commands still exit successfully because they
completed their calculation. Scientific infeasibility is a result, not a
software crash. Larger bounded searches can instead report
`no_feasible_bundle_found`; that is not a proof of impossibility.

Existing stage directories are immutable. Add `--reuse` to verify an unchanged
completed stage, or choose a new output directory. There is no comparative
`--overwrite` flag. `--output-dir` is an alias for `--output`. Every command also
has the equivalent prefix `python -m base_edit_safety.pipeline.comparison`.

## Mode declarations

```sh
bes-compare validate --config examples/comparison/workflow/mode-declarations/public-positive/study.json
bes-compare validate --config examples/comparison/workflow/mode-declarations/research-positive-artificial-declarations/study.json
python examples/comparison/workflow/validate_modes.py
```

These mock examples check input declarations and eligibility using
`synthetic_chr11`. They cover hypothetical specimens, confirmed specimen
declarations, permitted use, assay capability, sequence support and independent
groups. The script also checks rejection of unresolved specimen availability and
empirical credit in retained-callset simulation.

For real studies, follow the [comparison guide](../../../docs/comparison-workflow.md),
[source contracts](../../../docs/comparison-contracts.md) and
[selector benchmark guide](../../../docs/comparison-benchmark-command.md).
