# Validation

This page summarizes the validation work without publishing source tables or paper-derived records.

## What has been checked

- Coordinate and strand tests cover the fixed GRCh38 HBG1 and HBG2 promoter windows.
- Tests confirm that HBG1 and HBG2 remain separate loci and that a guide matching both receives the
  expected dual-locus and deletion warnings.
- Hand-written sample files exercise the analysis and reporting workflow from input parsing through
  final output.
- Tests check that missing predictions remain marked as unavailable instead of being treated as
  zero.
- The validation runner checks dataset roles, file checksums, code versions, and the one-use lock
  for reserved data.

The research workspace uses a written validation plan. Development data may be used to find and fix
problems; a separate reserved set is kept for the final check. Dataset roles and pass marks are
fixed before that final run.

## What is kept private

The public repository does not include publisher tables, source workbooks, paper-derived site
lists, accession-level sample maps, raw or normalized measurements, per-paper review notes, or
reserved-data results. Those records remain in the private research workspace.

## Reproducing the software checks

After installing the project with its test dependencies, run:

```sh
python -m pytest -q
```

The public test suite checks the exposed analysis workflow. The detailed dataset runner, source
records, and its dataset-specific checks remain in the private research workspace. Validation
against research datasets is still underway, so no final performance result is reported here.
