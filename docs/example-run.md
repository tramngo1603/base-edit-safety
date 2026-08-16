# Example workflow

This page describes the public workflow shape without publishing a current ranking, a measured
validation result, or a held-out outcome.

## 1. Configure external tools

Copy `.env.example` to the ignored `.env` file and set the local BE-DICT and CRISPOR paths. The
optional AlphaGenome credential should preferably be supplied through the process environment.

## 2. Enumerate candidate edits

```bash
python -m base_edit_safety.pipeline.enumerate_guides
```

The enumerator fetches the fixed GRCh38 HBG1/HBG2 windows, asks CRISPOR for guides and sequence-match
off-targets, asks BE-DICT for predicted allele spectra, and writes
`data/outputs/enumeration_HBG1_HBG2.json`.

Each emitted record carries the candidate identity, locus and strand, intended or bystander allele,
prediction status, separated off-target categories, HBG gene-copy flags, and source/tool details.
Malformed or unsupported predictions are represented as explicit exclusions or `not_evaluated`
states instead of being converted to zero without a record of the missing prediction.

## 3. Build review records

```bash
python -m base_edit_safety.pipeline.design_records
```

This writes `data/outputs/design_records_HBG1_HBG2.json`. A design record groups the relevant edit
records and presents evidence tier, registry relationship, locus-specific hazards, clonal-
hematopoiesis overlap flags, uncertainty, and the assays needed to resolve it.

Generated files remain ignored because they may contain machine paths or unpublished predictions.
Use the synthetic fixtures under `tests/fixtures/` for code review and automated tests.
