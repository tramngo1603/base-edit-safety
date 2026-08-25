# Workflow overview

The workflow starts with a guide, an editor, and fixed GRCh38 windows around the HBG1 and HBG2
promoters. It produces one review record for each guide/editor design.

## Main steps

1. `anchors.py` supplies the coordinate and allele definitions.
2. `pipeline/enumerate_guides.py` combines CRISPOR candidate sites with BE-DICT edit predictions.
3. `annotation.py` checks overlap with the documented clonal-hematopoiesis gene set.
4. `pipeline/design_records.py` groups the results into a compact record for each design.

## Checks built into the workflow

- HBG1 and HBG2 coordinates are stored explicitly in both genomic and promoter orientation.
- A guide that may engage both promoters is flagged for a locus-resolving assay.
- Missing predictions stay marked as `not_evaluated`; they are not converted to zero.
- Sequence-match candidates, editor-related concerns, and RNA effects remain separate fields.
- Each result records its reference build, source, tool version, and code version.

## Current limits

The current CRISPOR adapter searches the reference genome. Population variants, guide-independent
DNA deamination, RNA editing, and regulatory consequences need separate analyses.

See the [coordinate reference](coordinates.md), [local data layout](../data/README.md), and
[validation summary](validation.md).
