# Locked validation split

**Status:** locked on 2026-08-14 by Tram Ngo, protocol owner.

This document records which public datasets are reserved for method development and which are
reserved for the final evaluation. Terms are defined in the [glossary](../docs/glossary.md).

The role assignment is not permission to use a dataset. Every registry row still has
`QC status = stub` and remains excluded until its source files, reuse terms, checksums, sample and
guide mappings, reference build, and overlap with development inputs have been reviewed.

## Selection rule

The held-out set is reserved for the registry row that most closely matches the claim being tested:
HBG-promoter base editing in primary hematopoietic stem and progenitor cells (HSPCs). Datasets that
test a different nuclease, genomic location, cell system, or outcome remain in the tuning set.

| Dataset | Locked role | Reason | Current status |
| --- | --- | --- | --- |
| D1 | tuning | Base-editor off-target data can test off-target handling, but its sample, construct, guide, and reference-build mappings remain incomplete. | excluded pending QC |
| D2 | tuning | The registered data use an AsCas12a nuclease rather than a base editor, so they can test sequence-search and HBG1/HBG2 handling but not base-editing outcomes. | excluded pending QC |
| D3 | tuning | The data test a SpCas9 guide at the BCL11A enhancer rather than an HBG-promoter base editor. | excluded pending QC |
| D4 | tuning | The data provide cell-line base-editing outcomes, but not the primary-HSPC HBG context reserved for final evaluation. | excluded pending QC |
| D5 | held-out | Only the file-level subset matching an HBG-promoter assay, a base editor, and primary CD34+ HSPCs matches the final evaluation claim. HUDEP-2 and nuclease-comparator records are excluded. | excluded pending QC |

## Safeguards

- Every measurement and table join must include the edited allele, not only a gene, guide, or
  genomic position.
- D5 may be run once only after its admitted file list, protocol thresholds, and code version are
  locked.
- The D5 held-out subset must remain separate from rules or thresholds derived during development.
  Points that share the Han 2023 source are reported separately through
  [`HAN_OVERLAP_MAP.md`](HAN_OVERLAP_MAP.md).
- Split membership must not add patient risk, incidence, or clinical-outcome estimates.
- Changing any dataset role or the lock date requires a new protocol version. Version 1 is not
  edited in place.

## Items not yet verified

- No D1-D5 row has passed all registry admission checks; every row remains `stub`.
- D1 still needs a complete sample, editor, guide, reference-build, and validation-layer mapping.
- D2 still needs a source-file-to-assay and HBG1/HBG2 mapping.
- D3 still needs source confirmation of the product identity assigned to its guide.
- D4 still needs complete library, editor, sample, reference-build, and development-overlap records.
- D5 still needs file checksums, exact guide and primer mappings, reference-build review, reuse-term
  review, biological-replicate mapping, and development-overlap review.
- The hand-verified reference cases and held-out result files are not included in this public
  release. No validation metric is claimed.

Locked by: Tram Ngo | Date: 2026-08-14
