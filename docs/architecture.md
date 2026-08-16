# Architecture and review boundaries

The workflow is deliberately narrow: fixed GRCh38 windows around the HBG1 and HBG2 fetal-globin
promoters, ABE/CBE candidate designs, and records intended for research review before experimental
validation.

## Workflow

1. `anchors.py` supplies the fixed coordinate/allele definitions and strand conventions.
2. `pipeline/enumerate_guides.py` combines CRISPOR guide/off-target output with BE-DICT allele
   spectra and emits one record per predicted edit.
3. `pipeline/annotate_consequence.py` may add AlphaGenome output; without it, regulatory
   consequence remains explicitly deferred.
4. `annotation.py` adds bounded clonal-hematopoiesis driver-gene overlap flags.
5. `scoring/rubric.py` produces separate, provisional triage terms while retaining the evidence
   and missing-data state used to calculate each term.
6. `scoring/credibility.py` checks ranking behavior and records what remains unmeasured.
7. `pipeline/design_records.py` groups edit records into a file for human review.

## Checks on scientific claims

- HBG1/HBG2 alleles are stored in both genomic-plus and promoter-sense orientation.
- The gene copy is not inferred from spacing alone. An uncertain HBG1/HBG2 assignment remains
  flagged until an assay that distinguishes the two copies resolves it.
- A missing or out-of-window prediction is `not_evaluated`, not a fabricated zero.
- Sequence-match, base-editor-specific, guide-independent, and RNA off-target concerns remain
  separate categories.
- Clonal-hematopoiesis overlap is a locus/gene flag. The workflow does not emit patient risk,
  incidence, or clinical outcome estimates.
- Evidence labels cannot exceed what the accessible source supports. Tests enforce these evidence
  limits and preserve explicit `could_not_verify` states.
- Outputs carry code, tool, reference, and source details and remain ignored by Git until reviewed
  for a separate release.

## Current limits

- Reference-genome CRISPOR search only; population variation and ancestry-aware off-targets are not
  modeled here.
- No clinical interpretation, dose, delivery, engraftment, or long-term outcome model.
- Guide-independent DNA deamination and RNA off-targets are qualitative flags in v1.
- AlphaGenome consequence output is model-predicted and must not be presented as measured evidence.
- Public retrospective validation is still underway; no held-out metric is included in this
  release.

See [coordinates](coordinates.md), [validation safeguards](validation.md), and the
[public data policy](data-policy.md). Definitions of specialized terms are in the
[glossary](glossary.md).
