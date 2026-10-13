# How the parts fit together

BES has a local app for reference and sample searches, plus command-line workflows for
candidate review, design comparison and assay planning. They share guide and genome identities.

## Local app

`launcher.py` provides setup, input review, progress and saved-result browsing. `lab.py`
prepares guides, reference files and optional sample variants, then runs the search.
`lab_workflow.py` supplies the reference-only example and saves user-entered context.

`lab_report.py` writes the HTML report. `lab_results.py` reads completed analyses and
exports candidate tables, reports or analysis copies. Candidate exports can include all
rows, every row matching the applied filters, or selected rows. All and matching exports
use the complete results table, including rows beyond the app's display limit.

## Reference and sample search

The app and `bes-search` use `population/research.py`. It checks pinned inputs,
searches declared primary chromosomes, compares reference and sample candidates, and writes
coverage, exclusions and optional evidence links. `research_genotypes.py` reconstructs local
sample sequences from the VCF and declared ploidy. `research_native.py` runs Cas-OFFinder,
records execution limits and validates every returned row. See the [search guide](research-quickstart.md).

## Donor sequences and proposed assays

`population/` handles strict study contracts, donor-linked variants, local sequence paths,
bounded search, evidence applicability and assay actions. `pipeline/population.py` exposes
the file-based stages. `pipeline/donor_panel.py` resolves their source details before selecting
actions under the declared cost limits. The planner counts distinct biological donors and
keeps proposed coverage separate from measurements.

The core uses the Python standard library. Optional indexed public-data acquisition and
normalization use pysam and its bundled bcftools. The local scanner does not require an external
search tool. Read the [population guide](population.md) for contracts and interpretation limits.

## Guide searches and predictions

1. `anchors.py` defines the HBG coordinates and intended base changes.
2. `pipeline/enumerate_guides.py` calls CRISPOR for reference-genome candidates and BE-DICT
   for supported edit predictions.
3. `annotation.py` checks candidate overlap with the recorded blood-cell clonal expansion genes.
4. `pipeline/design_records.py` writes a record for each guide/editor design.

HBG1 and HBG2 are separate locations and both are intended targets here. The supplied target
list identifies them by coordinate. An identical sequence elsewhere stays a candidate.
The records flag the need to distinguish the two copies and to check for deletions or
rearrangements across their approximately 4.9-kb separation.

## Candidate lists and experiments

1. `pipeline/normalize_nominations.py` reads supported search or assay files and writes a common
   TSV format. A companion JSON file records the source command, versions, checksums and row counts.
2. `pipeline/candidate_panel.py` checks the locations and sequences, combines candidate lists,
   and checks each result against the DNA version it is declared to measure. Researcher choices
   about which sites to test are recorded separately.
3. `pipeline/candidate_report.py` reads the design records and candidate list to make a Markdown
   report of results, missing information and follow-up questions.

The readers support Cas-OFFinder, CRISPRme+, CHANGE-seq-BE identified tables and experimental
BED6 lists. They do not install or run the source tools. See the [candidate guide](candidate-panel.md).

## Three different kinds of information

- **Candidate nomination:** a search or experiment proposes a site for follow-up.
- **Prediction:** a sequence score or model estimates an outcome within its supported settings.
- **Measurement:** an experiment reports what happened under stated conditions.

A score cannot turn a candidate into a measured edit. A missing prediction cannot become zero.
A panel decision records what the researcher plans to test; it does not change the evidence.

A usable editing result must match the guide, editor, location, strand and tested DNA. It also
needs supplied allele evidence, cell context, controls and passed quality checks. A negative
needs a detection limit.

## Limits built into the checks

- HBG coordinates and alleles are stored in both genomic and promoter orientations.
- BE-DICT uses guide positions 3–10, counted from the PAM-distal end. Unsupported positions or
  PAMs receive no predicted rate. The candidate file has its own declared editor window.
- The ordinary CRISPOR path keeps sites with up to three mismatches and four-mismatch sites
  with CFD at least 0.10. It does not search gaps or population variants.
- Separate imports can cover broader search settings. This version accepts supported DNA-only
  or RNA-only gaps, not mixed DNA/RNA gaps. Unknown mappings and phase remain unresolved.
- Sequence matches, editor-specific effects, guide-independent DNA editing and RNA editing
  stay separate.

See the [coordinate definitions](coordinates.md), [local data rules](../data/README.md) and
[validation results](validation.md).
