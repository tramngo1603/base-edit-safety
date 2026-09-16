# How the parts fit together

I use two paths: one searches for guides and makes supported predictions; the other checks
candidate lists and experimental results supplied by the researcher. They share guide and
genome identities but run as separate commands.

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
needs a detection limit. These are checks of the supplied details, not an audit of the experiment.

## Limits built into the checks

- HBG coordinates and alleles are stored in both genomic and promoter orientations.
- BE-DICT uses guide positions 3–10, counted from the PAM-distal end. Unsupported positions or
  PAMs receive no predicted rate. The candidate file has its own declared editor window.
- The ordinary CRISPOR path keeps sites with up to three mismatches and four-mismatch sites
  with CFD at least 0.10. It does not search gaps or population variants.
- Separate imports can cover broader search settings. This version accepts supported DNA-only
  or RNA-only gaps, not mixed DNA/RNA gaps. Unknown mappings and phase remain unresolved.
- Sequence matches, editor-specific effects, guide-independent DNA editing and RNA editing
  stay separate. None of these checks establishes safety.

See the [coordinate definitions](coordinates.md), [local data rules](../data/README.md) and
[validation results](validation.md).
