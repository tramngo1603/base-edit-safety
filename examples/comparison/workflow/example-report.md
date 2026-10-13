# Comparative research report

## Research decision

Which design should be compared next in controlled, context-matched donor experiments?

Synthetic design A, Synthetic design B: **insufficient comparable evidence**. Run a matched comparison because current evidence cannot choose between these designs.

Recorded reasons: no_comparable_paired_source_measurements; no_universal_safety_or_efficacy_ranking.

## Designs and scope

Execution mode: **synthetic_test**; coverage basis: **synthetic_truth**.

Two artificial 23-base loci for a complete offline experiment

Declared reference assembly: **synthetic_GRCh38_geometry**.

The shared requested domain contains 2 declared intervals and 46 interval bases. Per-design evaluated regions, skips and retention are retained in search-audit.json and design-comparison.tsv.

This mock example uses the specimens, capabilities and prices declared in its inputs.

| Design | Editor | PAM | Mismatch / DNA gap / RNA gap caps | Scope | Donors evaluable / total | Unresolved donors | Complete loci | Complete sequences | Complete alignments |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Synthetic design A | synthetic-ABE-v1 | NGG | 1 / 0 / 0 | complete | 2 / 2 | 0 | 2 | 3 | 3 |
| Synthetic design B | synthetic-editor-B | NGG | 1 / 0 / 0 | complete | 2 / 2 | 0 | 2 | 3 | 3 |

Counts describe the declared sequence search. A missing complete count is unknown, including when some results were retained before a partial search.

## Observations, sequence findings and predictions

| Source | Design | Context | Endpoint | Raw source value | Normalized value | Controls / QC | Adjudicated outcome | Adjudication reason | Source outcome (verbatim) | Detection limit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Artificial observation example | Synthetic design A | Unrelated artificial historical cell line; synthetic electroporation; 48 hour; synthetic treated condition | intended_edit_fraction | 9.8000000000000007 percent | 0.098000000000000007 fraction | unknown / unknown | inconclusive | Legacy extraction map has no explicit outcome adjudication; conclusive credit withheld | measured_fraction | 0.01 fraction |

Adjudications are explicit four-state outcomes: observed positive, tested no signal, inconclusive and untested. Inconclusive and untested records withhold conclusive empirical credit even with a numeric value and passed QC. No-signal results include their stated detection limit. Censored or missing values retain their typed outcome and reason but withhold conclusive credit.

Each observation row is a source measurement. Technical repeats do not establish independent donors; missing source-unit metadata remains unknown.

| Source | Design | Endpoint | Measurements | Reported independent units | Measurements with unknown unit |
| --- | --- | --- | --- | --- | --- |
| Artificial observation example | Synthetic design A | intended_edit_fraction | 1 | 1 | 0 |

Applicability assessments: related design context: 2. Each assessment checks the donor identity, context, endpoint and assay quality separately.

Current empirical target-donor credits, counting only active observation revisions that pass the shared outcome and quality gate: 0. Current study permission must also enable empirical credit in research specimen planning mode. Stored credit bits in all retained applicability records: 0; a historical or non-conforming stored bit cannot grant current credit.

Intended loci were assessed independently of nomination. The protected intended-assessments.tsv retains every assessment, including lost PAMs, indels and unresolved sequences.

Recorded donor/reference sequence effects: changed: 2; unchanged: 8. These are sequence findings.

No prediction records were supplied.

Published evidence outside the admitted profile or with unresolved mapping remains in out-of-scope-evidence.tsv: 0 retained entries.

## Proposed experiments

Plan status: **heuristic_complete**. Completed experimental structures: **1**; frozen priority weight covered: **1**. Existing evidence weight: 0; newly planned weight: 1.

Planning and measurement status are shown separately below.

**Q1 — design by variant**. Evaluate design_by_variant at synthetic_chr11:200-223:+ under declared rule complete-design-by-variant

Measure the declared endpoint; genotype contrasts establish association, not variant causality Source-bound suggestion from variant_candidate; Each 2-unit treated panel measures T and L with its assay controls/QC; each 1-unit donor mock measures T and L with its assay controls/QC

Requirements: carrier: 1 independent donor(s); noncarrier: 1 independent donor(s); 1 technical measurement(s) per endpoint; within donor. Selected donors: Donor 1 in noncarrier; Donor 2 in carrier.

Planning: **fully_planned**; evidence: **no_measurements**; priority: 1; hard required: yes.

| Operation | Donor alias / arm | Design | Type | Readouts | Actual target DNA | QC / controls | Prerequisites |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Operation 1 | Donor 2; carrier | Synthetic design B | treatment | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG, ACGTGCGTACGTACGTACGTAGG | preparation only | none |
| Operation 2 | Donor 1; noncarrier | Synthetic design A | treatment | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG | preparation only | none |
| Operation 3 | Donor 2; carrier | Synthetic design A | measurement panel | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG, ACGTGCGTACGTACGTACGTAGG | assay_qc, intended_activity, assay_qc | Operation 4, Operation 10 |
| Operation 4 | Donor 2; carrier | design-neutral control | mock control | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG, ACGTGCGTACGTACGTACGTAGG | assay_qc, mock, assay_qc | none |
| Operation 5 | Donor 1; noncarrier | Synthetic design B | measurement panel | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG | assay_qc, intended_activity, assay_qc | Operation 8, Operation 9 |
| Operation 6 | Donor 1; noncarrier | Synthetic design A | measurement panel | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG | assay_qc, intended_activity, assay_qc | Operation 2, Operation 9 |
| Operation 7 | Donor 2; carrier | Synthetic design B | measurement panel | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG, ACGTGCGTACGTACGTACGTAGG | assay_qc, intended_activity, assay_qc | Operation 1, Operation 4 |
| Operation 8 | Donor 1; noncarrier | Synthetic design B | treatment | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG | preparation only | none |
| Operation 9 | Donor 1; noncarrier | design-neutral control | mock control | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG | assay_qc, mock, assay_qc | none |
| Operation 10 | Donor 2; carrier | Synthetic design A | treatment | unintended_edit_fraction at synthetic_chr11:200-223:+; intended_edit_fraction at synthetic-HBG1 | synthetic-HBG1: ACGTACGTACGTACGTACGTAGG; synthetic_chr11:200-223:+: ACGTACGTACGTACGTACGTAGG, ACGTGCGTACGTACGTACGTAGG | preparation only | none |

## Cost and unresolved requirements

| Resource package | Distinct instance subtotal |
| --- | --- |
| assay_library | 8 |
| donor_activation | 2 |
| mock_control | 2 |
| sequencing_batch | 0 |
| specimen_preparation | 0 |
| study_setup | 2 |
| treatment | 4 |

Total: **18 simulation_units** against budget **18** (declared cost scale 1). Each active named resource is charged once; resources.tsv lists its actual consumers and capacity demands.

Selector: complete_bundle_greedy; global optimum established: no; requested mandatory questions met: yes. Search limits and any infeasibility certificate remain in plan.json.

The selected plan has no recorded requirement gaps. Returned measurements and their controls still require reconciliation.

## Report format and human review records

Report format version: **2**; observation outcome contract: `comparison_observation_revision:v2`; accepted observation history versions: 1, 2. The machine-readable copy is [report-format.json](report-format.json).

Producing code identity captured when the study was resolved: package **0.5.0**, source SHA-256 `b131daeaa1be1f360f5cb1b0bc37ffa1e6e94ca3f2b94db1356b5c654b967f2c`, Python 3.12.12. It is not recomputed at render time.

Reviewer attestation: **not recorded**. No reviewer attestation is present in the report inputs.

Funding disclosure: **not recorded**. No funding disclosure is present in the report inputs.

Customer disclosure: **not recorded**. No customer disclosure is present in the report inputs.

Detailed donor, sequence and source records remain protected local research files. Donor aliases and aggregation do not guarantee anonymity.

Exports: [design comparison](design-comparison.tsv), [questions](questions.tsv), [experimental operations](experimental-plan.tsv), [resources](resources.tsv), [unmet requirements](unmet-requirements.tsv), [observations](observations.jsonl), [evidence applicability](evidence-applicability.tsv), [out-of-scope evidence](out-of-scope-evidence.tsv), [plan](plan.json), [report format](report-format.json). The stage manifest binds source hashes, policies, code and verified ancestors.
