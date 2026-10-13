# Candidate-panel workflow

Use this part of the tool to bring candidate lists and experiment results together before
choosing sites for follow-up. A candidate nomination proposes a site for testing; a measurement
reports what an experiment found. The tool keeps those records separate from your panel choices.

The JSON output is called the candidate ledger. It records the source of each site, the DNA
versions reported there, matching experiment results and what is missing. A site being in a list
does not show that it edits. Leaving it out of a panel does not show that it is safe.

## Checks before combining files

- Every guide needs a file listing its intended targets. For an HBG design, both HBG1 and HBG2
  must be present, 4,924 bp apart, on the same strand. The output keeps the need for assays that
  distinguish the two copies and check deletions or rearrangements between them.
- A perfect sequence match is intended only when its coordinates match that list. Any other
  zero-mismatch site remains a candidate.
- Standard chromosome aliases such as `1` and `chr1` use the same spelling. Exact coordinates merge;
  probable one-base convention disagreements remain separate and receive warnings. Conflicting
  known strands stop the run.
- Computational site sequences, mismatch counts, PAMs, and editable-window assertions are checked.
  For a real reference, the reported sequence and strand are also checked against an indexed FASTA.
- Reused source IDs and duplicate native rows are rejected. A `not_detected` measurement requires
  a positive numeric detection limit and unit.
- Each converted source has a companion JSON file, also called a sidecar. It records the original
  file's checksum, command, tool version or citation, reference checksum and output checksum.

These checks make conversion errors visible. They do not establish genome-wide sensitivity.

## Try the constructed example

After installing the package, run:

```bash
bes-candidate-panel \
  --profile examples/candidate-panel/search-profile.json \
  --guides examples/candidate-panel/guides.tsv \
  --intended-loci examples/candidate-panel/intended-loci.tsv \
  --nominations examples/candidate-panel/nominations.tsv \
  --source-provenance examples/candidate-panel/source-provenance.json \
  --measurements examples/candidate-panel/measurements.tsv \
  --decisions examples/candidate-panel/panel-decisions.tsv \
  --output /tmp/base-edit-safety-candidate-ledger.json \
  --panel-output /tmp/base-edit-safety-panel-review.tsv
```

Expected output:

```text
wrote 8 candidate sites to /tmp/base-edit-safety-candidate-ledger.json and /tmp/base-edit-safety-panel-review.tsv
```

Every sequence, coordinate, assay value, variant label and driver label in this example was
constructed for a software check. None is a paper's experimental data. No reference FASTA is needed because its profile has `reference_id=SYNTHETIC` and the
source sidecar states that the checksum is not applicable. A real panel must supply an indexed
reference FASTA.

The example has two intended HBG-like controls, sites proposed by more than one method, sites
without follow-up, and candidates with many mismatches or an alignment gap. It also includes
results with detection limits, a result without a matching nomination, and panel decisions.
All five example results lack tested-allele identity, so none establishes allele-specific editing.

## Normalize native candidate files

The package reads the supported formats below. It keeps valid rows regardless of what later
measurements say. The CRISPRme+ reader counts skipped rows and removes any claim of a complete
import if it skips a row.

First index the exact reference used for the source run:

```bash
samtools faidx /data/references/GRCh38.fa
```

### Cas-OFFinder

The Cas-OFFinder adapter accepts the official nine-column output: ID, bulge type, query, DNA,
chromosome, zero-based location, strand, mismatches, and bulge size. Native IDs should equal the
guide IDs; use `--guide-id` only for a single-guide file with numeric native IDs.
It does not invent CFD or MIT values. A profile that requires a score cutoff therefore stops unless
a separately validated scoring step supplies that score in a normalized source.

```bash
bes-normalize-nominations cas-offinder \
  --input /data/runs/hbg/cas-offinder.tsv \
  --profile search-profile.json \
  --guides guides.tsv \
  --intended-loci intended-loci.tsv \
  --source-version "Cas-OFFinder release or commit" \
  --source-command "cas-offinder input.txt G cas-offinder.tsv" \
  --reference-fasta /data/references/GRCh38.fa \
  --output cas-offinder.normalized.tsv \
  --provenance-output cas-offinder.source.json
```

The parser follows the [official Cas-OFFinder output definition](https://github.com/snugel/cas-offinder/blob/develop/README.md#running-cas-offinder).

### Cas-OFFinder 2.4.1

Cas-OFFinder 2.4.1 (upstream commit `9816b94c20c4cba2e79b039e1e2a6dee684b7b66`, tag
`2.4.1`, new BSD license) writes a different layout: tab-separated query, chromosome,
zero-based start, matched DNA with mismatch bases in lowercase, strand, and mismatch
count, optionally followed by a query ID. Declare it explicitly with `--format 2.4.1`
(six columns) or `--format 2.4.1-ids` (seven columns); the default `nine-column` path
is unchanged. Rows are never padded to imitate the nine-column producer.

Every row is checked before normalization: exact field count and no empty fields,
query symbols, query prefix against the bound guide's spacer, exact declared
spacer-plus-PAM geometry for both the query and the matched DNA, strand, integer
and non-negative coordinates and mismatch counts, matched-DNA bases, an
observed PAM inside the search profile, no exact duplicate rows, and the search
profile's own limits. Mismatch positions are recomputed from
the literal target and the IUPAC query under the pinned writer's rules; the
lowercase markers must equal those exact positions, and the declared count must
equal the recomputed spacer count used downstream by `score = -mismatches`. A
mismatch inside the query PAM mask is rejected rather than folded into that
count; a matching `NNN` or compatible PAM mask remains allowed. Target bases
outside A/C/G/T are unresolved ambiguity this bounded parser cannot validate
and are rejected; that is a documented limit, not a count. Rows map to guides
only through the seven-column query ID, `--guide-id` for a single declared
guide, or `--query-map`, a declared JSON object from query string to one guide
ID; duplicate serialized keys are rejected. A bare query never selects a guide
on its own, so two guides that share a spacer require one of these declarations.
The explicit 2.4.1 route also rejects a conflicting producer version:
`--source-version` must be `Cas-OFFinder 2.4.1`, optionally followed by a
parenthesized citation. The sidecar JSON records the declared format, the
binding rule, the query-map file and its checksum when used, complete row
accounting (lines read, blank or comment lines, rows normalized), and each
row's native line number and full raw-row SHA-256 beside the input checksums.

```bash
bes-normalize-nominations cas-offinder \
  --format 2.4.1-ids \
  --input tests/fixtures/cas-offinder-2.4.1/valid-seven-field.tsv \
  --profile examples/candidate-panel/search-profile.json \
  --guides tests/fixtures/cas-offinder-2.4.1/guides.tsv \
  --intended-loci tests/fixtures/cas-offinder-2.4.1/intended-loci.tsv \
  --source-version "Cas-OFFinder 2.4.1 (9816b94c20c4cba2e79b039e1e2a6dee684b7b66)" \
  --source-command "not-applicable: mock parser input" \
  --reference-fasta tests/fixtures/cas-offinder-2.4.1/reference.fa \
  --output cas-offinder.normalized.tsv \
  --provenance-output cas-offinder.source.json
```

This example imports mock Cas-OFFinder rows. The adapter checks their format,
coordinates and sequence against the supplied reference. Empty files are rejected
because the importer requires result rows.

### CRISPRme+

The importer reads CRISPRme+ 2.3.2 `*_integrated_results.tsv` files with a required JSON record
of the input checksum and search settings. A separate format check can read the permitted sheet
of the older deposited sg1617 workbook. Neither path retains sample values or private-target
sheet content. The record names the source version and its `CFD` scoring model.

The recorded HBG pilot checked a real 2.3.2 run: search and merge finished, plotting failed, and
a separate unchanged export step produced the imported tables. Different inputs and settings
need their own runtime check. See the [CRISPRme+ guide](crisprme-import.md) for commands,
required fields, coordinate rules and the optional local-workbook test.

### CHANGE-seq-BE identified table

The `change-seq-be` format reads the identified-matched header emitted by the official BE branch.
It retains the experimental nomination separately from measurements and carries source read,
control, overlap, replicate-file, and context fields in `source_metadata`.

```bash
bes-normalize-nominations experimental \
  --format change-seq-be \
  --input /data/runs/hbg/sample_identified_matched.txt \
  --guide-id HBG_guide_1 \
  --guides guides.tsv \
  --intended-loci intended-loci.tsv \
  --reference-id GRCh38 \
  --reference-fasta /data/references/GRCh38.fa \
  --source-method CHANGE-seq-BE \
  --source-version "CHANGE-seq BE branch commit or publication DOI" \
  --source-command "change_seq.py all --manifest manifest.yaml" \
  --output change-be.normalized.tsv \
  --provenance-output change-be.source.json
```

The expected header is documented by the [official CHANGE-seq BE branch](https://github.com/tsailabSJ/changeseq/tree/BE). This adapter imports an output that the researcher is authorized to use; it does not download or redistribute publisher tables.

### CIRCLE-seq or another experimental BED6 list

Use `--format bed6` for a zero-based, half-open BED6 nomination list. The BED score is retained as
`BED_score`; it is not assumed to be an edit rate. Converting another source format to BED6 remains
a declared upstream step and must not change coordinates or filter rows.

```bash
bes-normalize-nominations experimental \
  --format bed6 \
  --input /data/runs/hbg/circle-seq-sites.bed \
  --guide-id HBG_guide_1 \
  --guides guides.tsv \
  --intended-loci intended-loci.tsv \
  --reference-id GRCh38 \
  --reference-fasta /data/references/GRCh38.fa \
  --source-method CIRCLE-seq \
  --source-version "pipeline release, commit, or publication DOI" \
  --source-command "exact source command or an explicit acquisition record" \
  --output circle-seq.normalized.tsv \
  --provenance-output circle-seq.source.json
```

Repeat `--nominations` and `--source-provenance` when building the union. Each normalized basename
must have exactly one matching sidecar.

## Combine real candidate lists

```bash
bes-candidate-panel \
  --profile search-profile.json \
  --guides guides.tsv \
  --intended-loci intended-loci.tsv \
  --reference-fasta /data/references/GRCh38.fa \
  --nominations cas-offinder.normalized.tsv \
  --source-provenance cas-offinder.source.json \
  --nominations change-be.normalized.tsv \
  --source-provenance change-be.source.json \
  --measurements measurements.tsv \
  --decisions panel-decisions.tsv \
  --output candidate-ledger.json \
  --panel-output panel-review.tsv
```

## Input files

All TSV coordinates are zero-based and half-open: `start_0based` is included and
`end_0based_exclusive` is excluded.

### Guide manifest

`--guides` defines the spacer, editor, PAM class, one-based inclusive editing window, and structured
`locus_system`. Use `HBG1_HBG2` for the duplicated fetal-globin promoter or `other` for another
declared system. Guide IDs are stable keys shared by every input.

### Intended-locus manifest

`--intended-loci` requires `guide_id`, `genome_build`, `chrom`, `start_0based`,
`end_0based_exclusive`, `strand`, and `locus_name`. These coordinates—not mismatch count and not a
source label—control the intended/candidate classification.

If either `HBG1` or `HBG2` appears for a guide, both are required. The manifest also causes an
intended control to remain in the ledger if a source search failed to nominate it.

### Search profile

The JSON profile declares the reference, mismatch and bulge limits, IUPAC PAM patterns, variant
scope, intended use, and retention rule. A score cutoff is represented structurally by
`retention_score_name` and `minimum_score`; both are null when no cutoff is used. Computational rows
outside those limits are rejected.

The companion JSON records the source search settings and how many rows were
imported or skipped.

### Nomination table

Required fields identify the guide, nomination kind, source, build, interval, and stable
`source_record_id`. Computational rows also require PAM, mismatch count, search profile, strand,
and site sequence. Bulge fields default to zero. Ungapped rows have their mismatch and editable-
window assertions recomputed. A generic gapped row leaves `editable_base_in_window` blank unless
an alignment-aware adapter supplies the guide/target mapping. The CRISPRme adapter supplies that
mapping, and ledger v3 derives allele-specific editability without trusting a source boolean.

`site_role` is an optional assertion for readability. The workflow derives the role again from the
intended-locus manifest and stops if the assertion disagrees. `source_metadata_json` preserves
method-specific fields in the JSON ledger.

Population-variant and clonal-hematopoiesis labels are currently imported annotations. They are
placed in explicit review groups but are not independently recomputed by this workflow.

### Measurement table

Measurements require an assay, citation/version, biological context, stable source ID, and one of
`edited`, `not_detected`, or `inconclusive`. Values use separate numeric and unit columns. Detection
limits use `limit_of_detection_value` and `limit_of_detection_unit`; a `not_detected` row requires
both.

Usable editing non-detections require a limit in `percent`, `percentage`, `%`, `fraction` or
`proportion`, ignoring unit case and surrounding spaces. Percentage values must be greater
than zero and at most 100; fractions must be greater than zero and at most 1. Source values are
preserved. Unknown units or out-of-range limits remain in the ledger but earn no usable editing
coverage. Population studies can also declare a maximum acceptable limit in their assay policy;
see the [population guide](population.md#match-evidence-and-declare-assays).

Optional fields retain `replicate_id`, treated and control group labels, read depth, QC status,
statistical test, raw p-value, and adjusted p-value. `not_detected` means no detected editing within the stated detection limit in
that assay and context, never universal absence of editing.

The optional columns below identify the DNA that was tested. Blank values stay unknown; the
importer does not fill them from the guide file or from reference DNA.

| Field | Meaning |
| --- | --- |
| `editor` | The editor actually used; must match the guide manifest for editing coverage. |
| `measurement_kind` | `editing`, `consequence`, or `not_reported`; only editing results contribute to editing coverage. |
| `allele` | `ref`, `alt`, or `unknown` for the DNA version tested. |
| `allele_sequence` | Full target sequence, including PAM, in genomic-forward orientation. Use A/C/G/T/N without alignment gaps; unknown bases cannot establish an exact match. |
| `allele_resolution` | `resolved`, `mixed`, or `not_reported`. Resolved means the result and its detection limit apply to that DNA version, through a known template or an allele-resolving assay. |
| `allele_evidence` | Source evidence supporting the tested sequence and allele resolution. A candidate's predicted sequence alone is insufficient. |
| `configuration_id` | Optional exact source configuration to restrict the match; this identifier is not genotype evidence. |
| `haplotype_id`, `phase_status` | Non-identifying template label and `phased`, `unphased`, or `not_reported`. They do not supply missing phase for an upstream candidate. |
| `qc_note` | Source explanation for a quality result, when available. |

A usable editing result requires matching guide, editor, build, interval, strand and known DNA
sequence, explicit allele evidence and resolution, an assay and cell context, treated and control
groups, and passed QC. A negative also needs its detection limit. This checks declared source
metadata; it does not independently audit the experiment or infer genotypes from population data.

The original `status`, values, QC and other source fields stay in `measurements`.
`reported_measurement_statuses` lists the source outcomes. Each entry in `measurement_coverage`
identifies one reported DNA configuration, lists how each measurement applies or fails to apply,
and gives a plain-language follow-up question. A reference-only negative never covers ALT.
Missing identity information, failed/unreported QC, inconclusive results and unresolved candidate
mappings cannot produce usable negatives. Multi-change ALT sequences with unresolved upstream
phase remain unresolved even when a measurement supplies a haplotype label.

`measurement_summary` is now derived from those per-configuration checks. `partial_coverage`
means usable evidence exists for only some reported DNA versions; `coverage_unresolved` means
coverage cannot be established. `unusable` means matching records cannot support an editing
conclusion. `not_measured` means no matching measurement is established. An `edited` or
`not_detected` summary remains limited to the recorded assays and cells. Different usable outcomes
stay `mixed`; rates, detection limits and populations are never averaged. One valid record can
support a bounded result while another unusable record remains visible with its reason.

Observations labelled as consequences stay separate; this tool does not assess their
coverage or infer a consequence from editing. Design-level predictions and evidence tiers do not
supply candidate measurement coverage. The TSV and Markdown report show the original outcomes,
coverage gaps and follow-up questions alongside human panel decisions and HBG locus checks.
No allele is removed because it is rare or because a panel decision excludes it.

Keep individual genotypes and donor identifiers outside Git. Regenerate older ledgers to obtain
per-allele measurement fields; the report labels their old summaries as unassessed until then. See the
[validation summary](validation.md#bounded-hbg-pilot) for the bounded replay and checks.

### Panel decisions

Each decision is `include`, `exclude`, or `undecided`. The reason is one of `not_reviewed`,
`low_mm_high_cfd`, `bulge_site`, `variant_created_pam`, `variant_reduced_mm`,
`ch_driver_overlap`, `genic_context`, `lab_nominated`, `no_base_in_usual_window`,
`redundant_with`, or `other_note`. An optional note is limited to 140 characters.
`redundant_with` requires a valid `related_site_id`; that field is forbidden for other reasons.
Missing decisions default to `undecided` / `not_reviewed`. An unknown coordinate stops the run.

## Outputs and interpretation

The JSON ledger is the full record behind the report. It contains the retained source rows,
sequence checks, intended targets, source details, checksums, decisions and limits. It does
not add absolute paths for its own input files. The exact external source command is preserved
verbatim and may itself contain paths. Never put credentials, patient information, or donor/sample
identifiers in that command, and review the source details before sharing a resulting ledger.

The TSV is a shorter table for review. It keeps excluded and undecided sites, computational source
buckets, reconciliation state, allele-specific editability, per-source scores, detection limits,
identity warnings, locus requirements, required assays, replicate IDs, read depths, and
measurement QC.
Review groups divide large panels into intended controls, measured/reconciliation sites,
experimental nominations, annotation/bulge/variant flags, and broad mismatch bands. These groups
organize inspection. Measurement status and panel decisions remain separate fields.

Use `--panel-notes` for notes about the whole panel, such as sequencing depth. Do not include
sample identifiers. Use `--previous-ledger` to compare sites with an earlier JSON record. The
comparison lists added, removed and regrouped site IDs; it does not add removed sites to the
current candidates.

Make a Markdown report from a design-record v7 file and a candidate-ledger v3 file:

```bash
bes-candidate-report \
  --design-records design-records.json \
  --candidate-ledger candidate-ledger.json \
  --output candidate-review.md
```

The report records both input checksums and explains the limits of each section. The two input
formats remain separate. Keep reports containing paper-derived rows or individual measurements
private; public summaries contain aggregate findings and citations under the [data policy](data-policy.md).
