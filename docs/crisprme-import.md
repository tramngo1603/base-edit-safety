# CRISPRme+ candidate import

Base Edit Safety imports a CRISPRme+ candidate search; it does not install or execute CRISPRme+.
The importer expects the native `*_integrated_results.tsv` format from CRISPRme+
2.3.2. The deposited sg1617 XLSX layout from the original CRISPRme paper is supported only for an
older-format check. Both formats use the same row checks. No extra Python dependency is needed.

This version defines the accepted file format. The [bounded HBG pilot](validation.md#bounded-hbg-pilot)
checked an actual 2.3.2 run with recorded reference and variant inputs. Search/merge completed,
but plotting failed; a separate unchanged export step produced the imported tables. This does
not establish complete-command success or compatibility with other inputs and search settings.
Later changes address incorrect variant frequencies and missing candidates; review the
[upstream changelog](https://github.com/pinellolab/crisprme-plus/blob/main/CHANGELOG.md)
before choosing the runtime for another analysis.

The official
[`v2.3.2` release tag](https://github.com/pinellolab/crisprme-plus/releases/tag/v2.3.2) is at source
revision `f33bcba150c116bf294125d4b5c5a7d413661474`. The recorded Docker image is
`pinellolab/crisprme:v2.3.2`; its multi-platform image-index digest is
`sha256:17e1f20464d0b942a363f7aeeb5359b5cde6f227fa45e52369ceda2776641ca8`.
A separate check verified that source revision and observed runtime version 2.3.2 using the
Linux ARM64 image `pinellolab/crisprme@sha256:33ddaa089d449d458402c89a08b6e44a4ab413f0899af9b3e8dd725833354038`.
Record the resolved platform digest and exact invocation with the run's files. CRISPRme+ calls
its sequence-similarity score `CFD`; neither the release nor its CLI defines a `CFD++` score.

CRISPRme+ is distributed separately. Review its upstream installation instructions and license
before running it. The upstream project describes AGPL-3.0 academic use and requires a separate
license for commercial or for-profit use.

## Produce the native result

The exact reference, variant files, guide file, and PAM file are inputs to the external run.
Treat this command as a template. Its broad PAM setting and data paths must be reviewed for the
chosen guide/editor design before use.

```bash
docker run --rm \
  -v /data/crisprme:/DATA \
  -w /DATA \
  pinellolab/crisprme@sha256:17e1f20464d0b942a363f7aeeb5359b5cde6f227fa45e52369ceda2776641ca8 \
  crisprme.py complete-search \
    --genome /DATA/Genomes/hg38 \
    --vcf /DATA/vcf.config.txt \
    --samplesID /DATA/samplesIDs.config.txt \
    --guide /DATA/guides-crisprme/hbg-guide.txt \
    --pam /DATA/PAMs/20bp-NNN-NO-PAM.txt \
    --mm 6 \
    --bDNA 2 \
    --bRNA 2 \
    --max-total-edits 10 \
    --be-window 4,10 \
    --be-base A \
    --merge 0 \
    --output hbg-v2 \
    --thread 32
```

Record the installed version with `crisprme.py --version`; the importer currently requires
`v2.3.2` output. CRISPRme+ expects each guide entry to include the PAM-width placeholder, so a 20-base
spacer searched with the NNN PAM file must be copied to a run input as `SPACERNNN`. This does not
change the spacer. The PAM argument is a file path, not the literal string `NNN`.

CRISPRme+ writes results below `Results/<output>/`. Select the native
`*_integrated_results.tsv`; do not import `*_all_results_with_alternative_alignments.tsv` or a
spreadsheet exported from another application. `--max-total-edits` is explicit because the tagged
source and its help text have differed about the default. A cap of 10 includes the declared limit of six mismatches plus two DNA and two RNA gap bases.
This broad template can produce mixed-gap rows that this main-v2 importer does not support.
Those rows are counted as skipped; this template is not an example of a fully supported import. The base-editor flags annotate ABE
susceptibility in the 1-based positions 4--10 window; they do not nominate or retain candidates for
Base Edit Safety. Parameters in the declaration must describe the actual run.

## Earlier format checks

The format checks below used the stated test inputs. They do not establish search coverage
for the real HBG pilot inputs.

The v2.3.2 same-input CD7 smoke run produced two tab-separated files and no workbook:

- `*_integrated_results.tsv`: the primary import table;
- `*_all_results_with_alternative_alignments.tsv`: a separate diagnostic table that the importer
  rejects.

With ABE annotation enabled, the v2.1.14-tagged image and v2.3.2 smoke outputs had 86 columns in the
same order. The older tagged image self-reported v2.1.13, so it is retained only as format-comparison
evidence. All 21 fields in the Base Edit Safety required row format retained their names; no
required field was renamed, added, or removed. CRISPRme+ provides additive seed/non-seed, CRISTA,
and `Susceptible_to_ABE` fields that do not alter the logical row mapping. Native output is a TSV file, so it has no worksheets. Of the older workbook's integrated,
alternative-alignment and private-target sheets, only `1617_integrated_results` is read. Private
sheet names and sample identifiers are not retained.

The coordinate convention was checked with a synthetic CD7 sequence placed at reference offset 10:
v2.3.2 emitted start 10, confirming zero-based native starts. Base Edit Safety records the resulting
half-open interval after validating the declared convention.

`--merge 0` is required before a missing CRISPRme coordinate can be treated as a conditional QC
warning. A truncated, filtered, or merged import remains useful for candidate reconciliation, but
it cannot support that warning.

## Record the input and search settings

Place one JSON declaration beside the native result. Replace every example value and compute the
checksum from the exact input file.

```json
{
  "record_type": "crisprme_import_declaration",
  "schema_version": "1",
  "format_schema": "crisprme-2.3.2-integrated-results-v1",
  "native_input": {
    "file": "hbg.integrated_results.tsv",
    "sha256": "REPLACE_WITH_64_HEX_CHARACTERS"
  },
  "selected_sheet": null,
  "guide_id": "HBG_GUIDE_ID",
  "guide_protospacer": "REPLACE_WITH_20_UPPERCASE_DNA_BASES",
  "genome_build": "GRCh38",
  "coordinate_system": "zero_based_half_open",
  "crisprme_version": "2.3.2",
  "scoring_model": "CFD",
  "variant_databases": [
    {"name": "1000G", "release": "REPLACE_WITH_EXACT_RELEASE"},
    {"name": "HGDP", "release": "REPLACE_WITH_EXACT_RELEASE"}
  ],
  "pam_patterns": ["NNN"],
  "max_mismatches": 6,
  "max_dna_bulge": 2,
  "max_rna_bulge": 2,
  "retention": {
    "completeness": "complete",
    "top_n": null,
    "score_cutoff": null,
    "merge_window_bp": 0,
    "all_candidate_coordinates_retained": true
  },
  "source_command": "pinellolab/crisprme@sha256:17e1... crisprme.py complete-search ... exact command ..."
}
```

Compute the checksum without copying the native table into this repository:

```bash
shasum -a 256 /data/crisprme/Results/hbg-v2/hbg.integrated_results.tsv
```

Declaration fields are copied into the source record. Use non-identifying run labels and paths, and do
not put donor/sample identifiers, credentials, or patient information in filenames,
`source_command`, `guide_id`, or other declaration text.

For a legacy publication workbook, keep its original CRISPRme version or state that the version
is not reported. Do not assign it the CRISPRme+ runtime version. Both the normalized table and
the source record preserve that distinction.

`coordinate_system` is mandatory and applies to the whole file. The importer never guesses it
from a coordinate. Use `truncated_top_n` with a positive `top_n` for a top-N file, or `filtered`
with a structured `score_cutoff` when applicable. Incomplete inputs must set
`all_candidate_coordinates_retained` to `false`.

## Convert the candidate file

Index the same reference FASTA used for source verification, then run:

```bash
samtools faidx /data/references/GRCh38.fa

bes-normalize-nominations crisprme \
  --input /data/crisprme/Results/hbg-v2/hbg.integrated_results.tsv \
  --declaration /data/crisprme/Results/hbg-v2/hbg.crisprme-import.json \
  --profile search-profile.json \
  --guides guides.tsv \
  --intended-loci intended-loci.tsv \
  --reference-fasta /data/references/GRCh38.fa \
  --output crisprme.normalized.tsv \
  --provenance-output crisprme.source.json
```

The normalized TSV uses the existing nomination schema. Its generated source record is schema v2
and carries the verified declaration plus aggregate input/emitted/skipped counts. Recoverable bad
rows are counted by fixed reason; row contents are not written to diagnostics.

When rows are skipped, the importer sets `all_candidate_coordinates_retained` to `false` and
changes `complete` to `filtered`. The source record preserves the original declaration alongside the
one used for the import. The panel loader rejects an incomplete import that still claims
completeness. These counts describe imported rows, not exhaustive genome or allele coverage.

The importer drops `Variant_samples*` columns before constructing a logical row. It never reads a
private-target sheet through the validation adapter, and no sample-level values are copied into
the normalized TSV, source record, ledger, logs, or report.

## Compare with CRISPOR

Pass both normalized sources and both source records to the candidate-panel command:

```bash
bes-candidate-panel \
  --profile search-profile.json \
  --guides guides.tsv \
  --intended-loci intended-loci.tsv \
  --reference-fasta /data/references/GRCh38.fa \
  --nominations crispor.normalized.tsv \
  --source-provenance crispor.source.json \
  --nominations crisprme.normalized.tsv \
  --source-provenance crisprme.source.json \
  --output candidate-ledger.json \
  --panel-output panel-review.tsv
```

Ledger v3 groups sites as `both`, `crisprme_only` or `crispor_only`. Reference and alternative
DNA versions stay attached to their exact location and source row. Scores stay with their source;
they are not averaged. The tool keeps distinct reported DNA versions, uses the supplied editor
and window, and leaves unknown sequence unresolved. REF-supported substitutions receive
alignment-derived genomic-forward changes. ALT indels have no asserted reference interval;
ALT-only rows retain a one-base reported start anchor with an unresolved target span. This
anchor is not an inferred inserted-base coordinate. Alternative alignment axes and multi-change
ALT configurations without phase do not receive exact substrate positions.

[Upstream documents](https://github.com/pinellolab/crisprme-plus/blob/main/docs/INPUT_FORMATS.md#variant-driven-pam-creation-and-disruption)
that CRISPRme+ drops alternate alleles that break the required DNA motif (PAM). Native output
therefore cannot enumerate these destroyed-PAM configurations; their absence from the ledger
is not evidence that no such variants exist. The tests classify supplied REF/ALT sequences as
`destroyed` where appropriate; they do not establish upstream discovery of those configurations.

The ledger keeps the original `Variant_MAF` text with the declared source/release. It does not assign
that text to a particular allele, population, AC/AN, or callability, and does not interpret it
as ALT AF or patient prevalence. This importer does not read genotypes or phase.
The reported highest-CFD sequences do not establish exhaustive haplotypes. Rejected rows
remove the effective import's completeness claim; the source record preserves the original declaration.
See the [validation summary](validation.md#bounded-hbg-pilot) for tests, replay results, and remaining limits.

Regenerate older converted files and ledgers to include these allele and window checks. Source-row
IDs include the retained DNA versions; an unsupported ALT mapping may have no interval.
Ledger v3 also checks measurements for each reported DNA version and keeps unusable results with
the reason and a follow-up question. See the [measurement fields](candidate-panel.md#measurement-table)
and [validation summary](validation.md#bounded-hbg-pilot).

## Legacy deposited-workbook validation

The full Zenodo archive is not tracked. To run the opt-in local validation test:

```bash
BES_CRISPRME_GOLDEN_ZIP=/data/reference/File_S1_1617_search_results.zip \
  python -m pytest -q tests/test_crisprme_import.py
```

The test first checks the pinned ZIP SHA-256. It reads only `1617_integrated_results`; the
alternative-merge and private-target sheets are not imported. The CRISPRme+ v2.3.2 release has no
deposited sg1617 workbook asset. The old workbook therefore validates only the legacy parser and
source-row sequence and coordinate consistency; it does not validate the v2.3.2 runtime or its scores. A new reference result would require a recorded v2.3.2 sg1617 rerun with the stated population
files and settings. Keep the workbook and detailed imported rows outside the public repository;
share aggregate check results and citations under the [data policy](data-policy.md).
