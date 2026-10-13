# One bounded computational comparison

After installing the package, run this command from the source checkout:

```sh
bes-compare computational --config examples/comparison/computational/request.json --output /tmp/bes-computational-demo
```

Choose a new output directory each time. Open `report.md`; `candidates.tsv` preserves the
complete candidate comparison, and `warnings.tsv` preserves every warning. This mock example
runs offline without a native search tool or optional genomic dependency.

The wheel also installs the example under
`<environment>/share/base-edit-safety/examples/comparison/computational/`. Use that directory's
`request.json` as `--config` when working outside a source checkout. Input paths resolve from
the request file; the current working directory does not change them.

## What this command does

The command validates a checksum-bound reference window and small VCF, reuses the existing
local chromosome reconstruction, searches the reference and both modeled copies of each
sample, compares candidate identities against the reference, and writes all evidence together.
The supported geometry is one 20-base spacer, NGG, both strands, zero to six mismatches and no
alignment bulges. Simple genomic insertions and deletions remain coordinate paths; they are
not alignment bulges. Paths without an exact counterpart remain unpaired.

This route requires `purpose: computational_only`. Every successful run reports biological
`abstain` with its reasons: no assay evidence, detection limits or donor/specimen linkage.
Other commands retain their study and biological-domain checks.

All inputs use the retained-callset sequence model, including mock fixtures. Absent
VCF records supply conditional reference sequence, not evidence of a donor's reference
genotype. Evidence-supported carriage stays `unknown`. VCF sample IDs do not establish
physical donor, specimen or assay linkage. This command has no callability input.

The report calls a model with no local retained heterozygote
`no_heterozygote_in_retained_model`. The reused raw reconstruction/carriage records call that
local model state `homozygous`; neither label supplies confident reference-genotype evidence.

## Input declarations

Use the included `request.json` as the exact field example. Unknown fields are rejected.

- `format_version` is `"1"`; `purpose` must be `"computational_only"`.
- `input_kind` must be `artificial`, `public_genome` or `research_sample`. Mock inputs
  require assembly `artificial` and a contig starting `artificial_`. Real inputs require
  GRCh38 and a primary autosome (`1`–`22` or `chr1`–`chr22`).
- `region` is a single `{contig,start,end}` in 0-based, half-open reference coordinates.
- `guide` declares `id`, 20-base uppercase ACGT `spacer`, `pam: "NGG"` and `max_mismatches`.
  No editor or biological context is inferred from sequence geometry.
- `samples` lists every VCF sample ID exactly once. No additional undeclared samples are
  permitted in the VCF. IDs and labels use letters, digits, underscores, dots and hyphens.
- `reference` declares a local `.json` path and actual file SHA-256. It uses the existing
  `reference_windows` version 1 format, exactly one window equal to `region`, with its
  sequence checksum. The reference assembly and VCF `##reference=` must match the request.
- `variants` declares an uncompressed local `.vcf`, actual file SHA-256 and the actual
  `normalization` history. Keep original ALT indexes. Literal REF bases are checked
  against the pinned reference; supported indels must pass the existing normalization checks.
- `genotype_policy` explicitly declares `phase` (`explicit_ps` or `vcf_implicit_phase_set`),
  nullable nonnegative integer `min_gq`/`min_dp`, and lists of accepted site/sample filters
  drawn from `PASS` and `.`. Dot remains distinct from PASS. An unset threshold is not a
  passed quality check. An explicit threshold with missing quality leaves the call unresolved.

Missing/partial calls, unsupported ploidy, rejected filters, unsupported carried alleles and
overlapping variants preserve gaps. One heterozygote can resolve an unordered local pair;
multiple nearby heterozygotes require compatible phase. The implicit-VCF option is a declared
phase assumption, not measured linkage. Unknown reference bases are not searched. Unsupported
unbounded variant footprints can exclude the whole region. Out-of-scope VCF rows are counted
and reported as omitted; their REF bases are not validated. Do not infer their status from a
completed in-region search.

## Fixed limits and failure behavior

| Resource | Ceiling |
| --- | ---: |
| Guides / regions | 1 / 1 |
| Region length | 23–10,000 bases |
| VCF samples | 8 |
| VCF rows, including omitted/duplicate rows | 1,000 |
| Supported REF/ALT allele length | 50 bases |
| Each input, including the request | 1 MiB |
| Reconstructed bases / sequence materializations | 250,000 / 17 |
| Search hit rows | 10,000 |
| Complete output bundle | 64 MiB |
| Host time / process peak resident memory | 60 seconds / 512 MiB |

Time and resident memory are cooperative checks, not OS isolation. A bounded operation can
finish after a ceiling and then fail; no completed report is committed on such a failure.
Input/candidate caps bound the individual operations. Unsupported carried alleles are retained
as unresolved scope rather than silently dropped. A resource limit exits 4; invalid input
exits 2 with a plain JSON reason on stderr. Success exits 0. Output creation uses the existing
atomic bundle writer (private directory and files). Inputs are rehashed before committing.
An existing output or input/output overlap is rejected; this command has no overwrite option.

## Evidence and repeatability

`result.json` contains the report's producer facts. Existing JSONL records preserve
`normalized-variants`, `genotype-audit`, `sequence-configurations`, `candidate-configurations`,
`search-hits`, `sample-carriage` and `nomination-changes`. `search-audit.json` accounts for
evaluated/excluded windows and paired paths before filtering; `population-audit.json` accounts
for ingestion and reconstruction. The raw `donor_id` field is reused as a VCF sample identifier
only; it supplies no donor linkage. Editing windows are empty and editability is explicitly
`not_assessed_computational_only`.

Candidate identity includes reference identity, guide geometry, strand, exact allele sequence
and coordinate path. Inserted bases retain their reference boundary and insertion offset.
Reference and alternate configurations are not collapsed into one interval. The candidate TSV
has one row per configuration/sample, including absent and unresolved sample models.
`reference_target_sequence` is a pinned reference counterpart, not a donor genotype. Created
PAMs and changed nominations describe differences between the modeled sequences.

Repeat with a second fresh output directory and compare the manifest's output hashes. With
the same inputs and code, scientific output files and readable exports are byte-identical.
Only `run-manifest.json` has variable UTC timing and measured resource usage. It records
package/source identity and hashes for every other output file. Preserve the request and its
pinned files with the report. The example's independently hand-specified `expected.json`
lists the expected mock sequences and search results.

Protect donor-derived inputs and outputs as described in the [data policy](data-policy.md).
This command does not run native search tools, whole-genome searches or raw-read processing.
It accepts only the guide, region and sample counts listed above.
