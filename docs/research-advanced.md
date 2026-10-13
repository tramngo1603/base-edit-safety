# Advanced full-reference and sample sequence search

`bes-search` searches complete GRCh38 chromosomes 1–22, X and Y with Cas-OFFinder 2.4.1 plus the supplied compatibility patch. It uses 20-base spacers, an NGG PAM, both strands and at most six spacer mismatches, without alignment bulges. Optional sample reconstruction searches variant-overlapping windows even when no reference candidate exists there.

The public example searches all 24 primary reference chromosomes and processes HG003 across all 22 autosomes, retaining local exclusions where sequence or phase cannot be resolved. Sample X/Y calls and ploidy are not supplied. See [donor coverage](donor-coverage.md) for the current scope; the [historical retrospective benchmark](research-benchmark.md) retains the earlier results.

## Install

Use the [0.5.0 wheel quickstart](research-quickstart.md). It reuses the pinned runtime
and patched backend from the historical 0.4.1 handoff. That bundle's
[Start.command installation](../handoff/QUICKSTART.html) remains available with its
source, licenses and reproducible build commands. The following sections describe
the advanced request interface.

## Prepare public inputs

Choose a directory on researcher-approved protected storage outside Git. The command cannot determine whether a filesystem is encrypted. Keep downloads, VCFs, temporary files, receipts and search outputs there. Run one heavy genomic job at a time. Preparation checks free space before downloading: the configured output allowance plus 10 GiB for new preparation, or just the output allowance when reusing verified inputs. Keep search output on that filesystem or check its capacity separately; other jobs can consume space after the check.

```sh
python tools/prepare_public_example.py \
  --workdir /protected-storage/public-example \
  --backend /absolute/path/to/patched/cas-offinder \
  --device G0 \
  --guides examples/research/guides.tsv
```

The supplied panel contains 23 distinct published spacers. Preparation checks fixed SHA-256 values for the UCSC hg38 reference download and GIAB HG003 NIST v4.2.1 GRCh38 VCF. It extracts full primary chromosomes, retains their sequence in uppercase, and creates a FASTA index. It normalizes the VCF with the pinned pysam version's bundled `bcftools norm --check-ref e`, without splitting multiallelic rows or atomizing alleles. It preserves ALT indexes and source genotype/phase fields, checks row accounting and REF alleles, and verifies that the final header contains exactly one `##reference=GRCh38` entry. It never invents phase.

The result is `request.json`, with checksums for the reference, guides, backend and sample VCF. `sample-scope.json` records the public sample’s supported and excluded chromosomes. The request declares two copies across every autosome. X/Y sample ploidy remains unknown. The sample is HG003; the VCF does not supply experimental DNA-lot identity. Site FILTER values `PASS` and `.` are both admitted as supplied; the example includes 26,310 FILTER-dot rows. Dot is not relabelled as PASS. The explicit-PS policy leaves incompatible multi-site phase unresolved. Reference fill at absent VCF rows is conditional on the retained callset; it is not evidence that the donor has a measured reference genotype.

This public VCF contains 306 `*` spanning-deletion markers. Of these, 259 have an independent covering literal deletion; 47 lack that deletion record but have accepted flanking calls that conservatively bound the uncertainty. Reconstruction evaluates the original allele dosage and phase before resolving each relationship. Unsupported sequence and affected search windows remain excluded locally. All 4,000,097 source records stay in accounting, and all 22 autosomes retain their declared diploid ploidy. Excluded regions cannot produce a loss claim. The [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf) defines `*` as an allele absent because of an overlapping deletion.

Add `--reference-only` to omit the sample from the search request. Preparation still verifies the pinned public sample files. Add `--evidence /protected-storage/evidence.tsv` to link a separately held evidence ledger after the sequence searches. Evidence is not used to choose search windows and is not bundled in the example.

Preparation refuses an existing request or preparation output. Matching raw downloads can be reused; mismatched files and interrupted `.part` files are never overwritten. Failure logs and intermediate files remain for inspection. Choose a new work directory before retrying a failed normalization.

### Reuse verified prepared inputs

An existing directory containing the pinned `primary.fa`, `primary.fa.fai`, `final.vcf.gz` and `final.vcf.gz.tbi` can be used directly:

```sh
python tools/prepare_public_example.py \
  --workdir /protected-storage/replay-preparation \
  --backend /absolute/path/to/patched/cas-offinder \
  --device G0 \
  --guides examples/research/guides.tsv \
  --reuse-prepared /protected-storage/public-preparation/out
```

All four files must match the recipe's exact prepared-file hashes. Reuse neither downloads nor normalizes data and does not copy the large inputs. The request points to their absolute paths, so keep those files in place. Newly normalized files can differ in compressed bytes because normalization headers contain command history; their actual hashes are recorded in their own request and receipt.

## Input contracts

`guides.tsv` has exactly `guide_id` and `spacer_sequence` columns, with unique IDs and uppercase 20-base ACGT spacers. This preparation recipe requires the complete 23-guide panel. The search command accepts a checksum-pinned request with `format_version`, `input_kind`, `assembly`, `reference`, `guides` and `backend`; `contigs`, `max_mismatches`, `limits`, `sample` and `evidence` are optional. Each file specification contains `path` and `sha256`. The generated request supplies the full example, including all autosomal ploidy intervals.

For your own panel, copy this request to a new path and replace the guide path and checksum; do not pass a different panel to the fixed public preparation recipe. File paths resolve relative to the request, or may be absolute. `format_version` is the string `"1"`, `input_kind` is `public_genome` or `research_sample`, and `assembly` is `GRCh38`. The reference must be an uncompressed FASTA with its adjacent `.fai` index. Default search scope is all 24 primary chromosomes and up to 32 distinct guides. Omit `sample` for reference-only analysis. Never relabel a different assembly as GRCh38.

An optional evidence TSV requires `case_id`, `guide_id`, `contig`, `start0`, `end0`, `strand`, `target_sequence`, `stratum`, `outcome`, `citation`, `source_locator`, `eligibility` and `reason`. Case IDs must be unique. Coordinates are 0-based, half-open; `target_sequence` includes spacer and PAM in guide orientation. `eligibility=supported` enables sequence/locus matching; unsupported and unresolved rows remain in accounting. Optional `sample_id` distinguishes matched-donor linkage from other sample nominations. Keep this ledger outside Git.

For a different sample, replace `sample.id`, `variants.path` and `variants.sha256`, and state the actual preparation in `normalization`. The VCF must be sorted, REF-checked against the pinned reference and declare `##reference=GRCh38`; keep original ALT indexes and phase. `genotype_policy` declares `phase` (`explicit_ps` or `vcf_implicit_phase_set`), `min_gq` and `min_dp` (nonnegative integers or null), and `accepted_site_filters` / `accepted_sample_filters` (explicit lists of `PASS` and/or `.`). `ploidy_intervals` contains nonoverlapping `{contig, start, end, copies}` records, with 0-based half-open coordinates and explicit 0, 1 or 2 copies; undeclared regions remain unknown. Derive these declarations from the sample, rather than copying HG003 assumptions. `sample.ploidy_source` records the supporting source and is required for real X/Y declarations. Describe the PAR representation explicitly: for example, whether both PAR copies are represented on X and the corresponding Y interval is excluded with zero declared copies. Adjacent regions with different copy counts exclude crossing windows; windows wholly inside a supported region remain eligible. No copy count is inferred from chromosome or demographic labels.

Optional `sample.excluded_contigs` entries contain `contig` and a nonempty `reason`; they mask the entire sample chromosome without changing its declared ploidy. Compatible carried spanning-deletion markers are linked to their deletion. Otherwise, explicit event bounds or accepted flanking present-base calls can establish a local exclusion. Without defensible bounds the affected sample contig remains excluded. Breakend partner chromosomes are also excluded when admitted; unknown partners require a wider exclusion. Uncarried alternatives do not enlarge a complete supported call's footprint.

Sample reconstruction supports SNVs, equal-length substitutions and simple normalized insertions/deletions within the configured reconstruction limit (100,000 bases by default). Compatible overlaps are evaluated per carrying copy and changed bases. Complex replacements and unsupported symbolic alleles retain explicit gaps. Missingness, original ALT indexes, filters and phase are retained; a bounded footprint is excluded locally, while an unbounded carried or missing event retains the appropriate chromosome exclusion. The final preparation REF check writes to `/dev/null`, not to a captured in-memory VCF.

## Run and inspect

```sh
bes-search --config /protected-storage/public-example/request.json \
  --output /protected-storage/public-example-search
```

For the reuse example, use `replay-preparation/request.json`. The output directory must be new. The request limits search resident memory to 10 GiB, output files to 15 GiB and elapsed time to 21,600 seconds. Override these through preparation's `--max-resident-bytes`, `--max-output-bytes` and `--timeout-seconds` flags. These are search limits, not preparation limits; reserve additional protected storage for raw and normalized inputs. Exceeding a limit is a failed/incomplete run, not a negative search result. Failed native execution receipts and raw output remain available.

| Output | Meaning |
| --- | --- |
| `run-manifest.json` | Completion state, input hashes, runtime and memory accounting; check sample status separately |
| `request.json` | Exact request used for replay |
| `candidates.tsv` | Configuration and locus IDs, exact target sequences, strand, 0-based half-open coordinates, mismatches, coordinate paths and locally supported sample copies |
| `locus-changes.tsv` | Gained, lost, retained, unresolved or unpaired-indel status, with sequence changes separate from locus changes |
| `coverage.tsv`, `coverage.json` | Per-chromosome bases/windows, full denominators, variant dispositions and definitions |
| `coverage-regions.tsv`, `coverage-windows.tsv` | Complete partitions of reference bases and reference-window starts, with exclusion reasons |
| `exclusions.tsv` | Unsupported or unresolved footprints and reasons |
| `evidence.tsv` | Optional source ledger, exact sequence recovery and separate locus nominations |

A gained locus has an evaluable reference counterpart without a nomination. Loss requires an evaluable locus with no nominated supported sample copy. Inserted/deleted coordinate paths without an exact reference counterpart remain unpaired. Local copy labels do not establish long-range haplotypes. A completed reference search can coexist with failed sample scope; do not interpret partial sample accounting as complete.

Alternate contigs, mitochondrial sequence, bulges and structural rearrangements are outside this example. Missing/filtered calls, unsupported alleles, conflicting carried edits, uncertain phase, undeclared ploidy and unknown reference bases retain explicit limitations. Only complete 23-base windows are searched; fewer than 23 remaining bases at a contig end cannot define a candidate. The native reference input is a directory of full chromosome FASTAs, preserving chromosome boundaries.

Keep biochemical nominations, cellular measurements and nuclease indels separate when interpreting linked evidence. Untested candidates are not negatives, and candidate counts do not estimate biological precision.

See [donor coverage](donor-coverage.md) for exact denominator definitions and the difference between computational coverage and benchmark-region evidence.
