# Donor-genome coverage

## A finished run is not full genomic resolution

A completed run accounts for every input record and requested chromosome, including explicit exclusions. It does **not** mean the whole genome was resolved. Large regions may carry
no evaluable data at all, and they stay out of the comparable set.

## Fixed denominators

Two denominators are fixed and never shrink because of missing data:

- **Full reference bases** — every base of the reference contig.
- **All possible oriented 23-base reference windows** — every 23-mer start position on
  both strands of the reference, before any PAM or mismatch filtering.

Reference-coordinate coverage uses these fixed denominators. Copy-specific modeled
windows are reported separately, with their own model and nominal-copy counts.

`modeled_reference_bases` requires a supported state on every declared copy at a reference position. A resolved deletion counts as a known absence there. An unknown copy excludes that position from this total even if another copy can still contribute searchable windows.

## Reference-grid windows vs. copy-specific modeled windows

Reference-grid windows live on the reference coordinate system. Copy-specific modeled
windows are built per declared sample copy and **include indel paths** — insertions and
deletions shift coordinates, so a modeled window may have no reference-grid counterpart.
The two window sets are different objects; never compare their raw counts directly.

## Comparability

A window is **comparable** only when:

1. its reference window is exactly evaluable (no missing/ambiguous reference bases), and
2. **every** declared sample copy has an evaluable path with the same exact
   reference coordinates through that window.

Locus gain/loss calls are made only on comparable windows. Sequence change is reported
**only at retained comparable loci**. Insertion or deletion paths without an exact
reference counterpart are **unpaired** — they are not gains or losses.

## Benchmark membership vs. confident genotypes

Membership in the source benchmark BED regions is a *region* statement; it is not the
same as a base being confidently genotyped. When a variant row is absent from the VCF,
any modeled reference sequence there is a **conditional reference fill**, not evidence of a
confident homozygous-reference call.

## Deletion overlap (`*`)

A `*` allele marks overlap with a deletion. It can be resolved to a concrete state from
compatible called alleles. If an exact deletion remains unresolved, its span may be
**locally bounded** by flanking called-present bases on all copies; when neither explicit event bounds nor accepted flanking calls establish a boundary, the relevant chromosome remains excluded. Breakends also exclude an admitted partner chromosome; unknown partners retain a wider exclusion.

## Missingness and phase boundaries persist

Missing regions and phase-block boundaries are not smoothed over. They persist into every
downstream count and exclude affected loci from the comparable set.

## Sex chromosomes

X/Y analysis requires an **explicit ploidy source** and explicit interval declarations,
including how PAR (pseudoautosomal regions) is represented. A chromosome name ("X",
"chrY") or a demographic label is **never** a valid copy-count source.

## Supplied example input

The existing HG003 NIST v4.2.1 GRCh38 file covers **chromosomes 1–22 (autosomes only)**.
Any report built from it must name two gaps: there is **no X/Y callset**, and there is
**no sample-specific PAR/nonPAR copy declaration**.

## What changed

The historical 0.4.0 workflow excluded 15 autosomes containing 306 spanning-deletion markers. The original and normalized VCF both contain 259 markers with an independent covering literal deletion and 47 without one. Normalization did not create the orphan markers. All 47 have flanking accepted calls supporting present bases on both copies; the widest conservative envelope is 3,260 bases. Their exact deletion alleles and breakpoints remain unknown. These are conditional callset bounds, not independently measured structural-variant breakpoints.

The current reconstruction preserves original ALT indexes, GT and PS. It applies compatible overlapping edits to their carrying copies, resolves uniquely supported deletion-marker relationships, and keeps contradictory calls or incompatible phase local. Literal substitutions and simple normalized indels are supported within the declared reconstruction limit (100,000 bases by default). Complex replacements and other bounded unsupported events retain local gaps. Oversized neighborhoods retain a local computational exclusion; exceeding the input-record or global resource limit remains an incomplete run.

For a bounded orphan marker with a complete literal allele on the other copy, that supported copy remains searchable. Phase checks cover the full possible-deletion interval, so an unrelated heterozygote cannot be silently assigned to the recovered copy. A deletion that contradicts a source-called present base stays unresolved. A region can therefore contribute searchable copy windows while remaining excluded from comparisons that require every copy.

`coverage.tsv` provides one row per chromosome. `coverage-regions.tsv` partitions every reference base; `coverage-windows.tsv` partitions every reference window start. A start interval `[a,b)` describes all starts from `a` through `b-1` of 23-base windows, on both strands. It does not describe just the bases covered by those windows. Reasons may overlap, but excluded totals use their union.

Variant disposition columns count records and can overlap: one diploid record can be applied on one copy and unsupported on the other. Exclusions retain variant identities, local copy labels and boundary evidence. The compressed `variant-windows.jsonl.gz` contains concatenated gzip members, each with one original-record cluster; native aliases record the byte offset and member length for exact replay.

The [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf) defines allele indexes, phase sets and spanning-deletion markers. The supplied HG003 file declares VCF 4.2 and retains its original genotype and phase fields; the pipeline does not invent phased calls from the presence of PS alone. The [historical benchmark](research-benchmark.md) remains unchanged.

## Reproduce coverage

Follow the [quickstart](research-quickstart.md), using the same pinned reference, full HG003 VCF and 23-guide panel. Inspect `report.html`, `coverage.tsv`, `coverage.json` and `exclusions.tsv`. Optional `sample.confidence_bed` must have a pinned path/checksum and `sample.confidence_source` must describe what that mask establishes. BED membership is reported separately and never changes the computational denominator.

The handoff's RELEASE.json records the exact source revision, installed-package checks, before/after totals and local reproduction commands. Historical results are preserved; before coverage is derived from the completed baseline's retained materializations using its unchanged installed reconstruction.
