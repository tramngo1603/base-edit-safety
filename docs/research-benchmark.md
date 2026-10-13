# Retrospective primary-reference benchmark — 0.4.0

The completed 23-guide workflow searched all GRCh38 primary chromosomes and modeled public sample HG003 on seven explicitly supported autosomes. It recovered 19 of 29 published human ABE off-target positives exactly; three positives are outside the search geometry and seven remain unresolved. This retrospective benchmark used evidence reviewed during development.

The measured implementation was `2c34ee79a36cfbff4b250b32f640315b8ca24ffb`, a successor of preserved `a97ccc98da60428ce416edd659880bb800771b79`. The preserved 0.4.0 release is `8982108ac9622285e52a9a056b503252c1818ace`. Its release record notes the release revision and wheel checksum. This page reports those historical results. The [quickstart](research-quickstart.md) supplies installation, preparation and the single analysis command; [permissions](permissions.md) records the components and public data used.

The frozen comparison uses 20-base spacers, NGG, both strands, at most six spacer mismatches and no alignment bulges. Cas-OFFinder 2.4.1 uses the supplied `local-workgroup-id-tail-guard-aligned-upload-v3` patch. Both arms use identical rules. The reference arm searches full chr1–22, X and Y: 3,088,269,832 bases, including 150,630,719 unknown bases that cannot support exact candidates. Sample searches cover every admitted variant footprint, including places with no reference candidate. Published off-target coordinates enter only the subsequent evidence join and evaluation.

HG003 NIST v4.2.1 supplies 4,000,097 VCF rows. The complete input is retained and accounted for. Its 306 spanning-deletion `*` records lack consistently resolvable deletion extents; the predeclared recipe therefore excludes entire sample chromosomes containing them. Sample modeling covers **chr5, 8, 11, 12, 13, 14 and 21** (863,156,855 reference bases). The other 15 autosomes are excluded without changing their diploid declarations; X/Y sample ploidy is unknown. Of the source rows, 2,786,627 lie on excluded chromosomes and 1,213,470 on modeled chromosomes. Missingness, unsupported alleles, overlapping calls and incompatible phase remain exclusions. Absence from this VCF is conditional reference fill, not measured callability; the retained GIAB confidence BED is not used to claim otherwise.

The ledger preserves 512 source cases and 425 mapped guide/locus units at 423 physical loci. Its geometry accounting is 275 supported, 233 unsupported and four unresolved; these are not biological-positive denominators. Positive units are deduplicated by spacer, reference span, strand, endpoint and on/off-target role. Repeated donor rows and overlapping studies remain traceable without multiplying recovered loci.

| Cellular endpoint | Published positive units | Exact reference sequence recovered | Not nominated | Unsupported | Unresolved |
| --- | ---: | ---: | ---: | ---: | ---: |
| Human ABE, off-target | 29 | 19 (65.5%) | 0 | 3 | 7 |
| Human ABE, on-target | 20 | 19 (95.0%) | 0 | 0 | 1 |
| Rabbit ABE, on-target | 3 | 0 | 0 | 3 | 0 |
| Human Cas9 indels, off-target | 1 | 0 | 1 | 0 | 0 |

The three unsupported human off-target positives have bulged or non-23-base target alignments. The seven unresolved positives comprise four conflicting published PAM identities, two unidentified positives within a three-locus group, and one unlocated Liang target. The unresolved human on-target assay does not uniquely identify HBG1 versus HBG2. All four members of the separate Liu sequence group were recovered, satisfying the frozen all-four rule; individual editing-rate assignments remain unresolved. The single Cas9 cellular positive is an alternate-allele-created NGG target absent from the reference. Its absence is expected under reference-only input, and its experimental donor is not linked to HG003. The search recovered every resolved human ABE positive within the supported scope.

Biochemical results remain separate. Of 241 Kim nominations without a primary Table 4 cellular test, 104 exact reference sequences were nominated and 137 have unsupported geometry. These are not cellular negatives. The four Lazzarotto allelic contrasts give:

| Biochemical contrast | Reference locus nominated | Published enriched sequence recovered in reference | HG003 outcome comparison |
| --- | --- | --- | --- |
| OT1, HG001 | No; location unresolved | No | Unlinked donor |
| OT2, HG003 | Yes | No; alternate sequence differs | Unresolved sample scope on excluded chr6 |
| OT3, HG002 | No; location unresolved | No | Unlinked donor |
| OT4, HG004 | Yes | Yes | Unlinked donor |

HG003 linkage is to the same named donor only; the original experimental DNA lot and callset are unverified. No suitable donor-linked processed outcome matrix was obtained from the public assets checked within the one-working-day ceiling. That search is closed. There is no matched experimental basis for claiming improved ABE recovery, reduced biological burden or genotype benefit. Untested nominations, unknown detection limits and inconclusive measurements are not negatives; biological precision is not calculated.

| Computational burden or change | Reference-only | HG003 conditional model |
| --- | ---: | ---: |
| All 24 primary chromosomes | 706,588 guide/locus nominations | Not comparable: only seven autosomes modeled |
| Same seven modeled autosomes | 191,771 guide/locus nominations | 195,149 nominated coordinate paths |
| Retained loci on those seven autosomes | — | 188,767 |
| Gained / lost loci | — | 2,681 / 1,126 |
| Retained loci with changed target sequences | — | 2,164 |
| Indel paths without an exact reference counterpart | — | 3,701 |
| Unresolved reference loci within those seven autosomes | — | 1,878 |
| Reference loci outside modeled sample scope | — | 514,817 unresolved |

Sample nominations equal retained plus gained plus unpaired paths; unresolved reference loci do not count as losses. Sequence changes are distinct from gains/losses. The combined output contains 715,143 sequence configurations. Reference burden per guide ranges from 5,361 to 123,866 (median 14,746). Comparing all 24 reference chromosomes against only seven sample chromosomes would manufacture a reduction and is not done here.

| Measured stage | Wall time | Peak resident memory |
| --- | ---: | ---: |
| Complete installed workflow | 474.82 s (7 min 55 s) | 3.92 GiB, parent process |
| Reference native search, included above | 24.95 s | 463.5 MiB, native process |
| Variant-window native search, included above | 6.61 s | 422.1 MiB, native process |

Measurements are from one completed run on an Apple M3 Pro, 12 CPU cores, 18 GiB RAM, macOS 14.4.1, Apple OpenCL GPU `G0`, Python 3.12.12 and pysam 0.24.1. The parent and child peaks are separate measurements, not a simultaneous total. Output files occupy 5.46 GiB. Runtime excludes initial downloads/normalization, installation and independent verification; checksum reuse preparation took 14.97 s. The final run was the sole heavy genomic job, within the measured 10 GiB resident-memory, 14 GiB output and six-hour limits.

Independent enumeration exactly matched 58/58 reference hits in 24 fixed 10 kb regions across all guides and strands. All six predeclared real-data reconstruction checks passed: gained sites on both strands, a retained changed sequence, insertion and deletion paths, and a local coordinate shift. A separate complete-input audit verified all 4,000,097 rows, all 306 stars, original ploidy, and zero sample claims on masked chromosomes. No discrepancies remain. Focused backend checks included a 160,000,003-base mock input crossing the device-upload boundary; the repaired backend recovered all seven expected sites at mismatch ceilings zero and six. An independent build produced the same executable hash.

The one installed-package acceptance suite passed 1,208 tests with 14 optional skips at the initial integrated revision. After the sample-exclusion fix, the two affected modules were rechecked from a fresh installation: 40 passed, with no skips or failures; the other 85 runtime files were byte-identical. This fresh environment then completed the realistic workflow above. The release's documentation-only integration retained the tested runtime and preparation helpers. The skipped checks require opt-in curation assets, optional external tools or unavailable offline network fixtures; reserved outcomes were untouched.

Two incomplete attempts are retained, not scored as final benchmarks. The first exposed an Apple OpenCL upload-alignment defect and failed strict native-row validation. The second completed reference search but correctly withheld the sample arm at an unbounded deletion marker; its automatic partial evidence join was not used to tune outcomes. A complete input-class audit led to the conservative whole-contig exclusions before the final run. The retained final execution record binds unchanged guides, reference, VCF, ledger and scorer, plus the repaired backend and explicit sample scope.

Earlier conflicting summaries used different units: 106 versus 109 executions meant the same 105 Kim computations plus either one seven-donor comparison or four Cas9 site/input-arm runs. The 449 versus 444 rows meant seven donor rows versus two scenarios; 531 mapping entries contained 89 duplicates, yielding 442 loci. The later illustrative claim of seven gained and two lost sites described seven alternate-sequence additions and two reference-sequence absences at one retained locus across seven donors. Its outcome-selected windows were development examples, not unbiased discovery or evidence of biological genotype benefit. This full-reference comparison supersedes those interpretations without deleting their historical records.

Primary sources: [Kim 2019](https://www.nature.com/articles/s41587-019-0050-1), [Liu 2022](https://www.nature.com/articles/s41421-022-00384-4), [Liang 2020](https://www.mdpi.com/2073-4409/9/8/1786), [Cancellieri 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10272994/) and [Lazzarotto 2020, Figure 6](https://www.nature.com/articles/s41587-020-0555-7/figures/6). Individual source locators and every unsupported/unresolved case remain in the encrypted evidence outputs.

For a new installation, follow the [quickstart](research-quickstart.md). Exact historical replay requires the retained request, checksum-pinned inputs, evidence ledger, scoring plans and verification scripts. Their original paths and commands are retained outside the public repository. They are not distributed with the package. The results and check counts above belong to the stated revisions.
