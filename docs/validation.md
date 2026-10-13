# Validation and benchmark summary

This tool is checked in two ways: software tests with constructed inputs, and comparisons with
published results. This page reports aggregate findings from those checks. Paper workbooks,
individual data points, copied tables, detailed extracts and raw reads are not distributed here.

## How to read the results

- **Repeating a source check** asks whether the selected rows, labels, coordinates and guides
  match the source. It does not show how many real off-targets a method can find.
- **Candidate recovery within a measured panel** asks which methods proposed sites later found
  edited in that panel. A panel selected by those methods cannot give an independent estimate
  of their sensitivity.
- **Unmeasured candidates** were proposed but did not receive the relevant cellular follow-up.
  They are not confirmed off-target edits.

## CHANGE-seq-BE comparison

This review used the official supplementary workbook for
[Lazzarotto et al., CHANGE-seq-BE](https://doi.org/10.1038/s41587-025-02948-7) and the official
CRISPOR CFD and MIT implementations. The recorded run on 2026-08-29 verified workbook SHA-256
`9b9d7f4493eed9b5050009b412ed4ac93cd3e7d7511c817ea652a4d02f0f47b3`.

The primary comparison is the 4,291-row panel nominated by CHANGE-seq-BE or Digenome-seq and
measured by hybrid-capture sequencing:

| Check | Result |
| --- | ---: |
| Qualifying cellular off-target rows | 95 |
| Nominated by CHANGE-seq-BE | 88/95 |
| Nominated by the union (true by panel construction) | 95/95 |
| Admitted by a ≤4-mismatch canonical-NGG sequence rule | 20/95 |
| Recovered at CFD average rank ≤100 per guide | 35/95 |
| Recovered at MIT average rank ≤100 per guide | 27/95 |
| CFD ROC AUC on scoreable, evaluable rows | 0.825 |
| CFD average precision on scoreable, evaluable rows | 0.130 |
| Published regular-parameter misses reconciled | 12/12 |

The workbook yields 100 qualifying rows labelled as off-targets plus six intended controls, while
the article reports 96 bona fide off-target sites. Five qualifying workbook rows are outside the
primary union. Within that union, 96 rows have a significant flag; one fails the read-count flag,
leaving 95 qualifying rows. Whether this exactly explains the article's intended inclusion list
remains unresolved. There is no established contaminant explanation, and no rows were relabelled.

The 4,291 primary rows include six intended controls. Among the remaining 4,285 rows, 455 are
unevaluable and 632 are unscoreable, including 31 qualifying positives; those exclusions overlap.
ROC AUC and average precision pool 3,235 scoreable, evaluable rows across guides (64 positives and
3,171 nonpositives). Ranking includes all scoreable primary off-target rows, even unevaluable ones.
Ties receive average ranks, so the cutoff can select more or fewer than 100 rows per guide.
Unscoreable positives remain in recovery denominators.

This supplementary-table reanalysis imports significance and read-count flags. Checking the mean
editing difference does not recalculate the authors' statistical tests, FDR, or read-count QC.
Union recovery is true by selection and provides no independent sensitivity estimate or validation
of this tool. CFD/MIT rankings are retrospective extensions; the mismatch/PAM comparison checks
sequences already in the panel. Results cannot establish performance beyond that panel. This is
not independent experimental replication, raw-read reprocessing, or a new genome-wide search.
No threshold was changed against these outcomes.

## HBG public-source replays

These checks reconstruct source material within declared limits. The denominators are different
and must not be compared as if they were five sensitivity estimates.

| Source replay | Accounted for | Scope and interpretation |
| --- | ---: | --- |
| Xu | 38/38 assay labels | Coordinate replay after a documented generic alternate-PAM adapter correction. Because this source informed the correction, it is development evidence, not untouched validation. |
| Han | 81/81 labels | The check showed that a short CFD-only list can omit a site with reported editing. Individual paper rows and their values are not reproduced here. |
| Mayuranathan | 95/95 ordinary two- or three-mismatch sites | A standard search through four mismatches recovered two of four measured ABE8e off-targets. The other two require five-mismatch or gapped-site handling and were experimentally nominated. |
| Martin-Rufino | 57/57 guides | Enumeration is limited to published guides inside the fixed HBG window. |
| Papaioannou | 10/10 rows | Reconstruction covers the declared published top-table rows. |

## Mayuranathan nomination-to-follow-up reconciliation

Cross-referencing the public CIRCLE-seq nomination table from
[Mayuranathan et al.](https://doi.org/10.1038/s41588-023-01434-7) with its targeted-sequencing panel
identified 66 exact four-mismatch sites with an NGN PAM that were present in the experimental nomination data but
absent from targeted follow-up. Sixty-two contain an editable adenine in the configured window.

These candidates were already in the authors' public data. This check connects that candidate
list to the follow-up panel; it does not discover new editing signals. The 66 are unmeasured
candidates and must not be described as confirmed off-target edits. A previous top-ten subset was sensitive to tied read
counts and is not used as the headline result.

## Bounded HBG pilot

The native search and the published-panel replay answer different questions. The native check
used **CRISPRme+ 2.3.2**, verified at source revision
`f33bcba150c116bf294125d4b5c5a7d413661474`, with the published -175 guide and ABE7.10-NG/ABE8e-NG.
It searched GRCh38 chr11:[0,5,260,000), with 315 retained SNVs/small indels in
chr11:[5,248,000,5,257,000), an NGN PAM, at most three mismatches and no bulges. These intervals
are zero-based and half-open. Three symbolic structural variants were excluded and remain unsearched.
The two native rows are the intended reference HBG1/HBG2 sites. No real variant-created site
was demonstrated. Highest-CFD configurations do not establish complete haplotype coverage.

The upstream command completed search/merge but failed during plotting. A separate, unchanged
upstream export step succeeded and produced the imported tables. This supports the recorded
native-output path, not a claim that the complete upstream command succeeded. The separate
tiny synthetic example checks a created-PAM site against direct enumeration of constructed
haplotypes; it does not validate completeness of the real variant search.

The published-panel replay uses Supplementary Tables 9–11 from
[Mayuranathan et al., 2023](https://doi.org/10.1038/s41588-023-01434-7). It accounts for all
1,028 source rows: 683 CIRCLE-seq nominations, 122 Cas-OFFinder nominations and 223 follow-up
rows. Every original and completed interval matches GRCh38 sequence and strand; no source
data row is rejected or skipped.

| Reconciled quantity | Result |
| --- | ---: |
| Distinct nominated locations | 762 |
| Distinct amplified locations | 212: 211 off-targets and HBG2 |
| Locations absent from the follow-up table | 540 |
| Planned locations without a successful assay | 10 |
| Amplified source rows, including the repeated HBG2 control | 213 |
| Source-result rows across two editors | 426 |
| Distinct amplified site/editor entries | 424 |
| Amplified entries with mapped candidate DNA | 230 |
| Amplified entries whose candidate DNA this code does not model | 194 |
| Amplified entries with usable results tied to the tested DNA version | 0/424 |

For the 194 entries this code does not model, the source reference sequences and locations are
still verified; 152 are ordinary targets without gaps. That software limit is separate from the
missing tested-allele and per-site quality details. A reference match does not identify donor
DNA. All 424 entries remain unresolved for allele-specific results. This describes what the
imported records can establish; missing metadata remains unresolved.

The comparison preserves Table 11's shared Yes/No label beside an editor-specific observation.
The four off-target signals are attributed to ABE8e; the corresponding ABE7.10 entries record
the article's reported no-editing observation. Published No labels remain recoverable despite
normalization to `inconclusive` when numeric detection limits are unavailable. Missing-control
flags describe fields left blank in this conversion, not experiments proven to lack controls.
No editing rates, per-site P values, QC passes, donor alleles or phase were inferred.

The two HBG copies remain separate. HBG1 has no successful panel measurement; repeated HBG2
source rows are retained without treating them as independent replicates. Locus discrimination
and rearrangement checks remain required. The article's 121-versus-122 Cas-OFFinder count, exact
mismatch convention and individual figure-label mapping remain open source questions.

Separate computational reviews reproduced the published-panel replay counts and checked both
reporting corrections using fresh installed-package runs. No independent researcher review has
been performed. These software checks do not establish experimental validity. Detailed source tables
and local review records are not distributed.

## Software checks

The public tests cover coordinate and strand normalization, the HBG1/HBG2 4,924-bp spacing,
intended-paralog handling, output schemas, prediction source tracking, missing-value semantics, and
record construction from synthetic inputs. Candidate-panel tests also cover retention of
unclassified perfect matches, manifest-derived intended roles, chromosome aliases and one-base
warnings, duplicate evidence, numeric detection limits, reference FASTA validation, native
Cas-OFFinder and CHANGE-seq-BE/BED6 adapters, source-run checksums, and a 10,000-site union.

On 2026-09-11, two fresh environments passed **367 tests with 3 skips**: first an editable
installation, then the built wheel with source-path injection disabled. The skips were the
unconfigured external guide tools, the opt-in deposited workbook and the forced-offline UCSC
test. A separate live UCSC run passed all four coordinate tests. The 29 installed package
modules matched the source hashes. These checks used Python 3.14.4 on macOS ARM64; other Python
versions and operating systems were not newly checked.

Those are historical checks. The study-specific regression tests, native-output fixture and
replay scripts are excluded from the current checkout and its normal test run, so the suite is
smaller. The README's candidate example now produces the ledger and TSV; the design file used
for the earlier report example is also excluded from the checkout.

On 2026-09-16, the remaining offline suite passed **269 tests with 3 skips**. Those 98 regression
cases were excluded from this run. The same three optional integrations were skipped.

For the public-file checkout on 2026-09-16, a fresh editable installation passed **270 tests
with 2 skips**, including the live UCSC sequence check. A separate installed-wheel run, with
network checks disabled and source-package imports excluded, passed **269 tests with 3 skips**.
All 29 packaged modules matched the exported source. Both README examples passed: two design
records without predictions, and an eight-site candidate ledger with five unresolved results,
three unmeasured sites and zero usable allele-specific measurements. These checks used Python
3.14.4 on macOS ARM64. External prediction tools and the optional deposited workbook were not run.

In the earlier installation checks, the no-data and synthetic ledger/report commands passed.
The installed wheel also replayed the authorized HBG inputs: source accounting, normalized
nominations, measurements, comparison TSV and candidate content matched the corrected replay
record. No upstream search, predictor or new reference download was required. Tests of source
tracking, citation-review limits, coordinates, allele coverage and missing values passed.

```bash
python -m pytest -q
```

## Computational checks

Fresh Python 3.12 editable and isolated-wheel runs each completed with **459 passed, 2 skipped**.
The skips were the unconfigured external guide/prediction workflow and the optional deposited-workbook
check. The live UCSC coordinate check passed. A separate installed-wheel environment ran the
complete synthetic population workflow through console entry points, plus module validation,
without pysam or other genomic dependencies. The README examples and complete mock
ledger/linkage example also ran successfully. The previous baseline remains 270 passed and
2 skipped. There are 189 added cases and no removed cases; 23 were added for the latest review.

The new checks cover missing and partial GTs, original multiallelic indexes, cis/trans phase,
PS namespaces and incompatible blocks, quality filters, reference-callability requirements,
unsupported footprints, both strands, insertion-only targets, deletion junctions, stable
candidate identities, exact evidence context, assay resolution, source tracking and cost accounting.
Full-workflow regressions verify that a mixed literal/symbolic record cannot create an assay
on another contig, and that a valid source ledger with a different spacer under the same guide
ID cannot reduce the assay plan. Insertion-only targets remain evaluable beside an unrelated
heterozygous SNP and at reference boundaries. Fetch errors after the live-test connectivity
gate now fail instead of skipping. Several specimens from one person preserve one genotype
denominator, one activation charge and one donor credit per unit. Paired-window accounting
includes windows with no hit; fixed subset-frequency strata do not change selector plans.
Two mock PAM-created targets are discovered without supplying candidate coordinates.
The activation-cost counterexample produces objective 2 at cost 4 for greedy selection and
objective 3 at cost 6 for the independent exact solver.

The latest regressions prevent scientific settings from drifting through a preparation
template, use one exact-path eligibility rule for paired-window and nomination comparisons,
and reject empty or duplicate selector lists. An actual bcftools 1.24 check on mock repeat
sequence reproduced the cropped-reference normalization problem: adding genuine fixture
context changed normalized placement. Preparation now stops at that crop boundary and records
an unresolved normalization rather than emitting a prepared study.

Evidence checks normalize supported editing-frequency units, reject impossible percentage or
fraction limits from usable credit, and enforce an optional prespecified sensitivity threshold.
Full workflows retain three proposed assays for invalid negative evidence and two for the
matched 0.1-percent example. Unmatched guide names retain explicit diagnostics; definitions
matching under another name do not automatically receive credit. A protected assay review
table exposes target, specimen, assay and evidence requirements. Aggregate reports count
positive, non-detection, conflicting and excluded evidence separately.

The subsequent review fixes two localized issues. Inclusive sensitivity comparisons now use
decimal arithmetic: 0.07 percent and 0.0007 fraction yield identical credit and assay plans at
a maximum of 0.0007. Values immediately below and above that threshold are tested without a
tolerance. Carrier recounts reject incompatible assembly metadata, duplicate ALT entries and
excess sample FORMAT values. Their summaries distinguish expected protocol source from
verified acquisition and unverified local files; arbitrary local inputs receive no NYGC source
attribution. Missing assembly metadata requires the exact recorded acquisition contract or
direct acquisition from the documented source. A caller's checksum alone supplies no source
verification. These cases pass through installed-package tests and the full evidence workflow.

The frozen [public protocol](../examples/population/public-protocol.json) was executed using
the NYGC `20220422_3202_phased_SNV_INDEL_SV` release and the first eight sample IDs in sorted
VCF-header order. The source acquisition covered chr11:[5249500,5255400); the actual search
covered **chr11:[5249700,5255200)**, with zero-based, half-open coordinates. The declared HBG
−198 guide was `GTGGGGAAGGGGCCCCCAAG`, with ABE8e, NGG PAMs, at most four spacer mismatches,
both strands and no alignment bulges. Pysam 0.24.1 and its bcftools 1.24 normalized the bounded
input while retaining original GT allele indexes and unsupported records.
Fresh acquisition and preparation after the second review reproduced the scientific settings
of the earlier generated configuration, followed by all population stages and selector checks.

The acquisition contained 202 variant rows; 191 overlapped the search scope and 11 were outside
it. The scan evaluated 17 reference/donor-copy materializations containing 93,584 sequence
bases. It retained 34 hit rows representing **two intended configurations and zero additional
off-target configurations**. No variant-created off-target was demonstrated within this scope.
The resulting off-target planning universe was empty, so all four selectors produced empty
plans across the declared budgets and seeds. No selector advantage was established. Declared
study-control overhead remains separate from the zero selected-assay cost.

These results use `retained_callset_model`, not complete donor sequence support. The source
does not establish confident reference calls at every unlisted position, specimen availability,
assay feasibility, singleton completeness or independent donor generalization. Source release
notes and limits are available in the
[official release directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/).
Individual genotypes, memberships, source subsets and detailed run records remain outside Git.

The separate BCL11A/SpCas9 rs114518452 regression checks the published GRCh38 minus-strand
target, three spacer mismatches and TGC→TGG PAM change using a checksum-bound reference window
and mock genotypes. A separate query of the exact G>C allele in the same public release
found **97 ALT copies among 6,404 called alleles, in 96 donors**. Source FILTER was dot, recorded
as unfiltered rather than PASS; no new GQ/DP pass was inferred. This verifies public retained
genotype records, not callable whole-target DNA, available specimens or editing outcomes. It
remains a nuclease known-case check, separate from HBG base editing. The original study is
[Cancellieri et al.](https://doi.org/10.1038/s41588-022-01257-y).
The [carrier-count recipe](population.md#public-data-protocol) reproduced these counts both
from the checksum-verified preserved subset and through fresh public acquisition, including
after the assembly and attribution correction. The preserved subset lacks a reference-build header;
its attribution is bound to the documented recorded acquisition, not inferred from its coordinates.

The guarded experimental-reference representation repair passes checks on both strands and
preserves source measurements. Representation alone adds zero donor-measurement credits.
A complete source replay remains blocked
because the recorded publisher workbook and indexed FASTA are missing. Re-evaluating the
preserved 1,524-candidate ledger changes 1,060 representations, including 152 measured entries,
preserves every source measurement, and leaves usable measurement links at zero before and
after. The source ledger remains unchanged. This checks representation semantics without
renewing workbook or reference verification. Historical results above remain unchanged. No
new native external search backend or biological experiment was validated by these computational checks.

The eight-donor HBG example cannot resolve nonsingleton frequencies below 5% and contains no
competing off-target assays. The broader published-data computation below adds competing
assays. An independent researcher-led workflow pilot remains unperformed.

### Broader published-data computational pilot

A separately frozen study sampled 24 donors from the original 2,504-person panel, using
`random.Random(20260918).sample(sorted(eligible_ids),24)` after joining the NYGC release header.
All selected donors belong to distinct published pedigree components. Search regions came
from all HBG rows in supplementary tables 1, 9 and 15 of
[Liang et al.'s EndoV-seq study](https://www.nature.com/articles/s41467-018-07988-z), before
filtering by measured outcome or donor carriage. The source guide matches the existing HBG
−198 sequence. Three sex-chromosome rows were outside the scope of these checks; all 64 autosomal source rows
passed unique local DNA matching, contiguous hg19-to-hg38 conversion and exact hg38 sequence
verification. The source assembly is inferred by these checks, not explicitly named in the paper.

After deduplication, 58 source intervals were expanded by 250 bases per side and merged into
57 search windows. Search covered 29,799 reference bases; materialization covered 41,199.
The profile remained NGG, both strands, at most four mismatches and no alignment bulges.
Required reference context was added for normalization without changing the search scope.
Acquisition retained 2,038 source records, including 1,002 within materialized scope.

The completed search produced 1,205 hit rows representing two intended HBG configurations
and 24 off-target configurations. Twenty-three off-target configurations were present on
reference; one donor-dependent sequence configuration occurred at an already nominated locus.
No new off-target locus was demonstrated by donor variation. Independent reconstruction and
exhaustive search agreed for 2,771 complete materializations, 2,002,909 bases and coordinate
paths, and 2,775,678 oriented guide windows. Twenty-two partial materializations remained outside
that complete comparison because of unsupported calls or overlapping records. Their four
reported hits were independently reconstructed locally, confirming all 1,205 reported hits;
complete negative-search verification is not claimed for their remaining windows.

The frozen budget-10 simulation selected seven assays in one donor: seven assay units, one
activation unit and two overhead units. Seventeen configurations remained uncovered. All
412 selector/budget/seed results passed separate cost and coverage checks; a 12-action case
also matched independent exhaustive enumeration. The selected donor had unresolved sequence
at both intended HBG controls, so the plan is not ready for experimental use. All seven selected
configurations were shared across the panel; the donor-dependent configuration remained
uncovered. The selector comparison illustrates activation-cost sharing. Costs, specimen eligibility and assay capability
remain simulation assumptions, with zero empirical measurement credits.

A separate HG002 challenge used the paired
[GIAB v5.0q draft VCF and benchmark BED](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/)
and its exact named masked GRCh38 reference at chr20:[10000000,10010000). The interval contains
22 variants, 15 heterozygotes and three indels. Four reference-derived computational probes
exercise both strands and variant contexts; they have no demonstrated cellular activity.
Independent replacement and bcftools consensus agreed, as did this tool's three sequences and paths,
239,640 oriented guide windows and 15 hits. All 22 allele dosages also agreed with v4.2.1 inside
the intersected benchmark masks; those unphased calls do not independently confirm haplotype
phase. Direct Q100 assembly comparison was not performed.

This challenge exposed a normalization defect: a jointly normalized multiallelic deletion was
rejected because one ALT could move independently within its repeat. The correction validates
the complete literal multiallelic record and preserves ALT indexes. Four synthetic regressions
cover both phased GT copy orders and records that still require normalization. The corrected
installed wheel passed 463 tests with two optional integration skips; the live UCSC coordinate
test passed. The failed baseline and corrected rerun are recorded separately. Because the case
guided a code correction, the successful rerun is development validation, not a first-pass
held-out result.

Twenty source-specific HBG cellular measurement rows were audited separately. Four of six
source-reported positive off-target alignments lie outside the frozen search profile: two need
more mismatches and two require bulges. These ABE7.10/HEK293T measurements do not provide
matched evidence for NYGC donors in the ABE8e/HSPC scenario and were not imported as usable
donor evidence. Reported zeros remain distinct from detection-limit-qualified non-detections.
The source defined the search regions, so this comparison is not an independent sensitivity
estimate. Individual genotypes, donor mappings, detailed plans and source subsets remain
outside Git.

## Comparative workflow checks

The performance candidate after the public selector timeout has source identity
`e5e01ba1d50b141050346ef4af537ffcd972866682fe2b905221c962d9e9cce7`.
Both fresh installations passed 958 tests with 3 expected optional skips each.
Separate review checked preparation reuse, saved-plan policy enforcement, coverage
cache isolation and preservation of completed cases after interruption. A compatibility
fix keeps the original exception when progress cannot be saved and exception notes
are unavailable; that API absence was simulated on Python 3.12, not tested in a
Python 3.10 installation.

A bounded development profile used the unchanged full catalog of 701,541 records,
90 actions and 2,847 questions. All 181 coverage-helper comparisons matched. The
old helper took 43.598 seconds; cache construction and queries together took 0.199
seconds. These instrumented helper measurements include small concurrent checks
on the same host and do not establish an end-to-end selector speedup. The profile
stopped at its fixed 30-minute limit during baseline preparation; no full selector
plan finished. Its final scientific-input check was not reached.

The public benchmark gate remains incomplete. A new source-bound full execution
and its plan audits are still required. The retained results below describe their
stated earlier sources; the new installation passes do not complete that gate.

Two attempts to execute the repaired source stopped before planning. The first
failed its interpreter preflight because two path spellings of the same Python
installation compared unequal. The next passed validation, then returned a
partial search after 673.31 seconds. Three native launch wrappers rejected the
Docker container-absence check before native execution. The supervisor's direct
parent recorded exit 4. All 28 benchmark cases and the separate default plan
were not started; none can receive a passing plan audit.

The later attempt retained its resource samples and partial output files. Its
original inner Docker inspection responses were not saved. Subsequent read-only
probes reproduce a case-sensitive parsing defect with the pinned Docker binary,
but do not reconstruct the missing historical responses. A corrected execution
must use a separate freeze with the same scientific inputs and numerical limits.
These failures do not establish an end-to-end selector runtime or a negative
search finding.

## Outcome-contract integration checks (30 September 2026)

### Historical producer/report checks

The outcome repair landed at `f291a1dc9d02fcb90d52174855f9129ce19988ca`
(schema v2, outcome gate and migration). The producer/report
follow-up at `b347ea25e11a31b1b3da9f9f713abf301575e17f` (later integrated as
`cbcafc1`) added explicit `adjudicated_outcome` and `adjudication_reason` fields
to the published-source mappings (`-reviewed-outcomes-v2` map revision), and the
comparative report now shows adjudicated outcome, adjudication reason, verbatim
source outcome, detection limit and missing-value reason, and counts current
credit only from active revisions passing the shared outcome and quality gate.
Treated published cells import as `inconclusive` (no reviewed four-state
adjudication yet); GFP/UT control arms import as `untested`; neither pinned
workbook reports a detection limit, so no cell can be `tested_no_signal`.

Checks on the producer/report source `b347ea2`, before evaluator integration
and the current-permission correction, used an offline, fixtures-only environment
(Python 3.12.12, pytest 8.4.2, fresh wheel install): the shared full suite passed
998 tests with 14 skips and no failures. The skips are the standing missing-evidence items:
11 checks on real published files (`BES_COMPARISON_CURATION_RUNTIME` unset),
one unconfigured external guide-design tool, one CRISPRme deposited-workbook
check and one offline UCSC coordinate check. Focused source-import/report tests:
54 passed, 11 skipped. The documented mock workflow was rerun end to end:
budget 18 completes one prospective question with ten operations, twenty
requirement slots and cost 18; budget 17 remains `infeasible_proven`; the
historical mock measurement keeps its exact `9.8000000000000007` percent
value, imports as `inconclusive` under the documented legacy-map fallback and
earns zero credit. The example manifest and report from that source remain available
at `b347ea2`. The later `a529450` copies described below are preserved.
Their manifest excludes the path-containing run manifest. These checks used mock inputs.

### Combined-source installed-wheel checks at a529450

The combined evaluator/report revision
`a5294505388130136bd60c482e0d627cdb2d58a1` adds the evaluator decimal-boundary
rejection (exact fields accept integers or decimal strings only), the outcome
contract adoption above, and a report correction: current empirical credit now
also requires the current study's permission (`research_sample_planning` mode
and an evidence policy allowing empirical donor credit); stored credit stays
visible but cannot override current permission. Package source SHA-256 is
`c613d531a38824265153d961c568bb6fc1e8255c858e5160080567584b1f13a0` (80 files).
Documentation-only commit `e508b65` records these checks; it was not the
full-suite source.

Checks on that exact committed source (Python 3.12.12, pytest 8.4.2, fresh
offline environment, wheel
`15f354c06e03c1c27ca288b1bedc034f47e312bb50ca17ee31d6fd1393b543e5`): installed
imports resolved under the environment's site-packages with source shadowing
disabled, and the installed code fingerprint matched the frozen manifest. The
shared full suite ran once against the installed wheel with the repository
source-path setting overridden: 1031 passed, 14 skipped, no failures. The skips
are the standing missing-evidence items: 11 checks on real published files
(`BES_COMPARISON_CURATION_RUNTIME` unset), one unconfigured external
guide-design tool, one CRISPRme deposited-workbook check and one offline UCSC
coordinate check. Evaluator CLI checks reused the two mock decimal
protocols: the numeric target `0.500000000000000000000000000001` exited 2 with
no output; the quoted decimal string evaluated with three required positives,
`unreachable`. The documented mock workflow ran from the installed
resources twice; all 14 commands exited 0 per run, the 18 deterministic exports
were byte-identical between runs, and the plan and search-audit files matched
with per-run resource fields excluded. Budget 18 completes one prospective
question with ten operations, twenty requirement slots and cost 18; budget 17
remains `infeasible_proven`; current empirical donor credit is zero; the
historical mock measurement keeps its exact `9.8000000000000007` percent
value as `inconclusive`. The historical example manifest and report
were regenerated from these outputs using mock inputs.

### Retained checks before the performance repair
The comparative workflow adds separate design policies, exact donor DNA,
source-bound observations and complete experimental contrasts. Its scientific
package source SHA-256 is
`bba9701cca31230c0ed642b8d3ab4beced30fe062470f97ed38deca27f167f9d`.
This identifies the Python and JSON package files in the recorded dirty checkout;
it is not a Git commit or a hash of the documentation.

Fresh editable and isolated wheel installations of this package source each
passed **938 tests, with 3 optional skips**, using Python 3.12.3 on macOS ARM64. The skips were an
unconfigured external guide-design tool, an opt-in historical deposited-workbook
check, and an offline UCSC API check. Both publication-import tests ran against
the retained physical source files in each environment. Package imports and all
76 package file hashes were checked before, during and after each suite. No
required test skipped. The suites took 91.956 and 92.879 seconds, respectively;
these timings were measured while the public search was also running.

The installation check passed 47 command groups across fresh editable,
wheel and reinstalled wheel environments. It ran the documented legacy examples
and the complete mock comparison. The wheel contained 105 approved
resources in addition to the package and distribution metadata; acquired genomes,
source workbooks, individual source extracts and private validation files were
excluded. The documentation build for that release used the same tested package bytes,
tests, metadata, scientific resources and dependencies. It reused those two full
suites and repeated the focused installation, resource, example and reinstall
checks. Its retained package record identifies both the original tested wheel
and the wheel containing that documentation.

The preceding integrated source (`b9024686`) passed 928 tests with the same three
optional skips in each installation. Its results and early test-run failures
remain in the retained history. The current encoding, window-validation,
record-cache and CLI lifetime changes also have separate affected reviews of
scientific identities, sequence paths, alternative alignments and source tracking.
The full public scale gate remains incomplete after the selector timeout.

The mock contrast has ten required operations and twenty requirement
slots. A separate checker parsed its literal DNA, genotype and resource inputs
and inspected all 1,024 action subsets. Exactly one subset completed the declared
question, at a cost of 18 simulation units; a budget of 17 could not. This tests
paired designs, carrier/noncarrier DNA, controls, prerequisites and shared costs.
It establishes neither real assay prices nor an adequate biological sample size.
Returned-observation checks also distinguish a completed plan from measured
evidence, preserve historical plans, and require a failed mock to be repeated.

### Earlier public comparison and selector timeout

The public comparison has completed search, observation import, evidence
reconciliation, comparison, question construction, the default plan and its
report. Independent sequence and observation checks have passed. The fixed
selector comparison reached its original 12-hour cutoff on 24 September 2026
at 11:43:21 p.m. PDT (25 September at 06:43:21 UTC). Both original processes
exited 4 after the guard sent SIGINT, and the process group was verified absent.
No selector plans were saved. The independently checked default plan remains
available, but the required algorithm comparison and public scale gate are
incomplete. Passing installation checks does not close that gate.

The study uses 100 public samples, with one selected sample per declared source
pedigree group, from the original 2,504-member panel represented in the NYGC
20220422 high-coverage phased release. This metadata rule does not exclude
undocumented relatedness. The shared domain contains 762 merged GRCh38 intervals,
65,569 reference bases across 22 autosomes, and 1,765 retained source variant rows.
Published nominations helped define those intervals. Results therefore cannot
estimate genome-wide coverage, off-target sensitivity or patient prevalence.

The three verified active designs are the Mayuranathan −175 spacer with
ABE7.10-SpCas9-NG or ABE8e-SpCas9-NG, and the −113 spacer with ABE7.10-SpCas9.
They span two spacers but one primary publication. Liang source measurements
remain a separate evidence source; they do not supply an invented additional
design combination. Source checks verified all 255 imported quantitative values
and all 26 published alignment entries, including 13 source-positive labels.
The actual reconciliation produced 765 applicability assessments across the
three designs: 36 exact sequence contexts, 360 related design contexts and 369
unresolved assessments. None supplied empirical evidence for the target donors.
The independent observation check passed all 3,164 checks. Its 94 reported numeric
zeros lack usable detection limits and do not establish absence of editing.

The published-alignment replay retained all 26 entries: three exact matches,
17 outside the declared search policies and six unresolved entries. Among the
13 source-positive labels, three matched, eight were outside policy and two
remained unresolved. Excluded or unresolved positives remain in the evidence
record; they are not treated as negative experimental results.

Each active search admits at most four mismatches, one DNA-bulged base and one
RNA-bulged base, subject to four combined edits, on both strands and under that
design's PAM. All admitted alignment alternatives are retained within the
declared output limits. The pinned native backend is CRISPRme+ v2.5.5, source
revision `c005098e05ebbf83f049d9de5fd2598ceed65685`. It searches complete
reconstructed strings and projects results back to the original donors, copies
and loci. The retained-callset model does not establish reference truth where
genotypes or source support are unresolved.

The original attempt was deliberately interrupted after observed per-sequence
post-processing indicated that three jobs would exceed its one-hour host limit.
No actual timeout, memory cap or disk cap was observed. Its incomplete results,
commands and protocol are preserved. An independently reviewed amendment
increased only host, native and outer wall limits to 7,200, 3,600 and 7,200 seconds.
Native memory and CPU limits remain 4 GiB and two CPUs; all scientific policies,
sources, subjects and audit selections are unchanged. The amended protocol is
`e76ae5451401272e5ff5ed077b8fded290e355b521741c7916035a7b727c6318`.
This is an amended computational evaluation informed by observed runtime,
not untouched holdout performance evidence.

That cold invocation reached its 7,200-second outer search limit after all three
native jobs completed. It exited 4 after 7,204.381 seconds; the benchmark was
incomplete. Its measured peak child RSS was 6,406,258,688 bytes. The failed
attempt and complete native output files are preserved. A separately frozen warm
continuation uses the same code, source inputs and limits, with ordinary cache
validation and a guard that rejects any new native invocation. Its protocol is
`ecd570b06a1306a9cdaf7e5f4de488a23de4cdd43fc9d50477f760a8264dd77f`.
The warm continuation also reached the outer search limit. It exited 4 after
7,208.239 seconds, with peak child RSS of 6,789,300,224 bytes. Its incomplete
attempt is preserved, and all native input files remain unchanged. Neither run
passes the public scale gate. A separate bounded diagnostic was stopped after
739.026 seconds once it had recorded source replay, canonical hit sorting and
repeated window validation. No diagnostic resource cap fired; peak child RSS was
7,649,689,600 bytes. Its instrumentation adds overhead, and the sampled operations
do not establish the exact cause of the earlier failures.

The reviewed repairs avoid repeated validation of owned source windows, reuse a
small bounded cache of validated record templates, and preserve canonical output
while using standard JSON encoding where equivalent. A prespecified sample of
3,072 retained native rows produced identical record, hit, annotation and state
outputs. Mapping took 5.070 seconds with the repairs and 8.813 seconds previously;
these single sample timings do not predict full-run performance. The final CLI
also releases unused native result containers while preserving hits and evidence
bindings. Its separate effect on full-workload memory has not been measured.

A new cold development run is frozen under protocol
`a591876c59b3ccc76c46478ae56d278d76b97d0a5c4127f8a928777831978423`.
The host and outer wall limits increase to 14,400 seconds. Native limits, the
8 GiB host memory limit and all scientific policies remain unchanged. Previous
failures and development observations remain part of the record. This rerun uses
ordinary current-source identities and fresh native execution. The search stage
completed with exit 0 after 8,865.855 seconds, with 6,361,759,744 bytes reported
peak child RSS. The following six workflow stages through the default report
also completed. The final selector command reached its 14,400-second outer limit
and exited 4 after 14,402.560 monotonic seconds. No partial selector result was
saved. The complete failed attempt remains preserved.

A separately reviewed selector continuation uses protocol
`1e22b59a9264a91f2ede83996bb74ac28bbb6bbad6ba99d0640741ab58c6fb63`.
It retains all 28 cases, the same scientific inputs and package source, ordinary
source verification, the 8 GiB host RSS limit and the 32 GiB aggregate disk
limit. Its wall allowance is 43,200 monotonic seconds. It reuses the eight
completed stages without repeating native search. The continuation stopped
after 1,188.795 monotonic seconds at the memory limit. Sampled process-tree RSS
reached 8,987,115,520 bytes; the child process reported a 9,245,114,368-byte peak.
Both exceed the 8 GiB ceiling. The separately reported 16,350,857,536-byte memory
footprint is a different measure. The process group was stopped and verified
absent, and source and prior-output hashes matched before and after execution.
No selector outputs were saved, and the dependent plan audit did not run. This
failed attempt remains separate from the earlier time-limit failure.

A new selector-only attempt is frozen under protocol
`5ae28caeb82b5db7fb1cca11e590ff3ba4ca0a8758a02c538b06569629a7b1b4`.
It raises the selector RSS allowance to 12 GiB on the measured 18 GiB host, keeps
the 43,200-second allowance, and includes the failed selector directory under
the unchanged 32 GiB aggregate disk ceiling. Scientific settings, all 28 cases,
package source and preceding outputs are unchanged. Its separate prelaunch
review passed. This attempt reached the unchanged, dated cutoff without saving
selector outputs. The guard sent SIGINT at 06:43:21.146 UTC and observed both
processes exit 4 at 06:43:27 UTC. It verified the process group absent without
needing SIGTERM or SIGKILL. No benchmark was restarted during recovery.

All three native jobs and their six indexing/search commands completed
successfully. A separate audit verified all 8,766 retained output-file bindings.
The sequential native commands took 4,991.876 seconds in total; recorded native
CPU usage was 8,471.391 seconds, and the largest native container memory peak
was 794,411,008 bytes. These measurements cover the native jobs; the host search
measurement above has a different scope. For the final selector attempt,
`/usr/bin/time` reported 43,205.58 seconds real time, 42,248.07 seconds user CPU
and 583.03 seconds system CPU. Reported child peak RSS was 8,153,825,280 bytes;
the separately reported memory footprint was 16,329,279,808 bytes. Cumulative
sampled peaks were 7,888,257,024 bytes RSS and 11,194,436,435 bytes aggregate disk.
The memory footprint is not the RSS measure used for the 12 GiB ceiling.

The original monitor was interrupted while the original computation survived.
The replacement guard observed the same process and retained the fixed cutoff.
The interval from the last original emitted progress record to recovery timer
start was 612.841759 seconds; the interval between emitted records was
612.990473 seconds. Neither establishes the actual resource-sample gap.
Uninterrupted execution, continuous resource enforcement and the full original
monotonic duration cannot be recovered. The replacement guard ran for
7,133.357 seconds through cleanup. These measurements have different scopes
and are not added together. Sampled peaks cannot certify the unobserved gap.
The strict uninterrupted-monitoring condition remains unmet.

Original results, failed attempts, logs and partial files are preserved. The
missing original wrapper result was not recreated; terminal input verification
is separately attributed to the recovery review. Calendar and monotonic
durations remain distinct, and their observed difference has no established
cause.

### Completed sequence and evidence results

The observed search records have the following counts. They include intended
loci and all admitted alignment alternatives; they are not donor prevalence or
estimates of cellular editing.

| Design | Loci | DNA configurations | Alignments |
| --- | ---: | ---: | ---: |
| −175 ABE7.10-SpCas9-NG | 809 | 857 | 3,334 |
| −175 ABE8e-SpCas9-NG | 809 | 857 | 3,334 |
| −113 ABE7.10-SpCas9 | 6 | 6 | 316 |

Complete totals remain unknown because donor search scope is incomplete. Each
design has one donor with complete modeled search scope and 99 with partial
scope. Intended-target status is modeled for 98 donors and unresolved for two.
Across 200 donor/locus pairs per design, 196 are modeled and four unresolved;
these are not counts of chromosome copies. Retained-callset reference filling
does not establish confident reference sequence.

The two −175 designs tie on these DNA counts and on their variant-effect records.
The −113 design has no demonstrated created, changed or lost candidate in this
domain. Neither a smaller candidate count nor a null modeled effect establishes
a safer or more active editor. Unresolved source regions remain unknown.

The prespecified independent sequence check passed for 16 donor/interval cases,
32 independently reconstructed copies and 30 searchable windows. All 10,684
expected alignments matched both native output and normalized records, with
zero missing or extra alignments. This verifies the selected source-supported
windows; it does not establish completeness in unevaluated regions. The audit
took 490.825 seconds with 790,659,072 bytes peak child RSS, within its original
1,200-second and 2 GiB limits. Two earlier checker failures are preserved. Their
reviewed corrections changed image-reference comparison and lookup of inherited
intended loci; source inputs and expected alignment files remained unchanged.

### Selector comparison incomplete after resource failures

The fixed comparison includes complete-bundle greedy, donor-first, adapted
coverage, canonical bundles, and random bundles with seeds 7, 19 and 43, at
budgets 0, 24, 48 and 96. All 28 results, including nulls, ties and losses, must
remain visible. Action construction is limited to six questions and three
canonical donor alternatives per arm; all question suggestions and construction
exclusions remain visible. Beam and enumeration limits remain
explicit; no exact optimum is claimed for the public question universe.
The question output retains 2,847 suggestions. Construction limits leave 2,841
without constructed actions; this is distinct from biological ineligibility or
insufficient budget. The six questions with constructed actions are paired
design comparisons. This public action catalog does not demonstrate a selected
carrier/noncarrier contrast, evidence replication or mandatory-question seed;
those behaviors have separate mock acceptance checks. The output includes
207 design-by-variant and 207 variant-association suggestions, with no actions
constructed for them under the frozen limits.

The separately saved default plan passed an independent requirement and cost
check. At budget 96 it selects 50 actions using 61 distinct resources, costs
90 simulation units and fully plans five paired-design questions of total
weight five. Each completed question pairs designs within the same two donors
from distinct declared independence groups. The checker reconstructed all 60
action-locus donor copies from literal reference and genotype inputs, checked
all 2,847 reported question completion flags, and verified the action, resource
and unmet-requirement exports. All questions still have no measurements. The
plan is heuristic, with no proven optimum. Its reported planner call took
747.694 seconds; the enclosing plan command, including input verification and
other work, took 3,552.406 monotonic seconds.

This default plan is the 29th audit input, separate from the 28 benchmark cases.
The 28 cases comprise four budgets times seven selector/seed settings; the
separately saved default is not an additional benchmark case. Case 22 repeats
its budget-96 greedy, seed-7 setting. That equality check is unavailable because
case 22 was not saved.

Final accounting is one completed, independently audited default plan and 28
incomplete, unavailable benchmark plans. That tested implementation retained intermediate plans in
memory and publishes the result directory only when the schedule finishes.
There is no durable case-level progress record, so the missing cases cannot be
classified individually as not started, failed or timed out. The run timed out;
internal case progress remains unknown. Missing costs, question counts and
comparison metrics remain unavailable. No selector ranking or superiority
claim is supported.

The retained full exports are large: the readable comparison is 2,406,945 bytes,
the protected donor matrix is 93,839,074 bytes and the detailed search audit is
1,531,090,352 bytes. Use the aggregate validation summary and the mock
researcher packet as entry points. The full donor and sequence exports remain
local research files and are excluded from the package.

The researcher evaluation packet is prepared with equivalent manual and tool
tasks, a source-derived answer key and blank feedback forms. No supervised
researcher session or prospective laboratory experiment has been performed.

## What should and should not be concluded

The checks support using this tool to review candidate lists, keep source details attached and
show where experiment information is missing. In the checked panels, simple sequence rules and
short score-ranked lists omitted sites with reported editing.

Cellular accuracy of five-mismatch and gapped-site discovery remains unmeasured.
Unmeasured candidates retain their status. Outcome rates for non-NGG designs remain
unevaluated because the available BE-DICT models are NGG-trained.
