# Validation and benchmark summary

I check this tool in two ways: software tests with constructed inputs, and comparisons with
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

I used the official supplementary workbook for
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

These candidates were already in the authors' public data. My check connects that candidate
list to the follow-up panel; it does not discover new editing signals. The 66 are unmeasured
candidates and must not be described as confirmed off-target edits. A previous top-ten subset was sensitive to tied read
counts and is not used as the headline result.

## V2 bounded HBG pilot

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
imported records can establish, not failed experiments, search sensitivity, time saved or safety.

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
reporting corrections using fresh installed-package runs. Independent researcher review remains
pending. These software checks do not establish experimental validity. Detailed source tables
and local review records are not distributed.

## Software checks

The public tests cover coordinate and strand normalization, the HBG1/HBG2 4,924-bp spacing,
intended-paralog handling, output schemas, prediction provenance, missing-value semantics, and
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
replay scripts are now kept in a local archive. They are excluded from the current checkout
and its normal test run, so the suite is smaller. The README's candidate example now produces
the ledger and TSV; the design file used for the earlier report example is also archived.

On 2026-09-16, the remaining offline suite passed **269 tests with 3 skips**. The 98 archived
regression cases were excluded from this run. The same three optional integrations were skipped.

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
These checks are not a new full-text review of every cited publication.

```bash
python -m pytest -q
```

## What should and should not be concluded

The checks support using this tool to review candidate lists, keep source details attached and
show where experiment information is missing. In the checked panels, simple sequence rules and
short score-ranked lists omitted sites with reported editing.

They do not establish clinical safety, replace editor-specific nomination or targeted sequencing,
estimate patient risk, validate five-mismatch or gapped-site discovery, or show that an unmeasured
candidate is edited in cells. Outcome rates for non-NGG designs remain unevaluated because the
available BE-DICT models are NGG-trained.
