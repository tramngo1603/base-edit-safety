# Validation Protocol — HBG1/HBG2 Base-Edit Safety Triage (v1)

Status: DRAFT — pre-registration. Freeze this document **before** viewing any
prediction-measurement agreement result. Once frozen, changes to thresholds are not allowed; only the
change log in §10 may grow.

Terms used below are defined in the [glossary](../docs/glossary.md).

- Protocol version: v1 (bump to retire a frozen split; see the §3 lock and
  `tests/test_split_lock_immutability.py`)
- Protocol author: {{name}}
- Date frozen: {{YYYY-MM-DD}}
- Code release under test: `v{{x.y.z}}` (commit `{{sha}}`)
- Workflow reference build: {{GRCh38}}
- Strand convention: all coordinates stored as reference-forward; the β-globin
  cluster is **minus strand** — all guide/PAM/window logic must pass through
  `anchors.py` (single source of truth). Protospacer positions are 1-based,
  5'→3' from the PAM-distal end.

---

## 1. Purpose and scope

This protocol defines how to validate the v1 safety-triage records for
base-editing guides at the HBG1/HBG2 fetal-globin promoters, using **published,
measured datasets only**. Validation here means agreement between predictions and public
measurements, not new laboratory confirmation.

In scope (the information the record contains):

- L1: sequence-matched off-target candidates from a Cas-nuclease-style search
- L2: editor-specific off-target candidates within the supported editing window
- L3: on-target outcome + bystander prediction within the editing window
- L4: location-specific concern flags (possible targeting of both HBG1 and HBG2, deletion between
  the genes, and clonal-hematopoiesis driver-gene overlap)
- L5: evidence-tier labeling (`predicted_only` / `analog_supported` /
  `empirically_measured`) and assay recommendations

Out of scope for v1 validation: guide-independent deamination and RNA
off-targets (recorded as flags only — explicitly *not* numerically validated),
delivery-context effects, clinical interpretation.

## 2. Dataset registry

Every dataset used in validation must appear here **before freeze**. The tuning set may be inspected
during development. The held-out set is evaluated once after the code and thresholds are frozen.

| ID | Source (citation) | Accession / table | Assay | Cell type | Nuclease/editor | Build | Guide(s) | Validates layer | Role | QC status | Known issues |
|----|-------------------|-------------------|-------|-----------|-----------------|-------|----------|-----------------|------|-----------|--------------|
| D1 | CHANGE-seq-BE, Lazzarotto & Katta et al., Nat Biotechnol (2026) | {{GEO/SRA}} | CHANGE-seq-BE | {{cell-free/ genomic DNA}} | ABE/CBE variants | {{build}} | {{list}} | L1, L2 | {{tuning/held-out}} | {{pass/warn}} | {{notes}} |
| D2 | reni-cel (EDIT-301) off-target characterization, Editas publication(s) | {{acc}} | GUIDE-seq-family + targeted amplicon | HSPCs | AsCas12a (nuclease, not BE) | {{build}} | HBG1/2 guide | L1 + paralog logic (L4) — **not** L3 | {{role}} | {{qc}} | Cas12a 5'-PAM, window N/A |
| D3 | Casgevy-guide (BCL11A enhancer) CHANGE-seq | {{acc}} | CHANGE-seq | cell-free | SpCas9 nuclease | {{build}} | BCL11A-e | L1 stress test (off-locus) | tuning only | {{qc}} | Not HBG locus — layer test only |
| D4 | BE-Hive (Arbab et al., Cell 2020) and/or CRISPRon-ABE/CBE benchmark data | {{acc}} | amplicon NGS | HEK293T | ABE7.10/ABE8e/BE4 | {{build}} | many | L3 bystander/outcome | {{role}} | {{qc}} | Cell type ≠ HSPC |
| D5 | HBG-promoter editing papers (e.g., Wu et al. Nat Med 2019; Weber et al. Sci Adv 2020; {{others}}) | {{acc}} | amplicon NGS | HSPCs / HUDEP-2 | various | {{build}} | HBG promoter guides | L3 at-locus outcomes, L4 geometry | {{role}} | {{qc}} | Small n; formats vary |

Registry rules:
- Inclusion: dataset must contain measured per-site or per-allele values for a
  defined guide sequence, with sufficient metadata to populate every column.
- Exclusion: review-article re-tabulations without primary data; datasets
  whose guide sequence cannot be unambiguously reconstructed.
- Source record: every ingested value carries `{dataset_id, table/figure,
  accession, tool+version, ingestion_date}`.
- Row-count check: parsed row counts are spot-checked against the
  source publication (≥3 tables per dataset, by hand).
- True negatives: record per dataset whether "assayed but not detected" sites
  are available. If not, precision metrics are marked *not computable* for
  that dataset (see §4 M2).

## 3. Tuning / held-out split (pre-registered)

- Split locked: 2026-08-14 by Tram Ngo (protocol owner), protocol version v1.
  Immutable from this commit; a later split requires a separately versioned
  protocol, not an in-place reassignment (enforced by
  `tests/test_split_lock_immutability.py`).
- Tuning set: D1, D2, D3, D4
- Held-out set: D5 — the at-locus HBG × base-editor × primary-CD34+-HSPC slice,
  the registered dataset that most closely matches the claim under test; HBG-specific logic is
  the claim under test.
- Rule: no threshold, weight, rule, or constant may be adjusted after first
  inspection of held-out prediction-measurement agreement. Violations void the validation;
  log and restart the freeze.

## 4. Metrics and pre-registered thresholds

All thresholds below are **proposed defaults — adjust before freeze, then
lock**. Report every metric per dataset and per assay; report cell types separately
where metadata allows.

| ID | Metric (definition) | Layer | Dataset(s) | Threshold | Result |
|----|--------------------|-------|-----------|-----------|--------|
| M1 | **Off-target site recall**: of measured off-target sites (positive per the source publication's own cutoff), fraction present anywhere in the tool's candidate list | L1 | D1, D2, D3 | ≥ {{0.90}} (GUIDE-seq-class); ≥ {{0.80}} (biochemical assays incl. CHANGE-seq-BE, which detects weaker sites) | {{pass/fail, value}} |
| M2 | **Candidate precision / rank quality**: fraction of measured sites ranked within the tool's top {{100}} candidates. If true negatives unavailable for a dataset, report rank metric only and mark precision N/A | L1 | D1–D3 | ≥ {{0.80}} | {{}} |
| M3 | **Bystander agreement**: Spearman ρ between predicted and measured per-position editing frequencies within the window, per editor class | L3 | D4, D5 | ρ ≥ {{0.6}} per editor class | {{}} |
| M4 | **On-target outcome agreement** (if emitted in v1): Spearman ρ of predicted vs measured total edited-allele frequencies | L3 | D4, D5 | ρ ≥ {{0.6}} | {{}} |
| M5 | **Both-gene targeting flag sensitivity**: reference-set guides known to target both HBG1 and HBG2 must all be flagged | L4 | hand-verified reference set | 100% (no false negatives) | {{}} |
| M6 | **Intergenic-deletion warning**: present for every both-gene-targeting guide; specificity ≥ {{0.95}} on reference-set negatives | L4 | hand-verified reference set | 100% / ≥0.95 | {{}} |
| M7 | **CH-driver overlap flag**: reference-set guides with a predicted off-target in a CH driver must all be flagged; zero false negatives | L4 | hand-verified reference set | 100% | {{}} |
| M8 | **Evidence-level accuracy**: manual review of {{20}} sampled records against registry evidence | L5 | all | ≥ {{95%}} correct | {{}} |
| M9 | **Complete result accounting**: every reference-set and registry guide yields either a complete record or an explicit, logged error | all | all | 100% | {{}} |

Interpretation rules (write results section against these):
- A failed threshold is reported, not hidden; failures feed the failure
  taxonomy (§5) and the limitations section (§8).
- Assay detection limits (~0.1–1% depending on assay) mean some experimentally detected
  sites are below the candidate-search method's reach; misses
  are classified before being counted (see §5).
- Agreement on HEK293-derived data (D4) does not establish agreement in
  HSPCs; say so wherever D4 supports a claim.

## 5. Failure categories

Every miss/failure gets exactly one root-cause label:

| Code | Category | Examples |
|------|----------|----------|
| F1 | Strand/orientation error | minus-strand flip, PAM-side inversion |
| F2 | Paralog multi-mapping | HBG1↔HBG2 mis-assignment, under-called paralog copy |
| F3 | Coordinate build error | hg19/hg38 mixing, lift-over loss |
| F4 | Editor-window convention | 0/1-based drift, deaminase-specific window mis-set |
| F5 | Assay detection limit | measured site below the candidate-search method's plausible reach |
| F6 | Prediction error with verified inputs | prediction wrong with correct inputs and conventions |
| F7 | Curation error | parser drop, table misread (fix in curation, re-run) |

Policy: F1–F4 and F7 are bugs — fix, re-run **tuning set only**, log in
§10. F5/F6 are findings — report, do not engineer around them after freeze.

## 6. Freeze procedure

1. All registry QC statuses = pass/warn; thresholds locked in this document.
2. Tag code: `git tag v{{x.y.z}}-validation && git push --tags`.
3. Run the full pipeline on the tuning set, record metrics, and complete at most one tuning-only correction cycle.
4. Run held-out set **once**. Record outputs to `data/outputs/validation/`
   with a run record containing the commit, timestamps, tool versions, and checksums.
5. Results are published as observed, whether they pass or fail.

## 7. Golden set (hand-verified reference cases)

{{Table of 10–20 known-answer cases: guide sequence, expected on-target edit,
expected flags (dual-engagement, deletion warning, CHIP), source of truth
(paper + figure). A reviewer must verify every expected value directly against
the primary source and record that review.}}

## 8. Limitations (pre-written; expand only)

- All validation is retrospective and uses public data. Prospective confirmation is outside this
  release.
- Measured datasets cover limited assays, cell types, and guide sets rather than the full set of
  off-target sites.
- Guide-independent deamination and RNA off-targets are flagged, not
  validated, in v1.
- The tool triages; it does not clear any candidate for the clinic.
- Part of the held-out dataset (D5) shares a source paper with entries in the bystander
  table (Han 2023 aggregate findings). Aggregate-level findings informed rules frozen
  before the held-out split; the per-replicate held-out values were not inspected.
- The D5 held-out slice (Han 2023 Figure S3C) treats rep1–rep5 as independent biological
  experiments per the figure legend; donor-level independence is unverified, so
  independence-assuming statistics may overstate replicate independence where donors are
  shared.

## 9. Validation report skeleton

1. Summary: metrics table with pass/fail, produced by the evaluation command
2. Per-dataset prediction-measurement agreement: plots and tables produced by the evaluation command
3. Failure log: every failure, its F-code, and how it was handled, with a reviewer note
4. Held-out results: a separate section that states the one-run rule; agreement is reported for two
   groups—positions and editor contexts with no rule derived from the source paper (clean group),
   and positions with such a rule (source-overlap group, stored as `partially_entangled`).
5. Limitations (from §8)
6. Appendix A: dataset registry (from §2)
7. Appendix B: source and run details (accessions, versions, checksums)
8. Appendix C: environment (container/requirements lock, tool versions)
9. Human sign-off (below)

## 10. Change log (post-freeze only)

| Date | Commit | Change | Metric-affecting? | Reason |
|------|--------|--------|-------------------|--------|
| {{}} | {{}} | {{}} | {{yes/no}} | {{}} |

## 11. Human sign-off

- [ ] I hand-verified every golden-set expected value against primary sources.
- [ ] I verified the tuning/held-out split was respected.
- [ ] I reviewed every F1–F4/F7 fix and its re-run.
- [ ] I reviewed the held-out run record (single run, correct commit).
- [ ] The numbers in §4 are the numbers produced by the pipeline, unedited.

Signed: {{name}} — {{date}}
