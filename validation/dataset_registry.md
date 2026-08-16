# Validation dataset registry

This table records public datasets considered for validation. The tuning/held-out assignment was
locked on 2026-08-14 ([`SPLIT_PROPOSAL.md`](SPLIT_PROPOSAL.md)), so the `role` column records each
dataset's reserved use. Every row still has `QC status = stub`, which means it is excluded until its
files, checksums, reuse terms, samples, guides, editors, and reference build have been reviewed. A
stub is not validation evidence and cannot change a candidate's evidence level.

Terms and status labels are defined in the [glossary](../docs/glossary.md).

The versioned protocol is [`validation-protocol.md`](validation-protocol.md). Unknown fields remain
`TBD` rather than being inferred from a title, abstract, or neighboring database record.

| ID | source | accession | assay | cell type | nuclease/editor | build | guides | validates layer | role | QC status | known issues |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| D1 CHANGE-seq-BE | [D1-S1](#source-review-records), [D1-S2](#source-review-records) | GSE298535; PRJNA1269627 | CHANGE-seq-BE | TBD | ABE and CBE families; exact constructs TBD | TBD | TBD | TBD | tuning | stub | GSE298535 is an umbrella series; the sample/subseries file list and guide-to-editor mapping remain TBD. |
| D2 reni-cel / EDIT-301 | [D2-S1](#source-review-records), [D2-S2](#source-review-records), [D2-S3](#source-review-records), [D2-S4](#source-review-records) | PRJNA1226931 | Digenome-seq candidate-site search with targeted follow-up assays; exact dataset-to-assay mapping TBD | Manufactured CD34+ HSPC products; source-specific sample mapping TBD | Engineered AsCas12a nuclease; not a base editor | TBD | HBG1/HBG2 promoter guide; exact sequence and identifier TBD | L1 + paralog logic; not L3 | tuning | stub | Table S6's 12 Test/Control pairs have no SRR/SRX/SAMN mapping. The 513 AMPLICON records do not identify a rhAmpSeq site, pair, or HBG1/HBG2 run, and 11 of 12 Digenome WGS records lack a public file location in the reviewed ENA report. No mapping was inferred. Nuclease data cannot validate base-editor outcomes. |
| D3 Casgevy-guide CHANGE-seq | [D3-S1](#source-review-records) | PRJNA625995; GSE149295 | CHANGE-seq | Human primary T-cell genomic DNA | SpCas9 nuclease | hg38 | BCL11A guide labeled in the initial registry as the Casgevy guide; exact source-to-product mapping TBD | L1 stress set, off-locus | tuning | stub | Nuclease data from another genomic location do not validate an HBG-promoter base editor. The assigned product identity remains unverified. |
| D4 BE-Hive + CRISPRon-ABE/CBE | [D4-S1](#source-review-records), [D4-S2](#source-review-records) | BE-Hive: TBD; CRISPRon: CNP0001031, PRJNA1321836, Supplementary Data 1 | Genome-integrated and surrogate target-library amplicon sequencing; benchmark assembly TBD | BE-Hive: mES, HEK293T, U2OS; CRISPRon raw data: HEK293T | Multiple ABE/CBE systems; exact benchmark panel TBD | TBD | Large target libraries; exact file lists TBD | TBD | tuning | stub | Sources combine editors, cell types, and target-library designs. Dataset allocation and overlap with development data must be recorded before use. |
| D5 HBG amplicon data | [D5-S1](#source-review-records), [D5-S2](#source-review-records), [D5-S3](#source-review-records), [D5-S4](#source-review-records) | Mayuranathan: HBG amplicon SRA accession TBD; GSE228822 is not the amplicon deposit. Ravi: Dryad 10.5061/dryad.bzkh1897h; GSE192801 is transcriptome only. Han: PRJNA961909, PRJNA1013694 | Targeted HBG amplicon sequencing, including rhAmpSeq in Mayuranathan; D5-S4 maps a 20-run Han subset to AMPLICON/GENOMIC/PCR; remaining per-study assay mapping TBD | HUDEP-2 and primary CD34+ HSPCs; D5-S4 maps the 20-run Han subset to ENA `cell_line=HSPCs`; remaining per-study sample mapping TBD | Base editors and nuclease comparators; D5-S4 maps the 20-run Han subset to tBE or hA3A-BE3; remaining construct-to-sample mapping TBD | TBD | HBG promoter guides; the Han subset is the Figure S3C `sgHBG-1` context at separately amplified HBG1/HBG2 loci; exact sequence and remaining source-to-guide mapping TBD | TBD | held-out | stub | The opened Mayuranathan article does not name its HBG amplicon SRA accession. GSE228822 contains other assay types. Ravi's Dryad file-to-amplicon mapping remains TBD. The 20-run Han subset still requires local checksum review, reference-build confirmation, exact guide/primer sequences, reuse-term review, and development-data overlap review. |

## Source review records

- **D1-S1.** Lazzarotto, Katta et al., “Sensitive and unbiased genome-wide profiling of
  base-editor-induced off-target activity using CHANGE-seq-BE,”
  [Nature Biotechnology](https://www.nature.com/articles/s41587-025-02948-7).
  `access_level=full_text`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-12. Supports the assay
  name and the inclusion of ABE and CBE experiments.
- **D1-S2.** NCBI GEO
  [GSE298535](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298535) and linked NCBI
  BioProject [PRJNA1269627](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1269627).
  `access_level=database_record`; `audit_status=NOT_REVIEWED`; accessed 2026-08-12. Supports only
  the accessions and umbrella-series warning.
- **D2-S1.** “Nonclinical evaluation of renizgamglogene autogedtemcel for sickle cell disease and
  transfusion-dependent beta-thalassemia,”
  [full text at PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC12925820/).
  `access_level=full_text`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-15. Supports
  PRJNA1226931, the engineered AsCas12a HBG1/HBG2 nuclease context, manufactured CD34+ products,
  the stated candidate-site and follow-up assay family, and the publisher-table structure. It does
  not supply an accession-level mapping for the Table S6 samples.
- **D2-S2.** Official ENA run, experiment, sample, and file records for
  [PRJNA1226931](https://www.ebi.ac.uk/ena/browser/view/PRJNA1226931).
  `access_level=database_record`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-15. Supports the
  525-record inventory, the 513 AMPLICON / 12 WGS assay split, paired-library metadata, and the
  file-location review. It does not establish which
  amplicon runs correspond to the publisher's Table S6 Test/Control pairs, rhAmpSeq sites, or
  HBG1-versus-HBG2 on-target columns; no such mapping was inferred.
- **D2-S3.** Hanna et al., “CRISPR-Cas12a Gene Editing of HBG1 and HBG2 Promoters to Treat Sickle
  Cell Disease,” [New England Journal of Medicine](https://doi.org/10.1056/NEJMoa2415550).
  `access_level=abstract`; `audit_status=ABSTRACT_READ`; accessed 2026-08-15. Supports the RUBY
  reni-cel clinical context only. The abstract does not supply a PRJNA1226931/Table S6 crosswalk;
  the supplementary appendix returned HTTP 403 and was not reviewed.
- **D2-S4.** Frangoul et al., “CRISPR-Cas12a Gene Editing of HBG1 and HBG2 Promoters to Treat
  β-Thalassemia,” [New England Journal of Medicine](https://doi.org/10.1056/NEJMoa2501277).
  `access_level=abstract`; `audit_status=ABSTRACT_READ`; accessed 2026-08-15. Supports the EdiThal
  reni-cel clinical context only. The abstract does not supply a PRJNA1226931/Table S6 crosswalk;
  the supplementary appendix returned HTTP 403 and was not reviewed.
- **D3-S1.** Lazzarotto et al., “CHANGE-seq reveals genetic and epigenetic effects on CRISPR-Cas9
  genome-wide activity,” [Nature Biotechnology](https://www.nature.com/articles/s41587-020-0555-7).
  `access_level=full_text`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-12. Supports CHANGE-seq,
  SpCas9, primary T-cell DNA, hg38, PRJNA625995, and GSE149295. The recorded review does not confirm
  the source-to-Casgevy product identity mapping.
- **D4-S1.** Arbab et al. (Cell, 2020),
  [full text at PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC7384975/).
  `access_level=full_text`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-12. Supports the BE-Hive
  target-library assay, editor-family scope, and listed cell systems. No repository accession was
  verified in this pass.
- **D4-S2.** Yuan et al. (Nature Communications, 2025),
  [publisher article](https://www.nature.com/articles/s41467-025-65200-5).
  `access_level=full_text`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-12. Supports the
  CRISPRon-ABE/CBE benchmark context, HEK293T raw-data context, CNP0001031, PRJNA1321836, and
  Supplementary Data 1.
- **D5-S1.** Mayuranathan et al. (PMID 37400614),
  [full text at PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC10722557/).
  `access_level=full_text`;
  `audit_status=FULL_TEXT_READ`; accessed 2026-08-12. Supports targeted amplicon and rhAmpSeq
  methods, the cell contexts, and the statement that GSE228822 contains CUT&RUN, RNA-seq, and
  Micro-Capture-C rather than the HBG amplicon deposit.
- **D5-S2.** Ravi et al. (PMID 35147495),
  [full text at eLife](https://elifesciences.org/articles/65421). `access_level=full_text`;
  `audit_status=FULL_TEXT_READ`;
  accessed 2026-08-12. Supports targeted HBG amplicon sequencing, GSE192801 as transcriptome data,
  and the all-raw-data Dryad locator 10.5061/dryad.bzkh1897h. File-level amplicon mapping in Dryad
  was not reviewed.
- **D5-S3.** Han et al. (PMID 37989316),
  [PubMed record](https://pubmed.ncbi.nlm.nih.gov/37989316/).
  `access_level=full_text`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-12. Supports targeted
  HBG/candidate-site DNA amplicon sequencing and accessions PRJNA961909 and PRJNA1013694.
- **D5-S4.** Official ENA run, experiment, sample, and file records for
  [PRJNA1013694](https://www.ebi.ac.uk/ena/browser/view/PRJNA1013694).
  `access_level=database_record`; `audit_status=FULL_TEXT_READ`; accessed 2026-08-13. Supports only
  the 20-run Figure S3C candidate subset, ENA's `cell_line=HSPCs`, the deposited tBE/hA3A-BE3 and
  HBG1/HBG2 locus labels, AMPLICON/GENOMIC/PCR metadata, paired FASTQ locators, publisher MD5s, and
  byte counts in the public database records.
  Raw files, local checksum reproduction, exact guide/primer sequences, reference build, reuse
  conditions, biological-replicate mapping, and overlap with development inputs were not verified.

## Dataset admission checks

No row may leave `stub` until its raw files, checksums, reuse terms, sample list, guide/editor mapping,
reference build, and overlap with development inputs are recorded. The `role` assignment was locked
on 2026-08-14. A later assignment requires a new protocol version. Passing these checks admits a
dataset only to its already assigned role.
