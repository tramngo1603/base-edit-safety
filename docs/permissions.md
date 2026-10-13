# Permissions and version inventory

## Research workflow 0.4.0

The table records the components used for the 0.4.0 benchmark. The package contains code,
instructions and public-source preparation commands. Genomic inputs, individual sample
outputs and external binaries are stored separately. Historical build and review receipts
are retained outside the public repository.

| Component | Exact identity and use | Permission evidence |
| --- | --- | --- |
| Base Edit Safety | Version 0.4.0; Python standard-library search orchestration and existing reconstruction/normalization code | Apache-2.0, [LICENSE](../LICENSE). |
| Cas-OFFinder | Upstream 2.4.1, commit `9816b94c20c4cba2e79b039e1e2a6dee684b7b66`, plus `local-workgroup-id-tail-guard-aligned-upload-v3`. The patch fixes non-uniform GPU workgroup indexing, a terminal host-buffer overread, host upload alignment and two device-buffer access flags. Search geometry and scoring rules are unchanged. | [Pinned upstream BSD-3-Clause license](https://github.com/snugel/cas-offinder/blob/9816b94c20c4cba2e79b039e1e2a6dee684b7b66/LICENSE); retained with the [compatibility patch](../tools/cas-offinder-2.4.1-compat.patch) in [Cas-OFFinder-LICENSE.txt](../tools/Cas-OFFinder-LICENSE.txt). The Windows-only dirent header is not compiled on this Mac. No kseq dependency is present in this pinned build. |
| Native runtime and build | Apple M3 Pro G0; macOS 14.4.1; Apple OpenCL 1.2; Apple clang 15.0.0; CMake 4.0.3; optional OpenMP absent | Existing local Apple runtime/toolchain. [macOS Sonoma license](https://www.apple.com/legal/sla/docs/macOSSonoma.pdf). No Apple framework or native executable is bundled. |
| Public-data preparation | Python 3.12.3; separately installed pysam 0.24.1 and bundled samtools/bcftools/HTSlib 1.24 | Installed pysam wheel COPYING is MIT. [pysam notice](https://github.com/pysam-developers/pysam/blob/v0.24.1/COPYING), [HTSlib notices](https://github.com/samtools/htslib/blob/1.24/LICENSE), [bcftools license](https://github.com/samtools/bcftools/blob/1.24/LICENSE). Used locally for indexing and normalization; dependency binaries are not bundled. |
| GRCh38 reference | UCSC root `hg38/bigZips/hg38.fa.gz`, initial GRCh38; full chr1–22, X, Y extracted and uppercased; other sequences excluded | [UCSC source and terms](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/), underlying GenBank sequence. No UCSC annotation track is used. Preserve source attribution and third-party exceptions. |
| Public sample | GIAB HG003 / NA24149, NIST v4.2.1 GRCh38 autosomal benchmark VCF and confidence BED. Source calls retained; no phase, quality or missing call is invented. | [Official release](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/), [GIAB sample information](https://github.com/genome-in-a-bottle/about_GIAB), [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/). Public benchmark analysis only; no identification or private-data linkage. |
| Published evidence | Existing Kim 2019, Cancellieri 2023, Liu 2022, Liang 2020 and Lazzarotto 2020 records, with their original source details and exposure status | Use factual results and citations for local retrospective evaluation. Article/supplement access does not by itself grant an open redistribution license. Publisher files, individual donor records and private curation are excluded from the package. |

The reviewed Mac executable SHA-256 is
`d557d6b51f25027446c6dd0d05aab1e0664cc952c73925e2f6d7adf6e28983d9`;
the patch SHA-256 is
`69821dce797a3d57bf580d44f70e17929839861efd41aabbd2d591471c20cbfa`.
Its banner still reports upstream 2.4.1, so the patch label and executable checksum must
travel with the request. The original and earlier patched executables are preserved as failure evidence.
Only the combined v3 patch identifies the tested Mac release route. Each native row is validated; malformed rows cause
failure rather than silent removal or rewriting. Other devices require their own native checks.

The downloaded reference gzip SHA-256 is
`c1dd87068c254eb53d944f71e51d1311964fce8de24d6fc0effc9c61c01527d4`
(upstream MD5 `1c9dcaddfa41027f17cd8f7a82c7293b`). The prepared primary FASTA is
`e55dde29c6e54352f68be3973034e63c37b5efd6fbaa878558f6ab6fa8eade0e`.
The raw HG003 VCF is
`a4f8fefa826f8d2c9457eacf6f4b2f6bfa29daef429d570ce04254c5e7b1121e`.
Acquisition, build and installed-package receipts retain exact transformations and remaining
output-file hashes. The confidence BED is retained as input history; this workflow does not turn
VCF absence into a measured reference call or claim that the whole genome is callable.

CRISPRme, learned predictors and their result derivatives are not used by `bes-search`.
Their older private records remain historical evidence only.

## Optional integrations and source terms

The following license observations were recorded on 1 October 2026. These integrations are
separate from `bes-search`; installation and output use depend on their own terms.

- CRISPRme+ 2.3.2 (`f33bcba150c116bf294125d4b5c5a7d413661474`) and 2.5.5
  (`c005098e05ebbf83f049d9de5fd2598ceed65685`) have version-specific
  [2.3.2](https://github.com/pinellolab/crisprme-plus/blob/f33bcba150c116bf294125d4b5c5a7d413661474/LICENSE)
  and [2.5.5](https://github.com/pinellolab/crisprme-plus/blob/c005098e05ebbf83f049d9de5fd2598ceed65685/LICENSE)
  license texts describing academic AGPL-3.0 and a separate commercial license. The main-branch
  [license](https://github.com/pinellolab/crisprme-plus/blob/main/LICENSE) reviewed that day
  used MGB Open Access License 1.0, including report-use restrictions. Historical output and
  commercial-use rights remain unresolved; later terms do not establish retroactive changes.
- BE-DICT, CRISPRon-BE/ABE8e weights and CRISPOR are separately installed predictors and search
  tools. The recorded labels were MIT for BE-DICT, BSL 1.1 for CRISPRon-BE, and mixed component
  terms for CRISPOR. Exact retained weight, training-data and output rights remain unresolved.
  See [setup](setup.md) for supported configurations.
- The NYGC 1000 Genomes 20220422 panel requires resource attribution under the
  [IGSR terms](https://www.internationalgenome.org/IGSR_disclaimer/) and applicable
  [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/). Individual genotype exports
  remain outside the package.
- The HGDP `hgdp_wgs.20190516` [producer README](https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/README.data-access.hgdp_wgs.20190516.txt)
  permits sequence analysis and requests the Bergström citation. The
  [IGSR reuse statement](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/README_HGDP_datareuse_statement.md)
  and [CEPH specimen terms](https://cephb.fr/en/hgdp_panel.php) describe separate reuse and
  physical-material conditions.
- The [Liang article](https://www.nature.com/articles/s41467-018-07988-z) records
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) with third-party exceptions.
  The [Mayuranathan article](https://www.nature.com/articles/s41588-023-01434-7) did not have
  an established CC grant in the recorded review. The
  [CHANGE-seq-BE article](https://www.nature.com/articles/s41587-025-02948-7) records
  CC BY-NC-ND 4.0. Their source workbooks and derived individual tables are not distributed.
- NCBI accessions PRJNA506085/SRP169835, PRJNA733110/SRP321601 and PRJNA515331/SRP179672
  follow [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/), which retain
  submitter rights. Molecular-data access does not supply rights to publisher figures or tables.
- Clonal-hematopoiesis annotations use the boostDM-CH 12-gene list and pinned Ensembl gene
  spans. boostDM-CH per-mutation scores are not included; the recorded upstream terms are
  CC-BY-NC. See [coordinates](coordinates.md) and the source attribution in `annotation.py`.

Input rights, specimen access and publication rights are separate. See the
[data policy](data-policy.md) for storage and distribution boundaries.
