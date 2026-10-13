# Pinned public inputs

Total download: 1,130,871,852 bytes. Reserve 25 GiB free for downloads, preparation
and one full search. Preparation may temporarily retain normalized VCFs as well as
the final VCF. Additional runs need additional output space. The launcher checks
both cache and output filesystems and at least 4 GiB available RAM (16 GiB installed).

| Input | Bytes | SHA-256 |
| --- | ---: | --- |
| UCSC hg38 FASTA gzip | 983,659,424 | `c1dd87068c254eb53d944f71e51d1311964fce8de24d6fc0effc9c61c01527d4` |
| GIAB HG003 NIST v4.2.1 GRCh38 VCF gzip | 147,212,428 | `a4f8fefa826f8d2c9457eacf6f4b2f6bfa29daef429d570ce04254c5e7b1121e` |

Sources: [UCSC hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
and [GIAB HG003](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz).
The bundled 23-guide panel has SHA-256
`eaede910777b7f18aa1df2e541ee542f0b7b791e10f0803c9d950aa9caf43f8f`.

Use the launcher's Download and prepare action, or the public example command in
ADVANCED.md. Preparation extracts full chr1–22/X/Y, uppercases sequence and builds
the FASTA index. The resulting FASTA hash is
`e55dde29c6e54352f68be3973034e63c37b5efd6fbaa878558f6ab6fa8eade0e`,
and its index hash is
`aeaf6f0a845900acea315871a62982711cebc6da5d982f4085b78ad1ba03a0ef`.

pysam 0.24.1's bundled bcftools normalizes with `--check-ref e`, without splitting
multiallelic rows or atomizing alleles, sorts if necessary, and checks all final REF
alleles. It retains original genotype/phase/filter fields and writes one GRCh38
reference header. Compression and header metadata can differ between preparation
locations and times; their exact hashes are recorded in preparation-receipt.json.
Scientific row and coverage comparisons use expected-hg003.json.

Complete raw files in `data/ucsc/hg38.fa.gz` and
`data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` are verified before reuse.
An existing complete prepared cache is also hash-verified. Interrupted `.part` files
and mismatches are retained and never accepted as complete downloads.

No private data, private evidence ledger or genomic records are bundled. The only
network operations are downloads from the two public URLs. User genomic inputs
are processed locally. Public input terms and citations are retained in the source
documentation's permissions page.
