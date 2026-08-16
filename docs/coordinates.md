# HBG1/HBG2 promoter coordinates

This is the coordinate table used by every downstream step. `src/base_edit_safety/anchors.py` is a direct copy
of it, and nothing in the code re-derives a position. Checked against RefSeq and ClinVar in June 2026.

## Conventions

- **Build.** GRCh38 is NC_000011.10, GRCh37 is NC_000011.9. All positions here are 1-based.
- **Strand.** HBG1 and HBG2 are on the minus strand, so their promoters sit at *higher* coordinates
  than their coding sequence.
- **Promoter numbering.** Positions like -198 or -113 are counted back from the transcription cap
  site: +1 is the first transcribed base, -1 is the base immediately before it, and there is no zero.
  For these minus-strand genes that makes the arithmetic `genomic position = TSS + |offset|`.
- **Store both allele orientations.** The literature writes promoter alleles on the gene's own
  strand ("-198 T>C"). Genome browsers and VCFs write them on the plus strand ("-198 A>G"), which is
  the reverse complement. A single `ref_base` field would guarantee a strand error, so both are
  stored. The mapping is: sense T>C = plus A>G, sense A>G = plus T>C, sense G>A = plus C>T,
  sense C>T = plus G>A.
- **Identifiers.** The standard database form is genomic HGVS, for example `NC_000011.10:g.5250055A>G`. The
  current transcripts are NM_000559.3 (HBG1) and NM_000184.3 (HBG2); the older `.2` versions are
  aliases only. Do not use a transcript HGVS such as `c.-53-198T>C` as a coordinate: transcript
  numbering counts from the ATG through the 5' UTR, which is not the same arithmetic, and ClinVar
  has invalidated that particular expression anyway.

## Row format

One row per gene per position:

```
genome_build | chrom | pos_1based | genomic_plus_ref | genomic_plus_alt | gene | gene_strand |
promoter_offset | promoter_sense_ref | promoter_sense_alt | editor_class | variant_db_status |
mechanism_provenance
```

## Gene anchors

Use these transcription start sites, not the wider Ensembl gene spans and not older transcripts.

| Gene | NCBI ID | Transcript | GRCh38 span | GRCh38 TSS | GRCh37 TSS |
|---|---|---|---|---|---|
| HBG1 | 3047 | NM_000559.3 / NG_000007.3 | chr11:5,248,269-5,249,857 | 5,249,857 | 5,271,087 |
| HBG2 | 3048 | NM_000184.3 / NG_000007.3 | chr11:5,253,188-5,254,781 | 5,254,781 | 5,276,011 |

Sources: NCBI annotation release RS_2025_08 (records updated 3 June 2026) and Ensembl release 116
(9 June 2026). RefSeq and Ensembl agree on the GRCh38 TSS for these transcripts.

## Modelling window

-210 to -80 from the cap site:

- HBG1: chr11:5,249,937-5,250,067
- HBG2: chr11:5,254,861-5,254,991

The two promoters are nearly identical, so most guides hit both. Every edit is projected onto both.

## Positions in scope

| Offset | Editor | HBG1 (plus) | HBG2 (plus) | Sense | Mechanism |
|---|---|---|---|---|---|
| -198 | ABE | 5,250,055 A>G | 5,254,979 A>G | T>C | creates a KLF1 motif |
| -175 | ABE | 5,250,032 A>G | 5,254,956 A>G | T>C | TAL1 E-box, drives looping |
| -124 | ABE | 5,249,981 A>G | 5,254,905 A>G | T>C | synthetic screen hit; KLF1 |
| -123 | ABE | 5,249,980 A>G | 5,254,904 A>G | T>C | synthetic screen hit; KLF1 |
| -113 | ABE | 5,249,970 T>C | 5,254,894 T>C | A>G | creates a GATA1 motif |
| -122 | CBE | 5,249,979 C>T | 5,254,903 C>T | G>A | screen hit; CBE on-target for gRNA-10/-11 |
| -117 | CBE | 5,249,974 C>T | 5,254,898 C>T | G>A | distal BCL11A motif |
| -114 | CBE | 5,249,971 G>A | 5,254,895 G>A | C>T | distal BCL11A motif |

KLF1 binding at -123/-124 was shown for the double substitution, not for either position on its own
(Ravi 2022), which is why both are treated as intended edits.

## Variant database status

Store identifiers per gene because the records differ between HBG1 and HBG2. An rsID assigned to
one gene is not copied to the other.

**db_confirmed** (ClinVar accession and rsID):

| Position | rsID | ClinVar |
|---|---|---|
| HBG1 -198 | rs35710727 | RCV001814966 / VCV000015031 |
| HBG1 -117 | rs35378915 | RCV000016173 / VCV000015030 |
| HBG1 -114 | rs281860601 | RCV001814969 or RCV000016179 / VCV000015035 |
| HBG2 -175 | rs63750654 | RCV001814961 / VCV000014983 |
| HBG2 -114 | rs34809449 | RCV000016130 / VCV000014990 (multiallelic: store the G>A allele rather than only the rsID) |

**db_unverified**: shown experimentally in CD34+ HSPCs, with a coordinate derived from RefSeq, but
with **no ClinVar record** for that gene. The workflow may target it, but does not describe it as a
catalogued variant.

A dbSNP rsID may still exist, and for two of these it does. A dbSNP submission is not curation, so it
does not promote the status. The five rows are:

| | dbSNP | ClinVar |
|---|---|---|
| HBG1 -175 | rs34153902 (and rs34942711, an unmerged duplicate of the same allele) | none. Do **not** borrow HBG2's rs63750654 |
| HBG2 -198 | rs1261339811 | none |
| HBG2 -117 | none | none |
| HBG1 -113 | none | none |
| HBG2 -113 | none | none |

HBG2 -198 and HBG2 -117 were missing from this list while `anchors.py` marked them unverified. A
curator working from this page alone would have skipped both.

All five rows were re-queried by position on 2026-07-26, not by exact-HGVS string, because SPDI is
0-based and HGVS is 1-based and an off-by-one could have hidden records. It had not: the position
query returns the same two rsIDs. The search returned no dbSNP record for the other three rows on
that access date.

**synthetic**: -123 and -124. Screen hits, not natural HPFH alleles.

## Out of scope

Kept here as explicit negative cases:

- **-110 A>C (sense)**, plus-strand T>G. It hits an NF-Y motif, but it is a transversion.
  HBG1 ~5,249,967 / HBG2 ~5,254,891. Note this is a *different allele* from the -110 A>G ABE
  bystander, which is in scope. Both tables are keyed by (offset, ref, alt) so the two cannot merge.
- **-114 C>A and -114 C>G.** Transversions. C>A is Doerfler's cellular allele and C>G its EMSA
  oligo; neither is the repo's C>T.
- **-69 C>T and -70 C>T.** Substitutions first reported by Ravi 2022, outside the -210 to -80 modelling window.
  Recorded so that omitting them is a decision rather than an oversight. Widening the window to
  reach them is a separate call.
- **Delta-13**, roughly -114 to -102. A 13 nt deletion, modelled with a Cas9 nuclease. Base editors
  do not make it.

## Conversions this tool cannot predict

Ravi 2022 reports an observed **-117 G>A** CBE conversion *outside* the protospacer and outside
the editing window. BE-DICT is only asked about positions inside the window, so no prediction will
ever contain it. It is recorded in `anchors.OBSERVED_OUTSIDE_WINDOW` rather than beside the
within-window bystanders, because filing it there would imply a capability the tool does not have.
It is distinct from the -117 CBE anchor above, which is the intended within-window edit.

## Positions with the highest evidence level

Build and demo on these first:

- **-175**, **-113**, and the **-123/-124** pair have the highest exact-allele evidence level in HSPCs.
  HBG2 -175 has a database record and HBG1 -175 does not. The exact -113 A>G was installed with
  ABE8e in patient CD34+ HSCs and in an in-vivo human-locus model.
- **-117** and **-114** were covered by tiling screens and motif-disruption work, not by an
  experiment on the isolated exact natural allele. They are labelled accordingly.

## The spacing check

Matching HBG1 and HBG2 positions must differ by exactly 4,924 bp, which is the TSS separation
(5,254,781 - 5,249,857) on GRCh38 and the same on GRCh37. `anchors.check_paralog_spacing()` runs this
when the coordinate table loads and stops the run if it fails. It once caught an off-by-one that had put HBG2
-114 at HBG2 -113's coordinate.

The check catches shifted rows. It does not identify the gene copy, so do not assign HBG1 or HBG2 by
adding or subtracting 4,924. When processing reads, require a flank that can only
belong to one gene; otherwise label the read `HBG1/HBG2_AMBIGUOUS`.

## Citations

A clinical report establishes that a variant exists. A modern paper explains the mechanism. The
mechanism papers came decades later, so do not read them as the same claim.

- **-198**: clinical Tate/Wood/Weatherall 1986 PMID 2430647 | KLF1 mechanism Wienert 2017
  PMID 28659276 (plus the competition model, Doerfler 2021 PMID 34341563) | HSPC ABE Mayuranathan
  2023 PMID 37400614
- **-175**: clinical (HBG2) Surrey 1988 PMID 2449926; (HBG1, in cis with HBG2 -158) Coleman 1993
  PMID 7679879 | TAL1 E-box mechanism Wienert 2015 PMID 25971621 | HSPC ABE Mayuranathan 2023
  PMID 37400614
- **-123/-124**: Ravi 2022 PMID 35147495 (synthetic screen hit; KLF1 shown for the double mutant;
  HSPC-derived erythroid cells)
- **-117**: clinical Waber 1986 PMID 2417646 | distal CCAAT binding Superti-Furga 1988 PMID 3181130
  (plus Doerfler 2021 PMID 34341563) | region tiled by Ravi 2022 PMID 35147495; motif base editing
  Han 2023 PMID 37989316
- **-114**: clinical (HBG2) Fucharoen 1990 PMID 1698280; (HBG1) Oner 1991 PMID 1704803 | framework
  Doerfler 2021 PMID 34341563 | regional and motif-level HSPC base editing Ravi 2022 / Han 2023
- **-113**: clinical Amato 2014 PMID 23621512 | GATA1 mechanism Martyn 2019 PMID 30617196 (plus
  Doerfler 2021 PMID 34341563) | exact -113 A>G with ABE8e Li 2022 PMID 36006707; also Mayuranathan
  2023 PMID 37400614
