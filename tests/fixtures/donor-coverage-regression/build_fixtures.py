#!/usr/bin/env python3
"""Generate donor-coverage fixture inputs: synthetic references and case VCFs.

Only inputs are generated here. REF alleles are copied from the generated
reference so every record passes a strict REF check. Expected outcomes are
hand-authored separately in tests/test_donor_coverage_oracle.py and are not
derived by this script.

Sequences come from a tiny LCG so they are deterministic, non-repetitive and
easy to regenerate: state = (state * 1103515245 + 12345) mod 2**31;
base = ACGT[(state >> 16) mod 4].
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def generate(length, seed):
    state = seed
    out = []
    for _ in range(length):
        state = (state * 1103515245 + 12345) % 2**31
        out.append("ACGT"[(state >> 16) % 4])
    return "".join(out)


def flip(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


REFERENCES = {"artificial_1": generate(1200, 424242), "chrX": generate(800, 777), "chrY": generate(400, 888)}

VCF_HEADER = ("##fileformat=VCFv4.2\n##reference=artificial\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")


def row(contig, pos0, ref, alt, gt, fmt="GT", extra="", filt="PASS", info=".", name=None):
    value = gt if fmt == "GT" else ":".join([gt] + extra)
    return "\t".join([contig, str(pos0 + 1), name or f"v{pos0}", ref, alt, ".", filt, info, fmt, value])


def deletion_row(seq, contig, pos0, span, gt, ps=None, filt="PASS", name=None):
    """Pure deletion of seq[pos0+1:pos0+span+1]: REF = anchor + deleted bases."""
    ref = seq[pos0:pos0 + span + 1]
    fmt, extra = ("GT", []) if ps is None else ("GT:PS", [str(ps)])
    return row(contig, pos0, ref, ref[0], gt, fmt=fmt, extra=extra, filt=filt, name=name)


def snv_row(seq, contig, pos0, gt, ps=None, star=False, filt="PASS", name=None):
    ref = seq[pos0]
    alts = [flip(ref)] + (["*"] if star else [])
    fmt, extra = ("GT", []) if ps is None else ("GT:PS", [str(ps)])
    return row(contig, pos0, ref, ",".join(alts), gt, fmt=fmt, extra=extra, filt=filt, name=name)


def anchor_row(seq, contig, pos0):
    """Complete homozygous-reference SNV record: called-present bases on every copy."""
    return row(contig, pos0, seq[pos0], flip(seq[pos0]), "0/0", name=f"anchor{pos0}")


CASES = {}


def case(name, contig, rows, note):
    CASES.setdefault(name, {"contig": contig, "rows": [], "note": note})["rows"].extend(rows)


seq = REFERENCES["artificial_1"]

# C01: 10 bp deletion phased to copy 0; interior record carries the star on copy 0
# and a SNV on copy 1 (unphased GT, the dominant real shape).
case("c01_star_explained_same_position", "artificial_1", [
    deletion_row(seq, "artificial_1", 100, 10, "1|0", ps=1441),
    snv_row(seq, "artificial_1", 105, "2/1", star=True),
], "star carried on the same GT position as a covering literal deletion")

# C02: 120 bp deletion (>50 allele cap) phased to copy 0; interior star + SNV.
case("c02_star_explained_large_deletion", "artificial_1", [
    deletion_row(seq, "artificial_1", 200, 120, "1|0", ps=2442),
    snv_row(seq, "artificial_1", 250, "2/1", star=True),
], "supporting deletion exceeds the 50-base applied-allele cap; footprint still bounded")

# C03: single orphan star bounded by flanking called-present anchors (GT 1/2:
# literal allele on GT position 0, star on position 1). Envelope [351, 450).
case("c03_star_orphan_single", "artificial_1", [
    anchor_row(seq, "artificial_1", 350),
    snv_row(seq, "artificial_1", 400, "1/2", star=True),
    anchor_row(seq, "artificial_1", 450),
], "no covering deletion; span bounded by flanking bases called present on every copy")

# C04: three orphan stars across 17 bp with the star as ALT index 0, one shared
# envelope [451, 560), mimicking the chr6 orphan cluster span.
case("c04_star_orphan_cluster", "artificial_1", [
    anchor_row(seq, "artificial_1", 450),
    row("artificial_1", 500, seq[500], "*," + flip(seq[500]), "1/2"),
    row("artificial_1", 508, seq[508], "*," + flip(seq[508]), "1/2"),
    row("artificial_1", 517, seq[517], "*," + flip(seq[517]), "1/2"),
    anchor_row(seq, "artificial_1", 560),
], "consecutive orphan stars, star allele first; one shared bounded envelope")

# C05: multiallelic record [7 bp deletion, star]; same-record deletion must not bound the star.
case("c05_star_multiallelic_same_record", "artificial_1", [
    anchor_row(seq, "artificial_1", 550),
    row("artificial_1", 600, seq[600:608], seq[600] + ",*", "1/2"),
    anchor_row(seq, "artificial_1", 660),
], "copy 0 carries the record's own 7 bp deletion; copy 1's star is bounded by flanks")

# C06: the only covering deletion record has a rejected FILTER, so the star stays
# unresolved; envelope [651, 750) from the flanking anchors.
case("c06_star_filtered_support", "artificial_1", [
    anchor_row(seq, "artificial_1", 650),
    deletion_row(seq, "artificial_1", 700, 10, "1|0", ps=666, filt="LowQual"),
    snv_row(seq, "artificial_1", 705, "2/1", star=True),
    anchor_row(seq, "artificial_1", 750),
], "filtered supporting call is not a complete carried call; star remains unexplained")

# C14: orphan star without any flanking anchor: the span cannot be bounded.
case("c14_star_orphan_unbounded", "artificial_1", [
    snv_row(seq, "artificial_1", 700, "1/2", star=True),
], "no anchors and no covering deletion: whole-contig exclusion names the missing record")

# C07: overlapping deletion and SNV both homozygous (same copies): conflict gap preserved.
case("c07_overlapping_same_copies_conflict", "artificial_1", [
    deletion_row(seq, "artificial_1", 801, 19, "1/1"),
    snv_row(seq, "artificial_1", 810, "1/1"),
], "both copies carry both overlapping edits; overlap cannot be applied")

# C08a: deletion phased to GT position 0; the interior record is unphased with the
# SNV on GT position 0 and the star on position 1 — resolution must pair the deletion
# with the star across GT orders (parity-flipped constraint).
case("c08a_overlapping_split_copies_starred", "artificial_1", [
    deletion_row(seq, "artificial_1", 901, 19, "1|0", ps=88),
    snv_row(seq, "artificial_1", 910, "1/2", star=True),
], "starred split copies with opposite GT orders; resolution must pair deletion and star")

# C08b: same overlap without the star encoding; the records contradict each other
# (copy 0 is asserted both deleted and reference at 910) and must stay a conflict.
case("c08b_overlapping_split_copies_contradictory", "artificial_1", [
    deletion_row(seq, "artificial_1", 901, 19, "1|0", ps=88),
    snv_row(seq, "artificial_1", 910, "0|1", ps=88),
], "overlapping records without star encoding contradict; conflict gap preserved")

# C09/C10: missing and partially missing calls mask the record footprint.
case("c09_missing_call", "artificial_1", [
    snv_row(seq, "artificial_1", 1000, "./."),
], "missing call masks the record on every declared copy")
case("c10_partial_missing_call", "artificial_1", [
    snv_row(seq, "artificial_1", 1050, "./1"),
], "partially missing call masks the record on every declared copy")

# C11: two SNVs in phase set A and a third in set B; windows spanning both sets unresolved.
case("c11_phase_block_boundary", "artificial_1", [
    snv_row(seq, "artificial_1", 1100, "0|1", ps=31),
    snv_row(seq, "artificial_1", 1108, "1|0", ps=31),
    snv_row(seq, "artificial_1", 1116, "0|1", ps=47),
], "phase-set boundary at 1116; windows touching only one set stay evaluable")

# C12: insertion on copy 0 and nearby SNV on copy 1; exact paths and sequences retained.
case("c12_supported_indel_paths", "artificial_1", [
    row("artificial_1", 150, seq[150], seq[150] + "GG", "1|0", fmt="GT:PS", extra=["52"]),
    snv_row(seq, "artificial_1", 158, "0|1", ps=52),
], "2 bp insertion retains boundary/offset path; SNV on the other copy")

# C13: complex length-changing allele (neither pure indel nor substitution), carried
# heterozygous; unsupported as an edit but fully bounded by its record footprint.
case("c13_complex_allele_bounded", "artificial_1", [
    row("artificial_1", 350, seq[350:354], seq[350] + "TG", "0/1"),
], "complex allele: bounded footprint; not applicable as an edit")

# X/Y configurations. Contig names carry no ploidy meaning; intervals decide.
xseq = REFERENCES["chrX"]
yseq = REFERENCES["chrY"]

case("x01_xy_region_specific_ploidy", "chrX", [
    snv_row(xseq, "chrX", 200, "0/1"),          # PAR-style diploid interval
    snv_row(xseq, "chrX", 600, "1"),            # non-PAR haploid interval
], "diploid and haploid intervals on one contig, declared explicitly")

case("x02_xy_diploid_call_in_haploid_region", "chrX", [
    snv_row(xseq, "chrX", 600, "0/1"),
], "diploid syntax in a declared haploid interval stays unsupported")

case("x03_xy_undeclared_ploidy", "chrY", [
    snv_row(yseq, "chrY", 100, "0/1"),
], "no ploidy intervals: contig name must not imply copies")

case("x04_xy_zero_copies", "chrY", [
    snv_row(yseq, "chrY", 100, "0/1"),
], "zero declared copies: no materialization, explicit exclusion")

case("x05_xy_haploid_star_explained", "chrX", [
    deletion_row(xseq, "chrX", 300, 8, "1"),
    snv_row(xseq, "chrX", 304, "2", star=True),
], "haploid carrier: star explained by the same copy's covering deletion")

case("x06_xy_par_star_explained", "chrX", [
    deletion_row(xseq, "chrX", 210, 8, "1|0", ps=99),
    snv_row(xseq, "chrX", 214, "2|0", ps=99, star=True),
], "star carried on copy 0, explained by the same copy's deletion in a diploid interval")

PLOIDY = {
    "default": [{"contig": "artificial_1", "start": 0, "end": 1200, "copies": 2}],
    "x01_xy_region_specific_ploidy": [{"contig": "chrX", "start": 0, "end": 400, "copies": 2},
                                      {"contig": "chrX", "start": 400, "end": 800, "copies": 1}],
    "x02_xy_diploid_call_in_haploid_region": [{"contig": "chrX", "start": 0, "end": 800, "copies": 1}],
    "x03_xy_undeclared_ploidy": [],
    "x04_xy_zero_copies": [{"contig": "chrY", "start": 0, "end": 400, "copies": 0}],
    "x05_xy_haploid_star_explained": [{"contig": "chrX", "start": 0, "end": 800, "copies": 1}],
    "x06_xy_par_star_explained": [{"contig": "chrX", "start": 0, "end": 400, "copies": 2},
                                  {"contig": "chrX", "start": 400, "end": 800, "copies": 1}],
}


def main():
    (HERE / "references").mkdir(exist_ok=True)
    (HERE / "cases").mkdir(exist_ok=True)
    live = {name + ".vcf" for name in CASES}
    for stale in (HERE / "cases").glob("*.vcf"):
        if stale.name not in live:
            stale.unlink()
    manifest = {"generator": "LCG state=(state*1103515245+12345) mod 2**31; base=ACGT[(state>>16) mod 4]",
                "references": {}, "cases": {}}
    for name, sequence in REFERENCES.items():
        (HERE / "references" / (name + ".fa")).write_text(">" + name + "\n" + sequence + "\n")
        kmers = {sequence[i:i + 23] for i in range(len(sequence) - 22)}
        assert len(kmers) == len(sequence) - 22, name + " has repeated 23-mers"
        manifest["references"][name] = {"length": len(sequence), "unique_23mers_forward": len(kmers)}
    for name, spec in CASES.items():
        text = VCF_HEADER + "".join(r + "\n" for r in spec["rows"])
        (HERE / "cases" / (name + ".vcf")).write_text(text)
        manifest["cases"][name] = {"vcf": name + ".vcf", "note": spec["note"],
                                   "ploidy": PLOIDY.get(name, PLOIDY["default"])}
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    print(f"{len(REFERENCES)} references, {len(CASES)} cases written")


if __name__ == "__main__":
    main()
