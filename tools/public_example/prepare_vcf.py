#!/usr/bin/env python3
"""Normalize the raw GIAB HG003 v4.2.1 VCF against out/primary.fa and emit a
final VCF whose single ##reference header is exactly '##reference=GRCh38'.

Guarantees:
- Raw input is never modified; its SHA256 is re-verified against the pin in
  source-manifest.json before use.
- bcftools norm is run WITHOUT -m (no multiallelic splitting) and WITHOUT
  --atomize; genotype/phase/quality/filter fields are untouched.
- --check-ref e: bcftools exits loudly on any REF/reference mismatch; we never
  fix alleles to conceal a mismatch and never exclude records silently.
- Sortedness is verified; bcftools sort is applied only if needed.
- Only aggregate statistics are printed; no calls or sequences.
"""
import gzip
import os
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pysam
from pysam import bcftools

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--workdir", type=Path, required=True)
ROOT = parser.parse_args().workdir.resolve()
RAW = ROOT / "data" / "giab" / "HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
REF = ROOT / "out" / "primary.fa"
NORM = ROOT / "out" / "normalized.vcf.gz"
FINAL = ROOT / "out" / "final.vcf.gz"
MANIFEST = ROOT / "source-manifest.json"
RECEIPT = ROOT / "preparation-receipt.json"

NORM_CMD = (
    f"bcftools norm -f {REF.relative_to(ROOT)} --check-ref e -O z "
    f"-o {NORM.relative_to(ROOT)} {RAW.relative_to(ROOT)}"
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path, default):
    return json.loads(path.read_text()) if path.exists() else default


def save_receipt(r: dict) -> None:
    r["updated_utc"] = datetime.now(timezone.utc).isoformat()
    RECEIPT.write_text(json.dumps(r, indent=2, sort_keys=True) + "\n")


def fail(msg: str) -> "None":
    print(f"FAIL {msg}")
    sys.exit(1)


def fai_contigs() -> dict:
    out = {}
    for line in open(str(REF) + ".fai"):
        name, length = line.split("\t")[:2]
        out[name] = int(length)
    return out


def main() -> int:
    manifest = load_json(MANIFEST, {})
    pin = manifest.get("files", {}).get("giab_hg003_vcf_gz", {}).get("sha256")
    raw_hash = sha256_file(RAW)
    if pin is None or pin != raw_hash:
        fail(f"raw VCF sha256 mismatch vs manifest pin (on-disk={raw_hash})")
    receipt = load_json(RECEIPT, {})
    receipt.setdefault("tools", {
        "python": sys.version.split()[0],
        "pysam": pysam.__version__,
        "pysam_bundled_samtools": getattr(pysam, "__samtools_version__", "unknown"),
        "pysam_license": "MIT (wheel: pysam==0.24.1)",
    })
    receipt.setdefault("vcf_pipeline", {})["raw_vcf"] = {
        "path": str(RAW.relative_to(ROOT)), "sha256": raw_hash,
        "bytes": RAW.stat().st_size,
    }

    # --- header inspection (metadata only) ---
    vf = pysam.VariantFile(str(RAW))
    hdr = vf.header
    samples = list(hdr.samples)
    if len(samples) != 1:
        fail(f"expected exactly 1 sample, found {len(samples)}")
    sample = samples[0]
    formats = sorted(hdr.formats.keys())
    filters = sorted(hdr.filters.keys())
    contig_lines = {c: int(hdr.contigs[c].length or 0) for c in hdr.contigs}
    reference_lines = [l.value for l in hdr.records if l.key == "reference"]
    raw_rows = 0
    contigs_used = Counter()
    for rec in vf:
        raw_rows += 1
        contigs_used[rec.chrom] += 1
    vf.close()
    receipt["vcf_pipeline"]["raw_header"] = {
        "sample": sample,
        "format_fields": formats,
        "filter_values": filters,
        "reference_header_values": reference_lines,
        "contig_count": len(contig_lines),
    }
    receipt["vcf_pipeline"]["raw_row_count"] = raw_rows

    # --- chromosome alias check against primary.fa (contigs used by records) ---
    ref_contigs = fai_contigs()
    bad = {c for c in contigs_used if c not in ref_contigs}
    if bad:
        fail(f"record contigs absent from primary.fa: {sorted(bad)}")
    len_mismatch = {c: (contig_lines.get(c), ref_contigs[c])
                    for c in contigs_used
                    if contig_lines.get(c) and contig_lines[c] != ref_contigs[c]}
    if len_mismatch:
        fail(f"contig length mismatches vs primary.fa: {len_mismatch}")
    header_only = sorted(set(contig_lines) - set(contigs_used))
    receipt["vcf_pipeline"]["raw_header"]["contigs_with_records"] = dict(contigs_used)
    receipt["vcf_pipeline"]["raw_header"]["header_only_contigs_count"] = len(header_only)

    # --- normalization (loud REF checking, no splitting) ---
    try:
        bcftools.norm("-f", str(REF), "--check-ref", "e", "-O", "z",
                      "-o", str(NORM), str(RAW), catch_stdout=False)
    except pysam.utils.SamtoolsError as exc:
        fail(f"bcftools norm failed (likely REF mismatch or malformed record): {exc}")
    norm_hash = sha256_file(NORM)

    # Keep the normalization command and version recorded in the header.
    nvh = pysam.VariantFile(str(NORM)).header
    norm_version = [l.value for l in nvh.records if l.key == "bcftools_normVersion"]
    norm_command = [l.value for l in nvh.records if l.key == "bcftools_normCommand"]

    # --- sortedness check ---
    order = {name: i for i, name in enumerate(ref_contigs)}
    prev = (-1, -1)
    sorted_ok = True
    vf = pysam.VariantFile(str(NORM))
    for rec in vf:
        key = (order.get(rec.chrom, 10**9), rec.pos)
        if key < prev:
            sorted_ok = False
            break
        prev = key
    vf.close()
    sort_applied = False
    src_for_final = NORM
    if not sorted_ok:
        sorted_tmp = ROOT / "out" / "normalized.sorted.vcf.gz"
        bcftools.sort("-m", "2G", "-T", str(ROOT / "out"), "-O", "z",
                      "-o", str(sorted_tmp), str(NORM), catch_stdout=False)
        sorted_tmp.replace(NORM)
        norm_hash = sha256_file(NORM)
        sort_applied = True

    # --- final header rewrite: exactly one ##reference=GRCh38 ---
    ref_line_count = 0
    inserted = False
    with gzip.open(NORM, "rt", encoding="utf-8") as fin, open(ROOT / "out" / "final.vcf", "wt", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("##reference="):
                ref_line_count += 1
                if not inserted:
                    fout.write("##reference=GRCh38\n")
                    inserted = True
                continue  # drop any duplicate ##reference lines only after count check
            fout.write(line)
            if not inserted and line.startswith("##fileformat="):
                fout.write("##reference=GRCh38\n")
                inserted = True
    if ref_line_count > 1:
        fail(f"expected at most 1 ##reference line in normalized header, found {ref_line_count}")
    if not inserted:
        fail("could not place ##reference=GRCh38 (no ##fileformat line?)")
    pysam.tabix_compress(str(ROOT / "out" / "final.vcf"), str(FINAL), force=False)
    (ROOT / "out" / "final.vcf").unlink()
    pysam.tabix_index(str(FINAL), preset="vcf", force=False)
    final_hash = sha256_file(FINAL)

    # --- final validation: header, rows, aggregate stats ---
    vf = pysam.VariantFile(str(FINAL))
    fhdr = vf.header
    final_refs = [l.value for l in fhdr.records if l.key == "reference"]
    if final_refs != ["GRCh38"]:
        fail(f"final ##reference header wrong: {final_refs}")
    final_rows = 0
    class_counts = Counter()
    filter_counts = Counter()
    phased = 0
    multiallelic = 0
    for rec in vf:
        final_rows += 1
        n_filters = 0
        for f in rec.filter:
            filter_counts[f] += 1
            n_filters += 1
        if n_filters == 0:
            filter_counts["MISSING(.)"] += 1
        alts = rec.alts or ()
        if len(alts) > 1:
            multiallelic += 1
        allele_lens = {len(rec.ref)} | {len(a) for a in alts}
        if allele_lens == {1}:
            class_counts["SNP"] += 1
        elif len(allele_lens) == 1:
            class_counts["MNP_or_complex_equal_length"] += 1
        else:
            class_counts["INDEL_or_mixed"] += 1
        if rec.samples[sample].phased:
            phased += 1
    vf.close()
    if final_rows != raw_rows:
        fail(f"row count changed: raw={raw_rows} final={final_rows}")

    # final REF check against primary.fa (no realignment, loud)
    try:
        bcftools.norm("-N", "-f", str(REF), "--check-ref", "e",
                      "-o", "/dev/null", str(FINAL), catch_stdout=False)
    except pysam.utils.SamtoolsError as exc:
        fail(f"final REF verification failed: {exc}")

    receipt["vcf_pipeline"].update({
        "normalization": {
            "command": NORM_CMD,
            "options": "no -m (multiallelics preserved), no --atomize, --check-ref e (exit on REF mismatch)",
            "bcftools_norm_version_header": norm_version,
            "bcftools_norm_command_header": norm_command,
            "sortedness_check": "passed" if not sort_applied else "failed; bcftools sort -m 2G applied",
            "normalized_vcf": {"path": str(NORM.relative_to(ROOT)), "sha256": norm_hash,
                               "bytes": NORM.stat().st_size},
        },
        "final_vcf": {
            "path": str(FINAL.relative_to(ROOT)), "sha256": final_hash,
            "bytes": FINAL.stat().st_size,
            "tbi_sha256": sha256_file(Path(str(FINAL) + ".tbi")),
            "reference_header": final_refs,
            "original_reference_header_values_preserved": reference_lines,
        },
        "final_stats": {
            "total_rows": final_rows,
            "variant_class_counts": dict(class_counts),
            "multiallelic_records": multiallelic,
            "filter_counts": dict(filter_counts),
            "phased_records": phased,
            "sample": sample,
            "format_fields": formats,
        },
    })
    save_receipt(receipt)
    print(f"OK raw_rows={raw_rows} final_rows={final_rows} sample={sample} "
          f"class_counts={dict(class_counts)} multiallelic={multiallelic} "
          f"phased={phased} sort_applied={sort_applied}")
    print(f"OK normalized sha256={norm_hash}")
    print(f"OK final sha256={final_hash}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
