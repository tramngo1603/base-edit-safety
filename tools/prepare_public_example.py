#!/usr/bin/env python3
"""Prepare pinned public GRCh38/HG003 inputs and a complete search request.

Choose a work directory on protected storage. Existing mismatched inputs and
files from failed preparation are retained; preparation never starts a search.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import urllib.request

PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
SOURCES = {
    "ucsc_hg38_fa_gz": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz", "data/ucsc/hg38.fa.gz", "c1dd87068c254eb53d944f71e51d1311964fce8de24d6fc0effc9c61c01527d4"),
    "giab_hg003_vcf_gz": ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", "data/giab/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz", "a4f8fefa826f8d2c9457eacf6f4b2f6bfa29daef429d570ce04254c5e7b1121e"),
}
PREPARED = {"primary.fa": "e55dde29c6e54352f68be3973034e63c37b5efd6fbaa878558f6ab6fa8eade0e", "primary.fa.fai": "aeaf6f0a845900acea315871a62982711cebc6da5d982f4085b78ad1ba03a0ef", "final.vcf.gz": "964a8c7ea2b672ef06ba893b5a1a6594b6326697d335e1c6af9228bc1a92bdb9", "final.vcf.gz.tbi": "2ca4fccd341f78714aaaf344b3c89867e8c6057f9ebbc79eca31c730a9702a29"}


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def pin(path):
    return {"path": str(path.resolve()), "sha256": sha(path)}


SOURCE_BYTES = {"ucsc_hg38_fa_gz": 983659424, "giab_hg003_vcf_gz": 147212428}


def acquire(root):
    files = {}
    for name, (url, relative, expected) in SOURCES.items():
        target = root / relative
        print(f"Checking {name}: {SOURCE_BYTES[name]:,} download bytes", flush=True)
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            partial = target.with_name(target.name + ".part")
            with partial.open("xb") as out, urllib.request.urlopen(url, timeout=120) as response:
                downloaded = 0
                while block := response.read(1024 * 1024):
                    out.write(block)
                    downloaded += len(block)
                    if downloaded // (64 * 1024**2) != (downloaded - len(block)) // (64 * 1024**2):
                        print(f"Downloading {name}: {downloaded:,} / {SOURCE_BYTES[name]:,} bytes", flush=True)
            if sha(partial) != expected:
                raise ValueError("Downloaded checksum mismatch; partial file retained")
            partial.rename(target)
        print(f"Verifying {name} checksum", flush=True)
        if sha(target) != expected:
            raise ValueError("Existing raw input checksum mismatch; refusing overwrite")
        files[name] = {"url": url, "path": relative, "sha256": expected, "bytes": target.stat().st_size}
    manifest = root / "source-manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text()).get("files") != files:
            raise ValueError("Existing source manifest differs")
    else:
        save(manifest, {"files": files})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True, help="Researcher-chosen protected storage directory")
    parser.add_argument("--backend", type=Path, required=True, help="Locally built patched Cas-OFFinder 2.4.1")
    parser.add_argument("--guides", type=Path, required=True)
    parser.add_argument("--device", default="G0")
    parser.add_argument("--reuse-prepared", type=Path, help="Directory containing the four checksum-pinned prepared files")
    parser.add_argument("--reference-only", action="store_true")
    parser.add_argument("--confidence-bed", type=Path, help="Optional GIAB benchmark BED; reported separately from computational coverage")
    parser.add_argument("--evidence", type=Path, help="Optional private evaluation TSV, retained outside Git")
    parser.add_argument("--max-resident-bytes", type=int, default=10 * 1024**3)
    parser.add_argument("--max-output-bytes", type=int, default=15 * 1024**3)
    parser.add_argument("--timeout-seconds", type=int, default=21600)
    args = parser.parse_args()
    os.umask(0o077)
    root = args.workdir.expanduser().resolve()
    backend, guides = args.backend.expanduser().resolve(), args.guides.expanduser().resolve()
    if not os.access(backend, os.X_OK) or not re.fullmatch(r"[CG](?:[0-9]+(?:,[0-9]+)*)?", args.device):
        raise ValueError("Backend must be executable and device must be a Cas-OFFinder selector")
    with guides.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames != ["guide_id", "spacer_sequence"] or len(rows) != 23:
            raise ValueError("Require the complete 23-guide TSV panel")
    if len({r["guide_id"] for r in rows}) != 23 or len({r["spacer_sequence"] for r in rows}) != 23 or any(not re.fullmatch(r"[ACGT]{20}", r["spacer_sequence"]) or not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", r["guide_id"]) for r in rows):
        raise ValueError("Guide IDs and uppercase 20-base spacers must be valid and unique")
    if min(args.max_resident_bytes, args.max_output_bytes, args.timeout_seconds) <= 0:
        raise ValueError("Limits must be positive")
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    if (root / "request.json").exists():
        raise ValueError("Request exists; choose a new work directory")
    required_free = args.max_output_bytes + (0 if args.reuse_prepared else 10 * 1024**3)
    from base_edit_safety.lab import preflight
    preflight(root, large=True, required_disk_bytes=required_free)
    if shutil.disk_usage(root).free < required_free:
        raise ValueError("Insufficient free space on workdir filesystem; need " + str(required_free) + " bytes for this preparation and output budget")
    import pysam
    if pysam.__version__ != "0.24.1":
        raise ValueError("Install pysam==0.24.1 for this preparation recipe")
    if args.reuse_prepared:
        prepared = args.reuse_prepared.expanduser().resolve()
        for name, expected in PREPARED.items():
            if sha(prepared / name) != expected:
                raise ValueError("Prepared reuse checksum mismatch: " + name)
        provenance = "Checksum-verified reuse of pinned public GRCh38/HG003 prepared files; bcftools norm --check-ref e, no splitting or atomization; reference header GRCh38."
    else:
        prepared = root / "out"
        if prepared.exists() or (root / "preparation-receipt.json").exists():
            raise ValueError("Preparation output exists; retain it and use a new work directory")
        acquire(root)
        prepared.mkdir()
        environment = dict(os.environ, TMPDIR=str(root))
        for helper in ("make_primary.py", "prepare_vcf.py"):
            print("Preparing " + helper + "; this may take several minutes", flush=True)
            log = root / (helper + ".log")
            with log.open("x") as handle:
                subprocess.run([sys.executable, str(Path(__file__).parent / "public_example" / helper), "--workdir", str(root)], env=environment, stdout=handle, stderr=subprocess.STDOUT, check=True)
        provenance = "Pinned UCSC hg38 and GIAB HG003 NIST v4.2.1; pysam 0.24.1 bundled bcftools norm -f primary.fa --check-ref e; no multiallelic splitting, no atomization, no invented phase; see preparation-receipt.json."
    reference, vcf = prepared / "primary.fa", prepared / "final.vcf.gz"
    if sha(reference) != PREPARED["primary.fa"] or sha(Path(str(reference) + ".fai")) != PREPARED["primary.fa.fai"]:
        raise ValueError("Prepared primary reference differs from the pinned 24-contig reference")
    with pysam.FastaFile(str(reference)) as fa:
        lengths = dict(zip(fa.references, fa.lengths))
    if set(lengths) != set(PRIMARY):
        raise ValueError("Reference must contain all and only 24 primary contigs")
    star_counts = {}
    source_rows = 0
    with pysam.VariantFile(str(vcf)) as vf:
        if list(vf.header.samples) != ["HG003"] or [h.value for h in vf.header.records if h.key == "reference"] != ["GRCh38"]:
            raise ValueError("Final VCF sample or reference header is wrong")
        if not args.reference_only:
            for variant in vf:
                source_rows += 1
                if "*" in (variant.alts or ()):
                    star_counts[variant.contig] = star_counts.get(variant.contig, 0) + 1
            if source_rows != 4000097 or sum(star_counts.values()) != 306 or not set(star_counts) <= set(PRIMARY[:22]):
                raise ValueError("Prepared sample differs from the pinned variant-class accounting")
    excluded_contigs = []
    request = {"format_version": "1", "input_kind": "public_genome", "assembly": "GRCh38", "contigs": PRIMARY,
               "reference": pin(reference), "guides": pin(guides), "backend": {**pin(backend), "version": "2.4.1", "device": args.device},
               "max_mismatches": 6, "limits": {"max_resident_bytes": args.max_resident_bytes, "max_output_bytes": args.max_output_bytes, "timeout_seconds": args.timeout_seconds}}
    if not args.reference_only:
        request["sample"] = {"id": "HG003", "variants": pin(vcf), "normalization": provenance,
                             "genotype_policy": {"phase": "explicit_ps", "min_gq": None, "min_dp": None, "accepted_site_filters": ["PASS", "."], "accepted_sample_filters": ["PASS", "."]},
                             "ploidy_intervals": [{"contig": c, "start": 0, "end": lengths[c], "copies": 2} for c in PRIMARY[:22]],
                             "excluded_contigs": excluded_contigs,
                             "ploidy_source": "GIAB HG003 NISTv4.2.1 GRCh38 autosomal benchmark VCF; explicit diploid chr1-22 calls. No X/Y callset or region-specific copy declaration supplied."}
        save(root / "sample-scope.json", {
            "source_vcf": request["sample"]["variants"], "source_rows": source_rows,
            "spanning_deletion_rows": sum(star_counts.values()), "rows_by_contig": star_counts,
            "excluded_contigs": excluded_contigs,
            "retained_autosomes": PRIMARY[:22],
            "basis": "Evaluate all 22 autosomes. Resolve supported spanning deletions and local overlaps from the original calls. Bound unsupported deletion regions only by explicit event spans or flanking called bases; retain whole-contig exclusions if no bounds exist. X/Y calls and ploidy are not supplied. Exact coverage is determined during search, not assumed here."})
    if args.confidence_bed and not args.reference_only:
        request["sample"]["confidence_bed"] = pin(args.confidence_bed.expanduser().resolve())
        request["sample"]["confidence_source"] = "GIAB HG003 NISTv4.2.1 GRCh38 benchmark regions; small-variant benchmark membership, not per-base or phased genotype certification."
    if args.evidence:
        request["evidence"] = pin(args.evidence.expanduser().resolve())
    local_guides = root / "guides.tsv"
    shutil.copyfile(guides, local_guides)
    request["guides"] = pin(local_guides)
    for spec in (request["reference"], request["guides"], request.get("sample", {}).get("variants")):
        if spec and Path(spec["path"]).is_relative_to(root):
            spec["path"] = str(Path(spec["path"]).relative_to(root))
    save(root / "request.json", request)
    print(json.dumps({"status": "prepared", "request": str(root / "request.json"), "guides": 23, "primary_contigs": 24, "search_started": False}))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("Preparation failed; files retained: " + str(exc), file=sys.stderr)
        raise SystemExit(1)
