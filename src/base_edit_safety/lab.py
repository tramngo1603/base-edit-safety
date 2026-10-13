"""Local setup, input preparation and example checks for the lab launcher."""
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import sysconfig

from .population.common import checksum, save_json
from .population import research

GIB = 1024**3
PUBLIC_DOWNLOAD_BYTES = 1130871852
GRCH38_LENGTHS = {'chr1': 248956422, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr2': 242193529, 'chr20': 64444167, 'chr21': 46709983, 'chr22': 50818468, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chrX': 156040895, 'chrY': 57227415}


def files():
    return Path(sysconfig.get_path("data")) / "share/base-edit-safety"


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": checksum(path)}


def preflight(folder, *, large=False, required_disk_bytes=None):
    import psutil
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        raise ValueError("This handoff requires an Apple silicon Mac (arm64). See QUICKSTART.html.")
    if tuple(map(int, platform.mac_ver()[0].split(".")[:2])) < (14, 4):
        raise ValueError("This handoff requires macOS 14.4 or later; it was tested on 14.4.1.")
    if sys.version_info[:3] != (3, 12, 12):
        raise ValueError("Use Start.command so the pinned Python 3.12.12 runtime is selected.")
    import pysam
    if pysam.__version__ != "0.24.1":
        raise ValueError("Setup is incomplete. Run Start.command again to install pysam 0.24.1.")
    for command in ("/usr/bin/time", "/bin/ps"):
        if not os.access(command, os.X_OK):
            raise ValueError(f"Required macOS tool is missing: {command}")
    if not Path("/System/Library/Frameworks/OpenCL.framework").exists():
        raise ValueError("The macOS OpenCL framework is missing. Restore the OS runtime before running.")
    memory = psutil.virtual_memory()
    required_memory = 4 * GIB if large else 512 * 1024**2
    if large and memory.total < 16 * GIB:
        raise ValueError("Full-reference work requires at least 16 GiB RAM; this Mac has less.")
    if memory.available < required_memory:
        raise ValueError(f"Close other applications: need {required_memory / GIB:g} GiB available RAM; "
                         f"{memory.available / GIB:.1f} GiB is available.")
    folder = Path(folder)
    if not folder.is_dir() or not os.access(folder, os.W_OK):
        raise ValueError("Choose an existing writable output folder on local storage.")
    free = shutil.disk_usage(folder).free
    needed = required_disk_bytes if required_disk_bytes is not None else 25 * GIB if large else GIB
    if free < needed:
        raise ValueError(f"Free space in the selected folder is {free/GIB:.1f} GiB; need {needed/GIB:g} GiB.")
    return {"os": platform.mac_ver()[0], "architecture": platform.machine(),
            "total_memory_bytes": memory.total, "available_memory_bytes": memory.available,
            "free_disk_bytes": free}


def read_guides(text):
    text = text.strip()
    if text.startswith("guide_id\t"):
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        if reader.fieldnames != ["guide_id", "spacer_sequence"]:
            raise ValueError("Guide TSV columns must be guide_id and spacer_sequence, in that order.")
        rows = list(reader)
    else:
        rows = [{"guide_id": f"guide_{i+1}", "spacer_sequence": line.strip().upper()}
                for i, line in enumerate(text.splitlines()) if line.strip()]
    if not 1 <= len(rows) <= 32:
        raise ValueError("Enter 1–32 guides, one 20-base spacer per line, without the PAM.")
    for row in rows:
        if set(row) != {"guide_id", "spacer_sequence"} or not re.fullmatch(r"[ACGT]{20}", row["spacer_sequence"]):
            raise ValueError("Every spacer must have exactly 20 A/C/G/T bases, without a PAM or spaces.")
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", row["guide_id"]):
            raise ValueError("Guide IDs may contain letters, numbers, dots, underscores and hyphens.")
    if len({r["guide_id"] for r in rows}) != len(rows) or len({r["spacer_sequence"] for r in rows}) != len(rows):
        raise ValueError("Remove duplicate guide IDs or spacer sequences before starting.")
    return rows


def samples(path):
    import pysam
    with pysam.VariantFile(str(path)) as handle:
        return list(handle.header.samples)


def base_request(reference, guides, backend, contigs, *, small=False):
    return {"format_version": "1", "input_kind": "artificial" if small else "research_sample",
            "assembly": "artificial" if small else "GRCh38", "reference": pin(reference),
            "guides": pin(guides), "backend": {**pin(backend), "version": "2.4.1", "device": "G0"},
            "contigs": contigs, "max_mismatches": 2 if small else 6,
            "limits": {"max_resident_bytes": 10 * GIB, "max_output_bytes": 15 * GIB,
                       "timeout_seconds": 21600}}


def sample_config(vcf, sample, lengths, normalization, *, excluded=(), small=False):
    return {"id": sample, "variants": pin(vcf), "normalization": normalization,
            "genotype_policy": {"phase": "explicit_ps", "min_gq": None, "min_dp": None,
                                "accepted_site_filters": ["PASS", "."],
                                "accepted_sample_filters": ["PASS", "."]},
            "ploidy_intervals": [{"contig": c, "start": 0, "end": n, "copies": 2}
                                 for c, n in lengths.items() if small or c in research.PRIMARY[:22]],
            "excluded_contigs": [{"contig": c, "reason": "unresolved_spanning_deletion_extent"}
                                  for c in sorted(excluded)]}


def prepare_public(cache, backend, progress=print):
    cache = Path(cache).resolve()
    cache.parent.mkdir(parents=True, exist_ok=True)
    preflight(cache.parent, large=True)
    if (cache / "request.json").exists():
        progress("Verifying cached public inputs")
        config = json.loads((cache / "request.json").read_text())
        for spec in (config["reference"], config["guides"], config["sample"]["variants"]):
            if checksum(cache / spec["path"]) != spec["sha256"]:
                raise ValueError("Cached public input changed. Keep this cache for inspection and choose a new folder.")
        receipt_path = cache / "preparation-receipt.json"
        receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
        vcf_index_sha = receipt.get("vcf_pipeline", {}).get("final_vcf", {}).get(
            "tbi_sha256", "2ca4fccd341f78714aaaf344b3c89867e8c6057f9ebbc79eca31c730a9702a29")
        indexes = [(str(config["reference"]["path"]) + ".fai",
                    "aeaf6f0a845900acea315871a62982711cebc6da5d982f4085b78ad1ba03a0ef"),
                   (str(config["sample"]["variants"]["path"]) + ".tbi", vcf_index_sha)]
        for name, expected in indexes:
            if not (cache / name).is_file() or checksum(cache / name) != expected:
                raise ValueError("Cached public index changed or is missing: " + Path(name).name +
                                 ". Keep this cache for inspection and choose a new folder.")
        # Rebuild the small request when the bundle has moved; the genomic files stay local.
        # The pinned old cache used a blanket star exclusion; current reconstruction bounds these events.
        config["sample"].setdefault("ploidy_source", "GIAB HG003 NISTv4.2.1 GRCh38 autosomal benchmark VCF; explicit diploid chr1-22 calls. No X/Y callset or region-specific copy declaration supplied.")
        config["sample"]["excluded_contigs"] = [r for r in config["sample"].get("excluded_contigs", [])
            if r["reason"] != "unresolved_spanning_deletion_extent"]
        config["backend"] = {**pin(backend), "version": "2.4.1", "device": "G0"}
        for spec in (config["reference"], config["guides"], config["sample"]["variants"]):
            spec["path"] = str((cache / spec["path"]).resolve())
        return config
    progress("Download and prepare: 1.13 GB download; 25 GiB free disk required. Verified files are reused.")
    cmd = [sys.executable, "-I", str(files() / "tools/prepare_public_example.py"),
           "--workdir", str(cache), "--backend", str(backend),
           "--guides", str(files() / "examples/research/guides.tsv")]
    cache.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (cache / "setup.log").open("a") as log:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                env={**os.environ, "PYTHONUNBUFFERED": "1"})
        for line in proc.stdout:
            log.write(line)
            log.flush()
            progress(line.strip())
        if proc.wait():
            raise ValueError("Public preparation failed. Read setup.log in the cache folder. "
                             "Check the connection and disk space; keep partial files for inspection.")
    return json.loads((cache / "request.json").read_text())


def prepare_own(inputs, backend, *, guide_text, reference, vcf=None, sample=None,
                diploid=False, progress=print):
    import pysam
    from pysam import bcftools
    rows = read_guides(guide_text)
    reference = Path(reference).expanduser().resolve()
    if not reference.is_file() or reference.suffix == ".gz":
        raise ValueError("Choose an uncompressed GRCh38 FASTA file (.fa or .fasta).")
    if not Path(str(reference) + ".fai").is_file():
        progress("Indexing the selected reference FASTA")
        try:
            pysam.faidx(str(reference))
        except Exception as exc:
            raise ValueError("Cannot index the reference. Put the FASTA in a writable folder.") from exc
    with pysam.FastaFile(str(reference)) as fa:
        lengths = {c: n for c, n in zip(fa.references, fa.lengths) if c in research.PRIMARY}
    if lengths != GRCH38_LENGTHS:
        raise ValueError("Reference chromosome names and lengths must match full GRCh38 chr1–chr22, chrX and chrY. Use the public reference option.")
    guide_path = inputs / "guides.tsv"
    research.write_tsv(guide_path, rows, ["guide_id", "spacer_sequence"])
    request = base_request(reference, guide_path, backend, research.PRIMARY)
    if not vcf:
        return request
    if not diploid:
        raise ValueError("Confirm two copies of each autosome for this sample. Otherwise use the advanced ploidy workflow.")
    if not sample or not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", sample):
        raise ValueError("Select a VCF sample with an ID containing only letters, numbers, dots, underscores or hyphens.")
    vcf = Path(vcf).resolve()
    with pysam.VariantFile(str(vcf)) as handle:
        if sample not in handle.header.samples:
            raise ValueError("The selected sample is not in this VCF. Select it again.")
        if [h.value for h in handle.header.records if h.key == "reference"] != ["GRCh38"]:
            raise ValueError("The VCF must declare exactly ##reference=GRCh38. Ask its provider to confirm the assembly.")
    progress("Checking REF alleles and normalizing the local VCF; original calls are preserved")
    normalized, prepared = inputs / "normalized.vcf.gz", inputs / "sample.vcf.gz"
    try:
        bcftools.norm("--no-version", "-f", str(reference), "--check-ref", "e", "-O", "z",
                      "-o", str(normalized), str(vcf), catch_stdout=False)
        bcftools.sort("-m", "1G", "-T", str(inputs), "-O", "z", "-o", str(prepared),
                      str(normalized), catch_stdout=False)
        pysam.tabix_index(str(prepared), preset="vcf")
    except Exception as exc:
        raise ValueError("VCF preparation failed. Check GRCh38, chromosome names and REF alleles with the VCF provider.") from exc
    excluded = set()
    with pysam.VariantFile(str(prepared)) as handle:
        for row in handle:
            if row.contig not in lengths:
                raise ValueError("The VCF has records outside chr1–22/X/Y. Supply a primary-chromosome VCF.")
    request["sample"] = sample_config(prepared, sample, lengths,
        "pysam 0.24.1 bundled bcftools norm --check-ref e, no splitting or atomization; sorted; original phase retained",
        excluded=excluded)
    return request


def scientific_summary(output):
    """Compare every scientific row regardless of display order, retaining duplicates."""
    output = Path(output)
    manifest = json.loads((output / "run-manifest.json").read_text())
    keys = ("status", "genotype_status", "reference_loci", "sample_loci", "candidate_configurations",
            "changed_sequence_loci", "locus_changes")
    summary = {k: manifest[k] for k in keys}
    digests = {}
    for name in ("candidates.tsv", "locus-changes.tsv", "exclusions.tsv"):
        rows = []
        with (output / name).open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row.pop("evidence_links", None)
                for field in ("sample_configurations", "reference_configurations", "sample_copies"):
                    if field in row:
                        row[field] = sorted(json.loads(row[field]))
                # Coordinate paths stay ordered; only membership lists above are unordered.
                rows.append(hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":")).encode()).digest())
        digests[name] = hashlib.sha256(b"".join(sorted(rows))).hexdigest()
    summary["scientific_rows_sha256"] = digests
    coverage = json.loads((output / "coverage.json").read_text())
    # Runtime and command paths vary between runs, so compare coverage separately.
    coverage.get("reference_native_accounting", {}).pop("execution_receipt_sha256", None)
    summary["coverage_sha256"] = hashlib.sha256(
        json.dumps(coverage, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return summary


def verify_example(output, expected):
    actual = scientific_summary(output)
    expected = json.loads(Path(expected).read_text())
    differences = [key for key, value in expected.items() if actual.get(key) != value]
    if differences:
        raise ValueError("Example results differ from the accepted baseline: " + ", ".join(differences))
    return {"status": "matched", "fields_checked": list(expected)}


def record_verification(output, expected):
    """Keep a completed search readable even when its example comparison fails."""
    from .lab_report import write_report
    output = Path(output)
    state = json.loads((output / "run-manifest.json").read_text())
    try:
        state["example_verification"] = verify_example(output, expected)
    except ValueError as exc:
        state["example_verification"] = {"status": "mismatched", "message": str(exc)}
    state["example_verification"].update(verified_utc=datetime.now(timezone.utc).isoformat(),
                                         verifier_source=research.code_identity(),
                                         expected_sha256=checksum(expected))
    write_report(output, state, json.loads((output / "coverage.json").read_text()),
                 json.loads((output / "request.json").read_text()))
    report = output / "report.html"
    state["outputs"]["report.html"] = {"sha256": checksum(report), "bytes": report.stat().st_size}
    save_json(output / "run-manifest.json", state)
    return state


def run_job(mode, root, backend, *, cache=None, progress=print, **own):
    root, backend = Path(root).resolve(), Path(backend).resolve()
    preflight(root.parent, large=mode != "small")
    if root.exists():
        raise ValueError("This run folder already exists. Choose a new folder name; previous results are preserved.")
    root.mkdir(mode=0o700)
    inputs = root / "inputs"
    inputs.mkdir(mode=0o700)
    if mode == "small":
        source = files() / "examples/small"
        for name in ("reference.fa", "guides.tsv", "sample.vcf"):
            shutil.copyfile(source / name, inputs / name)
        request = base_request(inputs / "reference.fa", inputs / "guides.tsv", backend, ["artificial_1"], small=True)
        request["sample"] = sample_config(inputs / "sample.vcf", "SMALL", {"artificial_1": 166},
                                           "Bundled artificial SNVs", small=True)
    elif mode == "public":
        request = prepare_public(cache, backend, progress)
        # Resolve relative paths from the prepared cache before writing the run request.
        for spec in (request["reference"], request["guides"], request["sample"]["variants"]):
            spec["path"] = str((Path(cache) / spec["path"]).resolve())
    elif mode == "own":
        request = prepare_own(inputs, backend, progress=progress, **own)
    else:
        raise ValueError("Choose an example or My analysis.")
    config_path = inputs / "request.json"
    save_json(config_path, request)
    result = research.run(config_path, root / "results", progress=progress)
    if mode in {"small", "public"}:
        progress("Comparing scientific results with the accepted example")
        expected = files() / ("examples/small/expected.json" if mode == "small" else "handoff/expected-hg003.json")
        result = record_verification(root / "results", expected)
    return root / "results", result
