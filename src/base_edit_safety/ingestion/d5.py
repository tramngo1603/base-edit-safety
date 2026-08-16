"""Import D5 HBG amplicon reads and convert CRISPResso2 outcomes to common fields.

The module only processes measured reads.  It imports no prediction or scoring code and exposes no
comparison operation.  Paired 150-bp reads do not overlap across the 431/435-bp nested amplicons,
so they are first aligned as pairs with Bowtie2 and then supplied to CRISPResso2 through its BAM
interface.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any
from zipfile import ZipFile

from .common import sha256_file, write_tsv


DATASET_ID = "D5"
REFERENCE_BUILD = "GRCh38/hg38"
GUIDE_SEQUENCE = "CTTGACCAATAGCCTTGACA"
CRISPRESSO_IMAGE = (
    "pinellolab/crispresso2:v2.3.4@"
    "sha256:b7e7971a87769aa584d31d21e199a9abca5f641c4f889c2e807cb6324fc8ce0c"
)
CRISPRESSO_VERSION = "2.3.4"
ANALYSIS_WINDOW_LENGTH = 150
NORMALIZATION_SCHEMA_VERSION = "1.4"
MAPPED_END_GATE_SCHEMA_VERSION = "1.0"
MAPPED_END_REFERENCE_GEOMETRY = "staged_F1_R1_flank_inclusive"
QC_REPORT_DATE = "2026-08-16"

RUN_QC_NOTES = {
    "SRR25936811": "multicopy HBG1 F1-primer burden (technical)",
}
HELD_OUT_SENSITIVITY_RUNS = frozenset(RUN_QC_NOTES)

# These technical input checks are not prediction-performance thresholds. They do not use a
# predictor or biological outcome. Minimum read depth is checked at every guide-matched position,
# without reference to an expected edit rate.
MIN_MAPPED_END_FRACTION = 0.70
MIN_PROTOSPACER_DEPTH = 1_000


@dataclass(frozen=True)
class LocusSpec:
    locus: str
    cache_filename: str
    sequence_sha256: str
    cache_url: str
    outer_forward_primer: str
    outer_reverse_primer: str
    inner_forward_primer: str
    inner_reverse_primer: str


LOCI = {
    "HBG1": LocusSpec(
        locus="HBG1",
        cache_filename="hg38_chr11_5247857_5251857.json",
        sequence_sha256="02cf185405ae05c0b028ecccef148e5dbe1ba64b428123130ae1d6e054610544",
        cache_url=(
            "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr11;"
            "start=5247856;end=5251857"
        ),
        outer_forward_primer="TATTACTGCGCTGAAACTGTGGCT",
        outer_reverse_primer="ATTGATAACCTCAGACGTTCCAGA",
        inner_forward_primer="CCTGGACCTATGCCTAAAACA",
        inner_reverse_primer="ATTCTTCATCCCTAGCCAGCCGC",
    ),
    "HBG2": LocusSpec(
        locus="HBG2",
        cache_filename="hg38_chr11_5252781_5256781.json",
        sequence_sha256="bffd94832ea846ac333654f1c75b72afa87cd106f55aa4128f993556dc82dbc5",
        cache_url=(
            "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr11;"
            "start=5252780;end=5256781"
        ),
        outer_forward_primer="GCACTGAAACTGTTGCTTTATAGG",
        outer_reverse_primer="ATTGATAACCTCAGACGTTCCAGA",
        inner_forward_primer="CCTGGACCTATGCCTAAAACA",
        inner_reverse_primer="ATTCTTCATCCCTAGCCAGCCGC",
    ),
}


@dataclass(frozen=True)
class AmpliconReference:
    locus: str
    contig: str
    build: str
    chrom: str
    start_0based: int
    end_0based_exclusive: int
    sequence_plus: str
    guide_start_0based: int
    guide_end_0based_exclusive: int
    guide_strand: str
    source_cache: str
    source_url: str
    source_sha256: str
    source_sequence_sha256: str


@dataclass(frozen=True)
class RunSpec:
    run_accession: str
    study_accession: str
    sample_accession: str
    experiment_accession: str
    editor: str
    locus: str
    replicate: int
    cell_type: str
    assay: str
    fastq_r1: str
    fastq_r2: str
    md5_r1: str
    md5_r2: str


@dataclass(frozen=True)
class AnalysisWindow:
    name: str
    start_0based: int
    end_0based_exclusive: int
    sequence_plus: str
    contains_protospacer: bool


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _unique_index(sequence: str, query: str, label: str) -> int:
    first = sequence.find(query)
    if first < 0:
        raise ValueError(f"{label} was not found in the reference cache")
    if sequence.find(query, first + 1) >= 0:
        raise ValueError(f"{label} is not unique in the reference cache")
    return first


def _derive_primer_bounded_reference(
    spec: LocusSpec,
    cache_dir: str | Path,
    *,
    flank_inclusive: bool,
) -> AmpliconReference:
    """Build a primer-bounded reference after checking every inner and outer primer position."""
    cache_path = Path(cache_dir) / spec.cache_filename
    observed_cache_sha256 = sha256_file(cache_path)
    source = json.loads(cache_path.read_text(encoding="utf-8"))
    filename_match = re.fullmatch(r"([^_]+)_(chr[^_]+)_(\d+)_(\d+)\.json", spec.cache_filename)
    if filename_match is None:
        raise ValueError(f"invalid reference-cache filename: {spec.cache_filename}")
    expected_build, expected_chrom, expected_start, expected_end = filename_match.groups()
    if (
        source.get("genome_build") != expected_build
        or source.get("chrom") != expected_chrom
        or int(source.get("start_1based", -1)) != int(expected_start)
        or int(source.get("end_1based", -1)) != int(expected_end)
    ):
        raise ValueError(f"reference-cache metadata does not match filename: {cache_path}")
    if source["genome_build"] not in {"hg38", "GRCh38", "GRCh38/hg38"}:
        raise ValueError(f"unexpected reference build in {cache_path}: {source['genome_build']}")
    sequence = source["plus_sequence"].upper()
    observed_sequence_sha256 = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if observed_sequence_sha256 != spec.sequence_sha256:
        raise ValueError(
            f"reference-sequence checksum mismatch for {cache_path}: "
            f"{observed_sequence_sha256} != {spec.sequence_sha256}"
        )
    cache_start = int(source["start_1based"]) - 1

    outer_left = _unique_index(sequence, spec.outer_reverse_primer, "outer reverse primer")
    outer_right = _unique_index(
        sequence, reverse_complement(spec.outer_forward_primer), "outer forward primer"
    )
    inner_left = _unique_index(sequence, spec.inner_reverse_primer, "inner reverse primer")
    inner_right = _unique_index(
        sequence, reverse_complement(spec.inner_forward_primer), "inner forward primer"
    )
    if not (outer_left < inner_left < inner_right < outer_right):
        raise ValueError(f"primer nesting/order failed for {spec.locus}")

    if flank_inclusive:
        reference_left = outer_left
        reference_end = outer_right + len(spec.outer_forward_primer)
        contig = f"{spec.locus}_mapped_end_gate"
    else:
        reference_left = inner_left
        reference_end = inner_right + len(spec.inner_forward_primer)
        contig = f"{spec.locus}_amplicon"
    amplicon = sequence[reference_left:reference_end]
    guide_on_plus = reverse_complement(GUIDE_SEQUENCE)
    guide_start = _unique_index(amplicon, guide_on_plus, "reverse-complemented guide")
    return AmpliconReference(
        locus=spec.locus,
        contig=contig,
        build=REFERENCE_BUILD,
        chrom=source["chrom"],
        start_0based=cache_start + reference_left,
        end_0based_exclusive=cache_start + reference_end,
        sequence_plus=amplicon,
        guide_start_0based=guide_start,
        guide_end_0based_exclusive=guide_start + len(GUIDE_SEQUENCE),
        guide_strand="-",
        source_cache=str(cache_path),
        source_url=spec.cache_url,
        source_sha256=observed_cache_sha256,
        source_sequence_sha256=observed_sequence_sha256,
    )


def derive_amplicon(spec: LocusSpec, cache_dir: str | Path) -> AmpliconReference:
    """Derive the F2/R2 inner amplicon used to put outcome counts on one scale."""
    return _derive_primer_bounded_reference(spec, cache_dir, flank_inclusive=False)


def derive_mapped_end_reference(
    spec: LocusSpec, cache_dir: str | Path
) -> AmpliconReference:
    """Build the outer-primer F1/R1 reference used for the read-mapping check."""
    return _derive_primer_bounded_reference(spec, cache_dir, flank_inclusive=True)


def read_run_manifest(path: str | Path) -> list[RunSpec]:
    runs: list[RunSpec] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["dataset_id"] != DATASET_ID:
                continue
            runs.append(
                RunSpec(
                    run_accession=row["run_accession"],
                    study_accession=row["study_accession"],
                    sample_accession=row["sample_accession"],
                    experiment_accession=row["experiment_accession"],
                    editor=row["editor"],
                    locus=row["locus"],
                    replicate=int(row["replicate"]),
                    cell_type=row["cell_type"],
                    assay=row["assay"],
                    fastq_r1=Path(row["fastq_r1"]).name,
                    fastq_r2=Path(row["fastq_r2"]).name,
                    md5_r1=row["md5_r1"],
                    md5_r2=row["md5_r2"],
                )
            )
    accessions = [run.run_accession for run in runs]
    if len(accessions) != len(set(accessions)):
        raise ValueError("duplicate D5 run accession in manifest")
    for run in runs:
        if run.locus not in LOCI:
            raise ValueError(f"unsupported locus in D5 manifest: {run.locus}")
    return runs


def count_fastq_records(path: str | Path) -> int:
    """Count complete FASTQ records and reject truncated files."""
    lines = 0
    with gzip.open(path, "rb") as handle:
        for lines, _line in enumerate(handle, start=1):
            pass
    if lines % 4:
        raise ValueError(f"truncated FASTQ (line count not divisible by four): {path}")
    return lines // 4


def validate_fastq_pair(run: RunSpec, fastq_dir: str | Path) -> dict[str, Any]:
    paths = [Path(fastq_dir) / run.fastq_r1, Path(fastq_dir) / run.fastq_r2]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    observed_md5 = []
    for path in paths:
        digest = hashlib.md5()  # nosec B324 - archival checksum supplied by ENA
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        observed_md5.append(digest.hexdigest())
    if observed_md5 != [run.md5_r1, run.md5_r2]:
        raise ValueError(f"ENA MD5 mismatch for {run.run_accession}")
    reads = [count_fastq_records(path) for path in paths]
    if reads[0] != reads[1]:
        raise ValueError(f"unequal mate counts for {run.run_accession}: {reads}")
    return {"read_pairs": reads[0], "read_ends": reads[0] * 2, "md5_status": "pass"}


def write_reference_bundle(references: dict[str, AmpliconReference], output_dir: str | Path) -> None:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for locus, reference in references.items():
        (destination / f"{locus}.fa").write_text(
            f">{reference.contig}\n{reference.sequence_plus}\n", encoding="ascii"
        )
        row = asdict(reference)
        row["amplicon_length"] = len(reference.sequence_plus)
        rows.append(row)
    write_tsv(
        destination / "amplicon_references.tsv",
        [
            "locus",
            "contig",
            "build",
            "chrom",
            "start_0based",
            "end_0based_exclusive",
            "sequence_plus",
            "guide_start_0based",
            "guide_end_0based_exclusive",
            "guide_strand",
            "source_cache",
            "source_url",
            "source_sha256",
            "source_sequence_sha256",
            "amplicon_length",
        ],
        rows,
    )


def analysis_windows(reference: AmpliconReference) -> tuple[AnalysisWindow, AnalysisWindow]:
    length = len(reference.sequence_plus)
    if length <= 2 * ANALYSIS_WINDOW_LENGTH:
        raise ValueError("D5 amplicon is too short for two non-overlapping analysis windows")
    target = AnalysisWindow(
        name="target",
        start_0based=0,
        end_0based_exclusive=ANALYSIS_WINDOW_LENGTH,
        sequence_plus=reference.sequence_plus[:ANALYSIS_WINDOW_LENGTH],
        contains_protospacer=True,
    )
    distal_start = length - ANALYSIS_WINDOW_LENGTH
    distal = AnalysisWindow(
        name="distal",
        start_0based=distal_start,
        end_0based_exclusive=length,
        sequence_plus=reference.sequence_plus[distal_start:],
        contains_protospacer=False,
    )
    if not (
        target.start_0based
        <= reference.guide_start_0based
        < reference.guide_end_0based_exclusive
        <= target.end_0based_exclusive
    ):
        raise ValueError("protospacer is not fully contained in the target analysis window")
    return target, distal


def _relative_to_workspace(path: Path, workspace: Path) -> str:
    return "/work/" + path.resolve().relative_to(workspace.resolve()).as_posix()


def _docker(workspace: Path, arguments: list[str], *, capture: bool = False) -> str:
    command = [
        "docker",
        "run",
        "--rm",
        "--volume",
        f"{workspace.resolve()}:/work",
        "--workdir",
        "/work",
        "--entrypoint",
        "/bin/bash",
        CRISPRESSO_IMAGE,
        "-lc",
        " ".join(arguments),
    ]
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=capture,
    )
    return (completed.stdout + completed.stderr) if capture else ""


def build_bowtie2_index(workspace: Path, reference_fasta: Path) -> None:
    fasta = _relative_to_workspace(reference_fasta, workspace)
    prefix = _relative_to_workspace(reference_fasta.with_suffix(""), workspace)
    _docker(workspace, ["bowtie2-build", fasta, prefix], capture=True)


def align_run(
    workspace: Path,
    run: RunSpec,
    reference_fasta: Path,
    fastq_dir: Path,
    run_dir: Path,
    *,
    threads: int,
    artifact_label: str | None = None,
) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    prefix = _relative_to_workspace(reference_fasta.with_suffix(""), workspace)
    r1 = _relative_to_workspace(fastq_dir / run.fastq_r1, workspace)
    r2 = _relative_to_workspace(fastq_dir / run.fastq_r2, workspace)
    stem = run.run_accession if artifact_label is None else f"{run.run_accession}.{artifact_label}"
    bam_path = run_dir / f"{stem}.bam"
    bam = _relative_to_workspace(bam_path, workspace)
    log_name = "bowtie2.log" if artifact_label is None else f"bowtie2_{artifact_label}.log"
    log = _relative_to_workspace(run_dir / log_name, workspace)
    # All interpolated paths derive from Path.relative_to and all names come from a validated TSV.
    pipeline = [
        "set", "-o", "pipefail;",
        "bowtie2", "--sensitive-local", "--no-unal", "--threads", str(threads),
        "-x", prefix, "-1", r1, "-2", r2, "2>", log,
        "|", "samtools", "sort", "-@", str(threads), "-o", bam, "-;",
        "samtools", "index", "-@", str(threads), bam,
    ]
    _docker(workspace, pipeline)
    return bam_path


def run_crispresso(
    workspace: Path,
    run: RunSpec,
    reference: AmpliconReference,
    bam_path: Path,
    run_dir: Path,
    *,
    threads: int,
) -> dict[str, Path]:
    bam = _relative_to_workspace(bam_path, workspace)
    output = _relative_to_workspace(run_dir, workspace)
    results: dict[str, Path] = {}
    for window in analysis_windows(reference):
        name = f"{run.run_accession}_{window.name}"
        # samtools/CRISPResso region coordinates are 1-based inclusive.
        bam_region = (
            f"{reference.contig}:{window.start_0based + 1}-{window.end_0based_exclusive}"
        )
        arguments = [
            "CRISPResso",
            "--bam_input", bam,
            "--bam_chr_loc", bam_region,
            "--amplicon_seq", window.sequence_plus,
            "--amplicon_name", f"{reference.locus}_{window.name}_window",
            "--quantification_window_coordinates",
            (
                f"{reference.guide_start_0based - window.start_0based}-"
                f"{reference.guide_end_0based_exclusive - window.start_0based - 1}"
                if window.contains_protospacer
                else f"0-{len(window.sequence_plus) - 1}"
            ),
            "--exclude_bp_from_left", "0",
            "--exclude_bp_from_right", "0",
            "--suppress_plots",
            "--suppress_report",
            "--output_folder", output,
            "--name", name,
            "--n_processes", str(threads),
        ]
        if window.contains_protospacer:
            arguments.extend(
                [
                    "--guide_seq", GUIDE_SEQUENCE,
                    "--guide_name", "sgHBG-1",
                    "--base_editor_output",
                    "--conversion_nuc_from", "C",
                    "--conversion_nuc_to", "T",
                    "--write_detailed_allele_table",
                ]
            )
        log_text = _docker(workspace, arguments, capture=True)
        (run_dir / f"crispresso_{window.name}.log").write_text(log_text, encoding="utf-8")
        result = run_dir / f"CRISPResso_on_{name}"
        if not result.is_dir():
            raise RuntimeError(f"CRISPResso output directory was not created: {result}")
        results[window.name] = result
    return results


def bam_counts(workspace: Path, bam_path: Path) -> dict[str, int]:
    bam = _relative_to_workspace(bam_path, workspace)
    values: dict[str, int] = {}
    for name, flags in {
        "mapped_read_ends": ["-F", "4"],
        "primary_mapped_read_ends": ["-F", "2308"],
        "properly_paired_read_ends": ["-f", "2", "-F", "2308"],
    }.items():
        output = _docker(workspace, ["samtools", "view", "-c", *flags, bam], capture=True)
        values[name] = int(output.strip())
    return values


def _find_one(directory: Path, pattern: str) -> Path:
    matches = list(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"expected one {pattern!r} in {directory}, found {len(matches)}")
    return matches[0]


def crispresso_mapping_counts(crispresso_dir: Path) -> dict[str, int | float]:
    """Read CRISPResso's own input/alignment accounting for one analysis window."""
    source = crispresso_dir / "CRISPResso_mapping_statistics.txt"
    with source.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"expected one CRISPResso mapping-statistics row in {source}")
    input_reads = int(rows[0]["READS IN INPUTS"])
    preprocessed_reads = int(rows[0]["READS AFTER PREPROCESSING"])
    aligned_reads = int(rows[0]["READS ALIGNED"])
    if not (0 <= aligned_reads <= preprocessed_reads <= input_reads):
        raise ValueError(f"invalid CRISPResso read accounting in {source}")
    return {
        "reads_input": input_reads,
        "reads_after_preprocessing": preprocessed_reads,
        "reads_aligned": aligned_reads,
        "alignment_fraction": 0.0 if input_reads == 0 else aligned_reads / input_reads,
    }


def _open_crispresso_alleles(path: Path):
    if path.suffix == ".zip":
        archive = ZipFile(path)
        names = archive.namelist()
        if len(names) != 1:
            archive.close()
            raise ValueError(f"expected one member in {path}, found {names}")
        handle = archive.open(names[0], "r")
        return archive, handle
    return None, path.open("rb")


def normalize_alleles(
    run: RunSpec,
    crispresso_dirs: dict[str, Path],
    output_path: Path,
) -> tuple[int, dict[str, dict[str, int]]]:
    rows: list[dict[str, object]] = []
    stats: dict[str, dict[str, int]] = {}
    for window_name in ("target", "distal"):
        source = _find_one(crispresso_dirs[window_name], "*Alleles_frequency_table*.zip")
        archive, binary = _open_crispresso_alleles(source)
        window_rows = 0
        window_reads = 0
        try:
            text = (line.decode("utf-8") for line in binary)
            reader = csv.DictReader(text, delimiter="\t")
            for row in reader:
                aligned = row["Aligned_Sequence"]
                reference = row["Reference_Sequence"]
                reads = int(row["#Reads"])
                window_rows += 1
                window_reads += reads
                key_payload = (
                    f"{run.locus}\0{window_name}\0{aligned}\0{reference}".encode("ascii")
                )
                rows.append(
                    {
                        "dataset_id": DATASET_ID,
                        "run_accession": run.run_accession,
                        "editor": run.editor,
                        "locus": run.locus,
                        "replicate": run.replicate,
                        "analysis_window": window_name,
                        "allele_key": hashlib.sha256(key_payload).hexdigest(),
                        "aligned_sequence": aligned,
                        "reference_sequence": reference,
                        "read_count": reads,
                        "read_percentage": row["%Reads"],
                        "n_deleted": row["n_deleted"],
                        "n_inserted": row["n_inserted"],
                        "n_mutated": row["n_mutated"],
                    }
                )
        finally:
            binary.close()
            if archive is not None:
                archive.close()
        stats[window_name] = {"allele_rows": window_rows, "accounted_reads": window_reads}
    write_tsv(
        output_path,
        [
            "dataset_id", "run_accession", "editor", "locus", "replicate",
            "analysis_window", "allele_key", "aligned_sequence", "reference_sequence",
            "read_count", "read_percentage", "n_deleted", "n_inserted", "n_mutated",
        ],
        rows,
    )
    return len(rows), stats


def normalize_positions(
    run: RunSpec,
    reference: AmpliconReference,
    crispresso_dirs: dict[str, Path],
    output_path: Path,
) -> tuple[int, int, int]:
    windows = {window.name: window for window in analysis_windows(reference)}
    position_data: dict[int, tuple[str, dict[str, int], dict[str, str]]] = {}
    required = {"A", "C", "G", "T", "N"}
    for window_name, crispresso_dir in crispresso_dirs.items():
        window = windows[window_name]
        expected_aligned = int(crispresso_mapping_counts(crispresso_dir)["reads_aligned"])
        frequency_path = _find_one(crispresso_dir, "*Nucleotide_frequency_table.txt")
        percentage_path = _find_one(crispresso_dir, "*Nucleotide_percentage_table.txt")
        with frequency_path.open("r", encoding="utf-8", newline="") as count_handle:
            counts = list(csv.reader(count_handle, delimiter="\t"))
        with percentage_path.open("r", encoding="utf-8", newline="") as pct_handle:
            percentages = list(csv.reader(pct_handle, delimiter="\t"))
        if len(counts) != len(percentages) or len(counts) < 2:
            raise ValueError(f"invalid CRISPResso nucleotide tables in {crispresso_dir}")
        positions = counts[0][1:]
        if positions != percentages[0][1:] or len(positions) != len(window.sequence_plus):
            raise ValueError(
                f"CRISPResso window/reference length mismatch for {run.run_accession}:{window_name}"
            )
        count_by_base = {row[0]: row[1:] for row in counts[1:]}
        pct_by_base = {row[0]: row[1:] for row in percentages[1:]}
        if not required.issubset(count_by_base) or not required.issubset(pct_by_base):
            raise ValueError(
                f"missing nucleotide rows in CRISPResso tables: {required - count_by_base.keys()}"
            )
        for local_position in range(len(window.sequence_plus)):
            amplicon_position = window.start_0based + local_position
            nucleotide_counts = {
                base: int(float(count_by_base[base][local_position])) for base in required
            }
            nucleotide_percentages = {
                base: pct_by_base[base][local_position] for base in required
            }
            # CRISPResso represents a deletion at a reference position with the `-` row. Include it
            # in per-position depth and non-reference frequency rather than treating deletion-bearing
            # reads as missing coverage. Older fixed test tables without
            # that row remain readable and normalize to zero deletions.
            nucleotide_counts["-"] = (
                int(float(count_by_base["-"][local_position])) if "-" in count_by_base else 0
            )
            nucleotide_percentages["-"] = (
                pct_by_base["-"][local_position] if "-" in pct_by_base else "0"
            )
            if sum(nucleotide_counts.values()) != expected_aligned:
                raise ValueError(
                    "CRISPResso per-position read accounting mismatch for "
                    f"{run.run_accession}:{window_name}:{local_position}"
                )
            position_data[amplicon_position] = (
                window_name, nucleotide_counts, nucleotide_percentages
            )

    rows: list[dict[str, object]] = []
    min_guide_aligned_depth: int | None = None
    min_guide_base_call_depth: int | None = None
    for position, reference_base in enumerate(reference.sequence_plus):
        if position in position_data:
            window_name, nucleotide_counts, nucleotide_percentages = position_data[position]
            coverage_status = "observed"
        else:
            window_name = "insert_gap"
            nucleotide_counts = {base: 0 for base in required | {"-"}}
            nucleotide_percentages = {base: "" for base in required | {"-"}}
            coverage_status = "unsequenced_between_mates"
        base_call_depth = sum(nucleotide_counts[base] for base in required)
        depth = base_call_depth + nucleotide_counts["-"]
        in_guide = reference.guide_start_0based <= position < reference.guide_end_0based_exclusive
        if in_guide:
            min_guide_aligned_depth = (
                depth if min_guide_aligned_depth is None else min(min_guide_aligned_depth, depth)
            )
            min_guide_base_call_depth = (
                base_call_depth
                if min_guide_base_call_depth is None
                else min(min_guide_base_call_depth, base_call_depth)
            )
        # The amplicon is reference-forward while sgHBG-1 maps to the minus strand.  Therefore the
        # specified sense-strand C>T conversion is represented as G>A on this sequence.
        conversion_frequency = ""
        if reference_base == "G" and depth:
            conversion_frequency = nucleotide_counts["A"] / depth
        rows.append(
            {
                "dataset_id": DATASET_ID,
                "run_accession": run.run_accession,
                "editor": run.editor,
                "locus": run.locus,
                "replicate": run.replicate,
                "build": reference.build,
                "chrom": reference.chrom,
                "amplicon_position_0based": position,
                "genomic_position_1based": reference.start_0based + position + 1,
                "reference_base_plus": reference_base,
                "in_protospacer": str(in_guide).lower(),
                "analysis_window": window_name,
                "coverage_status": coverage_status,
                "depth": depth,
                "base_call_depth": base_call_depth,
                "a_count": nucleotide_counts["A"],
                "c_count": nucleotide_counts["C"],
                "g_count": nucleotide_counts["G"],
                "t_count": nucleotide_counts["T"],
                "n_count": nucleotide_counts["N"],
                "deletion_count": nucleotide_counts["-"],
                "a_frequency": nucleotide_percentages["A"],
                "c_frequency": nucleotide_percentages["C"],
                "g_frequency": nucleotide_percentages["G"],
                "t_frequency": nucleotide_percentages["T"],
                "n_frequency": nucleotide_percentages["N"],
                "deletion_frequency": nucleotide_percentages["-"],
                "non_reference_frequency": "" if not depth else (depth - nucleotide_counts.get(reference_base, 0)) / depth,
                "strand_normalized_c_to_t_frequency": conversion_frequency,
            }
        )
    if min_guide_aligned_depth is None or min_guide_base_call_depth is None:
        raise ValueError("protospacer is outside the amplicon")
    write_tsv(
        output_path,
        [
            "dataset_id", "run_accession", "editor", "locus", "replicate", "build", "chrom",
            "amplicon_position_0based", "genomic_position_1based", "reference_base_plus",
            "in_protospacer", "analysis_window", "coverage_status", "depth", "base_call_depth",
            "a_count",
            "c_count", "g_count", "t_count", "n_count", "deletion_count",
            "a_frequency", "c_frequency", "g_frequency", "t_frequency", "n_frequency",
            "deletion_frequency",
            "non_reference_frequency", "strand_normalized_c_to_t_frequency",
        ],
        rows,
    )
    return len(rows), min_guide_aligned_depth, min_guide_base_call_depth


def remove_generated_bams(run_dir: Path, run_accession: str) -> None:
    """Drop reproducible mapping intermediates after normalized tables have been written."""
    exact = [
        run_dir / f"{stem}{suffix}"
        for stem in (run_accession, f"{run_accession}.mapped_end_gate")
        for suffix in (".bam", ".bam.bai")
    ]
    for path in exact:
        path.unlink(missing_ok=True)
    for crispresso_dir in run_dir.glob("CRISPResso_on_*"):
        for filename in (
            "CRISPResso_output.bam",
            "CRISPResso_output.bam.bai",
            "CRISPResso_output.bam.sam",
        ):
            (crispresso_dir / filename).unlink(missing_ok=True)


def _run_qc_failures(
    *,
    mapped_fraction: float,
    minimum_base_call_depth: int,
    position_count: int,
    expected_position_count: int,
    target_allele_reads: int,
    distal_allele_reads: int,
    target_aligned_reads: int,
    distal_aligned_reads: int,
) -> list[str]:
    failures: list[str] = []
    if mapped_fraction < MIN_MAPPED_END_FRACTION:
        failures.append("mapped_end_fraction_below_minimum")
    if minimum_base_call_depth < MIN_PROTOSPACER_DEPTH:
        failures.append("protospacer_depth_below_minimum")
    if position_count != expected_position_count:
        failures.append("position_count_mismatch")
    if target_allele_reads != target_aligned_reads:
        failures.append("target_allele_read_accounting_mismatch")
    if distal_allele_reads != distal_aligned_reads:
        failures.append("distal_allele_read_accounting_mismatch")
    return failures


def _mapped_end_gate_fields(reference: AmpliconReference) -> dict[str, Any]:
    return {
        "mapped_end_gate_schema_version": MAPPED_END_GATE_SCHEMA_VERSION,
        "mapped_end_reference_geometry": MAPPED_END_REFERENCE_GEOMETRY,
        "mapped_end_reference_start_0based": reference.start_0based,
        "mapped_end_reference_end_0based_exclusive": reference.end_0based_exclusive,
        "mapped_end_reference_length": len(reference.sequence_plus),
        "mapped_end_reference_sequence_sha256": hashlib.sha256(
            reference.sequence_plus.encode("ascii")
        ).hexdigest(),
    }


def _run_qc_annotations(run_accession: str) -> dict[str, Any]:
    return {
        "qc_note": RUN_QC_NOTES.get(run_accession, ""),
        "held_out_sensitivity_flag": run_accession in HELD_OUT_SENSITIVITY_RUNS,
    }


def renormalize_existing_run(
    run: RunSpec,
    reference: AmpliconReference,
    run_dir: Path,
) -> dict[str, Any]:
    """Refresh the common result tables from fixed CRISPResso files without rereading FASTQs."""
    qc_path = run_dir / "run_qc.json"
    prior = json.loads(qc_path.read_text(encoding="utf-8"))
    if prior.get("run_accession") != run.run_accession:
        raise ValueError(f"run provenance mismatch in {qc_path}")
    crispresso_dirs = {
        window: run_dir / f"CRISPResso_on_{run.run_accession}_{window}"
        for window in ("target", "distal")
    }
    window_mapping: dict[str, int | float] = {}
    for window_name, crispresso_dir in crispresso_dirs.items():
        for field, value in crispresso_mapping_counts(crispresso_dir).items():
            window_mapping[f"{window_name}_crispresso_{field}"] = value
    position_count, min_aligned_depth, min_base_call_depth = normalize_positions(
        run, reference, crispresso_dirs, run_dir / "position_frequencies.tsv"
    )
    allele_count, allele_stats = normalize_alleles(run, crispresso_dirs, run_dir / "alleles.tsv")
    mapped_fraction = float(prior["primary_mapped_end_fraction"])
    failures = _run_qc_failures(
        mapped_fraction=mapped_fraction,
        minimum_base_call_depth=min_base_call_depth,
        position_count=position_count,
        expected_position_count=len(reference.sequence_plus),
        target_allele_reads=allele_stats["target"]["accounted_reads"],
        distal_allele_reads=allele_stats["distal"]["accounted_reads"],
        target_aligned_reads=int(window_mapping["target_crispresso_reads_aligned"]),
        distal_aligned_reads=int(window_mapping["distal_crispresso_reads_aligned"]),
    )
    prior.pop("minimum_protospacer_depth", None)
    prior.pop("allele_accounted_reads", None)
    prior.update(
        {
            **window_mapping,
            "minimum_protospacer_aligned_depth": min_aligned_depth,
            "minimum_protospacer_base_call_depth": min_base_call_depth,
            "position_rows": position_count,
            "allele_rows": allele_count,
            "target_allele_rows": allele_stats["target"]["allele_rows"],
            "target_allele_accounted_reads": allele_stats["target"]["accounted_reads"],
            "distal_allele_rows": allele_stats["distal"]["allele_rows"],
            "distal_allele_accounted_reads": allele_stats["distal"]["accounted_reads"],
            "analysis_window_length": ANALYSIS_WINDOW_LENGTH,
            "normalization_schema_version": NORMALIZATION_SCHEMA_VERSION,
            "qc_result": "pass" if not failures else "fail",
            "qc_failures": ";".join(failures),
        }
    )
    qc_path.write_text(json.dumps(prior, indent=2) + "\n", encoding="utf-8")
    return prior


def process_run(
    workspace: Path,
    run: RunSpec,
    reference: AmpliconReference,
    mapped_end_reference: AmpliconReference,
    fastq_dir: Path,
    output_dir: Path,
    *,
    threads: int,
) -> dict[str, Any]:
    run_dir = output_dir / "runs" / run.run_accession
    input_qc = validate_fastq_pair(run, fastq_dir)
    reference_fasta = output_dir / "references" / f"{run.locus}.fa"
    bam_path = align_run(workspace, run, reference_fasta, fastq_dir, run_dir, threads=threads)
    mapped_end_fasta = output_dir / "references" / "mapped_end_gate" / f"{run.locus}.fa"
    mapped_end_bam = align_run(
        workspace,
        run,
        mapped_end_fasta,
        fastq_dir,
        run_dir,
        threads=threads,
        artifact_label="mapped_end_gate",
    )
    counts = bam_counts(workspace, mapped_end_bam)
    crispresso_dirs = run_crispresso(
        workspace, run, reference, bam_path, run_dir, threads=threads
    )
    window_mapping: dict[str, int | float] = {}
    for window_name, crispresso_dir in crispresso_dirs.items():
        for field, value in crispresso_mapping_counts(crispresso_dir).items():
            window_mapping[f"{window_name}_crispresso_{field}"] = value
    position_count, min_aligned_depth, min_base_call_depth = normalize_positions(
        run, reference, crispresso_dirs, run_dir / "position_frequencies.tsv"
    )
    allele_count, allele_stats = normalize_alleles(run, crispresso_dirs, run_dir / "alleles.tsv")
    remove_generated_bams(run_dir, run.run_accession)
    mapped_fraction = counts["primary_mapped_read_ends"] / input_qc["read_ends"]
    failures = _run_qc_failures(
        mapped_fraction=mapped_fraction,
        minimum_base_call_depth=min_base_call_depth,
        position_count=position_count,
        expected_position_count=len(reference.sequence_plus),
        target_allele_reads=allele_stats["target"]["accounted_reads"],
        distal_allele_reads=allele_stats["distal"]["accounted_reads"],
        target_aligned_reads=int(window_mapping["target_crispresso_reads_aligned"]),
        distal_aligned_reads=int(window_mapping["distal_crispresso_reads_aligned"]),
    )
    result = {
        "dataset_id": DATASET_ID,
        **asdict(run),
        "reference_build": reference.build,
        "reference_chrom": reference.chrom,
        "amplicon_start_0based": reference.start_0based,
        "amplicon_end_0based_exclusive": reference.end_0based_exclusive,
        "amplicon_length": len(reference.sequence_plus),
        "guide_start_0based": reference.guide_start_0based,
        "guide_end_0based_exclusive": reference.guide_end_0based_exclusive,
        **_mapped_end_gate_fields(mapped_end_reference),
        **input_qc,
        **counts,
        **window_mapping,
        "primary_mapped_end_fraction": mapped_fraction,
        "minimum_protospacer_aligned_depth": min_aligned_depth,
        "minimum_protospacer_base_call_depth": min_base_call_depth,
        "position_rows": position_count,
        "allele_rows": allele_count,
        "target_allele_rows": allele_stats["target"]["allele_rows"],
        "target_allele_accounted_reads": allele_stats["target"]["accounted_reads"],
        "distal_allele_rows": allele_stats["distal"]["allele_rows"],
        "distal_allele_accounted_reads": allele_stats["distal"]["accounted_reads"],
        "crispresso_version": CRISPRESSO_VERSION,
        "crispresso_image": CRISPRESSO_IMAGE,
        "analysis_window_length": ANALYSIS_WINDOW_LENGTH,
        "normalization_schema_version": NORMALIZATION_SCHEMA_VERSION,
        "generated_bams_retained": False,
        "qc_result": "pass" if not failures else "fail",
        "qc_failures": ";".join(failures),
        **_run_qc_annotations(run.run_accession),
    }
    (run_dir / "run_qc.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


QC_FIELDS = [
    "dataset_id", "run_accession", "study_accession", "sample_accession",
    "experiment_accession", "editor", "locus", "replicate", "cell_type", "assay",
    "fastq_r1", "fastq_r2", "md5_r1", "md5_r2",
    "reference_build", "reference_chrom", "amplicon_start_0based",
    "amplicon_end_0based_exclusive", "amplicon_length", "guide_start_0based",
    "guide_end_0based_exclusive", "mapped_end_gate_schema_version",
    "mapped_end_reference_geometry", "mapped_end_reference_start_0based",
    "mapped_end_reference_end_0based_exclusive", "mapped_end_reference_length",
    "mapped_end_reference_sequence_sha256", "read_pairs", "read_ends", "md5_status",
    "mapped_read_ends", "primary_mapped_read_ends", "properly_paired_read_ends",
    "target_crispresso_reads_input", "target_crispresso_reads_after_preprocessing",
    "target_crispresso_reads_aligned", "target_crispresso_alignment_fraction",
    "distal_crispresso_reads_input", "distal_crispresso_reads_after_preprocessing",
    "distal_crispresso_reads_aligned", "distal_crispresso_alignment_fraction",
    "primary_mapped_end_fraction", "minimum_protospacer_aligned_depth",
    "minimum_protospacer_base_call_depth", "position_rows",
    "allele_rows", "target_allele_rows", "target_allele_accounted_reads",
    "distal_allele_rows", "distal_allele_accounted_reads",
    "crispresso_version", "crispresso_image",
    "analysis_window_length", "normalization_schema_version", "generated_bams_retained",
    "qc_result", "qc_failures", "qc_note", "held_out_sensitivity_flag",
]


def render_qc_summary(results: list[dict[str, Any]]) -> str:
    """Render the technical QC summary without including measured outcome tables."""
    pass_count = sum(row["qc_result"] == "pass" for row in results)
    lines = [
        "# D5 amplicon outcomes — QC summary",
        "",
        f"Generated: {QC_REPORT_DATE}",
        "",
        "Scope: measured-read ingestion and technical QC only; no tool predictions were loaded or compared.",
        "",
        f"- Runs processed: {len(results)} ({pass_count} pass, {len(results) - pass_count} fail)",
        "- Primary mapped-end reference: staged F1/R1 flank-inclusive geometry, applied "
        "uniformly to every run",
        f"- Minimum primary mapped-end fraction: `{MIN_MAPPED_END_FRACTION:.2f}`",
        f"- Minimum callable depth at every protospacer position: `{MIN_PROTOSPACER_DEPTH:,}`",
        f"- CRISPResso2: `{CRISPRESSO_VERSION}`; normalization schema: `{NORMALIZATION_SCHEMA_VERSION}`",
        "",
        "| Run | Editor | Locus | Rep | Read pairs | Primary mapped ends | CRISPResso target | CRISPResso distal | Min callable protospacer depth | Alleles | QC | Failure flags | QC note | Held-out sensitivity |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- |",
    ]
    for row in results:
        lines.append(
            "| {run_accession} | {editor} | {locus} | {replicate} | {read_pairs:,} | "
            "{mapped:.2%} | {target:.2%} | {distal:.2%} | {depth:,} | {alleles:,} | "
            "{qc_result} | {qc_failures} | {display_qc_note} | {sensitivity} |".format(
                **row,
                mapped=float(row["primary_mapped_end_fraction"]),
                target=float(row["target_crispresso_alignment_fraction"]),
                distal=float(row["distal_crispresso_alignment_fraction"]),
                depth=int(row["minimum_protospacer_base_call_depth"]),
                alleles=int(row["allele_rows"]),
                sensitivity=(
                    "flagged" if bool(row.get("held_out_sensitivity_flag")) else "—"
                ),
                display_qc_note=row.get("qc_note") or "—",
            )
        )
    lines.extend(
        [
            "",
            "Full normalized position, allele, log, and run-QC artifacts are gitignored below this directory.",
            "Registry QC status remains `stub`; changing it requires the protocol's documented review.",
            "",
        ]
    )
    return "\n".join(lines)


def _write_qc_artifacts(destination: Path, results: list[dict[str, Any]]) -> None:
    write_tsv(destination / "run_qc.tsv", QC_FIELDS, results)
    summary = {
        "dataset_id": DATASET_ID,
        "run_count": len(results),
        "pass_count": sum(row["qc_result"] == "pass" for row in results),
        "fail_count": sum(row["qc_result"] != "pass" for row in results),
        "thresholds": {
            "minimum_primary_mapped_end_fraction": MIN_MAPPED_END_FRACTION,
            "minimum_protospacer_base_call_depth": MIN_PROTOSPACER_DEPTH,
        },
        "mapped_end_gate": {
            "schema_version": MAPPED_END_GATE_SCHEMA_VERSION,
            "reference_geometry": MAPPED_END_REFERENCE_GEOMETRY,
            "uniform_run_count": len(results),
        },
        "held_out_sensitivity_flags": [
            {
                "run_accession": row["run_accession"],
                "qc_note": row["qc_note"],
            }
            for row in results
            if row.get("held_out_sensitivity_flag")
        ],
        "crispresso_version": CRISPRESSO_VERSION,
        "crispresso_image": CRISPRESSO_IMAGE,
        "normalization_schema_version": NORMALIZATION_SCHEMA_VERSION,
        "scope_guard": (
            "mapped-end rerun was outcome-blind; no editing-frequency or allele-sequence field "
            "was loaded"
        ),
    }
    (destination / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (destination / "QC_SUMMARY.md").write_text(
        render_qc_summary(results), encoding="utf-8"
    )


def _recheck_existing_qc(
    run: RunSpec,
    mapped_end_reference: AmpliconReference,
    input_qc: dict[str, Any],
    counts: dict[str, int],
    run_dir: Path,
) -> dict[str, Any]:
    """Repeat technical QC checks from summary counts without opening outcome files."""
    qc_path = run_dir / "run_qc.json"
    prior = json.loads(qc_path.read_text(encoding="utf-8"))
    if prior.get("run_accession") != run.run_accession:
        raise ValueError(f"run provenance mismatch in {qc_path}")
    for field in ("read_pairs", "read_ends"):
        if int(prior[field]) != int(input_qc[field]):
            raise ValueError(f"{field} changed for {run.run_accession}")

    mapped_fraction = counts["primary_mapped_read_ends"] / input_qc["read_ends"]
    failures = _run_qc_failures(
        mapped_fraction=mapped_fraction,
        minimum_base_call_depth=int(prior["minimum_protospacer_base_call_depth"]),
        position_count=int(prior["position_rows"]),
        expected_position_count=int(prior["amplicon_length"]),
        target_allele_reads=int(prior["target_allele_accounted_reads"]),
        distal_allele_reads=int(prior["distal_allele_accounted_reads"]),
        target_aligned_reads=int(prior["target_crispresso_reads_aligned"]),
        distal_aligned_reads=int(prior["distal_crispresso_reads_aligned"]),
    )
    prior.update(
        {
            **_mapped_end_gate_fields(mapped_end_reference),
            **input_qc,
            **counts,
            "primary_mapped_end_fraction": mapped_fraction,
            "qc_result": "pass" if not failures else "fail",
            "qc_failures": ";".join(failures),
            "mapped_end_gate_rerun_scope": "outcome_blind",
            **_run_qc_annotations(run.run_accession),
        }
    )
    qc_path.write_text(json.dumps(prior, indent=2) + "\n", encoding="utf-8")
    return prior


def rerun_mapped_end_qc(
    workspace: str | Path,
    *,
    manifest_path: str | Path,
    fastq_dir: str | Path,
    cache_dir: str | Path,
    output_dir: str | Path,
    threads: int = 4,
    only_runs: set[str] | None = None,
    docker_stage_parent: str | Path | None = None,
) -> list[dict[str, Any]]:
    """Outcome-blind mapped-end/QC rerun against staged F1/R1 flank-inclusive references."""
    root = Path(workspace).resolve()
    source_fastq = Path(fastq_dir).resolve()
    destination = Path(output_dir).resolve()
    references = {
        locus: derive_mapped_end_reference(spec, cache_dir) for locus, spec in LOCI.items()
    }
    runs = read_run_manifest(manifest_path)
    if only_runs:
        unknown = only_runs - {run.run_accession for run in runs}
        if unknown:
            raise ValueError(f"requested run(s) absent from manifest: {sorted(unknown)}")
        runs = [run for run in runs if run.run_accession in only_runs]
    if len(runs) != 20 and only_runs is None:
        raise ValueError(f"uniform D5 mapped-end rerun requires 20 runs, found {len(runs)}")

    results: list[dict[str, Any]] = []
    if docker_stage_parent is None:
        reference_dir = destination / "references" / "mapped_end_gate"
        write_reference_bundle(references, reference_dir)
        for locus in references:
            build_bowtie2_index(root, reference_dir / f"{locus}.fa")
        for run in runs:
            run_dir = destination / "runs" / run.run_accession
            input_qc = validate_fastq_pair(run, source_fastq)
            bam_path = align_run(
                root,
                run,
                reference_dir / f"{run.locus}.fa",
                source_fastq,
                run_dir,
                threads=threads,
                artifact_label="mapped_end_gate",
            )
            counts = bam_counts(root, bam_path)
            remove_generated_bams(run_dir, run.run_accession)
            results.append(
                _recheck_existing_qc(
                    run, references[run.locus], input_qc, counts, run_dir
                )
            )
    else:
        stage_parent = Path(docker_stage_parent).resolve()
        stage_parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix="base-edit-safety-d5-mapped-end-", dir=stage_parent
        ) as name:
            stage_root = Path(name)
            stage_fastq = stage_root / "fastq"
            reference_dir = stage_root / "out" / "references" / "mapped_end_gate"
            stage_fastq.mkdir()
            write_reference_bundle(references, reference_dir)
            for locus in references:
                build_bowtie2_index(stage_root, reference_dir / f"{locus}.fa")
            shutil.copytree(
                reference_dir,
                destination / "references" / "mapped_end_gate",
                dirs_exist_ok=True,
            )
            for run in runs:
                for filename in (run.fastq_r1, run.fastq_r2):
                    shutil.copy2(source_fastq / filename, stage_fastq / filename)
                input_qc = validate_fastq_pair(run, stage_fastq)
                stage_run = stage_root / "out" / "runs" / run.run_accession
                bam_path = align_run(
                    stage_root,
                    run,
                    reference_dir / f"{run.locus}.fa",
                    stage_fastq,
                    stage_run,
                    threads=threads,
                    artifact_label="mapped_end_gate",
                )
                counts = bam_counts(stage_root, bam_path)
                remove_generated_bams(stage_run, run.run_accession)
                published_run = destination / "runs" / run.run_accession
                published_run.mkdir(parents=True, exist_ok=True)
                shutil.copy2(
                    stage_run / "bowtie2_mapped_end_gate.log",
                    published_run / "bowtie2_mapped_end_gate.log",
                )
                results.append(
                    _recheck_existing_qc(
                        run,
                        references[run.locus],
                        input_qc,
                        counts,
                        published_run,
                    )
                )
                for filename in (run.fastq_r1, run.fastq_r2):
                    (stage_fastq / filename).unlink()
    _write_qc_artifacts(destination, results)
    return results


def run_pipeline(
    workspace: str | Path,
    *,
    manifest_path: str | Path,
    fastq_dir: str | Path,
    cache_dir: str | Path,
    output_dir: str | Path,
    threads: int = 4,
    only_runs: set[str] | None = None,
    docker_stage_parent: str | Path | None = None,
    resume: bool = False,
) -> list[dict[str, Any]]:
    root = Path(workspace).resolve()
    destination = Path(output_dir).resolve()
    references = {locus: derive_amplicon(spec, cache_dir) for locus, spec in LOCI.items()}
    mapped_end_references = {
        locus: derive_mapped_end_reference(spec, cache_dir) for locus, spec in LOCI.items()
    }
    runs = read_run_manifest(manifest_path)
    if only_runs:
        unknown = only_runs - {run.run_accession for run in runs}
        if unknown:
            raise ValueError(f"requested run(s) absent from manifest: {sorted(unknown)}")
        runs = [run for run in runs if run.run_accession in only_runs]
    if docker_stage_parent is None:
        write_reference_bundle(references, destination / "references")
        write_reference_bundle(
            mapped_end_references, destination / "references" / "mapped_end_gate"
        )
        for locus in references:
            build_bowtie2_index(root, destination / "references" / f"{locus}.fa")
            build_bowtie2_index(
                root, destination / "references" / "mapped_end_gate" / f"{locus}.fa"
            )
        results = [
            process_run(
                root,
                run,
                references[run.locus],
                mapped_end_references[run.locus],
                Path(fastq_dir),
                destination,
                threads=threads,
            )
            for run in runs
        ]
    else:
        # If Docker Desktop cannot mount a macOS Desktop checkout, stage one compressed pair at a
        # time under a Docker-visible temporary directory,
        # then publish each completed run into the requested ignored output tree.
        stage_parent = Path(docker_stage_parent).resolve()
        stage_parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="base-edit-safety-d5-", dir=stage_parent) as name:
            stage_root = Path(name)
            stage_fastq = stage_root / "fastq"
            stage_output = stage_root / "out"
            stage_fastq.mkdir()
            write_reference_bundle(references, stage_output / "references")
            write_reference_bundle(
                mapped_end_references, stage_output / "references" / "mapped_end_gate"
            )
            for locus in references:
                build_bowtie2_index(stage_root, stage_output / "references" / f"{locus}.fa")
                build_bowtie2_index(
                    stage_root,
                    stage_output / "references" / "mapped_end_gate" / f"{locus}.fa",
                )
            destination.mkdir(parents=True, exist_ok=True)
            shutil.copytree(
                stage_output / "references", destination / "references", dirs_exist_ok=True
            )
            results = []
            for run in runs:
                published_run = destination / "runs" / run.run_accession
                published_qc = published_run / "run_qc.json"
                if resume and published_qc.is_file():
                    prior = json.loads(published_qc.read_text(encoding="utf-8"))
                    if (
                        prior.get("crispresso_image") == CRISPRESSO_IMAGE
                        and prior.get("analysis_window_length") == ANALYSIS_WINDOW_LENGTH
                        and prior.get("normalization_schema_version")
                        != NORMALIZATION_SCHEMA_VERSION
                        and all(
                            (
                                published_run
                                / f"CRISPResso_on_{run.run_accession}_{window}"
                                / "CRISPResso_mapping_statistics.txt"
                            ).is_file()
                            for window in ("target", "distal")
                        )
                    ):
                        prior = renormalize_existing_run(
                            run, references[run.locus], published_run
                        )
                    required_outputs = [
                        published_run / "position_frequencies.tsv",
                        published_run / "alleles.tsv",
                    ]
                    if (
                        prior.get("crispresso_image") == CRISPRESSO_IMAGE
                        and prior.get("mapped_end_gate_schema_version")
                        == MAPPED_END_GATE_SCHEMA_VERSION
                        and prior.get("mapped_end_reference_geometry")
                        == MAPPED_END_REFERENCE_GEOMETRY
                        and prior.get("normalization_schema_version")
                        == NORMALIZATION_SCHEMA_VERSION
                        and prior.get("position_rows") == len(references[run.locus].sequence_plus)
                        and all(path.is_file() for path in required_outputs)
                    ):
                        results.append(prior)
                        continue
                for filename in (run.fastq_r1, run.fastq_r2):
                    shutil.copy2(Path(fastq_dir) / filename, stage_fastq / filename)
                result = process_run(
                    stage_root,
                    run,
                    references[run.locus],
                    mapped_end_references[run.locus],
                    stage_fastq,
                    stage_output,
                    threads=threads,
                )
                results.append(result)
                shutil.copytree(
                    stage_output / "runs" / run.run_accession,
                    destination / "runs" / run.run_accession,
                    dirs_exist_ok=True,
                )
                for filename in (run.fastq_r1, run.fastq_r2):
                    (stage_fastq / filename).unlink()
    _write_qc_artifacts(destination, results)
    return results


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--manifest", default="validation/D5_FILE_ELIGIBILITY_MANIFEST.tsv")
    parser.add_argument("--fastq-dir", default="data/inputs/d5_fastq")
    parser.add_argument("--cache-dir", default="data/inputs/reference_cache")
    parser.add_argument("--output-dir", default="validation/d5_outcomes")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--run", action="append", dest="runs")
    parser.add_argument(
        "--docker-stage-parent",
        help="Docker-visible parent for temporary staged reads (for example /private/tmp)",
    )
    parser.add_argument("--resume", action="store_true", help="reuse complete matching run outputs")
    parser.add_argument(
        "--mapped-end-qc-only",
        action="store_true",
        help=(
            "outcome-blind rerun of the flank-inclusive mapped-end gate and all QC checks from "
            "existing technical aggregates"
        ),
    )
    args = parser.parse_args(argv)
    common = {
        "manifest_path": args.manifest,
        "fastq_dir": args.fastq_dir,
        "cache_dir": args.cache_dir,
        "output_dir": args.output_dir,
        "threads": args.threads,
        "only_runs": None if not args.runs else set(args.runs),
        "docker_stage_parent": args.docker_stage_parent,
    }
    if args.mapped_end_qc_only:
        rerun_mapped_end_qc(args.workspace, **common)
    else:
        run_pipeline(args.workspace, resume=args.resume, **common)
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised through the CLI
    raise SystemExit(main())
