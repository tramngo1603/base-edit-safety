"""Acquire and normalize D1 CHANGE-seq-BE site tables.

Guide/editor labels are copied only when the publisher supplies them.  In particular, the generic
``ABE`` label and the CD40L editor remain unresolved; this module does not infer a construct.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import asdict
import json
from pathlib import Path
import re
import tarfile
from typing import Iterable

from .common import SourceFile, download_source, read_dict_rows, verify_source_file, write_tsv
from .xlsx import iter_rows


DATASET_ID = "D1"
BUILD = "GRCh38/hg38"
ACCESSED_ON = "2026-08-14"

REGISTRY_RECORDS = (
    {
        "accession": "GSE298535",
        "record_type": "GEO umbrella series",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298535",
        "accessed_on": ACCESSED_ON,
    },
    {
        "accession": "PRJNA1269627",
        "record_type": "BioProject registry record",
        "url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1269627",
        "accessed_on": ACCESSED_ON,
    },
)

SOURCES = (
    SourceFile(
        DATASET_ID,
        "GSE298465 (subseries of GSE298535; PRJNA1269500)",
        "GSE298465_RAW.tar",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298465/suppl/GSE298465_RAW.tar",
        "1a0eb89d701db709ec18b7be40ebabba510d0bc058735813c8bde6c995516713",
        ACCESSED_ON,
    ),
    SourceFile(
        DATASET_ID,
        "GSE298534 (subseries of GSE298535; PRJNA1269623)",
        "GSE298534_combined.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298534/suppl/"
        "GSE298534_5_hybridcapture_guideseq_changeseq_digenomeseq_withVCF_"
        "collapsed_withAnnotation_withRloop_version6.csv.gz",
        "4713b37ffbf571064669ce660ef709340a2ffe57479299301c45f83ac22d2684",
        ACCESSED_ON,
    ),
    SourceFile(
        DATASET_ID,
        "GSE307740 (subseries of GSE298535; PRJNA1327241)",
        "GSE307740_Supplementary_Table_1.xlsx",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307740/suppl/"
        "GSE307740_Supplementary_Table_1.xlsx",
        "adcfa3c49fbe33e75817384936ef80eda2ded9bd9e3f6b9c250756a819b60316",
        ACCESSED_ON,
    ),
    SourceFile(
        DATASET_ID,
        "GSE307740 (subseries of GSE298535; PRJNA1327241)",
        "GSE307740_Supplementary_Table_3.xlsx",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307740/suppl/"
        "GSE307740_Supplementary_Table_3.xlsx",
        "a6affaa902bf7618c8ec733dfa4c35d5a8079c1854efba592b0a71ca6a5e018e",
        ACCESSED_ON,
    ),
)


FIELDS = [
    "dataset_id", "registry_accession", "source_accession", "source_file", "source_row",
    "guide_id", "guide_sequence", "editor", "editor_crosswalk_status", "assay_signal",
    "build", "chrom", "start_0based", "end_0based_exclusive", "start_1based",
    "end_1based", "strand", "source_coordinate", "coordinate_basis", "offtarget_sequence",
    "replicate_1_signal", "replicate_2_signal", "detected_by_changeseq_be", "annotation",
]


def acquire_sources(destination: str | Path) -> list[Path]:
    return [download_source(source, destination) for source in SOURCES]


def extract_cd40l_tar(tar_path: str | Path, destination: str | Path) -> list[Path]:
    """Safely extract the two publisher gzip tables from the GEO tar archive."""
    root = Path(destination).resolve()
    root.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    with tarfile.open(tar_path) as archive:
        for member in archive.getmembers():
            filename = Path(member.name).name
            if not filename.endswith("_identified_matched_Annot.tsv.gz"):
                continue
            target = root / filename
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError(f"tar member has no content: {member.name}")
            with target.open("wb") as output:
                while chunk := stream.read(1024 * 1024):
                    output.write(chunk)
            extracted.append(target)
    if len(extracted) != 2:
        raise ValueError(f"expected two CD40L replicate tables, found {len(extracted)}")
    return sorted(extracted)


def _xlsx_dicts(path: Path, sheet: str) -> Iterable[tuple[int, dict[str, str]]]:
    rows = iter_rows(path, sheet)
    header = next(rows)
    for source_row, values in enumerate(rows, start=2):
        padded = values + [""] * (len(header) - len(values))
        yield source_row, dict(zip(header, padded))


def _interval(start: str, end: str, sequence: str) -> tuple[int, int, int, int, str]:
    start_0 = int(float(start))
    end_0 = int(float(end))
    if end_0 <= start_0:
        raise ValueError(f"invalid source interval: {start}-{end}")
    ungapped = sequence.replace("-", "").replace("N", "")
    basis = "source_start/end treated as 0-based half-open"
    if sequence and ungapped and end_0 - start_0 not in {len(sequence), len(ungapped)}:
        basis += "; interval/sequence lengths differ"
    return start_0, end_0, start_0 + 1, end_0, basis


def _present(value: str) -> bool:
    return value.strip() not in {"", "0", "0.0", "NA", "N/A"}


def parse_cd40l(paths: Iterable[Path]) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    peak_pattern = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")
    for replicate, path in enumerate(sorted(paths), start=1):
        for source_row, row in enumerate(read_dict_rows(path, delimiter="\t"), start=2):
            peak_column = next(key for key in row if key.startswith("PeakID"))
            match = peak_pattern.fullmatch(row[peak_column])
            if match is None:
                raise ValueError(f"invalid HOMER PeakID in {path}:{source_row}")
            chrom, start_text, end_text = match.groups()
            start_0, end_0 = int(start_text), int(end_text)
            # HOMER's displayed Start is one-based while PeakID carries the BED-like interval.
            if int(row["Start"]) != start_0 + 1 or int(row["End"]) != end_0:
                raise ValueError(f"CD40L PeakID/Start/End mismatch in {path}:{source_row}")
            normalized.append(
                {
                    "dataset_id": DATASET_ID,
                    "registry_accession": "GSE298535;PRJNA1269627",
                    "source_accession": "GSE298465;PRJNA1269500",
                    "source_file": path.name,
                    "source_row": source_row,
                    "guide_id": "CD40L",
                    "guide_sequence": "TBD",
                    "editor": "TBD",
                    "editor_crosswalk_status": "unresolved; source table does not name construct",
                    "assay_signal": "Peak Score",
                    "build": BUILD,
                    "chrom": chrom,
                    "start_0based": start_0,
                    "end_0based_exclusive": end_0,
                    "start_1based": start_0 + 1,
                    "end_1based": end_0,
                    "strand": row["Strand"],
                    "source_coordinate": row[peak_column],
                    "coordinate_basis": "PeakID is 0-based half-open; HOMER Start is 1-based",
                    "offtarget_sequence": "",
                    "replicate_1_signal": row["Peak Score"] if replicate == 1 else "",
                    "replicate_2_signal": row["Peak Score"] if replicate == 2 else "",
                    "detected_by_changeseq_be": "true",
                    "annotation": row["Annotation"],
                }
            )
    return normalized


def parse_cbe(path: Path) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for source_row, row in enumerate(read_dict_rows(path), start=2):
        start_0, end_0, start_1, end_1, basis = _interval(
            row["start"], row["end"], row["OffTarget_Sequence"]
        )
        normalized.append(
            {
                "dataset_id": DATASET_ID,
                "registry_accession": "GSE298535;PRJNA1269627",
                "source_accession": "GSE298534;PRJNA1269623",
                "source_file": path.name,
                "source_row": source_row,
                "guide_id": row["sample"],
                "guide_sequence": row["Target_Sequence"],
                "editor": "eA3A-BE3",
                "editor_crosswalk_status": "explicit in GEO series metadata",
                "assay_signal": "CHANGEseqBE reads",
                "build": BUILD,
                "chrom": row["chr"],
                "start_0based": start_0,
                "end_0based_exclusive": end_0,
                "start_1based": start_1,
                "end_1based": end_1,
                "strand": "",
                "source_coordinate": row["coordinate"],
                "coordinate_basis": basis + "; convention not explicitly declared in table",
                "offtarget_sequence": row["OffTarget_Sequence"],
                "replicate_1_signal": row["CHANGEseqBE_rep1"],
                "replicate_2_signal": row["CHANGEseqBE_rep2"],
                "detected_by_changeseq_be": str(
                    _present(row["CHANGEseqBE_rep1"]) or _present(row["CHANGEseqBE_rep2"])
                ).lower(),
                "annotation": row["Annotation"],
            }
        )
    return normalized


def parse_abe_table_1(path: Path) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for source_row, row in _xlsx_dicts(path, "MKSR"):
        start_0, end_0, start_1, end_1, basis = _interval(
            row["Start"], row["End"], row["Site_Sequence"]
        )
        normalized.append(
            {
                "dataset_id": DATASET_ID,
                "registry_accession": "GSE298535;PRJNA1269627",
                "source_accession": "GSE307740;PRJNA1327241",
                "source_file": path.name,
                "source_row": source_row,
                "guide_id": row["Target_site"],
                "guide_sequence": "TBD",
                "editor": row["Genome_Editor"],
                "editor_crosswalk_status": "explicit in source row",
                "assay_signal": "MKSR reads",
                "build": BUILD,
                "chrom": row["Chromosome"],
                "start_0based": start_0,
                "end_0based_exclusive": end_0,
                "start_1based": start_1,
                "end_1based": end_1,
                "strand": "",
                "source_coordinate": row["Genomic Coordinate"],
                "coordinate_basis": basis + "; convention not explicitly declared in table",
                "offtarget_sequence": row["Site_Sequence"],
                "replicate_1_signal": row["MKSR_rep1"],
                "replicate_2_signal": row["MKSR_rep2"],
                "detected_by_changeseq_be": str(
                    _present(row["MKSR_rep1"]) or _present(row["MKSR_rep2"])
                ).lower(),
                "annotation": "",
            }
        )
    return normalized


def parse_abe_table_3(path: Path) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for source_row, row in _xlsx_dicts(path, "ABE"):
        start_0, end_0, start_1, end_1, basis = _interval(
            row["Start"], row["End"], row["Off-target_Sequence"]
        )
        normalized.append(
            {
                "dataset_id": DATASET_ID,
                "registry_accession": "GSE298535;PRJNA1269627",
                "source_accession": "GSE307740;PRJNA1327241",
                "source_file": path.name,
                "source_row": source_row,
                "guide_id": row["Target"],
                "guide_sequence": row["Target_Sequence"],
                "editor": "TBD (source label: ABE)",
                "editor_crosswalk_status": "unresolved; generic ABE label retained",
                "assay_signal": "CHANGE-seq-BE reads",
                "build": BUILD,
                "chrom": row["Chromosome"],
                "start_0based": start_0,
                "end_0based_exclusive": end_0,
                "start_1based": start_1,
                "end_1based": end_1,
                "strand": row["Strand"],
                "source_coordinate": row["Genomic_Coordinate"],
                "coordinate_basis": basis + "; convention not explicitly declared in table",
                "offtarget_sequence": row["Off-target_Sequence"],
                "replicate_1_signal": row["CHANGEseq-BE_rep1_reads"],
                "replicate_2_signal": row["CHANGEseq-BE_rep2_reads"],
                "detected_by_changeseq_be": str(
                    _present(row["CHANGEseq-BE_rep1_reads"])
                    or _present(row["CHANGEseq-BE_rep2_reads"])
                ).lower(),
                "annotation": row["Annotation"],
            }
        )
    return normalized


def ingest(source_dir: str | Path, output_dir: str | Path) -> dict[str, object]:
    source = Path(source_dir)
    destination = Path(output_dir)
    ledger = [
        asdict(item) | {"local_sha256_verified": verify_source_file(item, source / item.filename)}
        for item in SOURCES
    ]
    # Extract derived tables on every run so an older or modified extraction cannot bypass the
    # checksum applied to the publisher archive.
    cd40l_paths = extract_cd40l_tar(source / "GSE298465_RAW.tar", source / "gse298465")
    groups = {
        "cd40l_replicate_rows": parse_cd40l(cd40l_paths),
        "cbe_combined_rows": parse_cbe(source / "GSE298534_combined.csv.gz"),
        "mksr_rows": parse_abe_table_1(source / "GSE307740_Supplementary_Table_1.xlsx"),
        "abe_combined_rows": parse_abe_table_3(source / "GSE307740_Supplementary_Table_3.xlsx"),
    }
    all_rows = [row for rows in groups.values() for row in rows]
    write_tsv(destination / "sites.tsv", FIELDS, all_rows)
    detected_rows = [row for row in all_rows if row["detected_by_changeseq_be"] == "true"]
    write_tsv(destination / "changeseq_be_detected_sites.tsv", FIELDS, detected_rows)

    cd40l_keys = {
        (row["chrom"], row["start_0based"], row["end_0based_exclusive"], row["strand"])
        for row in groups["cd40l_replicate_rows"]
    }
    counts_by_guide = Counter(str(row["guide_id"]) for row in all_rows)
    detected_by_source = Counter(str(row["source_file"]) for row in detected_rows)
    detected_by_guide_editor = Counter(
        f"{row['guide_id']} | {row['editor']}" for row in detected_rows
    )
    summary: dict[str, object] = {
        "dataset_id": DATASET_ID,
        "registry_accessions": ["GSE298535", "PRJNA1269627"],
        "normalized_rows": len(all_rows),
        "changeseq_be_detected_rows": len(detected_rows),
        "source_row_counts": {name: len(rows) for name, rows in groups.items()},
        "rows_by_guide": dict(sorted(counts_by_guide.items())),
        "detected_rows_by_source_file": dict(sorted(detected_by_source.items())),
        "detected_rows_by_guide_editor": dict(sorted(detected_by_guide_editor.items())),
        "cd40l_unique_sites": len(cd40l_keys),
        "cd40l_publication_reported_sites": 81,
        "cd40l_row_reconciliation": "pass" if len(cd40l_keys) == 81 else "fail",
        "coordinate_handling": {
            "build": BUILD,
            "cd40l": "PeakID normalized as BED-like 0-based half-open and cross-checked against HOMER Start/End",
            "combined_tables": "source start/end retained as 0-based half-open with explicit convention warning",
            "liftover": "not performed",
        },
        "guide_editor_crosswalk": "partial; CD40L and generic ABE entries remain TBD",
        "scope_guard": "measurement/source ingestion only; no prediction records loaded",
    }
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    accession_ledger = {
        "dataset_id": DATASET_ID,
        "registry_records": list(REGISTRY_RECORDS),
        "source_files": ledger,
    }
    (destination / "accession_ledger.json").write_text(
        json.dumps(accession_ledger, indent=2) + "\n"
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", default="data/inputs/d1_changeseq_be")
    parser.add_argument("--output-dir", default="data/outputs/ingestion/d1")
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args(argv)
    if args.acquire:
        acquire_sources(args.source_dir)
    ingest(args.source_dir, args.output_dir)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
