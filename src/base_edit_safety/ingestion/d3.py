"""Acquire and normalize the Lazzarotto et al. 2020 CHANGE-seq site table (D3)."""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path

from .common import SourceFile, download_source, verify_source_file, write_tsv
from .xlsx import iter_rows


DATASET_ID = "D3"
BUILD = "GRCh38/hg38"
SOURCE = SourceFile(
    dataset_id=DATASET_ID,
    accession="doi:10.1038/s41587-020-0555-7; PMC7652380; PRJNA625995",
    filename="Lazzarotto_2020_SuppTables_2_8.xlsx",
    url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-020-0555-7/"
    "MediaObjects/41587_2020_555_MOESM2_ESM.xlsx",
    sha256="f4a48a985b8616d4fa29798da153da65672170a434f1d4d845fd9163d0ee91fa",
    accessed_on="2026-08-14",
)

FIELDS = [
    "dataset_id", "registry_accession", "source_accession", "source_file", "source_sheet",
    "source_row", "guide_id", "source_target_sequence", "editor_or_nuclease", "build",
    "chrom", "start_0based", "end_0based_exclusive", "start_1based", "end_1based",
    "strand", "source_coordinate", "coordinate_basis", "offtarget_sequence",
    "changeseq_reads", "edit_distance", "product_identity_crosswalk",
]


def acquire_source(destination: str | Path) -> Path:
    return download_source(SOURCE, destination)


def _rows(path: Path, sheet: str):
    iterator = iter_rows(path, sheet)
    header = next(iterator)
    for source_row, values in enumerate(iterator, start=2):
        values += [""] * (len(header) - len(values))
        yield source_row, dict(zip(header, values))


def parse_sites(path: str | Path) -> list[dict[str, object]]:
    source_path = Path(path)
    normalized: list[dict[str, object]] = []
    for source_row, row in _rows(source_path, "CHANGE-seq_Supp_Table_3"):
        start_0 = int(row["chromStart"])
        end_0 = int(row["chromEnd"])
        if end_0 <= start_0:
            raise ValueError(f"invalid CHANGE-seq interval at workbook row {source_row}")
        expected_coordinate = f"{row['chrom']}:{start_0 + 1}-{end_0}:{row['strand']}"
        if row["chromStart:chromEnd"] != expected_coordinate:
            raise ValueError(
                f"BED/display coordinate mismatch at workbook row {source_row}: "
                f"{row['chromStart:chromEnd']} != {expected_coordinate}"
            )
        normalized.append(
            {
                "dataset_id": DATASET_ID,
                "registry_accession": "PRJNA625995;GSE149295",
                "source_accession": SOURCE.accession,
                "source_file": source_path.name,
                "source_sheet": "CHANGE-seq_Supp_Table_3",
                "source_row": source_row,
                "guide_id": row["name"],
                "source_target_sequence": row["target"],
                "editor_or_nuclease": "SpCas9",
                "build": BUILD,
                "chrom": row["chrom"],
                "start_0based": start_0,
                "end_0based_exclusive": end_0,
                "start_1based": start_0 + 1,
                "end_1based": end_0,
                "strand": row["strand"],
                "source_coordinate": row["chromStart:chromEnd"],
                "coordinate_basis": "chromStart/chromEnd are BED 0-based half-open; display field is 1-based inclusive",
                "offtarget_sequence": row["offtarget_sequence"],
                "changeseq_reads": row["CHANGEseq_reads"],
                "edit_distance": row["distance"],
                "product_identity_crosswalk": "unverified",
            }
        )
    return normalized


def ingest(source_path: str | Path, output_dir: str | Path) -> dict[str, object]:
    source = Path(source_path)
    destination = Path(output_dir)
    observed_sha256 = verify_source_file(SOURCE, source)
    sites = parse_sites(source)
    write_tsv(destination / "sites.tsv", FIELDS, sites)
    guides = Counter(str(row["guide_id"]) for row in sites)
    product_terms = {"BCL11A", "Casgevy", "exagamglogene", "exa-cel"}
    matching_terms = sorted(
        term
        for term in product_terms
        if any(term.casefold() in str(value).casefold() for row in sites for value in row.values())
    )
    summary: dict[str, object] = {
        "dataset_id": DATASET_ID,
        "normalized_site_rows": len(sites),
        "unique_guide_ids": len(guides),
        "publication_reported_site_rows": 202_043,
        "publication_reported_guides": 110,
        "row_reconciliation": (
            "pass" if len(sites) == 202_043 and len(guides) == 110 else "fail"
        ),
        "coordinate_handling": {
            "build": BUILD,
            "source": "workbook BED columns and paper hg38 methods",
            "normalization": "both 0-based half-open and 1-based inclusive coordinates emitted",
            "liftover": "not performed",
        },
        "registry_accession_caveat": (
            "GSE149295 is the ChIP-seq subseries; the normalized CHANGE-seq site list is the "
            "publisher workbook associated with PRJNA625995/PMC7652380"
        ),
        "product_identity_terms_found": matching_terms,
        "product_identity_crosswalk": (
            "unverified: the 110-guide table contains no BCL11A/Casgevy designation; no guide "
            "was relabeled or selected as a product guide"
        ),
        "scope_guard": "source measurements only; no prediction records loaded",
    }
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    ledger = asdict(SOURCE) | {"local_sha256": observed_sha256}
    (destination / "accession_ledger.json").write_text(json.dumps([ledger], indent=2) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", default="data/inputs/d3_changeseq/Lazzarotto_2020_SuppTables_2_8.xlsx"
    )
    parser.add_argument("--output-dir", default="data/outputs/ingestion/d3")
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args(argv)
    if args.acquire:
        acquire_source(Path(args.source).parent)
    ingest(args.source, args.output_dir)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
