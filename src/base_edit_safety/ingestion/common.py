"""Shared primitives for validation-data ingestion.

This package deliberately has no dependency on predictor adapters, scoring, or validation
metrics.  Its outputs are measurements and provenance records only.
"""
from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
import csv
import gzip
import hashlib
import io
import os
from pathlib import Path
import tempfile
from typing import BinaryIO, TextIO
from urllib.request import urlopen


@dataclass(frozen=True)
class SourceFile:
    """One immutable acquisition-ledger entry."""

    dataset_id: str
    accession: str
    filename: str
    url: str
    sha256: str
    accessed_on: str


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_source_file(source: SourceFile, path: str | Path) -> str:
    """Return the local digest only when it matches the pinned acquisition record."""
    observed = sha256_file(path)
    if observed != source.sha256:
        raise ValueError(
            f"checksum mismatch for {Path(path)}: {observed} != {source.sha256}"
        )
    return observed


@contextmanager
def open_text(path: str | Path) -> Iterator[TextIO]:
    """Open plain text or gzip-compressed text with stable newline handling."""
    source = Path(path)
    if source.suffix == ".gz":
        with gzip.open(source, "rt", encoding="utf-8-sig", newline="") as handle:
            yield handle
    else:
        with source.open("r", encoding="utf-8-sig", newline="") as handle:
            yield handle


def read_dict_rows(path: str | Path, *, delimiter: str = ",") -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"source has no header: {path}")
        for row in reader:
            yield {str(key): "" if value is None else value for key, value in row.items()}


def write_tsv(
    path: str | Path,
    fieldnames: Iterable[str],
    rows: Iterable[Mapping[str, object]],
) -> int:
    """Write a deterministic UTF-8 TSV and return the written data-row count."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    columns = tuple(fieldnames)
    count = 0
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def download_source(
    source: SourceFile,
    destination_dir: str | Path,
    *,
    opener=urlopen,
) -> Path:
    """Acquire and checksum one raw source.

    ``opener`` is injectable so tests use local byte fixtures and never access the network.
    An existing checksum-matching file is reused.  A partial download is never installed.
    """
    destination = Path(destination_dir) / source.filename
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        verify_source_file(source, destination)
        return destination

    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", prefix=f".{source.filename}.", dir=destination.parent, delete=False
        ) as output:
            temporary = Path(output.name)
            response: BinaryIO = opener(source.url)
            try:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
            finally:
                response.close()
        observed = sha256_file(temporary)
        if observed != source.sha256:
            raise ValueError(
                f"checksum mismatch for downloaded {source.filename}: {observed} != {source.sha256}"
            )
        os.replace(temporary, destination)
        temporary = None
        return destination
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def bytes_opener(payload: bytes):
    """Return a network-free opener useful in fixture tests."""

    def _open(_url: str) -> BinaryIO:
        return io.BytesIO(payload)

    return _open
