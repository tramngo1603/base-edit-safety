"""Fetches reference sequence for a genomic window, and caches it.

The pipeline needs the surrounding sequence to build guides and to set reference alleles. Sequence
comes from UCSC on the first run and lands in data/inputs/ (gitignored), so later runs stay offline.

UCSC hg38 is the GRCh38 primary assembly for chr11. The fixed positions in
docs/coordinates.md are primary-chromosome coordinates; this module does not use patch-contig
coordinates.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import urllib.request
from dataclasses import dataclass

_UCSC_SEQUENCE_API = "https://api.genome.ucsc.edu/getData/sequence"
_CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data",
    "inputs",
    "reference_cache",
)


@dataclass(frozen=True)
class ReferenceWindow:
    """Plus-strand reference sequence for chrom:[start_1based, end_1based] inclusive."""

    genome_build: str
    chrom: str
    start_1based: int
    end_1based: int
    plus_sequence: str
    fetched_utc: str | None = None


def ucsc_chromosome(chrom: str) -> str:
    """Return a chromosome name with the UCSC chr prefix."""
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def fetch_window(
    chrom: str,
    start_1based: int,
    end_1based: int,
    genome_build: str = "hg38",
    use_cache: bool = True,
    timeout: int = 30,
) -> ReferenceWindow:
    """Get the plus-strand sequence for an inclusive window, from the cache or from UCSC."""
    key = f"{genome_build}_{ucsc_chromosome(chrom)}_{start_1based}_{end_1based}.json"
    cache_path = os.path.join(_CACHE_DIR, key)
    if use_cache and os.path.isfile(cache_path):
        with open(cache_path) as handle:
            payload = json.load(handle)
        return ReferenceWindow(**payload)

    start0 = start_1based - 1  # UCSC API uses 0-based, half-open coordinates
    url = (
        f"{_UCSC_SEQUENCE_API}?genome={genome_build};chrom={ucsc_chromosome(chrom)}"
        f";start={start0};end={end_1based}"
    )
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = json.load(response)
    seq = data["dna"].upper()
    window = ReferenceWindow(
        genome_build=genome_build,
        chrom=ucsc_chromosome(chrom),
        start_1based=start_1based,
        end_1based=end_1based,
        plus_sequence=seq,
        fetched_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    )
    if use_cache:
        os.makedirs(_CACHE_DIR, exist_ok=True)
        with open(cache_path, "w") as handle:
            json.dump(window.__dict__, handle)
    return window


def sequence_checksum(seq: str) -> str:
    """Short sha256 of a sequence, so an output can say which sequence it used."""
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def window_provenance(chrom: str, start_1based: int, end_1based: int,
                      genome_build: str = "hg38") -> dict:
    """What is in the cache for a window, without going to the network.

    Returns the cache key, whether it is cached, and if so its fetch time, sequence checksum, and
    length. Legacy cache entries honestly return ``fetched_utc=None`` because their original fetch
    time was not recorded and cannot be reconstructed.
    """
    key = f"{genome_build}_{ucsc_chromosome(chrom)}_{start_1based}_{end_1based}.json"
    path = os.path.join(_CACHE_DIR, key)
    if not os.path.isfile(path):
        return {"cache_key": key, "cached": False}
    with open(path) as handle:
        payload = json.load(handle)
    seq = payload.get("plus_sequence", "")
    return {"cache_key": key, "cached": True,
            "fetched_utc": payload.get("fetched_utc"),
            "checksum_sha256_16": sequence_checksum(seq), "length": len(seq)}
