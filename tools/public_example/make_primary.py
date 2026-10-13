#!/usr/bin/env python3
"""Stream-extract chr1..chr22,chrX,chrY from hg38.fa.gz into out/primary.fa.

- Retains full contigs, uppercase sequence, no other contigs.
- Writes .fai index via pysam.faidx.
- Records SHA256 of input and outputs into preparation-receipt.json (merge).
- Never prints sequence content; only aggregate stats.
"""
import gzip
import os
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--workdir", type=Path, required=True)
ROOT = parser.parse_args().workdir.resolve()
SRC = ROOT / "data" / "ucsc" / "hg38.fa.gz"
DST = ROOT / "out" / "primary.fa"
RECEIPT = ROOT / "preparation-receipt.json"

WANTED = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
LINE_WIDTH = 60


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def update_receipt(key: str, value) -> None:
    r = json.loads(RECEIPT.read_text()) if RECEIPT.exists() else {}
    r[key] = value
    r["updated_utc"] = datetime.now(timezone.utc).isoformat()
    RECEIPT.write_text(json.dumps(r, indent=2, sort_keys=True) + "\n")


def main() -> int:
    if not SRC.exists():
        print(f"FAIL missing input {SRC}")
        return 1
    DST.parent.mkdir(parents=True, exist_ok=True)
    in_hash = sha256_file(SRC)

    kept, skipped = {}, {}
    out_h = hashlib.sha256()
    with gzip.open(SRC, "rt", encoding="ascii") as fin, open(DST, "wb") as fout:

        def emit(text: str) -> None:
            out_h.update(text.encode("ascii"))
            fout.write(text.encode("ascii"))

        keep = False
        name, bases = None, 0
        for line in fin:
            if line.startswith(">"):
                if name is not None:
                    (kept if keep else skipped)[name] = bases
                name = line[1:].split()[0]
                bases = 0
                keep = name in WANTED
                if keep:
                    emit(f">{name}\n")
            else:
                seq = line.strip()
                bases += len(seq)
                if keep:
                    emit(seq.upper() + "\n")
        if name is not None:
            (kept if keep else skipped)[name] = bases

    fai_dst = Path(str(DST) + ".fai")
    import pysam  # noqa: E402
    pysam.faidx(str(DST))

    result = {
        "transform": "gzip -dc hg38.fa.gz | keep chr1..chr22,chrX,chrY full contigs; uppercase; drop all other contigs",
        "input": {"path": str(SRC.relative_to(ROOT)), "sha256": in_hash},
        "output": {
            "path": str(DST.relative_to(ROOT)),
            "sha256": out_h.hexdigest(),
            "bytes": DST.stat().st_size,
            "fai_path": fai_dst.name,
            "fai_sha256": sha256_file(fai_dst),
        },
        "contigs_kept": len(kept),
        "total_bases_kept": sum(kept.values()),
        "contigs_dropped": len(skipped),
        "kept_names": sorted(kept, key=lambda c: (len(c), c)),
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    if set(kept) != WANTED:
        print(f"FAIL contig set mismatch: missing={sorted(WANTED - set(kept))} extra={sorted(set(kept) - WANTED)}")
        return 1
    update_receipt("primary_fa", result)
    print(f"OK primary.fa contigs={len(kept)} bases={result['total_bases_kept']} "
          f"sha256={out_h.hexdigest()} dropped_contigs={len(skipped)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
