"""Recount the separate rs114518452 public-genotype example from a bounded VCF."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from importlib.resources import files
from pathlib import Path
import re
import shlex
import sys
import tempfile

from .common import PopulationError, checksum, digest, json_text, parse_json, record, write_bundle


def protocol():
    return parse_json(files(__package__).joinpath("protocols", "rs114518452-nygc-20220422-v1.json").read_text())


def source_provenance(lines, source_sha256, specification, acquisition=None):
    """Separate assembly declarations from source identity established by acquisition."""
    declarations = []
    for line in lines:
        if line.startswith("##reference="):
            declarations.append(line.split("=", 1)[1].strip().strip('"'))
        elif line.startswith("##contig=<"):
            if not line.endswith(">"):
                raise PopulationError("Malformed carrier-source contig header")
            tokens = shlex.shlex(line[len("##contig=<"):-1], posix=True)
            tokens.whitespace, tokens.whitespace_split, tokens.commenters = ",", True, ""
            pairs = [token.strip().split("=", 1) for token in tokens]
            if any(len(pair) != 2 for pair in pairs) or len({pair[0] for pair in pairs}) != len(pairs):
                raise PopulationError("Malformed or duplicate carrier-source contig metadata")
            attributes = dict(pairs)
            if attributes.get("ID") == specification["contig"]:
                if "assembly" in attributes:
                    declarations.append(attributes["assembly"])
                if "length" in attributes and int(attributes["length"]) != specification["contig_length"]:
                    raise PopulationError("Carrier-source contig length contradicts the required assembly")
    for value in declarations:
        # Recognize assembly names in plain declarations and reference filenames/URLs.
        names = re.findall(r"(?<![A-Za-z0-9])(?:GRCh\d+(?:\.p\d+)?|hg\d+)(?![A-Za-z0-9])", value, re.I)
        if not names or any(name.lower() != "hg38" and re.fullmatch(r"grch38(?:\.p\d+)?", name.lower()) is None for name in names):
            raise PopulationError("Carrier-source reference is incompatible or unresolved; GRCh38 is required")
    recorded = source_sha256 == specification["recorded_acquisition"]["source_subset_sha256"]
    direct = acquisition is not None and acquisition.get("mode") == "direct_public_acquisition" \
        and acquisition.get("source_vcf_url") == specification["source_vcf_url"] \
        and acquisition.get("source_subset_sha256") == source_sha256
    verified = recorded or direct
    if not declarations and not verified:
        raise PopulationError("Carrier-source assembly is missing and no recorded acquisition contract matches these bytes; a caller-supplied checksum is not source verification")
    status = "direct_public_acquisition" if direct else "recorded_public_acquisition" if recorded else "unverified_local_source"
    return {"assembly": specification["assembly"],
        "assembly_basis": "source_documentation_and_acquisition" if verified else "vcf_declaration",
        "source_verification": {"status": status,
            "release": specification["release"] if verified else None,
            "source_vcf_url": specification["source_vcf_url"] if verified else None,
            "assembly_evidence_url": specification["assembly_evidence_url"] if verified else None},
        "expected_provenance": {key: specification[key] for key in ("assembly", "release", "source_vcf_url")}}


def recount(source, specification, *, acquisition=None):
    """Count complete diploid calls; preserve missing and filtered-call exclusions."""
    lines = Path(source).read_text().splitlines()
    provenance = source_provenance(lines, checksum(source), specification, acquisition)
    samples, matches = None, []
    for line in lines:
        if line.startswith("#CHROM\t"):
            if samples is not None:
                raise PopulationError("Duplicate carrier-source sample header")
            samples = line.split("\t")[9:]
        elif not line.startswith("#") and line:
            row = line.split("\t")
            if len(row) < 10:
                raise PopulationError("Carrier source requires sample genotype columns")
            if row[0] == specification["contig"] and int(row[1]) == specification["position_1based"] \
                    and row[3] == specification["ref"] and specification["alt"] in row[4].split(","):
                matches.append(row)
    if not samples or len(samples) != len(set(samples)) or len(matches) != 1:
        raise PopulationError("Carrier source needs unique sample IDs and exactly one matching REF/ALT record")
    row, = matches
    if len(row) != len(samples)+9 or row[6] not in specification["accepted_site_filters"]:
        raise PopulationError("Carrier source has inconsistent sample columns or an excluded site FILTER")
    alts, keys = row[4].split(","), row[8].split(":")
    if len(set(alts)) != len(alts) or any(not allele or allele == "." or allele == row[3] for allele in alts):
        raise PopulationError("Carrier source contains duplicate or invalid ALT entries")
    if keys.count("GT") != 1 or len(set(keys)) != len(keys) or any(not key for key in keys):
        raise PopulationError("Carrier source requires one unambiguous GT field")
    selected_alt = alts.index(specification["alt"])+1
    ac, an, carriers, excluded = 0, 0, [], Counter()
    for sample, value in zip(samples, row[9:]):
        values = value.split(":")
        if len(values) > len(keys):
            raise PopulationError("Carrier source contains excess sample FORMAT values")
        call = dict(zip(keys, values))
        if call.get("FT", ".") not in specification["accepted_sample_filters"]:
            excluded["sample_filter"] += 1
            continue
        gt = call.get("GT", ".")
        if not re.fullmatch(r"\d+[|/]\d+", gt):
            excluded["missing_partial_or_non_diploid_gt"] += 1
            continue
        indexes = [int(a) for a in re.split(r"[|/]", gt)]
        if any(i > len(alts) for i in indexes):
            raise PopulationError("Carrier source contains an out-of-range GT allele index")
        copies = indexes.count(selected_alt)
        ac += copies
        an += 2
        if copies:
            carriers.append({"source_sample_id": sample, "original_gt": gt, "alt_copies": copies})
    summary = {**provenance, "position_1based": specification["position_1based"], "contig": specification["contig"],
        "ref": specification["ref"], "alt": specification["alt"], "site_filter_raw": row[6],
        "source_samples": len(samples), "AC": ac, "AN": an, "carrier_donors": len(carriers),
        "excluded_calls": dict(sorted(excluded.items())), "quality_policy": specification["quality_policy"]}
    return summary, carriers


def acquire(output, specification):
    """Fetch only the exact allele and retain its complete original sample row."""
    try:
        import pysam
    except ImportError as exc:
        raise PopulationError("Install the genomics extra for explicit public acquisition") from exc
    from .public_data import retrieve
    index = output/"source.tbi"
    retrieve(specification["source_vcf_url"]+".tbi", index)
    target = output/"source.vcf"
    count = 0
    with pysam.VariantFile(specification["source_vcf_url"], index_filename=str(index)) as source:
        with pysam.VariantFile(str(target), "w", header=source.header) as dest:
            for row in source.fetch(specification["contig"], specification["position_1based"]-1, specification["position_1based"]):
                if row.pos == specification["position_1based"] and row.ref == specification["ref"] \
                        and specification["alt"] in (row.alts or ()):
                    dest.write(row)
                    count += 1
    if count != 1:
        raise PopulationError("Public query did not return exactly one matching REF/ALT record")
    return target, {"mode": "direct_public_acquisition", "source_vcf_url": specification["source_vcf_url"],
                    "source_subset_sha256": checksum(target), "pysam_version": pysam.__version__,
                    "source_index_sha256": checksum(index)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-vcf", type=Path, help="Use a previously saved bounded, plain VCF without network access")
    parser.add_argument("--expected-sha256", help="Optional checksum required for the supplied VCF")
    args = parser.parse_args(argv)
    try:
        specification = protocol()
        if args.expected_sha256 and not args.source_vcf:
            raise PopulationError("An expected checksum requires --source-vcf")
        with tempfile.TemporaryDirectory(prefix="bes-known-case-") as tmp:
            source, acquisition = (args.source_vcf, {"mode": "local_source"}) if args.source_vcf else acquire(Path(tmp), specification)
            source_sha = checksum(source)
            if args.expected_sha256 and source_sha != args.expected_sha256:
                raise PopulationError("Carrier source checksum differs from the expected subset")
            summary, carriers = recount(source, specification, acquisition=acquisition)
            rid = "known_case:" + digest({"protocol": specification, "source_sha256": source_sha,
                                         "source_verification": summary["source_verification"]})
            outputs = {"source.vcf": source.read_bytes().decode("utf-8"), "protocol.json": json_text(specification),
                "summary.json": json_text(record("known_case_genotype_check", rid, **summary,
                    source_subset_sha256=source_sha, protocol_sha256=hashlib.sha256(json_text(specification).encode()).hexdigest(), acquisition=acquisition,
                    interpretation=specification["interpretation"])),
                "protected-carriers.json": json_text(record("known_case_carriers", rid, rows=carriers))}
            write_bundle(args.output_dir, outputs, inputs=[source])
        print(f"Input genotypes: {summary['AC']} ALT copies / {summary['AN']} called alleles; {summary['carrier_donors']} carriers. Source: {summary['source_verification']['status']}. No editing outcome established.")
        return 0
    except (PopulationError, OSError, ValueError, KeyError) as exc:
        print("Carrier recount failed: " + (str(exc) if isinstance(exc, PopulationError) else type(exc).__name__), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
