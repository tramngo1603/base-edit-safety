"""Bounded reference access and donor-linked normalized VCF ingestion."""
from __future__ import annotations

from collections import Counter
import gzip
import hashlib
from pathlib import Path
import re

from .common import (PopulationError, ResourceLimit, checksum, fields, identity, integer,
                     nonempty, read_json, read_tsv, record, unique)

DONOR_COLUMNS = "donor_id genotype_sample_id source_kind family_id specimen_id specimen_context_id availability use_permitted activation_cost simulation_eligible sampling_stratum".split()
CALLABLE_COLUMNS = "donor_id contig start_0based end_0based_exclusive evidence_type source_id".split()


def boolean(value):
    if value not in {"true", "false"}:
        raise PopulationError("TSV boolean must be true or false")
    return value == "true"


def canonical_donors(rows):
    """Resolve one genotype record per person while retaining separate specimen rows."""
    shared = ("genotype_sample_id", "source_kind", "family_id", "activation_cost", "sampling_stratum")
    people = {}
    for row in rows:
        previous = people.setdefault(row["donor_id"], row)
        if any(previous.get(k) != row.get(k) for k in shared):
            raise PopulationError("Specimens for one donor disagree on genotype mapping or donor metadata")
    result = [people[k] for k in sorted(people)]
    unique(result, "genotype_sample_id")
    return result


def donors_from_tsv(path, basis):
    rows = read_tsv(path, DONOR_COLUMNS)
    for row in rows:
        for key in ("donor_id", "genotype_sample_id", "specimen_id", "specimen_context_id"):
            nonempty(row[key], key)
        for key in ("family_id", "sampling_stratum"):
            row[key] = row[key] or None
        if row["availability"] not in {"verified_available", "unavailable", "unverified"}:
            raise PopulationError("Unknown specimen availability")
        if row["source_kind"] not in {"synthetic", "public_genome", "research_sample"}:
            raise PopulationError("Unknown donor source kind")
        if (row["source_kind"] == "synthetic") != (basis == "synthetic_truth"):
            raise PopulationError("Synthetic and real donor sources cannot share a study")
        row["use_permitted"] = boolean(row["use_permitted"])
        row["simulation_eligible"] = boolean(row["simulation_eligible"])
        try:
            row["activation_cost"] = integer(int(row["activation_cost"]), "activation cost")
        except ValueError as exc:
            raise PopulationError("Activation cost must be a nonnegative integer") from exc
    canonical_donors(rows)
    unique(rows, "specimen_id")
    return sorted(rows, key=lambda r: (r["donor_id"], r["specimen_id"]))


class Reference:
    """Fetch global coordinates from an indexed FASTA or checksum-bound windows."""

    def __init__(self, path, policy, synthetic=False):
        self.path = Path(path)
        self.sha256 = checksum(path)
        self.aliases = policy["contig_alias_map"]
        self.windows = {}
        self.index = {}
        self.index_sha256 = None
        self.build = policy["genome_build"]
        if self.path.suffix == ".json":
            data = read_json(path)
            fields(data, {"record_type", "schema_version", "genome_build", "windows"}, label="reference windows")
            if data["record_type"] != "reference_windows" or data["schema_version"] != "1" or data["genome_build"] != self.build:
                raise PopulationError("Reference-window build or schema mismatch")
            for row in data["windows"]:
                fields(row, {"contig", "start", "end", "sequence", "sequence_sha256"}, label="reference window")
                integer(row["start"], "reference start")
                integer(row["end"], "reference end", row["start"] + 1)
                seq = row["sequence"]
                if not isinstance(seq, str) or len(seq) != row["end"] - row["start"] or set(seq) - set("ACGTN") \
                        or hashlib.sha256(seq.encode("ascii")).hexdigest() != row["sequence_sha256"]:
                    raise PopulationError("Reference-window length, alphabet or checksum mismatch")
                name = row["contig"]
                for start, end, _ in self.windows.get(name, []):
                    if start < row["end"] and row["start"] < end:
                        raise PopulationError("Reference windows overlap")
                self.windows.setdefault(name, []).append((row["start"], row["end"], seq))
        elif Path(str(path) + ".fai").is_file():
            index_path = Path(str(path) + ".fai")
            self.index_sha256 = checksum(index_path)
            for line in index_path.read_text().splitlines():
                parts = line.split("\t")
                if len(parts) != 5:
                    raise PopulationError("Unsupported FASTA index row")
                name = parts[0]
                if name in self.index:
                    raise PopulationError("Duplicate FASTA contig")
                length, offset, bases, width = map(int, parts[1:])
                if length < 1 or offset < 0 or bases < 1 or width < bases:
                    raise PopulationError("Invalid FASTA index geometry")
                self.index[name] = (length, offset, bases, width)
        elif synthetic and self.path.stat().st_size <= 1000000:
            name, parts = None, []
            for line in self.path.read_text().splitlines() + [">"]:
                if line.startswith(">"):
                    if name is not None:
                        seq = "".join(parts).upper()
                        if name in self.windows or set(seq) - set("ACGTN"):
                            raise PopulationError("Invalid synthetic FASTA")
                        self.windows[name] = [(0, len(seq), seq)]
                    name, parts = (line[1:].split() or [None])[0], []
                else:
                    if name is None:
                        raise PopulationError("FASTA sequence precedes a header")
                    parts.append(line.strip())
        else:
            raise PopulationError("Reference needs a FASTA .fai index or explicit global reference windows")
        available = set(self.windows) | set(self.index)
        if not set(policy["permitted_contigs"]) <= available:
            raise PopulationError("Permitted contigs do not resolve in the reference")
        for alias, name in self.aliases.items():
            if alias != name and alias in available and name in available:
                raise PopulationError("Reference contigs collide under the alias map")

    def fetch(self, contig, start, end):
        if start < 0 or end <= start:
            raise PopulationError("Invalid reference interval")
        if contig in self.windows:
            for lo, hi, seq in self.windows[contig]:
                if lo <= start and end <= hi:
                    return seq[start-lo:end-lo]
            raise PopulationError("Reference window lacks required global coordinates")
        if contig not in self.index:
            raise PopulationError("Reference contig does not resolve")
        length, offset, bases, width = self.index[contig]
        if end > length:
            raise PopulationError("Reference interval exceeds contig")
        first = offset + start // bases * width + start % bases
        last = offset + (end-1) // bases * width + (end-1) % bases
        with self.path.open("rb") as handle:
            handle.seek(first)
            seq = handle.read(last-first+1).replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(seq) != end-start or set(seq) - set("ACGTN"):
            raise PopulationError("FASTA and its index disagree")
        return seq

    def bounds(self, contig):
        if contig in self.windows:
            return [[lo, hi] for lo, hi, _ in sorted(self.windows[contig])]
        return [[0, self.index[contig][0]]]


def callable_intervals(path, donors, reference, source):
    if path is None:
        return []
    rows = read_tsv(path, CALLABLE_COLUMNS)
    for row in rows:
        if row["donor_id"] not in {d["donor_id"] for d in donors}:
            raise PopulationError("Callability donor does not resolve")
        row["start_0based"] = integer(int(row["start_0based"]), "callability start")
        row["end_0based_exclusive"] = integer(int(row["end_0based_exclusive"]), "callability end", row["start_0based"]+1)
        if row["evidence_type"] not in {"explicit_reference_calls", "callable_mask_with_source_contract"} \
                or row["source_id"] != source["source_id"]:
            raise PopulationError("Unresolved reference-callability source contract")
        reference.fetch(row["contig"], row["start_0based"], row["end_0based_exclusive"])
    rows.sort(key=lambda r: (r["donor_id"], r["contig"], r["start_0based"]))
    for a, b in zip(rows, rows[1:]):
        if (a["donor_id"], a["contig"]) == (b["donor_id"], b["contig"]) and a["end_0based_exclusive"] > b["start_0based"]:
            raise PopulationError("Overlapping callability rows require upstream reconciliation")
    return rows


def parse_call(text, formats, site_filter, policy, source, donor, contig, *, expected_ploidy=2):
    values = text.split(":")
    if len(values) > len(formats):
        raise PopulationError("VCF sample has excess FORMAT fields")
    if any(not value for value in values):
        raise PopulationError("VCF sample has an empty FORMAT value; use '.' for missing data")
    raw = dict(zip(formats, values))
    gt = raw.get("GT", ".")
    alleles = re.split(r"[/|]", gt)
    status, indexes = "complete", None
    if expected_ploidy not in (1, 2):
        raise PopulationError("Genotype parsing requires explicit haploid or diploid ploidy")
    if len(alleles) != expected_ploidy:
        status = "missing" if gt == "." else "unsupported_ploidy"
    elif all(a == "." for a in alleles):
        status = "missing"
    elif "." in alleles:
        status = "partial_missing"
    elif any(not re.fullmatch(r"[0-9]+", a) for a in alleles) or ("/" in gt and "|" in gt):
        status = "invalid"
    else:
        indexes = [int(a) for a in alleles]
    reasons = []
    if site_filter not in policy["accepted_site_filters"]:
        reasons.append("site_filter_rejected")
    if raw.get("FT", ".") not in policy["accepted_sample_filters"]:
        reasons.append("sample_filter_rejected")
    for field, key in (("GQ", "min_gq"), ("DP", "min_dp")):
        if policy[key] is not None:
            value = raw.get(field, ".")
            if not re.fullmatch(r"[0-9]+", value) or int(value) < policy[key]:
                reasons.append(field.lower() + "_missing_or_failed")
    if reasons:
        status = "filtered"
    phase_key, phase_status = None, "unsupported"
    if status == "complete":
        if len(indexes) == 1:
            phase_status = "haploid"
        elif indexes[0] == indexes[1]:
            phase_status = "homozygous"
        elif "|" in gt and raw.get("PS", ".") != ".":
            phase_status = "phased_explicit_ps"
            phase_key = [source, donor, contig, "explicit", raw["PS"]]
        elif "|" in gt and policy["phase"] == "vcf_implicit_phase_set":
            phase_status = "phased_implicit_ps"
            phase_key = [source, donor, contig, "implicit"]
        else:
            phase_status = "unresolved_unphased"
    return {"gt": gt, "selected_indexes": indexes, "separator": "|" if "|" in gt else "/" if "/" in gt else None,
            "ps": raw.get("PS"), "quality": {k: raw.get(k) for k in ("GQ", "DP", "FT")},
            "site_filter": site_filter, "call_status": status, "phase_status": phase_status,
            "phase_key": phase_key, "reasons": reasons if reasons else ([] if status == "complete" else [status])}


def allele_edit(pos, ref, alt, limit):
    """Map a minimal literal allele without splitting its original GT index."""
    if not alt or set(alt) - set("ACGT") or max(len(ref), len(alt)) > limit:
        return None
    if len(ref) == len(alt):
        return {"start": pos, "end": pos+len(ref), "sequence": alt, "kind": "substitution"}
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    r, a = ref[prefix:], alt[prefix:]
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    if r and a:
        return None
    return {"start": pos+prefix, "end": pos+prefix+len(r), "sequence": a,
            "kind": "insertion" if not r else "deletion"}


def load_vcf(path, config, reference, donors, run_id):
    donors = canonical_donors(donors)
    source = config["input_sources"]["variants"]["source_id"]
    policy = config["genotype_policy"]
    scope = config["scope"]["materialized_intervals"]
    aliases = config["reference_policy"]["contig_alias_map"]
    records, audit, samples, duplicates = [], [], None, {}
    seen_contigs, raw_count, omitted, build_seen = {}, 0, 0, False
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##reference="):
                if line.strip().split("=", 1)[1] != reference.build:
                    raise PopulationError("VCF reference declaration must match the configured assembly")
                build_seen = True
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM\t"):
                header = line.rstrip("\r\n").split("\t")
                if header[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
                    raise PopulationError("Invalid VCF header")
                samples = header[9:]
                if len(samples) != len(set(samples)) or not {d["genotype_sample_id"] for d in donors} <= set(samples):
                    raise PopulationError("VCF sample mapping is missing or ambiguous")
                continue
            if not line.strip():
                continue
            if samples is None:
                raise PopulationError("VCF header is missing")
            raw_count += 1
            if raw_count > config["resource_limits"]["max_variant_rows"]:
                raise ResourceLimit("max_variant_rows exceeded before complete ingestion")
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != 9 + len(samples):
                raise PopulationError("VCF row width does not match header")
            chrom, pos1, name, ref, alt_text, qual, filt, info_text, fmt = parts[:9]
            if chrom not in aliases:
                raise PopulationError("VCF contig is absent from the explicit alias map")
            contig = aliases[chrom]
            if contig in seen_contigs and seen_contigs[contig] != chrom:
                raise PopulationError("VCF contig alias collision")
            seen_contigs[contig] = chrom
            pos = integer(int(pos1), "VCF POS", 1) - 1
            if not ref or set(ref) - set("ACGTN"):
                raise PopulationError("VCF REF is not a literal reference allele")
            alts = alt_text.split(",")
            if any(not alt for alt in alts):
                raise PopulationError("Empty VCF ALT index; literal alleles require a reference anchor")
            if len(alts) != len(set(alts)):
                raise PopulationError("Duplicate original ALT indexes")
            formats = fmt.split(":")
            if any(not key for key in formats) or len(set(formats)) != len(formats):
                raise PopulationError("Empty or duplicate VCF FORMAT key")
            info = {}
            for item in info_text.split(";"):
                if item == ".":
                    continue
                key, value = item.split("=", 1) if "=" in item else (item, True)
                if not key or key in info:
                    raise PopulationError("Empty or duplicate VCF INFO key")
                info[key] = value
            if "END" in info and (not isinstance(info["END"], str) or not re.fullmatch(r"[0-9]+", info["END"])
                                  or int(info["END"]) < pos+len(ref)):
                raise PopulationError("VCF END must bound the complete reference allele")
            symbolic = any(set(a) - set("ACGT") for a in alts)
            unbounded = symbolic and "END" not in info
            end = max(pos + len(ref), int(info.get("END", pos + len(ref))))
            relevant = list(scope) if unbounded else [s for s in scope if s["contig"] == contig and s["end"] > pos and s["start"] < end]
            if not relevant:
                omitted += 1
                continue
            if reference.fetch(contig, pos, pos+len(ref)) != ref:
                raise PopulationError("VCF REF disagrees with pinned reference")
            edits = [None] + [allele_edit(pos, ref, a, config["variant_policy"]["max_allele_length"]) for a in alts]
            if len(alts) > 1 and not symbolic:
                # Normalization preserves this record's ALT indexes. One ALT can
                # shift alone even when the complete multiallelic record cannot.
                alleles = [ref] + alts
                has_indel = any(e and e["kind"] in {"insertion", "deletion"} for e in edits[1:])
                common_suffix = len({a[-1] for a in alleles}) == 1
                common_prefix = len({a[0] for a in alleles}) == 1
                can_trim = min(map(len, alleles)) > 1
                if has_indel and ((common_suffix and (pos > 0 or can_trim)) or (common_prefix and can_trim)):
                    raise PopulationError("Indel record is not jointly normalized; normalize against the pinned reference")
            else:
                for edit in edits[1:]:
                    if edit and edit["kind"] in {"insertion", "deletion"}:
                        changed = edit["sequence"] if edit["kind"] == "insertion" else reference.fetch(contig, edit["start"], edit["end"])
                        previous = reference.fetch(contig, edit["start"]-1, edit["start"]) if edit["start"] > pos else None
                        if changed and previous == changed[-1]:
                            raise PopulationError("Indel is not left normalized; normalize against the pinned reference")
            vid = identity("var", [reference.sha256, contig, pos, ref, alts])
            calls = {}
            for donor in donors:
                call = parse_call(parts[9+samples.index(donor["genotype_sample_id"])] , formats, filt,
                                  policy, source, donor["donor_id"], contig)
                if call["selected_indexes"] is not None and any(i > len(alts) for i in call["selected_indexes"]):
                    call.update(call_status="invalid", reasons=["allele_index_out_of_range"])
                call["call_id"] = identity("call", [source, vid, donor["donor_id"], call])
                calls[donor["donor_id"]] = call
            if vid in duplicates:
                if duplicates[vid] != calls:
                    raise PopulationError("Discordant duplicate VCF record")
                continue
            duplicates[vid] = calls
            ac = Counter(i for call in calls.values() if call["call_status"] == "complete" for i in call["selected_indexes"])
            an = sum(ac.values())
            rec = record("normalized_variant", run_id, variant_id=vid, source_id=source,
                         source_record_index=raw_count, original_id=name, contig=contig, pos0=pos, end=end,
                         unbounded_footprint=unbounded, ref=ref, original_alts=alts, edits=edits, calls=calls,
                         normalization=config["variant_policy"]["normalizer"],
                         frequency={"AC": [ac[i] for i in range(1, len(alts)+1)], "AN": an,
                                    "AF": [ac[i]/an if an else None for i in range(1, len(alts)+1)],
                                    "evaluable_donors": an//2, "subset": "declared_donor_catalog"})
            records.append(rec)
            for donor, call in calls.items():
                audit.append(record("genotype_audit", run_id, donor_id=donor, source_id=source,
                                    source_record_index=raw_count, variant_id=vid, original_alleles=[ref]+alts, **call))
    if samples is None or not build_seen:
        raise PopulationError("VCF requires an explicit matching reference and sample header")
    return sorted(records, key=lambda r: (r["contig"], r["pos0"], r["variant_id"])), audit, {
        "input_rows": raw_count, "unique_retained_rows": len(records), "outside_scope_rows": omitted,
        "duplicate_rows": raw_count - omitted - len(records)}
