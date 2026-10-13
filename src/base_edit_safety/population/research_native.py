"""Run and validate pinned Cas-OFFinder jobs with retained completion evidence."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import resource
import signal
import subprocess
import sys
import time

from .common import PopulationError, ResourceLimit, checksum, save_json
from ..pipeline.candidate_panel import SearchProfile
from ..pipeline.normalize_nominations import normalize_cas_offinder_2_4_1


def _rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)


def _check_limits(limits, check_limits, start):
    check_limits()
    if time.monotonic()-start > limits["timeout_seconds"]:
        raise ResourceLimit("Native job or parsing exceeded timeout_seconds")
    if _rss() > limits["max_resident_bytes"]:
        raise ResourceLimit("Parent memory exceeded max_resident_bytes")


def _process_group_rss(group):
    """Measure the timed process and its children without reading command lines."""
    result = subprocess.run(["/bin/ps", "-axo", "pgid=,rss="], capture_output=True,
                            text=True, timeout=5)
    if result.returncode:
        raise PopulationError("Cannot monitor native process memory with ps")
    total = 0
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 2:
            raise PopulationError("Unexpected process memory monitor output")
        pgid, rss = map(int, fields)
        if pgid == group:
            total += rss*1024
    return total


def _tail(path, size=65536):
    if not path.exists():
        return ""
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size-size))
        return handle.read(size).decode("utf-8", errors="replace")


def _version_in(path):
    if not path.exists():
        return False
    with path.open("rb") as handle:
        previous = b""
        while chunk := handle.read(65536):
            text = previous + chunk
            if re.search(rb"Cas-OFFinder v2\.4\.1(?:[ \r\n(]|$)", text):
                return True
            previous = text[-64:]
    return False


def _stop(proc):
    if proc is None or proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=3)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait(timeout=3)


def run_native(fasta, directory, guides, config, binary, limits, *, check_limits=lambda: None):
    """Retain every attempt; completion binds the exact raw output to a successful job."""
    directory = Path(directory)
    directory.mkdir(mode=0o700)
    raw = directory / "raw.tsv"
    request = directory / "input.txt"
    command = [str(binary), str(request), config["backend"]["device"], str(raw)]
    start = time.monotonic()
    proc = None
    audit = {"status": "failed", "command": command, "exit_code": None,
             "version": None, "version_probe_exit_code": None, "peak_resident_bytes": None,
             "sampled_process_group_peak_bytes": 0, "memory_sampling_seconds": 0.2,
             "raw_sha256": None, "raw_bytes": None}

    def check():
        _check_limits(limits, check_limits, start)
        size = sum(p.stat().st_size for p in directory.iterdir() if p.is_file())
        if size > limits["max_output_bytes"]:
            raise ResourceLimit("Native artifacts exceeded max_output_bytes")

    def execute(argv, stdout, stderr):
        nonlocal proc
        with stdout.open("wb") as out, stderr.open("wb") as err:
            proc = subprocess.Popen(argv, stdout=out, stderr=err, start_new_session=True)
            while proc.poll() is None:
                check()
                rss = _process_group_rss(proc.pid)
                audit["sampled_process_group_peak_bytes"] = max(audit["sampled_process_group_peak_bytes"], rss)
                if rss > limits["max_resident_bytes"]:
                    raise ResourceLimit("Running native process memory exceeded max_resident_bytes")
                time.sleep(0.2)
            check()
            return proc.returncode

    try:
        check()
        request.write_text(str(Path(fasta).resolve()) + "\n" + "N"*21 + "GG\n" +
                           "".join(g["spacer_sequence"] + "NNN " + str(config.get("max_mismatches", 6)) + " " + g["guide_id"] + "\n" for g in guides))
        version_out, version_err = directory/"version-stdout.log", directory/"version-stderr.log"
        audit["version_probe_command"] = [str(binary)]
        audit["version_probe_exit_code"] = execute([str(binary)], version_out, version_err)
        # Usage can fail while listing unavailable devices after printing the version.
        # Keep that diagnostic; the requested-device search must still succeed.
        if not (_version_in(version_out) or _version_in(version_err)):
            raise PopulationError("Native executable did not report Cas-OFFinder v2.4.1 in its usage probe")
        audit["version"] = "2.4.1"
        timed = ["/usr/bin/time", "-l" if sys.platform == "darwin" else "-v", *command]
        audit["timed_command"] = timed
        audit["exit_code"] = execute(timed, directory/"stdout.log", directory/"stderr.log")
        stderr = _tail(directory/"stderr.log")
        match = re.search(r"(\d+)\s+maximum resident set size", stderr) if sys.platform == "darwin" else re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
        audit["peak_resident_bytes"] = int(match[1])*(1 if sys.platform == "darwin" else 1024) if match else None
        if audit["exit_code"] or not raw.is_file():
            raise PopulationError("Native search failed or produced no output file")
        if audit["peak_resident_bytes"] is None:
            raise PopulationError("Native OS peak-memory measurement is missing")
        if audit["peak_resident_bytes"] > limits["max_resident_bytes"]:
            raise ResourceLimit("Native OS peak memory exceeded max_resident_bytes")
        audit["raw_sha256"] = checksum(raw, check_limits=check)
        audit["raw_bytes"] = raw.stat().st_size
        check()
        audit["status"] = "complete"
    except BaseException as exc:
        audit["error"] = str(exc)
        raise
    finally:
        try:
            _stop(proc)
        except (OSError, subprocess.SubprocessError) as cleanup_error:
            audit["cleanup_error"] = str(cleanup_error)
            audit["status"] = "failed"
        if "timed_command" in audit:
            audit["exit_code"] = proc.returncode if proc is not None else None
        audit["wall_seconds"] = time.monotonic()-start
        if raw.is_file():
            audit["raw_bytes"] = raw.stat().st_size
        stderr = _tail(directory/"stderr.log")
        match = re.search(r"(\d+)\s+maximum resident set size", stderr) if sys.platform == "darwin" else re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
        if match:
            audit["peak_resident_bytes"] = int(match[1])*(1 if sys.platform == "darwin" else 1024)
        save_json(directory/"execution.json", audit)
    return raw, audit


def _validate_native_row(parts, by_guide, mismatch_limit):
    """Validate seven native fields against the submitted spacer-plus-NNN queries."""
    if len(parts) != 7 or any(not field for field in parts):
        raise PopulationError("Expected seven nonempty native fields")
    query, chrom, pos, dna, strand, mm, guide_id = parts
    guide = by_guide.get(guide_id)
    if guide is None or query != guide["spacer_sequence"] + "NNN":
        raise PopulationError("Native query or guide ID differs from submitted guide")
    if any(c.isspace() for c in chrom) or strand not in {"+", "-"}:
        raise PopulationError("Invalid native chromosome or strand")
    if not re.fullmatch(r"[0-9]+", pos) or not re.fullmatch(r"[0-9]+", mm):
        raise PopulationError("Native coordinates and mismatch counts must be nonnegative integers")
    if len(dna) != 23 or set(dna.upper())-set("ACGTN"):
        raise PopulationError("Native target must contain exactly 23 ACGTN bases")
    sequence = dna.upper()
    mismatches = [i for i in range(20) if sequence[i] != query[i]]
    # The pinned writer lowercases N against a concrete spacer base, but the
    # query's wildcard PAM positions retain their original uppercase letters.
    if [i for i, base in enumerate(dna) if base.islower()] != mismatches:
        raise PopulationError("Native lowercase mismatch markers disagree with the pinned writer")
    if int(mm) != len(mismatches) or int(mm) > mismatch_limit:
        raise PopulationError("Native mismatch count disagrees with the query or exceeds its threshold")
    if sequence[21:] != "GG":
        raise PopulationError("Native target violates the submitted NGG pattern")


def parsed_rows(raw, guides, config, limits, *, check_limits=lambda: None):
    """Validate every raw row; an empty result requires a checksum-bound completion."""
    raw = Path(raw)
    start = time.monotonic()
    check = lambda: _check_limits(limits, check_limits, start)
    check()
    receipt_path = raw.with_name("execution.json")
    if not receipt_path.is_file() or receipt_path.stat().st_size > 1024*1024:
        raise PopulationError("Native parsing requires a retained execution receipt")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("status") != "complete" or receipt.get("exit_code") != 0 or receipt.get("version") != "2.4.1" \
            or receipt.get("raw_bytes") != raw.stat().st_size or receipt.get("raw_sha256") != checksum(raw, check_limits=check):
        raise PopulationError("Native parsing requires successful completion bound to the unchanged raw output")
    if raw.stat().st_size > limits["max_output_bytes"]:
        raise ResourceLimit("Native raw output exceeded max_output_bytes")
    clean = raw.with_name("concrete.tsv")
    retained_bytes = sum(p.stat().st_size for p in raw.parent.iterdir() if p.is_file() and p != clean)
    if retained_bytes > limits["max_output_bytes"]:
        raise ResourceLimit("Native artifacts exceeded max_output_bytes before parsing")
    unknown, unknown_seen = [], set()
    count, concrete_count = 0, 0
    by_guide = {g["guide_id"]: g for g in guides}
    with raw.open() as source, clean.open("w") as target:
        while True:
            check()
            line = source.readline(65537)
            if not line:
                break
            if len(line) > 65536:
                raise PopulationError("Native row exceeds the supported bounded line length")
            count += 1
            if count > limits["max_native_rows"]:
                raise ResourceLimit("Native row count exceeded max_native_rows; no truncation accepted")
            if retained_bytes + target.tell() + len(line.encode()) > limits["max_output_bytes"]:
                raise ResourceLimit("Native parsing artifacts exceeded max_output_bytes")
            parts = line.rstrip("\r\n").split("\t")
            _validate_native_row(parts, by_guide, config.get("max_mismatches", 6))
            if "N" in parts[3].upper():
                key = hashlib.sha256(line.rstrip("\r\n").encode()).digest()
                if key in unknown_seen:
                    raise PopulationError("Duplicate native unknown-sequence row")
                unknown_seen.add(key)
                unknown.append(parts)
                target.write("\n")  # Keep concrete-row line numbers aligned with raw output.
            else:
                concrete_count += 1
                target.write(line)
    profile = SearchProfile("research-ngg-ungapped", config.get("max_mismatches", 6), 0, 0,
                            ("NGG",), False, "computational_only", config["assembly"],
                            config["reference"]["sha256"], None, "all_within_bounds", None, None)
    if concrete_count:
        rows, accounting = normalize_cas_offinder_2_4_1(clean, guides, [], profile,
            "Cas-OFFinder 2.4.1", with_ids=True, check_limits=check)
    else:
        rows, accounting = [], {"rows_normalized": 0}
    check()
    accounting.update(raw_rows=count, unknown_sequence_rows=len(unknown), successful_process_required=True,
                      execution_receipt_sha256=checksum(receipt_path, check_limits=check))
    return rows, unknown, accounting
