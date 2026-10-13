"""Run pinned CRISPRme+ on retained, exact donor-copy sequence windows.

Installation is separate. Execution never pulls an image, opens the network,
or mounts the source repository. Incomplete attempts remain available for audit.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import uuid
import resource

from .adapters import crisprme_plus as native
from .common import PopulationError, ResourceLimit, canonical, checksum, digest, integer, json_text, read_json
from .comparison import HostBudget, host_limits_from_study, _configured_limits


_CODE_AT_IMPORT = {Path(p).name: checksum(p) for p in (__file__, native.__file__)}
SEQUENCE_PROJECTION = "identical-full-window-sequence-v1"


def _no_limit_check():
    pass


def _window_source_digest(windows, check_limits=_no_limit_check):
    value = hashlib.sha256()
    for window in windows:
        check_limits()
        value.update(canonical(window))
        value.update(b"\n")
    return value.hexdigest()


def _write_window_sources(path, windows, expected, check_limits=_no_limit_check):
    value = hashlib.sha256()
    with Path(path).open("wb") as handle:
        for window in windows:
            check_limits()
            line = canonical(window) + b"\n"
            value.update(line)
            handle.write(line)
    if value.hexdigest() != expected:
        raise PopulationError("Submitted native windows changed while retaining their source")


def _sequence_aliases(windows, guard_length, check_limits=_no_limit_check):
    """Search each exact string once while retaining every source window."""
    groups = {}
    for window in windows:
        check_limits()
        group = groups.setdefault(window["sequence"], [])
        group.append(window)
    aliases, projection = {}, []
    for index, group in enumerate(sorted(groups.values(), key=lambda rows: rows[0]["window_id"])):
        check_limits()
        representative = group[0]
        aliases[f"chrW{index:08d}"] = {**representative, "native_guard": {
            "version": native.WINDOW_GUARD_TRANSFORM,
            "prefix": "A" * guard_length, "suffix": "A" * guard_length}}
        projection.append({"representative_window_id": representative["window_id"],
                           "sequence_sha256": hashlib.sha256(representative["sequence"].encode()).hexdigest(),
                           "window_ids": [row["window_id"] for row in group]})
    return aliases, projection


def _project_alignments(parsed, projection, profile, check_limits=_no_limit_check):
    """Associate native string coordinates with every retained source path."""
    by_representative = {row["representative_window_id"]: row for row in projection}
    result = {key: value for key, value in parsed.items() if key != "hits"}
    result["hits"] = []
    result["reason_codes"] = list(parsed["reason_codes"])
    result["sequence_projection_transform"] = SEQUENCE_PROJECTION
    for original in parsed["hits"]:
        group = by_representative[original["window_id"]]
        for identifier in group["window_ids"]:
            check_limits()
            if len(result["hits"]) >= profile["max_rows"]:
                result["retention_status"] = "truncated"
                result["reason_codes"].append("projected_alignment_cap")
                result["projected_alignment_count"] = len(result["hits"])
                return result
            hit = deepcopy(original)
            hit["window_id"] = identifier
            hit["sequence_projection"] = {"version": SEQUENCE_PROJECTION,
                "representative_window_id": group["representative_window_id"],
                "target_window_id": identifier, "sequence_sha256": group["sequence_sha256"]}
            for ref in hit["native_row_refs"]:
                ref["note"] = ref.get("note", "Verbatim native alignment row") + (
                    "; associated with the source path by exact full-window sequence identity; "
                    "see retained sequence-projection.json")
            result["hits"].append(hit)
    result["projected_alignment_count"] = len(result["hits"])
    return result


def _verify_window_source(directory, payload):
    if payload.get("adapter_contract") != "materialized-donor-windows-v2":
        return
    source = payload.get("window_source", {})
    if source.get("filename") != "source-windows.jsonl" or type(source.get("count")) is not int \
            or source["count"] < 0 or checksum(Path(directory) / source["filename"]) != source.get("sha256"):
        raise PopulationError("Native window source differs from its effective input identity")
    projection = payload.get("sequence_projection", {})
    if projection.get("version") != SEQUENCE_PROJECTION \
            or digest(read_json(Path(directory) / "sequence-projection.json")) != projection.get("sha256"):
        raise PopulationError("Native sequence projection differs from its effective input identity")


@dataclass(frozen=True)
class NativeLimits:
    cpus: int = 2
    memory_bytes: int = 4 * 1024 ** 3
    disk_bytes: int = 1024 ** 3
    timeout_seconds: int = 600
    pids: int = 128
    max_rows: int | None = None

    def validate(self):
        integer(self.cpus, "native CPUs", 1, 4)
        integer(self.memory_bytes, "native memory bytes", 128 * 1024 ** 2, 6 * 1024 ** 3)
        integer(self.disk_bytes, "native disk bytes", 1024 ** 2)
        integer(self.timeout_seconds, "native timeout seconds", 1)
        integer(self.pids, "native process ceiling", 16, 512)
        if self.max_rows is not None:
            integer(self.max_rows, "native retained rows", 1)
        return self


def native_limits_from_study(study):
    """Resolve configured native execution ceilings; row caps only narrow retention."""
    try:
        return NativeLimits(**_configured_limits(study, "native")).validate()
    except TypeError as exc:
        raise PopulationError("Invalid search.execution_limits.native fields") from exc


# The driver is part of the cache key. It records the actual binaries before
# running either index construction or search, then records cgroup measurements.
_DRIVER = r'''
import hashlib, json, os, pathlib, resource, subprocess, sys, time
root = pathlib.Path('/DATA')
plan = json.loads((root / 'plan.json').read_text())
audit = {'commands': [], 'executable_hashes': {}, 'version': None, 'complete': False}
def save():
    audit['peak_child_rss_kib'] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    for key, path in [('cgroup_memory_peak_bytes', '/sys/fs/cgroup/memory.peak'),
                      ('cgroup_cpu_stat', '/sys/fs/cgroup/cpu.stat')]:
        try:
            value = pathlib.Path(path).read_text().strip()
            audit[key] = int(value) if value.isdigit() else value
        except OSError:
            audit[key] = None
    (root / 'native-audit.json').write_text(json.dumps(audit, indent=2, sort_keys=True)+'\n')
exit_code = 1
try:
    for path, expected in plan['executable_hashes'].items():
        actual = hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
        audit['executable_hashes'][path] = actual
        if actual != expected:
            raise RuntimeError('Unknown native executable hash: ' + path)
    version = subprocess.run(['crisprme.py', '--version'], text=True, capture_output=True, timeout=30)
    audit['version'] = version.stdout.strip()
    audit['version_stderr'] = version.stderr
    if version.returncode or audit['version'] != plan['version']:
        raise RuntimeError('Unknown native version: ' + str(audit['version']))
    for i, command in enumerate(plan['commands']):
        started = time.monotonic()
        with (root/'logs'/('%02d.stdout.txt'%i)).open('w') as out, (root/'logs'/('%02d.stderr.txt'%i)).open('w') as err:
            result = subprocess.run(command, cwd=root, stdout=out, stderr=err)
        audit['commands'].append({'argv': command, 'exit_code': result.returncode,
                                   'wall_seconds': time.monotonic()-started})
        save()
        if result.returncode:
            raise RuntimeError('Native command exited nonzero: ' + str(result.returncode))
    audit['complete'] = True
    exit_code = 0
except Exception as exc:
    audit['error'] = str(exc)
finally:
    save()
sys.exit(exit_code)
'''


def _files(root):
    return {str(p.relative_to(root)): {"sha256": checksum(p), "bytes": p.stat().st_size}
            for p in sorted(Path(root).rglob("*")) if p.is_file() and not p.is_symlink()}


def _disk_bytes(root):
    total = 0
    for p in Path(root).rglob("*"):
        try:
            if p.is_file():
                total += p.stat().st_size
        except FileNotFoundError:
            pass
    return total


def _retain_raw(scratch, attempt):
    """Copy regular output files; record container-only index links as metadata."""
    raw = attempt / "raw"
    raw.mkdir(exist_ok=True)
    links = {}
    def ignore_links(directory, names):
        ignored = []
        for name in names:
            path = Path(directory) / name
            if path.is_symlink():
                links[str(path.relative_to(scratch))] = os.readlink(path)
                ignored.append(name)
        return ignored
    for member in ("Results", "logs", "Genomes", "PAMs", "window-map.json", "guides.txt",
                   "plan.json", "driver.py", "native-audit.json"):
        source = scratch / member
        if source.is_dir():
            shutil.copytree(source, raw / member, ignore=ignore_links, dirs_exist_ok=True)
        elif source.is_file():
            shutil.copy2(source, raw / member)
    (raw / "container-links.json").write_text(json_text(links))
    return raw


def _validate_windows(windows, check_limits=_no_limit_check):
    ids = set()
    for window in windows:
        check_limits()
        name, sequence, path = window["window_id"], window["sequence"], window["path"]
        if name in ids:
            raise PopulationError("Duplicate submitted native window identity")
        ids.add(name)
        if not sequence or set(sequence) - set("ACGT") or len(path) != len(sequence):
            raise PopulationError("Native windows require exact ACGT sequence and one path point per base")
        if any(point.get("base") != base or point.get("kind") not in {"reference", "insertion"}
               for point, base in zip(path, sequence)):
            raise PopulationError("Native input path does not match its retained sequence")
        if window.get("copy_index") not in (0, 1) \
                or window.get("phase_status") not in {"confirmed", "unphased", "homozygous"}:
            raise PopulationError("Native windows require explicit copy and source phase status")


def _validate_equivalence(equivalence, windows, profile, check_limits=_no_limit_check):
    status = equivalence.get("status")
    if status not in {"verified", "candidate_only", "not_applicable_synthetic"}:
        raise PopulationError("Invalid source-equivalence status")
    if status == "not_applicable_synthetic" and any(w.get("coverage_basis") != "synthetic_truth" for w in windows):
        raise PopulationError("Real donor windows cannot use synthetic source equivalence")
    if status != "verified":
        return
    audit = equivalence.get("reconstruction", {})
    required = {"verified", "source_panel_ids", "donor_sample_membership", "input_sha256",
                "genotype_policy", "variant_policy", "ordered_alt_rule", "phase_rule", "missingness_rule",
                "coverage_basis", "material_audit"}
    if not required <= audit.keys() or audit.get("verified") is not True \
            or profile["source_panel_id"] not in audit["source_panel_ids"]:
        raise PopulationError("Verified donor search requires the retained source reconstruction audit")
    membership = {row["donor_id"]: row["sample_ids"] for row in audit["donor_sample_membership"]}
    from .comparison_inputs import validate_prepared_window
    for window in windows:
        check_limits()
        validate_prepared_window(window)
        if window["source_version"] not in audit["input_sha256"].values():
            raise PopulationError("Native window source version is absent from its equivalence audit")
        if window["donor_id"] is not None and window["sample_id"] not in membership.get(window["donor_id"], []):
            raise PopulationError("Native window donor/sample is absent from its frozen panel")


def _commands(design, profile, pam, cpus):
    pam_path = f"/DATA/PAMs/20bp-{pam}-SpCas9.txt"
    common = ["--genome", "/DATA/Genomes/windows", "--pam", pam_path,
              "--bDNA", str(profile["max_dna_bulges"]), "--bRNA", str(profile["max_rna_bulges"]),
              "--thread", str(cpus)]
    commands = []
    if profile["max_dna_bulges"] or profile["max_rna_bulges"]:
        commands.append(["crisprme.py", "build-index-only", *common, "--path", "/DATA"])
    commands.append(["crisprme.py", "complete-search", *common, "--guide", "/DATA/guides.txt",
                     "--mm", str(profile["max_mismatches"]), "--max-total-edits",
                     str(profile["max_total_edits"]), "--merge", "0", "--output", "search",
                     "--ci-cd-test"])
    if profile["max_dna_bulges"] or profile["max_rna_bulges"]:
        commands[-1].extend(["--index-path", "/DATA/genome_library"])
    return commands


def _cached(directory, key):
    path = directory / "result.json"
    if not path.exists():
        return None
    seal = directory / "COMPLETE.json"
    if not seal.is_file() or read_json(seal) != {"stage_key": key, "result_sha256": checksum(path)}:
        raise PopulationError("Native cache completion checksum is missing or changed")
    result = read_json(path)
    if result.get("stage_key") != key or result.get("execution_status") != "complete" \
            or result.get("retention_status") != "complete":
        raise PopulationError("Native cache is not a complete matching job")
    for name, metadata in result["artifact_manifest"].items():
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise PopulationError("Invalid native cache artifact path")
        member = directory / name
        if member.is_symlink() or not member.is_file() or checksum(member) != metadata["sha256"]:
            raise PopulationError("Native cache artifact changed: " + name)
    if "input.json" not in result["artifact_manifest"]:
        raise PopulationError("Native cache is missing its effective input artifact")
    payload = read_json(directory / "input.json")
    if digest(payload) != key:
        raise PopulationError("Native cache effective input identity changed")
    _verify_window_source(directory, payload)
    result["reused"] = True
    result["artifact_path"] = str(directory)
    return result


def native_job_binding(result):
    """Bind native evidence into the enclosing comparison-stage manifest."""
    directory = Path(result["artifact_path"]).resolve()
    completion = directory / "COMPLETE.json"
    return {"path": str(directory), "stage_key": result["stage_key"],
            "result_sha256": checksum(directory / "result.json"),
            "completion_sha256": checksum(completion) if completion.exists() else None}


def verify_native_job(binding):
    """Verify externally bound complete or incomplete native evidence offline.

    The caller stores this binding inside its own checksummed stage manifest.
    A newly edited result and COMPLETE seal cannot replace that bound evidence.
    Partial evidence is verified for audit; it is never a reusable complete job.
    """
    required = {"path", "stage_key", "result_sha256", "completion_sha256"}
    if not isinstance(binding, dict) or set(binding) != required:
        raise PopulationError("Native evidence requires an enclosing-stage checksum binding")
    directory = Path(binding["path"])
    if directory.is_symlink() or not directory.is_dir():
        raise PopulationError("Native evidence directory is missing or is a symbolic link")
    directory = directory.resolve()
    path = directory / "result.json"
    if path.is_symlink() or not path.is_file() or checksum(path) != binding["result_sha256"]:
        raise PopulationError("Native result changed from its enclosing-stage checksum")
    result = read_json(path)
    if result.get("stage_key") != binding["stage_key"]:
        raise PopulationError("Native evidence stage identity changed")
    completion = directory / "COMPLETE.json"
    complete = result.get("execution_status") == "complete" and result.get("retention_status") == "complete"
    if complete:
        if completion.is_symlink() or not completion.is_file() \
                or checksum(completion) != binding["completion_sha256"] \
                or read_json(completion) != {"stage_key": binding["stage_key"], "result_sha256": binding["result_sha256"]}:
            raise PopulationError("Native completion seal changed from its enclosing-stage checksum")
    elif binding["completion_sha256"] is not None or completion.exists():
        raise PopulationError("Incomplete native evidence cannot carry a complete-job seal")
    manifest = result.get("artifact_manifest")
    if not isinstance(manifest, dict) or "input.json" not in manifest:
        raise PopulationError("Native evidence is missing its input/artifact manifest")
    for name, metadata in manifest.items():
        member = directory / name
        if Path(name).is_absolute() or ".." in Path(name).parts or member.is_symlink() \
                or not member.resolve().is_relative_to(directory) or not member.is_file() \
                or checksum(member) != metadata["sha256"]:
            raise PopulationError("Retained native artifact changed: " + name)
    payload = read_json(directory / "input.json")
    if digest(payload) != binding["stage_key"]:
        raise PopulationError("Native effective input/argument/code key does not match its evidence")
    _verify_window_source(directory, payload)
    return result


def run_native(windows, design, editor, profile, *, output_dir, pam=None, limits=None,
               cancelled=None, source_equivalence=None, scratch_root=None, host_budget=None):
    """Run one exact PAM policy and return hits, source details and execution scope.

    ``source_equivalence`` supplies the source reconstruction audit. Calls without
    that audit return candidate-only results.
    """
    limits = limits or NativeLimits()
    limits.validate()
    host = host_budget or HostBudget()
    host.check("native_input_validation")
    native.validate_policy(design, editor, profile)
    retained_windows = []
    for window in windows:
        host.check()
        retained_windows.append(window)
    windows = sorted(retained_windows, key=lambda row: row["window_id"])
    _validate_windows(windows, host.check)
    pam = pam or profile["pam_patterns"][0]
    if pam not in profile["pam_patterns"]:
        raise PopulationError("Native retrieval PAM must be a declared design-specific PAM")
    equivalence = source_equivalence or {"status": "candidate_only", "reason": "source_reconstruction_audit_absent"}
    _validate_equivalence(equivalence, windows, profile, host.check)
    commands = _commands(design, profile, pam, limits.cpus)
    code = {Path(p).name: checksum(p) for p in (__file__, native.__file__)}
    if code != _CODE_AT_IMPORT:
        raise PopulationError("Native adapter source changed after import; restart the process before executing another job")
    guard_length = len(design["guide_sequence"]) + len(pam) + max(profile["max_dna_bulges"], profile["max_rna_bulges"]) + 2
    skipped = [{"window_id": w["window_id"], "interval": w["interval"], "donor_id": w["donor_id"],
                "copy_index": w["copy_index"], "reason": "max_window_variants"} for w in windows
               if len(w.get("source_genotypes", w.get("source_variants", []))) > profile["max_window_variants"]]
    skipped_ids = {s["window_id"] for s in skipped}
    selected = [w for w in windows if w["window_id"] not in skipped_ids]
    host.check("native_sequence_projection_preparation")
    aliases, projection = _sequence_aliases(selected, guard_length, host.check)
    window_source = {"filename": "source-windows.jsonl", "sha256": _window_source_digest(windows, host.check), "count": len(windows)}
    row_cap = min(profile["max_rows"], limits.max_rows or profile["max_rows"])
    retention_profile = {**profile, "max_rows": row_cap}
    payload = {"adapter_contract": "materialized-donor-windows-v2", "window_source": window_source,
               "sequence_projection": {"version": SEQUENCE_PROJECTION, "sha256": digest(projection)},
               "design": design, "editor": editor, "profile": profile, "pam": pam,
               "source_equivalence": equivalence, "image": native.IMAGE,
               "source_revision": native.SOURCE_REVISION, "executable_hashes": native.EXECUTABLE_HASHES,
               "commands": commands, "driver": _DRIVER, "code": code, "limits": asdict(limits),
               "host_limits": asdict(host.limits), "effective_row_cap": row_cap,
               "window_guard": {"version": native.WINDOW_GUARD_TRANSFORM, "length_each_side": guard_length,
                                "base": "A", "admission": "all target and PAM bases must lie inside the original window"},
               "environment": {"MPLCONFIGDIR": "/tmp/matplotlib", "PYTHONHASHSEED": "0", "TZ": "UTC"}}
    key = digest(payload)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    complete = output / key
    cached = _cached(complete, key)
    host.check("native_cache_verification")
    if cached is not None:
        host.retain_native_binding(native_job_binding(cached))
        return cached
    attempt = output / (key + ".attempt-" + uuid.uuid4().hex[:12])
    attempt.mkdir(mode=0o700)
    (attempt / "input.json").write_text(json_text(payload))
    _write_window_sources(attempt / window_source["filename"], windows, window_source["sha256"], host.check)
    (attempt / "sequence-projection.json").write_text(json_text(projection))
    # Desktop bind mounts can be refused by macOS privacy controls. A dedicated
    # temporary directory mounts only the submitted sequences and tool inputs.
    scratch_parent = scratch_root or ("/private/tmp" if os.uname().sysname == "Darwin" else None)
    scratch = Path(tempfile.mkdtemp(prefix="bes-crisprme-", dir=scratch_parent))
    result = {"stage_key": key, "image": native.IMAGE, "version": native.VERSION,
              "source_revision": native.SOURCE_REVISION, "design_id": design["design_id"],
              "search_profile_id": profile["search_profile_id"], "limits": asdict(limits),
              "effective_limits": {"native": asdict(limits), "host": asdict(host.limits)}, "effective_row_cap": row_cap,
              "timeout_seconds": limits.timeout_seconds, "source_equivalence": equivalence["status"],
              "evaluated_window_ids": [], "skipped_windows": skipped, "hits": [],
              "execution_status": "failed", "retention_status": "unknown", "exit_code": None,
              "reason_codes": [], "reused": False,
              "source_window_count": len(windows), "eligible_source_window_count": len(selected),
              "native_sequence_count": len(aliases)}
    name = "bes-" + key[:12] + "-" + uuid.uuid4().hex[:8]
    process = None
    started = time.monotonic()
    peak_disk = 0
    try:
        host.check("native_input_serialization")
        for directory in ("Genomes/windows", "PAMs", "logs"):
            (scratch / directory).mkdir(parents=True)
        for alias, window in aliases.items():
            host.check()
            (scratch / "Genomes/windows" / (alias + ".fa")).write_text(
                ">" + alias + "\n" + native.native_sequence(window) + "\n")
        (scratch / "window-map.json").write_text(json_text(aliases))
        (scratch / "guides.txt").write_text(design["guide_sequence"] + "N" * len(pam) + "\n")
        (scratch / "PAMs" / f"20bp-{pam}-SpCas9.txt").write_text("N" * 20 + pam + " " + str(len(pam)) + "\n")
        (scratch / "driver.py").write_text(_DRIVER)
        (scratch / "plan.json").write_text(json_text({"version": native.VERSION,
            "executable_hashes": native.EXECUTABLE_HASHES, "commands": commands}))
        argv = ["docker", "run", "--name", name, "--pull=never", "--network=none", "--platform", "linux/arm64",
                "--cpus=" + str(limits.cpus), "--memory=" + str(limits.memory_bytes),
                "--memory-swap=" + str(limits.memory_bytes), "--pids-limit=" + str(limits.pids),
                "--mount", "type=bind,source=" + str(scratch) + ",target=/DATA", "-w", "/DATA"]
        for env_key, value in payload["environment"].items():
            argv.extend(["-e", env_key + "=" + value])
        argv.extend([native.IMAGE, "python3", "/DATA/driver.py"])
        result["argv"] = argv
        if not selected:
            result["reason_codes"].append("no_evaluable_windows")
        else:
            with (attempt / "stdout.txt").open("wb") as out, (attempt / "stderr.txt").open("wb") as err:
                process = subprocess.Popen(argv, stdin=subprocess.DEVNULL, stdout=out, stderr=err)
                while process.poll() is None:
                    host.check("native_execution")
                    peak_disk = max(peak_disk, _disk_bytes(scratch) + _disk_bytes(attempt))
                    reason = ("cancelled" if cancelled and cancelled() else "native_timeout"
                              if time.monotonic() - started > limits.timeout_seconds else "native_disk_cap"
                              if peak_disk > limits.disk_bytes else None)
                    if reason:
                        result["reason_codes"].append(reason)
                        subprocess.run(["docker", "kill", name], capture_output=True, timeout=15)
                        try:
                            process.wait(timeout=15)
                        except subprocess.TimeoutExpired:
                            process.kill()
                            process.wait(timeout=5)
                        break
                    time.sleep(0.1)
                result["exit_code"] = process.wait()
            inspected = subprocess.run(["docker", "inspect", "--format", "{{json .State}}", name],
                                       capture_output=True, text=True, timeout=15)
            if inspected.returncode == 0:
                result["container_state"] = json.loads(inspected.stdout)
                if result["container_state"].get("OOMKilled"):
                    result["reason_codes"].append("native_memory_cap")
            if result["exit_code"] != 0:
                result["reason_codes"].append("native_nonzero_exit")
        # Retain all native results, logs and exact inputs before opening a row.
        raw = _retain_raw(scratch, attempt)
        host.check("native_raw_retention")
        result["index_manifest"] = _files(scratch / "genome_library")
        result["index_id"] = "sha256:" + digest({"files": result["index_manifest"],
            "fasta": _files(raw / "Genomes"), "pam": pam, "image": native.IMAGE})
        expected_index = profile.get("index_id", {}).get("value")
        if expected_index is not None and expected_index != result["index_id"]:
            raise PopulationError("Constructed native index does not match the configured index identity")
        result["native_audit"] = read_json(raw / "native-audit.json") if (raw / "native-audit.json").exists() else None
        target_dir = raw / "Results/search/crispritz_targets"
        suffix = (f"20bp-{pam}-SpCas9.txt_guides.txt_{profile['max_mismatches']}_"
                  f"{profile['max_dna_bulges']}_{profile['max_rna_bulges']}.targets.txt")
        primary = target_dir / ("windows_" + suffix)
        secondary = target_dir / ("windows+__" + suffix)
        targets = [p for p in (primary, secondary) if p.is_file()]
        if set(target_dir.glob("*.targets.txt")) - {primary, secondary}:
            raise PopulationError("Unexpected native raw targets filename")
        log = raw / "Results/search/log_verbose.txt"
        log_text = log.read_text(errors="replace") if log.is_file() else ""
        zero_hit = ("No off-targets found for the provided guide(s) with the given parameters." in log_text
                    and "The search completed successfully; the result set is empty." in log_text)
        finished = "JOB END" in log_text or zero_hit
        if selected and not result["reason_codes"] and result["native_audit"] \
                and result["native_audit"].get("complete") and primary.is_file() and finished:
            host.check("native_raw_parsing")
            parsed = native.parse_targets(targets, windows=aliases, design=design, editor=editor,
                                          profile=retention_profile, check_limits=host.check)
            if zero_hit and parsed["raw_rows"]:
                raise PopulationError("Native zero-hit completion disagrees with its retained raw rows")
            result["native_zero_hit"] = zero_hit
            (attempt / "native-alignments.json").write_text(json_text(parsed))
            host.check("native_alignment_derivation")
            parsed = native.extend_leading_dna(parsed, windows=aliases, profile=retention_profile, check_limits=host.check)
            host.check("native_alignment_projection")
            parsed = _project_alignments(parsed, projection, retention_profile, host.check)
            for hit in parsed["hits"]:
                host.check()
                hit["design_id"] = design["design_id"]
                hit["search_profile_id"] = profile["search_profile_id"]
            result.update(parsed)
            result["evaluated_window_ids"] = [w["window_id"] for w in selected]
            result["execution_status"] = "complete" if not skipped and parsed["retention_status"] == "complete" else "partial"
            if equivalence["status"] == "candidate_only":
                result["execution_status"] = "partial"
                result["reason_codes"].append("source_panel_equivalence_unverified")
        elif selected and not result["reason_codes"]:
            result["reason_codes"].append("native_completion_or_raw_outputs_missing")
        if skipped:
            result["reason_codes"].append("high_density_windows_skipped")
            if result["execution_status"] == "failed" and not selected:
                result["execution_status"] = "partial"
    except ResourceLimit as exc:
        result["execution_status"] = "partial"
        result["retention_status"] = "truncated"
        result["error"] = str(exc)
        result["host_execution_audit"] = getattr(exc, "execution_audit", {})
        result["reason_codes"].extend(result["host_execution_audit"].get("reason_codes", ["host_resource_cap"]))
    except (OSError, PopulationError, subprocess.SubprocessError, ValueError) as exc:
        result["error"] = str(exc)
        result["reason_codes"].append("native_execution_or_parse_error")
    except KeyboardInterrupt:
        result["reason_codes"].append("interrupted")
    finally:
        if process is not None and process.poll() is None:
            try:
                subprocess.run(["docker", "kill", name], capture_output=True, timeout=15)
            finally:
                process.kill()
                process.wait(timeout=5)
        if process is not None and result["exit_code"] is None:
            result["exit_code"] = process.returncode
        if process is not None:
            try:
                subprocess.run(["docker", "rm", "-f", name], capture_output=True, timeout=15)
            except (OSError, subprocess.SubprocessError):
                result["reason_codes"].append("container_cleanup_failed")
                result["execution_status"] = "failed"
        result["wall_seconds"] = time.monotonic() - started
        result["host_resource_usage"] = host.metrics()
        if not (attempt / "raw").exists():
            try:
                _retain_raw(scratch, attempt)
            except OSError as exc:
                result["raw_retention_error"] = str(exc)
                result["reason_codes"].append("raw_retention_failed")
                result["execution_status"] = "failed"
        result["peak_disk_bytes"] = max(peak_disk, _disk_bytes(scratch) + _disk_bytes(attempt))
        if result["peak_disk_bytes"] > limits.disk_bytes:
            result["reason_codes"].append("native_disk_cap")
            result["execution_status"] = "partial" if result["hits"] else "failed"
        shutil.rmtree(scratch)
    result["artifact_manifest"] = _files(attempt)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result["adapter_peak_rss_bytes"] = rss if os.uname().sysname == "Darwin" else rss * 1024
    result_text = json_text(result)
    try:
        host.check("native_result_serialization")
    except ResourceLimit as exc:
        result["hits"] = []
        result["execution_status"] = "partial"
        result["retention_status"] = "truncated"
        result["host_execution_audit"] = getattr(exc, "execution_audit", {})
        result["reason_codes"].extend(result["host_execution_audit"].get("reason_codes", ["host_resource_cap"]))
        result_text = json_text(result)
    result["host_resource_usage"] = host.metrics()
    retained_size = _disk_bytes(attempt) + len(result_text.encode("utf-8"))
    if retained_size > limits.disk_bytes:
        result["hits"] = []
        result["execution_status"] = "failed"
        result["retention_status"] = "truncated"
        result["reason_codes"].append("normalized_result_disk_cap")
        result_text = json_text(result)
    result["peak_disk_bytes"] = max(result["peak_disk_bytes"], retained_size)
    (attempt / "result.json").write_text(json_text(result))
    try:
        host.check("native_result_retention")
    except ResourceLimit as exc:
        result["hits"] = []
        result["execution_status"] = "partial"
        result["retention_status"] = "truncated"
        result["host_execution_audit"] = getattr(exc, "execution_audit", {})
        result["reason_codes"].extend(result["host_execution_audit"].get("reason_codes", ["host_resource_cap"]))
        result["host_resource_usage"] = host.metrics()
        (attempt / "result.json").write_text(json_text(result))
    destination = attempt
    if result["execution_status"] == "complete" and result["retention_status"] == "complete":
        (attempt / "COMPLETE.json").write_text(json_text({"stage_key": key,
            "result_sha256": checksum(attempt / "result.json")}))
        os.replace(attempt, complete)
        destination = complete
    result["artifact_path"] = str(destination)
    host.retain_native_binding(native_job_binding(result))
    return result


def search_prepared(prepared, *, output_dir, limits=None, cancelled=None, host_budget=None):
    """Search every effective design policy over the same prepared source panel."""
    from .comparison_provenance import verify_sources
    host = host_budget or HostBudget(host_limits_from_study(prepared.study))
    limits = limits or native_limits_from_study(prepared.study)
    limits.validate()
    host.check("native_source_verification")
    verify_sources(prepared.study)
    windows = []
    host.check("native_source_reconstruction")
    for window in prepared.iter_windows():
        host.check()
        windows.append(window)
    host.check()
    verify_sources(prepared.study)
    host.check()
    audit = getattr(prepared, "source_equivalence_audit", {})
    equivalence = {"status": "verified" if audit.get("verified") is True else "candidate_only",
                   "reconstruction": audit, "reconstruction_code": prepared.study["code"]}
    if audit.get("coverage_basis") == "synthetic_truth":
        equivalence["status"] = "not_applicable_synthetic"
    hits, executions, jobs = [], [], []
    for item in prepared.profiles:
        design, editor, profile = item["design"], item["editor"], item["profile"]
        profile_jobs = []
        row_cap = min(profile["max_rows"], limits.max_rows or profile["max_rows"])
        profile_rows, profile_truncated = 0, False
        for pam in sorted(set(profile["pam_patterns"])):
            host.check("native_profile_execution")
            job = run_native(windows, design, editor, profile, output_dir=output_dir, pam=pam,
                             limits=limits, cancelled=cancelled, source_equivalence=equivalence, host_budget=host)
            jobs.append(job)
            profile_jobs.append(job)
            for hit in job["hits"]:
                host.check("native_result_materialization")
                if profile_rows >= row_cap:
                    profile_truncated = True
                    break
                hits.append(hit)
                profile_rows += 1
            host.check("native_profile_complete")
        statuses = {job["execution_status"] for job in profile_jobs}
        retention = {job["retention_status"] for job in profile_jobs}
        evaluated_ids = set.intersection(*(set(job["evaluated_window_ids"]) for job in profile_jobs))
        skipped_by_window = {row["window_id"]: row for job in profile_jobs for row in job["skipped_windows"]}
        for window in windows:
            if window["window_id"] not in evaluated_ids and window["window_id"] not in skipped_by_window:
                skipped_by_window[window["window_id"]] = {"window_id": window["window_id"],
                    "interval": window["interval"], "donor_id": window["donor_id"], "copy_index": window["copy_index"],
                    "reason": "native_pam_job_incomplete"}
        executions.append({"design_id": design["design_id"], "search_profile_id": profile["search_profile_id"],
            "execution_status": "partial" if profile_truncated else "complete" if statuses == {"complete"} else "failed" if statuses == {"failed"} else "partial",
            "retention_status": "truncated" if profile_truncated else "complete" if retention == {"complete"} else "truncated" if "truncated" in retention else "unknown",
            "exit_code": 0 if all(job["exit_code"] == 0 for job in profile_jobs) else None,
            "timeout_seconds": limits.timeout_seconds,
            "effective_limits": {"native": asdict(limits), "host": asdict(host.limits)},
            "effective_row_cap": row_cap, "host_resource_usage": host.metrics(),
            "reason_codes": sorted({r for job in profile_jobs for r in job["reason_codes"]} | ({"native_profile_row_cap"} if profile_truncated else set())),
            "evaluated_window_ids": sorted(evaluated_ids),
            "skipped_windows": [skipped_by_window[key] for key in sorted(skipped_by_window)],
            "source_equivalence": equivalence["status"] if equivalence else "candidate_only",
            "native_manifest": [native_job_binding(job) for job in profile_jobs]})
    return {"hits": hits, "execution": executions, "windows": windows, "jobs": jobs}
