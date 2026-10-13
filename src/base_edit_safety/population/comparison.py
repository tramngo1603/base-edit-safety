"""Compare exact donor sequences under each active design's own search policy."""
from __future__ import annotations

from collections import Counter, OrderedDict, defaultdict
from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import resource
import sys
import time

from .common import PopulationError, ResourceLimit, known, unknown
from .comparison_contracts import make_record, validate_catalog
from .comparison_identity import canonical_bytes
from .comparison_inputs import (contains, oriented_path, path_bounds, prepare_inputs,
                                support_bounds, validate_prepared_window)
from .comparison_provenance import verify_sources
from .identity import matches_pam, reverse_complement

LOCAL_BACKEND_NAME = "comparative_local"
LOCAL_BACKEND_VERSION = "v1"
LOCAL_BACKEND_DIGEST = hashlib.sha256(b"comparative_local_v1").hexdigest()


@dataclass(frozen=True)
class LocalLimits:
    timeout_seconds: int = 600
    max_resident_bytes: int = 2 * 1024 ** 3

    def validate(self):
        if type(self.timeout_seconds) is not int or self.timeout_seconds < 1 or type(self.max_resident_bytes) is not int or self.max_resident_bytes < 1:
            raise PopulationError("Local runtime and memory ceilings must be positive integers")
        return self


class _LocalBudget:
    def __init__(self, limits, *, started=None, cpu_started=None):
        self.limits = limits.validate()
        self.started = time.monotonic() if started is None else started
        self.cpu_started = time.process_time() if cpu_started is None else cpu_started

    def metrics(self):
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return {"wall_time_seconds": time.monotonic() - self.started,
                "cpu_time_seconds": time.process_time() - self.cpu_started,
                "peak_resident_bytes": int(peak if sys.platform == "darwin" else peak * 1024),
                "timeout_seconds": self.limits.timeout_seconds,
                "max_resident_bytes": self.limits.max_resident_bytes}

    def check(self):
        metrics = self.metrics()
        if metrics["wall_time_seconds"] > self.limits.timeout_seconds:
            self.fail("Local comparison runtime ceiling exceeded; no complete result retained", "local_timeout", metrics)
        if metrics["peak_resident_bytes"] > self.limits.max_resident_bytes:
            self.fail("Local comparison memory ceiling exceeded; no complete result retained", "local_memory_cap", metrics)

    def fail(self, message, reason, metrics):
        error = ResourceLimit(message)
        error.execution_audit = {"execution_status": "partial", "retention_status": "incomplete",
            "reason_codes": [reason], "phase": getattr(self, "phase", "local_search"),
            "effective_limits": {"host" if isinstance(self, HostBudget) else "local": asdict(self.limits)},
            "resource_usage": metrics, "native_manifest": deepcopy(getattr(self, "native_manifest", []))}
        raise error


def _configured_limits(study, name):
    config = study.get("config", study)
    return config["policies"]["search"].get("execution_limits", {}).get(name, {})


def local_limits_from_study(study):
    """Resolve the declared local ceilings without changing the search profile."""
    try:
        return LocalLimits(**_configured_limits(study, "local")).validate()
    except TypeError as exc:
        raise PopulationError("Invalid search.execution_limits.local fields") from exc


def host_limits_from_study(study):
    """Host limits cover total workflow wall time and process peak resident memory."""
    try:
        return LocalLimits(**_configured_limits(study, "host")).validate()
    except TypeError as exc:
        raise PopulationError("Invalid search.execution_limits.host fields") from exc


class HostBudget(_LocalBudget):
    """Share one host deadline across source preparation, native work and mapping.

    Checks occur before and after bulk operations and during window/row loops.
    Peak memory belongs to this Python process, including earlier operations.
    Native waiting counts toward the same overall wall-time ceiling.
    """

    def __init__(self, limits=None, *, started=None, cpu_started=None):
        super().__init__(limits or LocalLimits(), started=started, cpu_started=cpu_started)
        self.phase = "source_preparation"
        self.native_manifest = []

    def check(self, phase=None):
        if phase is not None:
            self.phase = phase
        metrics = self.metrics()
        if metrics["wall_time_seconds"] > self.limits.timeout_seconds:
            self.fail("Host reconstruction/normalization runtime ceiling exceeded; no complete result retained", "host_timeout", metrics)
        if metrics["peak_resident_bytes"] > self.limits.max_resident_bytes:
            self.fail("Host reconstruction/normalization memory ceiling exceeded; no complete result retained", "host_memory_cap", metrics)

    def retain_native_binding(self, binding):
        if binding not in self.native_manifest:
            self.native_manifest.append(deepcopy(binding))


def _no_limit_check():
    pass


def _alignment_counts(guide_alignment, target_alignment):
    guide_pos, target_pos = 0, 0
    mismatches, dna, rna, columns = [], [], [], []
    events, previous_gap = 0, None
    for index, (guide, target) in enumerate(zip(guide_alignment, target_alignment)):
        gap = "dna" if guide == "-" else "rna" if target == "-" else None
        if guide == target == "-":
            raise PopulationError("Alignment cannot contain two gaps in one column")
        if gap and gap != previous_gap:
            events += 1
        previous_gap = gap
        if guide != "-":
            guide_pos += 1
        if guide == "-":
            dna.append(index)
        elif target == "-":
            rna.append(guide_pos)
        elif guide != target:
            mismatches.append(guide_pos)
        columns.append({"guide": None if guide == "-" else guide,
                        "target": None if target == "-" else target,
                        "guide_position_1based": None if guide == "-" else guide_pos,
                        "target_path_index": None if target == "-" else target_pos})
        if target != "-":
            target_pos += 1
    return columns, mismatches, dna, rna, events


def _admitted(guide_alignment, target_alignment, pam, profile):
    _, mismatches, dna, rna, events = _alignment_counts(guide_alignment, target_alignment)
    total = len(mismatches) + (events if profile["total_edit_accounting"] == "mismatch_plus_bulge_events" else len(dna) + len(rna))
    return (matches_pam(pam, profile["pam_patterns"]) and len(mismatches) <= profile["max_mismatches"]
            and len(dna) <= profile["max_dna_bulges"] and len(rna) <= profile["max_rna_bulges"]
            and total <= profile["max_total_edits"])


def enumerate_alignments(guide, target, profile, *, check_limits=_no_limit_check):
    """Enumerate complete global columns, retaining distinct admitted alignments.

    DNA/RNA bulge limits count gapped bases. Total-edit accounting may instead
    count contiguous gap events, as explicitly declared by the profile.
    """
    states = 0
    def visit(i, j, mismatches, dna, rna, gap_events, previous_gap, ga, ta):
        nonlocal states
        states += 1
        if states % 128 == 1:
            check_limits()
        if states > profile["max_rows"] * 100:
            raise ResourceLimit("Local alignment-state ceiling exceeded; no complete result retained")
        total = mismatches + (gap_events if profile["total_edit_accounting"] == "mismatch_plus_bulge_events" else dna + rna)
        if mismatches > profile["max_mismatches"] or dna > profile["max_dna_bulges"] or rna > profile["max_rna_bulges"] or total > profile["max_total_edits"]:
            return
        if len(target) - j - (len(guide) - i) > profile["max_dna_bulges"] - dna or len(guide) - i - (len(target) - j) > profile["max_rna_bulges"] - rna:
            return
        if i == len(guide) and j == len(target):
            yield ga, ta
            return
        if i < len(guide) and j < len(target):
            yield from visit(i + 1, j + 1, mismatches + (guide[i] != target[j]), dna, rna, gap_events, None, ga + guide[i], ta + target[j])
        if j < len(target) and dna < profile["max_dna_bulges"]:
            yield from visit(i, j + 1, mismatches, dna + 1, rna, gap_events + (previous_gap != "dna"), "dna", ga + "-", ta + target[j])
        if i < len(guide) and rna < profile["max_rna_bulges"]:
            yield from visit(i + 1, j, mismatches, dna, rna + 1, gap_events + (previous_gap != "rna"), "rna", ga + guide[i], ta + "-")
    yield from visit(0, 0, 0, 0, 0, 0, None, "", "")


def local_hits(window, item, *, check_limits=_no_limit_check):
    design, editor, profile = item["design"], item["editor"], item["profile"]
    geometry = editor["guide_geometry"]
    spacer_length, pam_length = geometry["spacer_length"], geometry["pam_length"]
    for length in range(max(1, spacer_length - profile["max_rna_bulges"]), spacer_length + profile["max_dna_bulges"] + 1):
        width = length + pam_length
        for start in range(len(window["sequence"]) - width + 1):
            check_limits()
            for strand in profile["strands"]:
                indexes = list(range(start, start + width))
                if strand == "-":
                    indexes.reverse()
                bases = "".join(window["sequence"][index] for index in indexes)
                if strand == "-":
                    bases = reverse_complement(bases[::-1])
                if geometry["pam_side"] == "3prime":
                    target, pam = bases[:length], bases[length:]
                    target_indexes, pam_indexes = indexes[:length], indexes[length:]
                else:
                    pam, target = bases[:pam_length], bases[pam_length:]
                    pam_indexes, target_indexes = indexes[:pam_length], indexes[pam_length:]
                if not matches_pam(pam, profile["pam_patterns"]):
                    continue
                for ga, ta in enumerate_alignments(design["guide_sequence"], target, profile, check_limits=check_limits):
                    yield {"window_id": window["window_id"], "design_id": design["record_id"],
                           "start": start, "end": start + width, "strand": strand,
                           "guide_alignment": ga, "target_alignment": ta,
                           "target_indexes": target_indexes, "pam_sequence": pam, "pam_indexes": pam_indexes,
                           "native_row_refs": []}


def _union_intervals(rows):
    result = []
    for row in sorted(rows, key=lambda r: (r["contig"], r["start"], r["end"])):
        if result and result[-1]["contig"] == row["contig"] and row["start"] <= result[-1]["end"]:
            result[-1]["end"] = max(result[-1]["end"], row["end"])
        else:
            result.append({key: row[key] for key in ("contig", "start", "end")})
    return result


class ComparisonBuilder:
    """Accumulate canonical hit records without retaining all reconstructed windows."""

    def __init__(self, prepared, *, check_limits=_no_limit_check):
        self.prepared = prepared
        self.check_limits = check_limits
        self.records = deepcopy(prepared.catalog)
        self.original_ids = set(self.records)
        self.profiles = {item["design"]["record_id"]: item for item in prepared.profiles}
        self.nominees = defaultdict(set)
        self.hit_alignments = defaultdict(set)
        self.hit_donors = defaultdict(set)
        self.hit_rows = []
        self.admission_exclusions = Counter()
        self.annotations = {}
        self._record_templates = OrderedDict()
        for row in prepared.source_records:
            self.add(row)

    @staticmethod
    def _same_values(left, right):
        if type(left) is not type(right):
            return False
        if isinstance(left, dict):
            return left.keys() == right.keys() and all(
                ComparisonBuilder._same_values(value, right[key]) for key, value in left.items())
        if isinstance(left, list):
            return len(left) == len(right) and all(
                ComparisonBuilder._same_values(a, b) for a, b in zip(left, right))
        if type(left) is float:
            return repr(left) == repr(right)
        return left == right

    def _record(self, kind, **values):
        """Reuse a bounded validated template only for exactly equal input values."""
        encoded = canonical_bytes(values)
        key = (kind, encoded)
        template = self._record_templates.get(key)
        if template is not None and all(
                field in template and self._same_values(value, template[field])
                for field, value in values.items()):
            row = deepcopy(template)
            self._record_templates.move_to_end(key)
        else:
            row = make_record(kind, **values)
            if len(encoded) <= 8192:
                self._record_templates[key] = deepcopy(row)
                self._record_templates.move_to_end(key)
                if len(self._record_templates) > 256:
                    self._record_templates.popitem(last=False)
        return self.add(row)

    def add(self, row):
        previous = self.records.get(row["record_id"])
        if previous is not None and previous != row:
            kind = row["record_type"]
            if kind == "comparison_sequence_configuration":
                previous["variants"] = sorted({canonical_bytes(value): value for value in previous["variants"] + row["variants"]}.values(), key=lambda value: (value["contig"], value["pos0"], value["source_row"]))
                if row["support_status"] == "supported":
                    previous["support_status"], previous["support_reason"] = "supported", None
                return previous
            if kind in {"comparison_sequence_occurrence", "comparison_alignment"}:
                field = "source_refs" if kind == "comparison_sequence_occurrence" else "native_row_refs"
                merged = sorted({canonical_bytes(value): value for value in previous[field] + row[field]}.values(), key=canonical_bytes)
                left, right = deepcopy(previous), deepcopy(row)
                left[field] = right[field] = merged
                if left != right:
                    raise PopulationError("Conflicting scientific records share an identity")
                previous[field] = merged
                return previous
            if kind == "comparison_design_configuration" and previous["role"] != row["role"]:
                if row["role"] == "intended":
                    previous["role"] = "intended"
                return previous
            raise PopulationError("Conflicting comparative records share an identity")
        self.records[row["record_id"]] = row
        return row

    def locus(self, contig, path, strand, design):
        lo, hi = path_bounds(path)
        for identifier in design["intended_locus_ids"]:
            row = self.records[identifier]
            if (row["contig"], row["start"], row["end"], row["strand"]) == (contig, lo, hi, strand):
                return row, "intended"
        return self._record("locus", reference_bundle_id=self.prepared.reference_bundle_id,
            contig=contig, start=lo, end=hi, strand=strand, locus_name=f"{contig}:{lo}-{hi}:{strand}",
            gene=unknown("no gene annotation supplied for this nominated locus"), path=[]), "off_target"

    def sequence(self, design, locus, sequence, path, variants, status, role):
        seq = self._record("sequence_configuration", reference_bundle_id=self.prepared.reference_bundle_id,
            locus_id=locus["record_id"], contig=locus["contig"], strand=locus["strand"], sequence=sequence,
            path=path, variants=variants, support_status=status,
            support_reason=None if status == "supported" else "retained-callset modeled genome" if status == "modeled" else "unresolved source sequence")
        dc = self._record("design_configuration", design_id=design["record_id"], sequence_configuration_id=seq["record_id"], role=role)
        return seq, dc

    def occurrence(self, window, sequence):
        if window["donor_id"] is None:
            return None
        return self._record("sequence_occurrence", sequence_configuration_id=sequence["record_id"],
            donor_id=window["donor_id"], sample_id=window["sample_id"], copy_index=window["copy_index"],
            phase_block=known(window["phase_block"]) if window["phase_block"] is not None else unknown("homozygous or single locally unordered heterozygous site"),
            phase_status="unphased" if window["phase_status"] == "unphased" else "confirmed",
            source_version=window["source_version"], source_refs=window["source_refs"])

    def source_refs(self, raw_refs, derivation=None):
        result = []
        for raw in raw_refs:
            if "source_artifact_id" in raw:
                result.append(raw)
                continue
            sha = raw.get("sha256", raw.get("raw_sha256"))
            filename = raw.get("file", raw.get("filename", "native-targets.txt"))
            source = self._record("source_artifact", sha256=sha,
                source_id=raw.get("source_id", "native-search-output:" + sha), release=raw.get("version", "native-search-v1"),
                uri="urn:bes:native-search-output:" + sha, media_type="text/tab-separated-values", scope="complete_file",
                parent_artifact_id=unknown("raw output captured by the pinned native execution"),
                locator={"kind": "file", "member": filename}, access="synthetic" if self.prepared.study["config"]["mode"] == "synthetic_test" else "restricted" if self.prepared.study["config"]["mode"] == "research_sample_planning" else "public",
                use_terms="Native search result; backend terms remain in the execution manifest")
            note = raw.get("note", "Verbatim native alignment row")
            if derivation:
                note += "; derived alignment: " + canonical_bytes(derivation).decode("utf-8")
            result.append({"source_artifact_id": source["record_id"],
                           "locator": {"kind": "text", "page": 1, "line_start": raw["line"], "line_end": raw["line"]},
                           "note": note})
        return result

    def alignment(self, dc, ga, ta, pam, pam_path, pam_side, refs):
        columns, mismatches, dna, rna, _ = _alignment_counts(ga, ta)
        if pam_side == "5prime":
            for column in columns:
                if column["target_path_index"] is not None:
                    column["target_path_index"] += len(pam)
        return self._record("alignment", design_configuration_id=dc["record_id"], columns=columns,
            mismatch_positions=mismatches, dna_bulge_positions=dna, rna_bulge_positions=rna,
            pam_sequence=pam, pam_path=pam_path, native_row_refs=refs)

    def add_hit(self, window, hit):
        validate_prepared_window(window)
        return self._map_hit(window, hit)

    def _map_hit(self, window, hit):
        """Map one hit after the caller has validated its source window."""
        if hit["window_id"] != window["window_id"] or hit["design_id"] not in self.profiles:
            raise PopulationError("Search hit does not resolve its exact window and active design")
        item = self.profiles[hit["design_id"]]
        design, geometry, profile = item["design"], item["editor"]["guide_geometry"], item["profile"]
        strand, ga, ta = hit["strand"], hit["guide_alignment"], hit["target_alignment"]
        if strand not in {"+", "-"} or len(ga) != len(ta) or ga.replace("-", "") != design["guide_sequence"]:
            raise PopulationError("Native alignment guide/strand does not match the active design")
        target_indexes, pam_indexes = hit["target_indexes"], hit["pam_indexes"]
        indexes = target_indexes + pam_indexes if geometry["pam_side"] == "3prime" else pam_indexes + target_indexes
        if any(type(index) is not int or not 0 <= index < len(window["path"]) for index in indexes):
            raise PopulationError("Native alignment indexes lie outside the submitted sequence")
        expected_order = list(range(min(indexes), max(indexes) + 1))
        if strand == "-":
            expected_order.reverse()
        if indexes != expected_order or (hit["start"], hit["end"]) != (min(indexes), max(indexes) + 1):
            raise PopulationError("Native alignment span/path indexes disagree")
        path = [{**window["path"][index], "base": window["sequence"][index] if strand == "+" else reverse_complement(window["sequence"][index])} for index in indexes]
        pam_path = path[-len(pam_indexes):] if geometry["pam_side"] == "3prime" else path[:len(pam_indexes)]
        target = "".join(window["sequence"][index] if strand == "+" else reverse_complement(window["sequence"][index]) for index in target_indexes)
        pam = "".join(point["base"] for point in pam_path)
        if ta.replace("-", "") != target or hit["pam_sequence"] != pam or len(pam) != geometry["pam_length"]:
            raise PopulationError("Native alignment bases do not match the exact submitted donor path")
        lo, hi = path_bounds(path)
        span = {"contig": window["contig"], "start": lo, "end": hi}
        support_lo, support_hi = support_bounds(path, window["source_variants"])
        support_span = {"contig": window["contig"], "start": support_lo, "end": support_hi}
        if not any(contains(interval, support_span) for interval in self.prepared.domain["intervals"]):
            self.admission_exclusions["outside_shared_domain"] += 1
            return None
        if strand not in profile["strands"] or not _admitted(ga, ta, pam, profile):
            self.admission_exclusions["outside_exact_design_profile"] += 1
            return None
        local_variants = [variant for variant in window["source_variants"] if variant["contig"] == span["contig"] and variant["pos0"] < hi and lo < variant["pos0"] + len(variant["reference"])]
        locus, role = self.locus(window["contig"], path, strand, design)
        sequence, dc = self.sequence(design, locus, "".join(point["base"] for point in path), path,
                                     local_variants, window["support_status"], role)
        self.occurrence(window, sequence)
        alignment = self.alignment(dc, ga, ta, pam, pam_path, geometry["pam_side"], self.source_refs(hit.get("native_row_refs", []), hit.get("derivation")))
        self.nominees[design["record_id"]].add(dc["record_id"])
        self.hit_alignments[design["record_id"]].add(alignment["record_id"])
        self.hit_donors[(design["record_id"], window["donor_id"])].add(dc["record_id"])
        self.hit_rows.append({"design_id": design["record_id"], "donor_id": window["donor_id"], "copy_index": window["copy_index"],
                             "design_configuration_id": dc["record_id"], "alignment_id": alignment["record_id"],
                             "window_id": window["window_id"], "derivation": hit.get("derivation")})
        if "sequence_projection" in hit:
            self.hit_rows[-1]["sequence_projection"] = deepcopy(hit["sequence_projection"])
        self.annotate(item, dc, alignment)
        return alignment

    def annotate(self, item, dc, alignment):
        editor = item["editor"]
        positions = [column["guide_position_1based"] for column in alignment["columns"]
                     if column["guide_position_1based"] is not None and column["target"] == editor["conversion"]["substrate_base"]]
        window = editor["editing_window"]
        within = unknown(window["reason"]) if "value" in window else known([position for position in positions if window["start_1based"] <= position <= window["end_1based"]])
        sequence = self.records[dc["sequence_configuration_id"]]
        intended_positions = {edit["pos0"] for edit in item["design"]["desired_edit_definitions"] if edit["locus_id"] == sequence["locus_id"]}
        bystanders = []
        for column in alignment["columns"]:
            index = column["target_path_index"]
            if index is not None and column["guide_position_1based"] in positions and sequence["path"][index].get("pos0") not in intended_positions:
                bystanders.append(column["guide_position_1based"])
        self.annotations[alignment["record_id"]] = {"alignment_id": alignment["record_id"], "design_configuration_id": dc["record_id"],
            "substrate_positions_1based": positions, "within_declared_window": within,
            "potential_bystander_positions_1based": bystanders, "interpretation": "Sequence annotation; no editing probability or assumed outcome"}

    def fixed_target_alignments(self, item, dc, target, pam_path, refs):
        """Retain permitted gap alternatives and an available ungapped diagnostic."""
        guide = item["design"]["guide_sequence"]
        pam = "".join(point["base"] for point in pam_path)
        side = item["editor"]["guide_geometry"]["pam_side"]
        profile = item["profile"]
        rows = {}
        alternatives = enumerate_alignments(guide, target, profile, check_limits=self.check_limits)
        for ga, ta in alternatives:
            alignment = self.alignment(dc, ga, ta, pam, pam_path, side, refs)
            rows[alignment["record_id"]] = (alignment, _admitted(ga, ta, pam, profile))
        if len(guide) == len(target):
            alignment = self.alignment(dc, guide, target, pam, pam_path, side, refs)
            rows[alignment["record_id"]] = (alignment, _admitted(guide, target, pam, profile))
        return sorted(rows.values(), key=lambda pair: (not pair[1], pair[0]["record_id"]))

    def intended_assessments(self):
        for item in self.prepared.profiles:
            design, geometry = item["design"], item["editor"]["guide_geometry"]
            for locus_id in design["intended_locus_ids"]:
                locus = self.records[locus_id]
                interval = {key: locus[key] for key in ("contig", "start", "end")}
                for donor in sorted(self.prepared.donor_records.values(), key=lambda row: row["record_id"]):
                    for copy_index in (0, 1):
                        self.check_limits()
                        material = self.prepared.reconstruct(interval, donor["record_id"], copy_index)
                        status, reason = self.prepared.sequence_support(material, interval)
                        variants, _, refs = self.prepared.variant_metadata(interval, donor["record_id"], copy_index)
                        alignments, seq_id = [], unknown(reason)
                        pam_result = unknown(reason)
                        reasons = [] if status == "supported" else [reason]
                        if status != "unresolved" and material["sequence_path"]:
                            path = oriented_path(material["sequence_path"], locus["strand"])
                            sequence, dc = self.sequence(design, locus, "".join(point["base"] for point in path), path, variants, status, "intended")
                            seq_id = known(sequence["record_id"])
                            local_window = self.prepared._window(material, material["sequence_path"], donor["record_id"], copy_index, interval)
                            self.occurrence(local_window, sequence)
                            pam_len = geometry["pam_length"]
                            original_path = locus["path"] or oriented_path([{"kind": "reference", "pos0": locus["start"] + i, "base": base} for i, base in enumerate(self.prepared.reference.fetch(locus["contig"], locus["start"], locus["end"]))], locus["strand"])
                            original_pam = original_path[-pam_len:] if geometry["pam_side"] == "3prime" else original_path[:pam_len]
                            pam_positions = {point["pos0"] for point in original_pam}
                            observed_pam = [point for point in path if point["kind"] == "reference" and point["pos0"] in pam_positions or point["kind"] == "insertion" and min(pam_positions) < point["boundary0"] <= max(pam_positions)]
                            pam = "".join(point["base"] for point in observed_pam)
                            compatible = len(pam) == pam_len and matches_pam(pam, item["profile"]["pam_patterns"])
                            pam_result = known(compatible)
                            if not compatible:
                                reasons.append("intended_pam_incompatible")
                            if len(path) > pam_len and len(observed_pam) == pam_len:
                                target = sequence["sequence"][:-pam_len] if geometry["pam_side"] == "3prime" else sequence["sequence"][pam_len:]
                                alternatives = self.fixed_target_alignments(item, dc, target, observed_pam, refs)
                                for alignment, admitted in alternatives:
                                    alignments.append(alignment["record_id"])
                                    self.annotate(item, dc, alignment)
                                if alternatives and not any(admitted for _, admitted in alternatives):
                                    reasons.append("intended_alignment_outside_nomination_profile")
                                elif not alternatives:
                                    reasons.append("intended_alignment_not_representable_with_declared_bulge_limits")
                            else:
                                reasons.append("intended_pam_geometry_changed_no_fixed_alignment")
                        self._record("intended_assessment", donor_id=donor["record_id"], design_id=design["record_id"], locus_id=locus_id,
                            sequence_configuration_id=seq_id, support_status=status, pam_compatible=pam_result,
                            alignment_ids=alignments, reference_differences=variants, reason_codes=sorted(set(reasons)))

    def scopes(self, execution):
        by_design = {row["design_id"]: row for row in execution}
        if set(by_design) != set(self.profiles) or len(execution) != len(by_design):
            raise PopulationError("Search execution must account for every active design")
        if not self.prepared.source_equivalence_audit.get("verified") or not self.prepared.material_audit:
            raise PopulationError("Search scope requires a completed source reconstruction audit")
        all_windows = {identifier for audit in self.prepared.material_audit for identifier in audit["window_intervals"]}
        for design_id, item in sorted(self.profiles.items()):
            run = by_design[design_id]
            if run["search_profile_id"] != item["profile"]["record_id"]:
                raise PopulationError("Search execution claims a different effective profile")
            if "evaluated_window_ids" not in run or set(run["evaluated_window_ids"]) - all_windows:
                raise PopulationError("Backend execution must resolve its exact evaluated source windows")
            evaluated_ids = set(run["evaluated_window_ids"])
            skipped_windows = {row["window_id"]: row for row in run.get("skipped_windows", [])}
            if skipped_windows.keys() - all_windows or skipped_windows.keys() & evaluated_ids:
                raise PopulationError("Backend skipped and evaluated windows contradict each other")
            for donor in sorted(self.prepared.donor_records.values(), key=lambda row: row["record_id"]):
                relevant = [audit for audit in self.prepared.material_audit if audit["donor_id"] == donor["record_id"]]
                skipped, evaluated = [], []
                def skip(span, reason):
                    for requested in self.prepared.domain["intervals"]:
                        if requested["contig"] == span["contig"] and requested["start"] < span["end"] and span["start"] < requested["end"]:
                            skipped.append({"interval": {"contig": span["contig"], "start": max(span["start"], requested["start"]), "end": min(span["end"], requested["end"])}, "reason": reason})
                for audit in relevant:
                    for gap in audit["skipped_intervals"]:
                        span = {"contig": audit["interval"]["contig"], "start": gap["start"], "end": gap["end"]}
                        skip(span, gap["reason"])
                    for identifier, span in audit["window_intervals"].items():
                        if identifier in evaluated_ids:
                            for requested in self.prepared.domain["intervals"]:
                                if requested["contig"] == span["contig"] and requested["start"] < span["end"] and span["start"] < requested["end"]:
                                    evaluated.append({"contig": span["contig"], "start": max(span["start"], requested["start"]), "end": min(span["end"], requested["end"])})
                        else:
                            native_skip = skipped_windows.get(identifier)
                            if native_skip and (native_skip["donor_id"], native_skip["copy_index"], native_skip["interval"]) != (donor["record_id"], audit["copy_index"], span):
                                raise PopulationError("Backend skip metadata differs from its immutable source window")
                            skip(span, native_skip["reason"] if native_skip else "backend_window_not_evaluated")
                reasons = set(run["reason_codes"]) | {row["reason"] for row in skipped}
                verified = run["source_equivalence"] == ("not_applicable_synthetic" if self.prepared.study["config"]["mode"] == "synthetic_test" else "verified")
                if not verified:
                    reasons.add("source_panel_equivalence_unverified")
                complete = run["execution_status"] == "complete" and run["retention_status"] == "complete" and run["exit_code"] == 0 and not skipped and verified
                if not complete and not reasons:
                    reasons.add("backend_execution_incomplete")
                self._record("search_scope", design_id=design_id, donor_id=known(donor["record_id"]), source_panel_id=donor["source_panel_id"],
                    search_profile_id=item["profile"]["record_id"], requested_intervals=self.prepared.domain["intervals"],
                    evaluated_intervals=_union_intervals(evaluated), skipped_intervals=list({canonical_bytes(row): row for row in skipped}.values()),
                    execution_status="complete" if complete else "failed" if run["execution_status"] == "failed" else "partial",
                    retention_status=run["retention_status"], backend_limits={"max_rows": min(item["profile"]["max_rows"], run.get("effective_row_cap", item["profile"]["max_rows"])), "max_window_variants": item["profile"]["max_window_variants"], "timeout_seconds": max(1, int(run["timeout_seconds"]))},
                    exit_code=known(run["exit_code"]) if run["exit_code"] is not None else unknown("backend exit unavailable"), reason_codes=sorted(reasons))

    def finish(self, execution):
        self.intended_assessments()
        self.scopes(execution)
        self.variant_effects_and_carriage()
        self.summaries()
        new_rows = sorted((row for identifier, row in self.records.items() if identifier not in self.original_ids), key=lambda row: (row["record_type"], row["record_id"]))
        validate_catalog([*self.prepared.study["config"]["catalogs"], *new_rows])
        self.check_limits()
        result = {"comparison-records.jsonl": new_rows,
            "search-audit.json": {"domain_id": self.prepared.domain["record_id"], "execution": execution,
                "source_equivalence": self.prepared.source_equivalence_audit,
                "admission_exclusions": dict(sorted(self.admission_exclusions.items())),
                "hits": sorted({canonical_bytes(row): row for row in self.hit_rows}.values(), key=canonical_bytes)},
            "editability-annotations.json": sorted(self.annotations.values(), key=lambda row: row["alignment_id"])}
        from .comparison_domain import assess_domain
        domain = assess_domain(self.prepared.study["config"], result)
        if domain is not None:
            result["search-audit.json"]["supported_domain"] = domain
        return result

    def variant_effects_and_carriage(self):
        """Assess counterparts independently of whether either sequence was nominated."""
        cache = {}
        def counterpart(item, locus, donor, copy_index):
            self.check_limits()
            key = item["design"]["record_id"], locus["record_id"], donor["record_id"] if donor else None, copy_index
            if key in cache:
                return cache[key]
            interval = {field: locus[field] for field in ("contig", "start", "end")}
            donor_id = donor["record_id"] if donor else None
            material = self.prepared.reconstruct(interval, donor_id, copy_index)
            status, reason = self.prepared.sequence_support(material, interval)
            path = oriented_path(material["sequence_path"], locus["strand"])
            variants, _, refs = self.prepared.variant_metadata(interval, donor_id, copy_index)
            if not path:
                cache[key] = None
                return None
            sequence, dc = self.sequence(item["design"], locus, "".join(point["base"] for point in path), path, variants,
                status, "intended" if locus["record_id"] in item["design"]["intended_locus_ids"] else "contextual")
            occurrence = None
            if donor:
                phase_valid, phase = self.prepared.sequence_support(material, interval)
                phase_keys = [site["phase_key"] for site in material["heterozygous_sites"]]
                block = known(canonical_bytes(phase_keys[0]).decode()) if phase_keys and all(k is not None and k == phase_keys[0] for k in phase_keys) else unknown("homozygous or single locally unordered heterozygous site" if phase_valid != "unresolved" else "unresolved phase or source sequence")
                occurrence = self._record("sequence_occurrence", sequence_configuration_id=sequence["record_id"], donor_id=donor_id,
                    sample_id=self.prepared.sample_id(donor_id), copy_index=copy_index, phase_block=block,
                    phase_status="unresolved" if phase_valid == "unresolved" else "unphased" if phase_keys and any(k is None for k in phase_keys) else "confirmed",
                    source_version=next(source["sha256"] for source in self.prepared.study["sources"].values() if source["role"] == "variants"), source_refs=refs)
            geometry = item["editor"]["guide_geometry"]
            alignment, admitted = None, False
            if status != "unresolved" and len(path) > geometry["pam_length"]:
                length = geometry["pam_length"]
                pam_path = path[-length:] if geometry["pam_side"] == "3prime" else path[:length]
                pam = "".join(point["base"] for point in pam_path)
                target = sequence["sequence"][:-length] if geometry["pam_side"] == "3prime" else sequence["sequence"][length:]
                alternatives = self.fixed_target_alignments(item, dc, target, pam_path, refs)
                if alternatives:
                    alignment, admitted = alternatives[0]
            cache[key] = sequence, dc, occurrence, status, reason, variants, refs, alignment, admitted
            return cache[key]

        for design_id, nominee_ids in sorted(self.nominees.items()):
            item = self.profiles[design_id]
            for nominee_id in sorted(nominee_ids):
                nominee = self.records[nominee_id]
                nominated_sequence = self.records[nominee["sequence_configuration_id"]]
                locus = self.records[nominated_sequence["locus_id"]]
                reference = counterpart(item, locus, None, 0)
                for donor in sorted(self.prepared.donor_records.values(), key=lambda row: row["record_id"]):
                    for copy_index in (0, 1):
                        current = counterpart(item, locus, donor, copy_index)
                        if current is None:
                            self._record("variant_effect", donor_id=donor["record_id"], donor_configuration_id=unknown("locus has no retained reconstructable bases"),
                                reference_configuration_id=known(reference[1]["record_id"]) if reference else unknown("reference unavailable"),
                                donor_alignment_id=unknown("deleted or unresolved target"), reference_alignment_id=known(reference[7]["record_id"]) if reference and reference[7] else unknown("no representable reference alignment"),
                                status="unknown", reason="No complete donor locus sequence could be reconstructed", source_variants=[])
                            continue
                        sequence, dc, occurrence, status, reason, variants, refs, alignment, admitted = current
                        carrier = sequence["record_id"] == nominated_sequence["record_id"]
                        carrier_status = "unresolved" if status == "unresolved" else status + ("_carrier" if carrier else "_noncarrier")
                        calls = [call for variant in self.prepared.variants_for({key: locus[key] for key in ("contig", "start", "end")})
                                 for label, call in variant["calls"].items() if label == self.prepared._legacy_by_canonical[donor["record_id"]]]
                        genotype = "retained_callset_model" if status == "modeled" else "filtered" if any(call["call_status"] == "filtered" for call in calls) else "missing" if any(call["call_status"] != "complete" for call in calls) else "called" if variants else "callable_reference"
                        self._record("donor_configuration", donor_id=donor["record_id"], sample_id=self.prepared.sample_id(donor["record_id"]),
                            sequence_occurrence_id=occurrence["record_id"], design_configuration_id=dc["record_id"], copy_index=copy_index,
                            carrier_status=carrier_status, carrier_of_sequence_id=nominated_sequence["record_id"], genotype_support=genotype,
                            phase_support=occurrence["phase_status"], reason_codes=[reason] if status in {"unresolved", "modeled"} else [], source_refs=refs)
                        if status == "unresolved" or reference is None or reference[3] == "unresolved":
                            effect, explanation = "unknown", reason
                        elif [{key: value for key, value in point.items() if key != "base"} for point in sequence["path"]] != [{key: value for key, value in point.items() if key != "base"} for point in reference[0]["path"]]:
                            effect, explanation = "unpaired", "Genomic insertion/deletion path has no exact reference-path counterpart"
                        elif sequence["sequence"] == reference[0]["sequence"]:
                            effect, explanation = "unchanged", "Same supported local DNA and path as the reference"
                        elif admitted and not reference[8]:
                            effect, explanation = "created", "Donor sequence admitted under this design profile; reference sequence excluded"
                        elif reference[8] and not admitted:
                            effect, explanation = "lost", "Reference sequence admitted under this design profile; donor sequence excluded"
                        else:
                            effect, explanation = "changed", "Supported local sequence differs under the same declared design profile"
                        self._record("variant_effect", donor_id=donor["record_id"], donor_configuration_id=known(dc["record_id"]),
                            reference_configuration_id=known(reference[1]["record_id"]) if reference else unknown("reference unavailable"),
                            donor_alignment_id=known(alignment["record_id"]) if alignment else unknown("no alignment within the declared bulge limits for this actual sequence"),
                            reference_alignment_id=known(reference[7]["record_id"]) if reference and reference[7] else unknown("no reference alignment within the declared bulge limits"),
                            status=effect, reason=explanation, source_variants=variants)

    def summaries(self):
        donors = sorted(self.prepared.donor_records.values(), key=lambda row: row["record_id"])
        for design_id in sorted(self.profiles):
            scopes = [row for row in self.records.values() if row["record_type"] == "comparison_search_scope" and row["design_id"] == design_id]
            assessments = [row for row in self.records.values() if row["record_type"] == "comparison_intended_assessment" and row["design_id"] == design_id]
            all_evaluable, all_modeled = set(), set()
            for donor in donors:
                local_assessments = [row for row in assessments if row["donor_id"] == donor["record_id"]]
                local_scopes = [row for row in scopes if row["donor_id"]["value"] == donor["record_id"]]
                reasons = sorted({reason for row in local_assessments + local_scopes for reason in row["reason_codes"]})
                intended_status = "unresolved" if any(row["support_status"] in {"unresolved", "not_evaluated"} for row in local_assessments) else "modeled" if any(row["support_status"] == "modeled" for row in local_assessments) else "supported"
                if intended_status != "unresolved":
                    all_evaluable.add(donor["record_id"])
                    if intended_status == "modeled":
                        all_modeled.add(donor["record_id"])
                complete = all(row["execution_status"] == "complete" for row in local_scopes)
                dc_ids = sorted(self.hit_donors[(design_id, donor["record_id"])])
                alignment_ids = sorted({row["alignment_id"] for row in self.hit_rows if row["design_id"] == design_id and row["donor_id"] == donor["record_id"]})
                self._record("donor_design_summary", donor_id=donor["record_id"], design_id=design_id,
                    domain_id=self.prepared.domain["record_id"], intended_assessment_ids=sorted(row["record_id"] for row in local_assessments),
                    search_scope_ids=sorted(row["record_id"] for row in local_scopes), intended_evaluability=intended_status,
                    search_status="complete" if complete else "partial", nominated_configuration_ids=dc_ids,
                    nominated_alignment_ids=alignment_ids, observed_configuration_count=len(dc_ids), observed_alignment_count=len(alignment_ids),
                    complete_configuration_count=known(len(dc_ids)) if complete else unknown("incomplete donor search scope"),
                    reason_codes=reasons, evidence_observation_ids=[], evidence_availability="not_assessed")
            nominees = [self.records[identifier] for identifier in self.nominees[design_id]]
            sequences = [self.records[row["sequence_configuration_id"]] for row in nominees]
            effects = [row for row in self.records.values() if row["record_type"] == "comparison_variant_effect" and any(row[field]["value"] is not None and self.records[row[field]["value"]]["design_id"] == design_id for field in ("donor_configuration_id", "reference_configuration_id"))]
            complete = all(row["execution_status"] == "complete" for row in scopes)
            eligible = {row["donor_id"]["value"] for row in scopes if row["execution_status"] == "complete"}
            supported_carriers, modeled_carriers = set(), set()
            for row in self.records.values():
                if row["record_type"] != "comparison_donor_configuration" or self.records[row["design_configuration_id"]]["design_id"] != design_id:
                    continue
                if row["carrier_status"] == "supported_carrier":
                    supported_carriers.add(row["donor_id"])
                elif row["carrier_status"] == "modeled_carrier":
                    modeled_carriers.add(row["donor_id"])
            locus_count, seq_count, alignment_count = len({row["locus_id"] for row in sequences}), len({row["record_id"] for row in sequences}), len(self.hit_alignments[design_id])
            self._record("design_summary", design_id=design_id, domain_id=self.prepared.domain["record_id"],
                search_scope_ids=sorted(row["record_id"] for row in scopes), donor_count=len(donors),
                intended_evaluable_donors=len(all_evaluable), intended_unresolved_donors=len(donors) - len(all_evaluable),
                search_evaluable_donors=len(eligible), search_unresolved_donors=len(donors) - len(eligible),
                observed_locus_count=locus_count, observed_sequence_configuration_count=seq_count, observed_alignment_count=alignment_count,
                complete_locus_count=known(locus_count) if complete else unknown("incomplete donor search scope"),
                complete_sequence_configuration_count=known(seq_count) if complete else unknown("incomplete donor search scope"),
                complete_alignment_count=known(alignment_count) if complete else unknown("incomplete donor search scope"),
                supported_carrier_donors=len(supported_carriers), modeled_carrier_donors=len(modeled_carriers),
                variant_effect_ids=sorted(row["record_id"] for row in effects), evidence_observation_ids=[],
                reason_codes=sorted({reason for row in scopes for reason in row["reason_codes"]} | {"empirical_effectiveness_not_assessed"}))
        changes = [row for row in self.records.values() if row["record_type"] == "comparison_variant_effect" and row["status"] in {"created", "lost", "changed"}]
        self._record("decision", design_ids=sorted(self.profiles), domain_id=self.prepared.domain["record_id"],
            state="sequence_based_hypothesis" if changes else "insufficient_comparable_evidence", observation_ids=[], prediction_ids=[],
            reason_codes=["sequence_changes_require_context_matched_experiments" if changes else "no_comparable_experimental_measurements", "nomination_counts_do_not_establish_a_preferred_design"],
            question="Which design should be tested next in matched donor experiments under the declared scope?")


def normalize_search(prepared, hits, execution, *, windows=None, host_budget=None):
    budget = host_budget or HostBudget(host_limits_from_study(prepared.study))
    for row in execution:
        for binding in row.get("native_manifest") or []:
            budget.retain_native_binding(binding)
    budget.check("normalization_source_verification")
    verify_sources(prepared.study)
    retained = []
    for hit in hits:
        budget.check("normalization_hit_materialization")
        retained.append(hit)
    hits = retained
    budget.check("normalization_hit_materialization")
    needed = {hit["window_id"] for hit in hits}
    supplied = windows if isinstance(windows, dict) else {} if windows is not None else None
    if windows is not None and not isinstance(windows, dict):
        for window in windows:
            budget.check("normalization_supplied_windows")
            supplied[window["window_id"]] = window
    windows = {}
    budget.check("normalization_source_reconstruction")
    for window in prepared.iter_windows():
        budget.check("normalization_source_reconstruction")
        if window["window_id"] in needed:
            windows[window["window_id"]] = window
    budget.check("normalization_source_reconstruction")
    if supplied is not None and any(identifier not in supplied or supplied[identifier] != window for identifier, window in windows.items()):
        raise PopulationError("Submitted native windows differ from their pinned source reconstruction")
    if needed - windows.keys():
        raise PopulationError("Native hit refers to a window absent from exact source reconstruction")
    # Own and check the windows after exhausting reconstruction. Retained source
    # aliases cannot then change the validated input while hits are mapped.
    for identifier, window in windows.items():
        budget.check("normalization_window_validation")
        windows[identifier] = validate_prepared_window(deepcopy(window))
    builder = ComparisonBuilder(prepared, check_limits=budget.check)
    budget.check("normalization_alignment_mapping")
    for hit in sorted(hits, key=canonical_bytes):
        budget.check()
        builder._map_hit(windows[hit["window_id"]], hit)
    budget.check("normalization_intended_and_counterpart_assessments")
    result = builder.finish(execution)
    budget.check("normalization_complete")
    for row in execution:
        row.setdefault("effective_limits", {})["host"] = asdict(budget.limits)
        row["host_resource_usage"] = budget.metrics()
    return result


def local_search(study, *, limits=None, host_budget=None):
    budget = _LocalBudget(limits or local_limits_from_study(study))
    host = host_budget or HostBudget(host_limits_from_study(study), started=budget.started, cpu_started=budget.cpu_started)
    def check():
        budget.check()
        host.check()
    check()
    prepared = prepare_inputs(study)
    host.phase = "local_source_reconstruction_and_search"
    check()
    for item in prepared.profiles:
        profile = item["profile"]
        if (profile["backend_name"], profile["backend_version"], profile["backend_digest"]) != (LOCAL_BACKEND_NAME, LOCAL_BACKEND_VERSION, LOCAL_BACKEND_DIGEST):
            raise PopulationError("Configured search profile does not identify this local backend contract")
    builder = ComparisonBuilder(prepared, check_limits=check)
    counts = Counter()
    evaluated, skipped = defaultdict(list), defaultdict(list)
    for window in prepared.iter_windows():
        check()
        for item in prepared.profiles:
            design_id = item["design"]["record_id"]
            if len(window["source_genotypes"]) > item["profile"]["max_window_variants"]:
                skipped[design_id].append({key: window[key] for key in ("window_id", "interval", "donor_id", "copy_index")} | {"reason": "max_window_variants_exceeded"})
                continue
            evaluated[design_id].append(window["window_id"])
            for hit in local_hits(window, item, check_limits=check):
                counts[item["design"]["record_id"]] += 1
                if counts[item["design"]["record_id"]] > item["profile"]["max_rows"]:
                    budget.fail("Local search row ceiling exceeded; complete results were not retained", "local_row_cap", budget.metrics())
                builder.add_hit(window, hit)
    execution = [{"design_id": item["design"]["record_id"], "search_profile_id": item["profile"]["record_id"],
                  "execution_status": "partial" if skipped[item["design"]["record_id"]] else "complete", "retention_status": "complete", "exit_code": 0,
                  "timeout_seconds": budget.limits.timeout_seconds,
                  "effective_limits": {"local": asdict(budget.limits), "host": asdict(host.limits)},
                  "evaluated_window_ids": sorted(evaluated[item["design"]["record_id"]]), "skipped_windows": skipped[item["design"]["record_id"]],
                  "reason_codes": ["max_window_variants_exceeded"] if skipped[item["design"]["record_id"]] else [], "source_equivalence": "not_applicable_synthetic" if study["config"]["mode"] == "synthetic_test" else "verified",
                  "native_manifest": None} for item in prepared.profiles]
    host.phase = "local_intended_and_counterpart_assessments"
    check()
    result = builder.finish(execution)
    check()
    for row in execution:
        row["local_resource_usage"] = budget.metrics()
        row["host_resource_usage"] = host.metrics()
    return result
