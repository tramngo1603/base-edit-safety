"""Exact local admissions, independent intended assessments and unknown scope."""
from collections import defaultdict
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from base_edit_safety.population.common import PopulationError, ResourceLimit
from base_edit_safety.population.comparison import (
    ComparisonBuilder, HostBudget, LocalLimits, enumerate_alignments, local_hits, local_search, normalize_search,
)
from base_edit_safety.population.comparison_contracts import (
    make_record, record_identity, record_registry, validate_catalog, validate_study,
)
from base_edit_safety.population.comparison_inputs import (
    _merge_ranges, _phase_segments, _validate_regional_bundle, prepare_inputs, support_bounds,
)
from base_edit_safety.population.comparison_provenance import resolve_study
from base_edit_safety.population.identity import reverse_complement

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/comparison/synthetic"
EXPECTED = json.loads((EXAMPLE / "expected.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace(value, mapping):
    if isinstance(value, str):
        return mapping.get(value, value)
    if isinstance(value, list):
        return [replace(item, mapping) for item in value]
    if isinstance(value, dict):
        return {key: replace(item, mapping) for key, item in value.items()}
    return value


def rekey(config):
    config = deepcopy(config)
    for _ in range(32):
        changed = {row["record_id"]: record_identity(row) for row in config["catalogs"]
                   if row["record_id"] != record_identity(row)}
        if not changed:
            return validate_study(config)
        config = replace(config, changed)
        for row in config["catalogs"]:
            if row["record_type"] == "comparison_design":
                row["definition_hash"] = row["record_id"].rsplit(":", 1)[1]
    raise AssertionError("Test identity graph did not converge")


def case(tmp_path, *, change_config=None, change_vcf=None, change_reference=None, reverse_donors=False):
    directory = tmp_path / "case"
    shutil.copytree(EXAMPLE, directory)
    config = json.loads((directory / "study.json").read_text())
    if change_config:
        change_config(config)
    if change_vcf:
        path = directory / "variants.vcf"
        path.write_text(change_vcf(path.read_text()))
    if change_reference:
        path = directory / "reference.fa"
        previous = sha(path)
        path.write_text(change_reference(path.read_text()))
        current = sha(path)
        for row in config["catalogs"]:
            if row["record_type"] == "comparison_reference_bundle":
                row["complete_reference_sha256"]["value"] = current
                row["source_checksums"] = [current if value == previous else value for value in row["source_checksums"]]
            elif row["record_type"] == "comparison_source_artifact" and row["sha256"] == previous:
                row["sha256"] = current
    if reverse_donors:
        path = directory / "donors.tsv"
        rows = path.read_text().splitlines()
        path.write_text("\n".join([rows[0], *reversed(rows[1:])]) + "\n")
    legacy_path = directory / "genome-source.json"
    legacy = json.loads(legacy_path.read_text())
    for key, path in legacy["inputs"].items():
        if path:
            legacy["input_sources"][key]["sha256"] = sha(directory / path)
    legacy_path.write_text(json.dumps(legacy))
    for source in config["input_refs"]:
        source["sha256"] = sha(directory / source["path"])
    config = rekey(config)
    (directory / "study.json").write_text(json.dumps(config))
    return resolve_study(directory / "study.json")


@pytest.fixture(scope="module")
def result():
    study = resolve_study(EXAMPLE / "study.json")
    output = local_search(study)
    return study, output, validate_catalog(study["config"]["catalogs"] + output["comparison-records.jsonl"])


def records(output, kind):
    return [row for row in output["comparison-records.jsonl"] if row["record_type"] == "comparison_" + kind]


def nominated(study, output, index):
    labels = {row["record_id"]: row["canonical_label"] for row in study["config"]["catalogs"] if row["record_type"] == "comparison_donor"}
    designs = {row["record_id"]: row["display_name"].rsplit(" ", 1)[1] for row in study["config"]["catalogs"] if row["record_type"] == "comparison_design"}
    actual = defaultdict(set)
    for hit in output["search-audit.json"]["hits"]:
        dc = index[hit["design_configuration_id"]]
        sequence = index[dc["sequence_configuration_id"]]
        actual[labels.get(hit["donor_id"], "reference") + ":" + designs[hit["design_id"]]].add(sequence["contig"] + ":" + sequence["strand"])
    return {key: sorted(value) for key, value in actual.items()}


def test_each_design_uses_its_own_pam_and_matches_independent_expected_table(result):
    study, output, index = result
    assert nominated(study, output, index) == EXPECTED["nominated_locus_strand_by_donor_design"]
    assert len({row["guide_sequence"] for row in study["config"]["catalogs"] if row["record_type"] == "comparison_design"}) == 1
    assert len(records(output, "design_summary")) == 2
    for hit in output["search-audit.json"]["hits"]:
        design = index[hit["design_id"]]
        aln = index[hit["alignment_id"]]
        assert aln["pam_sequence"][-2:] == "GG" if design["display_name"].endswith("A") else aln["pam_sequence"][-2] == "G"


def test_aliases_are_traceable_without_duplicate_design_execution(tmp_path, result):
    def change(config):
        design = next(row for row in config["catalogs"] if row["record_type"] == "comparison_design")
        source = next(row for row in config["catalogs"] if row["record_type"] == "comparison_source_artifact")
        for label in ("Source spelling one", "Source spelling two"):
            config["catalogs"].append(make_record("design_alias", source_name=label, design_id=design["record_id"],
                source_location={"source_artifact_id": source["record_id"], "locator": {"kind": "file", "member": "artificial-aliases"}, "note": "Synthetic alias fixture"}, justification="The complete design definition is unchanged"))
    study = case(tmp_path, change_config=change)
    output = local_search(study)
    assert len(output["search-audit.json"]["execution"]) == 2
    assert {row["record_id"] for row in records(output, "alignment")} == {row["record_id"] for row in records(result[1], "alignment")}


def test_identical_hbg_dna_keeps_two_loci_and_true_reverse_paths(result):
    study, output, index = result
    label = next(row["record_id"] for row in study["config"]["catalogs"] if row.get("canonical_label") == "N")
    intended = [row for row in records(output, "intended_assessment") if row["donor_id"] == label]
    sequences = [index[row["sequence_configuration_id"]["value"]] for row in intended]
    assert len({row["sequence"] for row in sequences}) == 1
    assert len({row["locus_id"] for row in sequences}) == 2
    reverse = next(row for row in sequences if row["strand"] == "-")
    assert [point["pos0"] for point in reverse["path"]] == list(range(22, -1, -1))
    assert reverse["sequence"] == EXPECTED["guide"] + "TGG"


def test_intended_pam_loss_indel_and_unresolved_phase_survive_no_hit(result):
    study, output, index = result
    donor = {row["canonical_label"]: row["record_id"] for row in study["config"]["catalogs"] if row["record_type"] == "comparison_donor"}
    relevant = [row for row in records(output, "intended_assessment") if row["donor_id"] == donor["C"]]
    loss = next(row for row in relevant if index[row["locus_id"]]["gene"] == "HBG1" and index[row["design_id"]]["display_name"].endswith("A"))
    assert loss["support_status"] == "supported" and loss["pam_compatible"]["value"] is False
    assert loss["alignment_ids"] and "intended_alignment_outside_nomination_profile" in loss["reason_codes"]
    indels = [row for row in relevant if index[row["locus_id"]]["gene"] == "HBG2"]
    assert all(len(index[row["sequence_configuration_id"]["value"]]["sequence"]) == 24 and row["pam_compatible"]["value"] is True and not row["alignment_ids"] for row in indels)
    unresolved = [row for row in records(output, "intended_assessment") if row["donor_id"] == donor["U"] and index[row["locus_id"]]["gene"] == "HBG1"]
    assert unresolved and all(row["support_status"] == "unresolved" and row["sequence_configuration_id"]["value"] is None for row in unresolved)


def gapped_gain_case(tmp_path):
    def profiles(config):
        for row in config["catalogs"]:
            if row["record_type"] == "comparison_search_profile":
                row.update(max_dna_bulges=1, max_rna_bulges=1, max_total_edits=2)
    return case(tmp_path, change_config=profiles,
        change_reference=lambda text: text.replace(">chrSYNGAIN\nACGTT", ">chrSYNGAIN\nTTGTT"),
        change_vcf=lambda text: text.replace("chrSYNGAIN\t23", "chrSYNGAIN\t1\tgap_created\tT\tA\t.\tPASS\t.\tGT:PS\t1|1:1\t0|0:1\t0|0:1\nchrSYNGAIN\t23"))


def test_gap_created_nomination_uses_an_admitted_alignment(tmp_path):
    study = gapped_gain_case(tmp_path)
    output = local_search(study)
    index = validate_catalog(study["config"]["catalogs"] + output["comparison-records.jsonl"])
    admitted = {row["alignment_id"] for row in output["search-audit.json"]["hits"]}
    effects = []
    for effect in records(output, "variant_effect"):
        key = effect["donor_configuration_id"]["value"]
        if key is None or index[effect["donor_id"]]["canonical_label"] != "C":
            continue
        sequence = index[index[key]["sequence_configuration_id"]]
        locus = index[sequence["locus_id"]]
        if (sequence["contig"], sequence["strand"], locus["start"], locus["end"]) == ("chrSYNGAIN", "+", 0, 23):
            effects.append(effect)
    assert len(effects) == 2
    assert all(row["status"] == "created" for row in effects)
    for effect in effects:
        key = effect["donor_alignment_id"]["value"]
        assert key in admitted
        assert len(index[key]["dna_bulge_positions"]) == len(index[key]["rna_bulge_positions"]) == 1


def test_intended_inserted_sequence_retains_admitted_gap_alternatives(tmp_path):
    study = gapped_gain_case(tmp_path)
    output = local_search(study)
    index = validate_catalog(study["config"]["catalogs"] + output["comparison-records.jsonl"])
    intended = [row for row in records(output, "intended_assessment")
                if index[row["donor_id"]]["canonical_label"] == "C" and index[row["locus_id"]]["gene"] == "HBG2"]
    assert len(intended) == 2
    for assessment in intended:
        expected = {row["alignment_id"] for row in output["search-audit.json"]["hits"]
                    if row["donor_id"] == assessment["donor_id"] and row["design_id"] == assessment["design_id"]
                    and index[row["design_configuration_id"]]["sequence_configuration_id"] == assessment["sequence_configuration_id"]["value"]}
        assert expected and expected <= set(assessment["alignment_ids"])
        assert all(len(index[key]["dna_bulge_positions"]) == 1 for key in expected)


def test_joint_multiallelic_indexes_and_each_copy_are_preserved():
    prepared = prepare_inputs(resolve_study(EXAMPLE / "study.json"))
    windows = list(prepared.iter_windows())
    donor = {row["canonical_label"]: row["record_id"] for row in prepared.donor_records.values()}
    for label in ("C", "N"):
        actual = sorted((row for row in windows if row["donor_id"] == donor[label] and row["contig"] == "chrSYNMULT"), key=lambda row: row["copy_index"])
        assert [row["sequence"] for row in actual] == EXPECTED["joint_multiallelic"][label]
        for row in actual:
            call = row["source_genotypes"][0]
            assert call["original_alleles"] == ["TTTGTTTG", "T", "TTTG"]
            assert call["selected_indexes"] == EXPECTED["joint_multiallelic"]["ordered_alt_indexes"][label]
            assert row["source_variants"][0]["selected_alt_index"] == call["selected_indexes"][row["copy_index"]]
    assert prepared.source_equivalence_audit["verified"]


def test_missing_and_filtered_calls_preserve_unknown_not_reference(tmp_path, result):
    study, output, index = result
    unknown_donor = next(row["record_id"] for row in study["config"]["catalogs"] if row.get("canonical_label") == "U")
    summaries = [row for row in records(output, "donor_design_summary") if row["donor_id"] == unknown_donor]
    assert all(row["intended_evaluability"] == "unresolved" and row["search_status"] == "partial" and row["complete_configuration_count"]["value"] is None for row in summaries)
    def filter_pam(text):
        return "\n".join(line.replace("\tPASS\t", "\tLowQual\t") if "pam_loss" in line else line for line in text.splitlines()) + "\n"
    changed = case(tmp_path, change_vcf=filter_pam)
    filtered = local_search(changed)
    assessments = records(filtered, "intended_assessment")
    assert any(row["support_status"] == "unresolved" and "unresolved_variant_footprint" in row["reason_codes"] for row in assessments)


def test_unknown_window_is_an_annotation_unknown_and_never_a_search_filter(tmp_path, result):
    def change(config):
        for editor in config["catalogs"]:
            if editor["record_type"] == "comparison_editor_definition":
                editor["editing_window"] = {"value": None, "reason": "Biological window is not reported in this fixture"}
    study = case(tmp_path, change_config=change)
    output = local_search(study)
    index = validate_catalog(study["config"]["catalogs"] + output["comparison-records.jsonl"])
    assert nominated(study, output, index) == nominated(*result)
    assert all(row["within_declared_window"]["value"] is None for row in output["editability-annotations.json"])


def test_global_alignment_alternatives_are_retained_and_gap_events_are_explicit():
    profile = {"max_mismatches": 0, "max_dna_bulges": 1, "max_rna_bulges": 0, "max_total_edits": 1,
               "total_edit_accounting": "mismatch_plus_bulged_bases", "max_rows": 10000}
    assert set(enumerate_alignments("AAA", "AAAA", profile)) == {(guide, "AAAA") for guide in ("-AAA", "A-AA", "AA-A", "AAA-")}
    profile.update(max_dna_bulges=2, max_total_edits=1, total_edit_accounting="mismatch_plus_bulge_events")
    assert set(enumerate_alignments("AC", "ATTC", profile)) == {("A--C", "ATTC")}
    profile.update(max_dna_bulges=1, max_rna_bulges=1, max_total_edits=2, total_edit_accounting="mismatch_plus_bulged_bases")
    assert set(enumerate_alignments("AC", "CA", profile)) == {("-AC", "CA-"), ("AC-", "-CA")}


@pytest.mark.parametrize("pam_side,pam,strand", [("3prime", "TG", "+"), ("3prime", "TG", "-"), ("5prime", "TTC", "+"), ("5prime", "TTC", "-")])
def test_native_two_base_and_five_prime_geometry_maps_exact_indexes(pam_side, pam, strand):
    guide = "AGCA"
    oriented = guide + pam if pam_side == "3prime" else pam + guide
    sequence = oriented if strand == "+" else reverse_complement(oriented)
    item = {"design": {"record_id": "design", "guide_sequence": guide},
            "editor": {"guide_geometry": {"spacer_length": 4, "pam_length": len(pam), "pam_side": pam_side}},
            "profile": {"pam_patterns": [pam], "strands": [strand], "max_mismatches": 0, "max_dna_bulges": 0, "max_rna_bulges": 0,
                        "max_total_edits": 0, "max_rows": 100, "total_edit_accounting": "mismatch_plus_bulged_bases"}}
    hits = list(local_hits({"window_id": "window", "sequence": sequence}, item))
    assert len(hits) == 1 and hits[0]["guide_alignment"] == guide and hits[0]["pam_sequence"] == pam
    indexes = list(range(len(sequence))) if strand == "+" else list(range(len(sequence) - 1, -1, -1))
    assert hits[0]["target_indexes"] == (indexes[:4] if pam_side == "3prime" else indexes[len(pam):])


def test_skipped_native_windows_do_not_appear_evaluated_and_cannot_be_zero(result):
    study, output, _ = result
    prepared = prepare_inputs(study)
    windows = list(prepared.iter_windows())
    execution = deepcopy(output["search-audit.json"]["execution"])
    target = next(row for row in windows if row["donor_id"] is not None and row["contig"] == "chrSYNGAIN")
    for row in execution:
        row["evaluated_window_ids"] = [identifier for identifier in row["evaluated_window_ids"] if identifier != target["window_id"]]
        row["skipped_windows"] = [{key: target[key] for key in ("window_id", "interval", "donor_id", "copy_index")} | {"reason": "independent_test_backend_skip"}]
    mapped = normalize_search(prepared, [], execution, windows=windows)
    scopes = [row for row in records(mapped, "search_scope") if row["donor_id"]["value"] == target["donor_id"]]
    assert all(row["execution_status"] == "partial" and "independent_test_backend_skip" in row["reason_codes"] for row in scopes)
    summaries = [row for row in records(mapped, "donor_design_summary") if row["donor_id"] == target["donor_id"]]
    assert all(row["observed_configuration_count"] == 0 and row["complete_configuration_count"]["value"] is None for row in summaries)
    # Both copies must be skipped before their genomic span disappears entirely.
    for row in execution:
        row["evaluated_window_ids"] = [identifier for identifier in row["evaluated_window_ids"] if not any(w["window_id"] == identifier and w["donor_id"] == target["donor_id"] and w["contig"] == target["contig"] for w in windows)]
    mapped = normalize_search(prepared, [], execution)
    assert all(not any(span["contig"] == target["contig"] for span in row["evaluated_intervals"]) for row in records(mapped, "search_scope") if row["donor_id"]["value"] == target["donor_id"])


@pytest.mark.parametrize("change", [{"retention_status": "truncated"}, {"exit_code": 2}, {"source_equivalence": "candidate_only"}])
def test_truncation_nonzero_exit_and_nonequivalent_source_prevent_complete_claim(change, result):
    study, output, _ = result
    execution = deepcopy(output["search-audit.json"]["execution"])
    for row in execution:
        row.update(change)
    mapped = normalize_search(prepare_inputs(study), [], execution)
    assert all(row["execution_status"] != "complete" for row in records(mapped, "search_scope"))
    assert all(row["complete_configuration_count"]["value"] is None for row in records(mapped, "donor_design_summary"))


def test_high_density_cap_is_an_explicit_skip(tmp_path):
    def change(config):
        for row in config["catalogs"]:
            if row["record_type"] == "comparison_search_profile":
                row["max_window_variants"] = 1
    output = local_search(case(tmp_path, change_config=change))
    assert all(row["execution_status"] == "partial" for row in output["search-audit.json"]["execution"])
    assert any("max_window_variants_exceeded" in row["reason_codes"] for row in records(output, "search_scope"))
    assert records(output, "intended_assessment")


def test_reordering_source_rows_preserves_biological_results(tmp_path, result):
    def reverse_variants(text):
        rows = text.splitlines()
        return "\n".join([line for line in rows if line.startswith("#")] + list(reversed([line for line in rows if not line.startswith("#")]))) + "\n"
    study = case(tmp_path, change_vcf=reverse_variants, reverse_donors=True)
    output = local_search(study)
    for kind in ("sequence_configuration", "design_configuration", "alignment"):
        assert {row["record_id"] for row in records(output, kind)} == {row["record_id"] for row in records(result[1], kind)}
    index = validate_catalog(study["config"]["catalogs"] + output["comparison-records.jsonl"])
    assert nominated(study, output, index) == nominated(*result)


def test_nested_paths_and_wrong_backend_cannot_bypass_declared_sources(tmp_path):
    study = case(tmp_path)
    path = tmp_path / "case/genome-source.json"
    value = json.loads(path.read_text())
    shutil.copy(path.parent / "variants.vcf", path.parent / "unbound.vcf")
    value["inputs"]["variants"] = "unbound.vcf"
    path.write_text(json.dumps(value))
    config = study["config"]
    next(row for row in config["input_refs"] if row["input_id"] == "genome_source")["sha256"] = sha(path)
    (path.parent / "study.json").write_text(json.dumps(config))
    with pytest.raises(PopulationError, match="not explicitly checksum-bound"):
        prepare_inputs(resolve_study(path.parent / "study.json"))
    original = json.loads((EXAMPLE / "study.json").read_text())
    for row in original["catalogs"]:
        if row["record_type"] == "comparison_search_profile":
            row["backend_name"] = "another_engine"
    (path.parent / "study.json").write_text(json.dumps(rekey(original)))
    shutil.copy(EXAMPLE / "genome-source.json", path)
    with pytest.raises(PopulationError, match="local backend"):
        local_search(resolve_study(path.parent / "study.json"))


def test_runtime_and_memory_are_enforced_and_measured(monkeypatch, result):
    study = result[0]
    usage = result[1]["search-audit.json"]["execution"][0]["local_resource_usage"]
    assert usage["timeout_seconds"] == 600 and usage["wall_time_seconds"] > 0 and usage["cpu_time_seconds"] > 0 and usage["peak_resident_bytes"] > 0
    with pytest.raises(ResourceLimit, match="memory ceiling"):
        local_search(study, limits=LocalLimits(max_resident_bytes=1))
    ticks = iter([0.0, 2.0])
    monkeypatch.setattr("base_edit_safety.population.comparison.time.monotonic", lambda: next(ticks))
    with pytest.raises(ResourceLimit, match="runtime ceiling"):
        local_search(study, limits=LocalLimits(timeout_seconds=1))


def test_inserted_bases_require_their_reference_anchor_and_adjacent_masks_merge():
    path = [{"kind": "insertion", "boundary0": 9, "offset0": 0, "base": "A"}, {"kind": "reference", "pos0": 9, "base": "C"}]
    variant = {"pos0": 8, "reference": "G", "alternates": ["GA"], "selected_alt_index": 1}
    assert support_bounds(path, [variant]) == (8, 10)
    assert _merge_ranges([(0, 10), (10, 23), (30, 35)]) == [(0, 23), (30, 35)]


def test_reference_bundle_replays_exact_source_sequence_not_only_matching_intervals(tmp_path):
    source = tmp_path / "source.json"
    source.write_text(json.dumps({"genome": "hg38", "chrom": "chr11", "start": 10, "end": 13, "dna": "acg"}))
    aggregate = tmp_path / "aggregate.json"
    aggregate.write_text(json.dumps({"windows": [{"contig": "chr11", "start": 10, "end": 13, "sequence": "ACG"}]}))
    study = {"sources": {"window": {"sha256": sha(source), "role": "reference", "path": str(source), "format": "json"}}}
    bundle = {"identity_kind": "indexed_windows", "assembly": "GRCh38", "windows": [{"interval": {"contig": "chr11", "start": 10, "end": 13}, "source_sha256": sha(source)}]}
    _validate_regional_bundle(study, aggregate, bundle)
    aggregate.write_text(aggregate.read_text().replace("ACG", "ACT"))
    with pytest.raises(PopulationError, match="pinned source bytes"):
        _validate_regional_bundle(study, aggregate, bundle)


def test_sequence_summary_never_infers_a_design_winner_or_empirical_zero(result):
    _, output, _ = result
    assert all(row["evidence_observation_ids"] == [] and row["evidence_availability"] == "not_assessed" for row in records(output, "donor_design_summary"))
    decisions = records(output, "decision")
    assert len(decisions) == 1 and decisions[0]["state"] == "sequence_based_hypothesis"
    assert decisions[0]["observation_ids"] == [] and decisions[0]["prediction_ids"] == []
    assert "nomination_counts_do_not_establish_a_preferred_design" in decisions[0]["reason_codes"]


def _mapping_case():
    prepared = prepare_inputs(resolve_study(EXAMPLE / "study.json"))
    windows = list(prepared.iter_windows())
    for window in windows:
        if window["donor_id"] is None:
            continue
        for item in prepared.profiles:
            hit = next(iter(local_hits(window, item)), None)
            if hit is not None:
                return prepared, windows, window, hit
    raise AssertionError("The fixture must contain a nominated hit")


def test_batch_mapping_rejects_window_changed_after_it_was_yielded(monkeypatch):
    prepared, windows, target, hit = _mapping_case()

    def changed_windows():
        yield from windows
        target["source_version"] = "e" * 64

    monkeypatch.setattr(prepared, "iter_windows", changed_windows)
    with pytest.raises(PopulationError, match="immutable content identity"):
        normalize_search(prepared, [hit], [])


@pytest.mark.parametrize("field,reason", [
    ("guide_alignment", "guide/strand"),
    ("pam_sequence", "bases do not match"),
    ("target_indexes", "indexes lie outside"),
])
def test_batch_mapping_checks_each_hit_against_its_window(field, reason):
    prepared, _, _, hit = _mapping_case()
    changed = deepcopy(hit)
    if field == "target_indexes":
        changed[field][0] = True
    else:
        changed[field] = "N" * len(changed[field])
    with pytest.raises(PopulationError, match=reason):
        normalize_search(prepared, [hit, changed], [])


def test_individual_hit_mapping_rechecks_changed_window():
    prepared, _, window, hit = _mapping_case()
    builder = ComparisonBuilder(prepared)
    builder.add_hit(window, hit)
    window["source_version"] = "e" * 64
    with pytest.raises(PopulationError, match="immutable content identity"):
        builder.add_hit(window, hit)


def test_batch_mapping_owns_validated_windows_during_resource_callbacks(monkeypatch, result):
    prepared, windows, target, hit = _mapping_case()
    original_version = target["source_version"]
    changed_version = "e" * 64
    assert original_version != changed_version
    monkeypatch.setattr(prepared, "iter_windows", lambda: iter(windows))

    class RetainedWindowBudget(HostBudget):
        def check(self, phase=None):
            super().check(phase)
            if phase == "normalization_alignment_mapping":
                target["source_version"] = changed_version

    execution = deepcopy(result[1]["search-audit.json"]["execution"])
    output = normalize_search(prepared, [hit], execution, host_budget=RetainedWindowBudget())
    occurrences = records(output, "sequence_occurrence")
    assert target["source_version"] == changed_version
    assert any(row["source_version"] == original_version for row in occurrences)
    assert all(row["source_version"] != changed_version for row in occurrences)


def test_repeated_record_construction_keeps_integer_contracts():
    prepared, _, window, hit = _mapping_case()
    builder = ComparisonBuilder(prepared)
    builder.add_hit(window, hit)
    occurrence = next(row for row in builder.records.values()
                      if row["record_type"] == "comparison_sequence_occurrence")
    values = {key: deepcopy(value) for key, value in occurrence.items()
              if key not in {"record_id", "record_type", "schema_version", "sequence_occurrence_id"}}
    builder._record("sequence_occurrence", **values)
    values["copy_index"] = float(values["copy_index"])
    with pytest.raises(PopulationError, match="incorrect type"):
        builder._record("sequence_occurrence", **values)


def test_repeated_record_construction_keeps_literal_nonidentity_text():
    prepared, _, _, _ = _mapping_case()
    builder = ComparisonBuilder(prepared)
    values = {"reference_bundle_id": prepared.reference_bundle_id, "contig": "chrSYNGAIN",
              "start": 10000, "end": 10023, "strand": "+", "locus_name": "café",
              "gene": {"value": None, "reason": "No gene annotation supplied"}, "path": []}
    builder._record("locus", **values)
    values["locus_name"] = "cafe\u0301"
    with pytest.raises(PopulationError, match="Conflicting comparative records"):
        builder._record("locus", **values)


def test_returned_record_mutation_does_not_change_validated_template():
    prepared, _, _, _ = _mapping_case()
    builder = ComparisonBuilder(prepared)
    values = {"reference_bundle_id": prepared.reference_bundle_id, "contig": "chrSYNGAIN",
              "start": 10000, "end": 10023, "strand": "+", "locus_name": "synthetic locus",
              "gene": {"value": None, "reason": "No gene annotation supplied"}, "path": []}
    row = builder._record("locus", **values)
    row["record_id"] = "comparison_locus:v1:" + "f" * 64
    with pytest.raises(PopulationError, match="Conflicting comparative records"):
        builder._record("locus", **values)
