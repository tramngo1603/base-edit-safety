"""Candidate-panel v1 keeps nomination, measurement, and panel choice distinct."""
from dataclasses import replace
from collections import Counter
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from base_edit_safety.pipeline.candidate_panel import (
    CandidatePanelError,
    SearchProfile,
    build_ledger,
    compare_ledgers,
    load_measurements,
    load_nominations,
    main,
    merge_evidence,
    panel_rows,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "candidate-panel"
PROFILE = EXAMPLE / "search-profile.json"
GUIDES = EXAMPLE / "guides.tsv"
INTENDED = EXAMPLE / "intended-loci.tsv"
NOMINATIONS = EXAMPLE / "nominations.tsv"
MEASUREMENTS = EXAMPLE / "measurements.tsv"
DECISIONS = EXAMPLE / "panel-decisions.tsv"
PROVENANCE = EXAMPLE / "source-provenance.json"


def example_ledger() -> dict:
    return build_ledger(
        PROFILE, GUIDES, INTENDED, [NOMINATIONS], MEASUREMENTS, DECISIONS, [PROVENANCE]
    )


def rewrite_tsv(source: Path, destination: Path, mutate) -> Path:
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)
    assert fieldnames is not None
    mutate(rows)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return destination


def provenance_for(
    nomination_path: Path,
    destination: Path,
    intended_path: Path = INTENDED,
) -> Path:
    payload = json.loads(PROVENANCE.read_text())
    checksum = hashlib.sha256(nomination_path.read_bytes()).hexdigest()
    with nomination_path.open(encoding="utf-8", newline="") as handle:
        row_count = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
    payload["native_input"] = {"file": nomination_path.name, "sha256": checksum}
    payload["normalized_output"] = {
        "file": nomination_path.name,
        "sha256": checksum,
        "row_count": row_count,
    }
    payload["intended_locus_manifest"] = {
        "file": intended_path.name,
        "sha256": hashlib.sha256(intended_path.read_bytes()).hexdigest(),
    }
    destination.write_text(json.dumps(payload), encoding="utf-8")
    return destination


def candidate(ledger: dict, chrom: str, start: int) -> dict:
    return next(
        item
        for item in ledger["candidates"]
        if item["chrom"] == chrom and item["start_0based"] == start
    )


def test_example_merges_evidence_without_collapsing_evidence_types():
    ledger = example_ledger()

    assert ledger["summary"]["candidate_sites"] == 8
    assert ledger["summary"]["nomination_rows"] == 8
    assert ledger["summary"]["measurement_rows"] == 5
    assert ledger["summary"]["by_review_group"]["measurement_reconciliation"] == 3

    joined = candidate(ledger, "chrSYN2", 200)
    assert {item["nomination_kind"] for item in joined["nominations"]} == {
        "computational",
        "experimental",
    }
    assert joined["reported_measurement_statuses"] == ["edited"]
    assert joined["measurement_summary"] == "coverage_unresolved"
    assert joined["evidence_state"] == "measurement_coverage_unresolved"
    assert joined["measurements"][0]["editor"] is None
    assert joined["measurements"][0]["allele"] is None
    assert joined["panel_decision"]["status"] == "include"


def test_both_synthetic_hbg_copies_remain_intended_controls():
    ledger = example_ledger()
    intended = [item for item in ledger["candidates"] if item["site_role"] == "intended"]

    assert [(item["chrom"], item["start_0based"]) for item in intended] == [
        ("chrSYN1", 100),
        ("chrSYN1", 5024),
    ]
    assert all(item["review_group"] == "intended_control" for item in intended)
    requirement = ledger["locus_requirements"][0]
    assert requirement["paralog_spacing_bp"] == 4924
    assert {item["validates"] for item in requirement["required_confirmation"]} == {
        "hbg_paralog_rearrangement",
        "hbg_paralog_cotargeting",
    }


def test_unmeasured_experimental_candidate_is_not_called_an_off_target_finding():
    site = candidate(example_ledger(), "chrSYN5", 500)

    assert site["evidence_state"] == "experimentally_nominated_unmeasured"
    assert site["measurement_summary"] == "not_measured"
    assert site["review_group"] == "experimental_nomination"


def test_measurement_without_nomination_is_exposed_for_reconciliation():
    site = candidate(example_ledger(), "chrSYN7", 700)

    assert site["evidence_state"] == "measurement_without_nomination"
    assert site["reported_measurement_statuses"] == ["inconclusive"]
    assert site["measurement_summary"] == "coverage_unresolved"
    assert site["nominations"] == []


def test_panel_table_keeps_excluded_and_undecided_sites_with_reasons():
    rows = panel_rows(example_ledger()["candidates"])
    decisions = {row["panel_decision"] for row in rows}

    assert decisions == {"include", "exclude", "undecided"}
    assert all(row["decision_reason"] for row in rows)
    excluded = next(row for row in rows if row["panel_decision"] == "exclude")
    assert excluded["reported_measurement_statuses"] == "not_detected"
    assert excluded["measurement_summary"] == "coverage_unresolved"


def test_profile_rejects_rows_beyond_predeclared_limits():
    profile = SearchProfile.load(PROFILE)
    base = {
        "search_profile": profile.name,
        "genome_build": profile.reference_id,
        "mismatch_count": 3,
        "dna_bulge": 0,
        "rna_bulge": 0,
        "pam": "AGG",
        "population_variant_id": None,
    }

    with pytest.raises(CandidatePanelError, match="mismatch_count 7"):
        profile.validate_computational_nomination({**base, "mismatch_count": 7}, "row")
    with pytest.raises(CandidatePanelError, match="dna_bulge 2"):
        profile.validate_computational_nomination({**base, "dna_bulge": 2}, "row")
    with pytest.raises(CandidatePanelError, match="does not match"):
        profile.validate_computational_nomination({**base, "pam": "AAA"}, "row")


def test_nonvariant_profile_rejects_variant_annotated_computational_row():
    profile = replace(SearchProfile.load(PROFILE), variant_aware=False)
    nomination = {
        "search_profile": profile.name,
        "genome_build": profile.reference_id,
        "mismatch_count": 3,
        "dna_bulge": 0,
        "rna_bulge": 0,
        "pam": "TGG",
        "population_variant_id": "rsSYNTHETIC",
    }

    with pytest.raises(CandidatePanelError, match="not variant-aware"):
        profile.validate_computational_nomination(nomination, "row")


def test_profile_accepts_iupac_pam_and_declared_broad_candidate():
    profile = SearchProfile.load(PROFILE)
    nomination = {
        "search_profile": profile.name,
        "genome_build": profile.reference_id,
        "mismatch_count": 6,
        "dna_bulge": 1,
        "rna_bulge": 1,
        "pam": "TAG",
        "population_variant_id": "rsSYNTHETIC",
    }

    profile.validate_computational_nomination(nomination, "row")


def test_structured_retention_threshold_is_enforced():
    profile = replace(
        SearchProfile.load(PROFILE), retention_score_name="CFD", minimum_score=0.10
    )
    nomination = {
        "search_profile": profile.name,
        "genome_build": profile.reference_id,
        "mismatch_count": 4,
        "dna_bulge": 0,
        "rna_bulge": 0,
        "pam": "AGG",
        "population_variant_id": None,
        "score_name": "CFD",
        "score_value": 0.09,
    }

    with pytest.raises(CandidatePanelError, match="score_value must be at least 0.1"):
        profile.validate_computational_nomination(nomination, "row")


def test_unknown_panel_decision_cannot_be_silently_added():
    decision = {
        ("synthetic_HBG_guide", "SYNTHETIC", "chrUNKNOWN", 1, 24): {
            "status": "include",
            "reason": "unknown",
            "source_file": "decisions.tsv",
        }
    }

    with pytest.raises(CandidatePanelError, match="does not match a nomination or measurement"):
        merge_evidence([], [], decision)


def test_panel_decision_reason_and_note_contract_is_fixed(tmp_path):
    invalid_reason = rewrite_tsv(
        DECISIONS,
        tmp_path / "invalid-reason.tsv",
        lambda rows: rows[0].update(reason="invented_reason"),
    )
    with pytest.raises(CandidatePanelError, match="decision reason must be one of"):
        build_ledger(
            PROFILE,
            GUIDES,
            INTENDED,
            [NOMINATIONS],
            decision_path=invalid_reason,
            source_provenance_paths=[PROVENANCE],
        )

    long_note = rewrite_tsv(
        DECISIONS,
        tmp_path / "long-note.tsv",
        lambda rows: rows[0].update(note="x" * 141),
    )
    with pytest.raises(CandidatePanelError, match="at most 140 characters"):
        build_ledger(
            PROFILE,
            GUIDES,
            INTENDED,
            [NOMINATIONS],
            decision_path=long_note,
            source_provenance_paths=[PROVENANCE],
        )

    missing_related = rewrite_tsv(
        DECISIONS,
        tmp_path / "missing-related.tsv",
        lambda rows: rows[0].update(reason="redundant_with", related_site_id=""),
    )
    with pytest.raises(CandidatePanelError, match="requires related_site_id"):
        build_ledger(
            PROFILE,
            GUIDES,
            INTENDED,
            [NOMINATIONS],
            decision_path=missing_related,
            source_provenance_paths=[PROVENANCE],
        )


def test_provenance_uses_basenames_and_checksums_not_absolute_input_paths():
    ledger = example_ledger()
    rendered = json.dumps(ledger)

    assert str(ROOT) not in rendered
    assert [item["file"] for item in ledger["inputs"]] == [
        "search-profile.json",
        "guides.tsv",
        "intended-loci.tsv",
        "nominations.tsv",
        "measurements.tsv",
        "panel-decisions.tsv",
        "source-provenance.json",
    ]
    expected = hashlib.sha256(NOMINATIONS.read_bytes()).hexdigest()
    assert ledger["inputs"][3]["sha256"] == expected


def test_cli_writes_candidate_ledger_and_review_table(tmp_path, capsys):
    ledger_path = tmp_path / "candidate-ledger.json"
    panel_path = tmp_path / "panel-review.tsv"

    status = main(
        [
            "--profile",
            str(PROFILE),
            "--guides",
            str(GUIDES),
            "--intended-loci",
            str(INTENDED),
            "--nominations",
            str(NOMINATIONS),
            "--measurements",
            str(MEASUREMENTS),
            "--decisions",
            str(DECISIONS),
            "--source-provenance",
            str(PROVENANCE),
            "--output",
            str(ledger_path),
            "--panel-output",
            str(panel_path),
            "--panel-notes",
            "Depth targets are assay-specific.",
        ]
    )

    assert status == 0
    written = json.loads(ledger_path.read_text())
    assert written["summary"]["candidate_sites"] == 8
    assert written["export_metadata"]["panel_notes"] == "Depth targets are assay-specific."
    assert len(panel_path.read_text().splitlines()) == 9
    assert "wrote 8 candidate sites" in capsys.readouterr().out


def test_v3_diff_keeps_removed_sites_at_ledger_level():
    previous = example_ledger()
    current = deepcopy(previous)
    removed = current["candidates"].pop(0)
    added = deepcopy(current["candidates"][0])
    added["site_id"] = "site-added000000000"
    added["start_0based"] = 9999
    added["end_0based_exclusive"] = 10022
    current["candidates"].append(added)
    shared = current["candidates"][0]
    old_shared = next(
        site for site in previous["candidates"] if site["site_id"] == shared["site_id"]
    )
    old_shared["source_bucket"] = "crispor_only"
    shared["source_bucket"] = "both"

    diff = compare_ledgers(previous, current)

    assert diff["removed_site_ids"] == [removed["site_id"]]
    assert diff["added_site_ids"] == ["site-added000000000"]
    assert diff["rebucketed"] == [
        {
            "site_id": shared["site_id"],
            "old_source_bucket": "crispor_only",
            "new_source_bucket": "both",
        }
    ]
    assert all("diff_status" not in site for site in current["candidates"])


def test_intended_role_is_derived_from_manifest_not_trusted_from_source(tmp_path):
    bad = rewrite_tsv(
        NOMINATIONS,
        tmp_path / "mislabelled.tsv",
        lambda rows: rows[1].update(site_role="candidate"),
    )

    with pytest.raises(CandidatePanelError, match="must be 'intended'.*HBG2"):
        build_ledger(
            PROFILE,
            GUIDES,
            INTENDED,
            [bad],
            MEASUREMENTS,
            DECISIONS,
            [provenance_for(bad, tmp_path / "mislabelled-source.json")],
        )


def test_blank_site_role_is_derived_from_intended_manifest(tmp_path):
    blank = rewrite_tsv(
        NOMINATIONS,
        tmp_path / "blank-role.tsv",
        lambda rows: rows[1].update(site_role=""),
    )
    ledger = build_ledger(
        PROFILE,
        GUIDES,
        INTENDED,
        [blank],
        source_provenance_paths=[
            provenance_for(blank, tmp_path / "blank-role-source.json")
        ],
    )

    assert candidate(ledger, "chrSYN1", 5024)["site_role"] == "intended"


def test_hbg_manifest_requires_both_copies(tmp_path):
    incomplete = rewrite_tsv(
        INTENDED,
        tmp_path / "one-copy.tsv",
        lambda rows: rows.pop(),
    )

    with pytest.raises(CandidatePanelError, match="must declare both HBG1 and HBG2"):
        build_ledger(
            PROFILE,
            GUIDES,
            incomplete,
            [NOMINATIONS],
            MEASUREMENTS,
            None,
            [provenance_for(NOMINATIONS, tmp_path / "one-copy-source.json", incomplete)],
        )


def test_missing_intended_nomination_remains_visible_as_a_control(tmp_path):
    without_hbg2 = rewrite_tsv(
        NOMINATIONS,
        tmp_path / "without-hbg2.tsv",
        lambda rows: rows.pop(1),
    )
    ledger = build_ledger(
        PROFILE,
        GUIDES,
        INTENDED,
        [without_hbg2],
        source_provenance_paths=[
            provenance_for(without_hbg2, tmp_path / "without-hbg2-source.json")
        ],
    )
    hbg2 = candidate(ledger, "chrSYN1", 5024)

    assert hbg2["site_role"] == "intended"
    assert hbg2["evidence_state"] == "intended_control_not_nominated"
    assert hbg2["nominations"] == []


def test_standard_chromosome_aliases_merge_before_identity_matching(tmp_path):
    def nomination_aliases(rows):
        rows[2]["chrom"] = "1"
        rows[3]["chrom"] = "chr1"

    def measurement_alias(rows):
        rows[2]["chrom"] = "CHR1"

    nominations = rewrite_tsv(NOMINATIONS, tmp_path / "aliases.tsv", nomination_aliases)
    measurements = rewrite_tsv(MEASUREMENTS, tmp_path / "measure-alias.tsv", measurement_alias)
    ledger = build_ledger(
        PROFILE,
        GUIDES,
        INTENDED,
        [nominations],
        measurements,
        source_provenance_paths=[
            provenance_for(nominations, tmp_path / "aliases-source.json")
        ],
    )
    site = candidate(ledger, "chr1", 200)

    assert len(site["nominations"]) == 2
    assert len(site["measurements"]) == 1
    assert {item.get("source_chrom") for item in site["nominations"]} == {None, "1"}


def test_conflicting_known_strands_stop_the_merge(tmp_path):
    bad = rewrite_tsv(
        MEASUREMENTS,
        tmp_path / "strand-conflict.tsv",
        lambda rows: rows[2].update(strand="+"),
    )

    with pytest.raises(CandidatePanelError, match="conflicting strands"):
        build_ledger(
            PROFILE,
            GUIDES,
            INTENDED,
            [NOMINATIONS],
            bad,
            source_provenance_paths=[PROVENANCE],
        )


def test_one_base_coordinate_disagreements_are_flagged_not_merged():
    base = {
        "guide_id": "g",
        "genome_build": "GRCh38",
        "chrom": "chr1",
        "end_0based_exclusive": 123,
        "strand": "+",
        "site_role": "candidate",
        "nomination_kind": "computational",
        "source_method": "search",
        "source_record_id": "one",
        "mismatch_count": 2,
        "population_variant_id": None,
        "dna_bulge": 0,
        "rna_bulge": 0,
    }
    adjacent = {**base, "start_0based": 101, "end_0based_exclusive": 124,
                "source_record_id": "two"}
    sites = merge_evidence([{**base, "start_0based": 100}, adjacent], [], {})

    assert len(sites) == 2
    assert all(site["identity_warnings"] for site in sites)


def test_duplicate_nomination_source_ids_are_rejected(tmp_path):
    duplicate = rewrite_tsv(
        NOMINATIONS,
        tmp_path / "duplicate-id.tsv",
        lambda rows: rows[1].update(source_record_id=rows[0]["source_record_id"]),
    )

    with pytest.raises(CandidatePanelError, match="duplicate nomination source_record_id"):
        load_nominations([duplicate], SearchProfile.load(PROFILE))


def test_duplicate_measurement_source_ids_are_rejected(tmp_path):
    duplicate = rewrite_tsv(
        MEASUREMENTS,
        tmp_path / "duplicate-measurement-id.tsv",
        lambda rows: rows[1].update(source_record_id=rows[0]["source_record_id"]),
    )

    with pytest.raises(CandidatePanelError, match="duplicate measurement source_record_id"):
        load_measurements(duplicate)


def test_nonnumeric_detection_limit_is_rejected(tmp_path):
    bad = rewrite_tsv(
        MEASUREMENTS,
        tmp_path / "bad-lod.tsv",
        lambda rows: rows[3].update(limit_of_detection_value="unknown"),
    )

    with pytest.raises(CandidatePanelError, match="limit_of_detection_value must be numeric"):
        load_measurements(bad)


def test_sequence_derived_mismatch_assertion_is_checked(tmp_path):
    bad = rewrite_tsv(
        NOMINATIONS,
        tmp_path / "bad-mismatch.tsv",
        lambda rows: rows[4].update(mismatch_count="4"),
    )

    with pytest.raises(CandidatePanelError, match="sequence comparison gives 5"):
        build_ledger(
            PROFILE,
            GUIDES,
            INTENDED,
            [bad],
            source_provenance_paths=[
                provenance_for(bad, tmp_path / "bad-mismatch-source.json")
            ],
        )


def test_panel_rows_put_limitations_and_locus_assays_in_front_of_reviewers():
    rows = panel_rows(example_ledger()["candidates"])
    hbg1 = next(row for row in rows if row["site_role"] == "intended")

    assert "locus-resolving" in hbg1["required_confirmation"]
    assert hbg1["locus_requirements"]
    assert hbg1["interpretation_warning"] == "Panel exclusion is not evidence that a site is safe."
    assert all(
        not row["locus_requirements"] for row in rows if row["site_role"] == "candidate"
    )
    flagged = next(row for row in rows if row["ch_driver_genes"])
    assert "not independently verified" in flagged["annotation_warnings"]
    measured = next(row for row in rows if row["replicate_ids"])
    assert measured["measurement_qc_statuses"] in {"passed", "not_reported"}
    assert measured["read_depths"]


def test_missing_or_tampered_source_provenance_stops_the_run(tmp_path):
    with pytest.raises(CandidatePanelError, match="source-provenance"):
        build_ledger(PROFILE, GUIDES, INTENDED, [NOMINATIONS])

    payload = json.loads(PROVENANCE.read_text())
    payload["normalized_output"]["sha256"] = "0" * 64
    tampered = tmp_path / "tampered-source.json"
    tampered.write_text(json.dumps(payload))
    with pytest.raises(CandidatePanelError, match="checksum does not match"):
        build_ledger(
            PROFILE,
            GUIDES,
            INTENDED,
            [NOMINATIONS],
            source_provenance_paths=[tampered],
        )


def test_large_candidate_union_is_split_into_reviewable_mismatch_bands():
    nominations = []
    for index in range(10_000):
        mismatches = index % 6
        nominations.append(
            {
                "guide_id": "g",
                "genome_build": "GRCh38",
                "chrom": f"chr{1 + index % 22}",
                "start_0based": 1000 + index * 30,
                "end_0based_exclusive": 1023 + index * 30,
                "strand": "+",
                "site_role": "candidate",
                "nomination_kind": "computational",
                "source_method": "search",
                "source_record_id": f"row-{index}",
                "mismatch_count": mismatches,
                "population_variant_id": None,
                "ch_driver_gene": None,
                "dna_bulge": 0,
                "rna_bulge": 0,
            }
        )

    sites = merge_evidence(nominations, [], {})
    groups = Counter(site["review_group"] for site in sites)

    assert len(sites) == 10_000
    assert set(groups) == {
        "sequence_match_0_to_2_mismatches",
        "sequence_match_3_to_4_mismatches",
        "sequence_match_5_plus_mismatches",
    }
