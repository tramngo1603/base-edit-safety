"""Nomination readers preserve rows, identity and source-run details."""
import json
from pathlib import Path

import pytest

from base_edit_safety.pipeline.candidate_panel import (
    CandidatePanelError,
    SearchProfile,
    build_ledger,
    load_guides,
    load_intended_loci,
    validate_reference_sequences,
)
from base_edit_safety.pipeline.normalize_nominations import (
    load_query_map,
    main,
    normalize_cas_offinder,
    normalize_cas_offinder_2_4_1,
    normalize_experimental,
    write_normalized,
    write_provenance,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "candidate-panel"
PROFILE_PATH = EXAMPLE / "search-profile.json"
GUIDES_PATH = EXAMPLE / "guides.tsv"
INTENDED_PATH = EXAMPLE / "intended-loci.tsv"


def manifests():
    return (
        load_guides(GUIDES_PATH),
        load_intended_loci(INTENDED_PATH),
        SearchProfile.load(PROFILE_PATH),
    )


def test_cas_offinder_keeps_unintended_perfect_matches(tmp_path):
    native = tmp_path / "cas-offinder.txt"
    query = "ACGTACGTACGTACGTACGTNNN"
    sequence = "ACGTACGTACGTACGTACGTAGG"
    native.write_text(
        "\n".join(
            [
                f"synthetic_HBG_guide\tX\t{query}\t{sequence}\tchrSYN1\t100\t+\t0\t0",
                f"synthetic_HBG_guide\tX\t{query}\t{sequence}\tchr3\t900\t+\t0\t0",
            ]
        )
        + "\n"
    )
    guides, intended, profile = manifests()

    rows = normalize_cas_offinder(native, guides, intended, profile, "Cas-OFFinder test")

    assert [row["site_role"] for row in rows] == ["intended", "candidate"]
    assert rows[1]["mismatch_count"] == 0
    assert rows[1]["end_0based_exclusive"] == 923
    assert rows[1]["source_record_id"].startswith("casoff:synthetic_HBG_guide:")


def test_cas_offinder_adapter_rejects_exact_duplicate_native_rows(tmp_path):
    native = tmp_path / "duplicate.txt"
    row = (
        "synthetic_HBG_guide\tX\tACGTACGTACGTACGTACGTNNN\t"
        "ACGTACGTACGTACGTACGTAGG\tchr3\t900\t+\t0\t0"
    )
    native.write_text(row + "\n" + row + "\n")
    guides, intended, profile = manifests()

    with pytest.raises(CandidatePanelError, match="exact duplicate Cas-OFFinder row"):
        normalize_cas_offinder(native, guides, intended, profile, "Cas-OFFinder test")


def test_change_seq_be_adapter_preserves_read_and_control_metadata(tmp_path):
    native = tmp_path / "change-be.txt"
    native.write_text(
        "\t".join(
            [
                "#Chromosome",
                "Start",
                "End",
                "Genomic Coordinate",
                "Nuclease_Read_Count",
                "Strand",
                "Control_Read_Count",
                "Site_Sequence",
                "Site_Substitution_Number",
                "Nuclease_overlap_bp_list",
                "Control_overlap_bp_list",
                "File_Name",
                "Cell",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "chrSYN1",
                "5024",
                "5047",
                "synthetic-site",
                "42",
                "+",
                "2",
                "ACGTACGTACGTACGTACGTAGG",
                "0",
                "4,3",
                "0,0",
                "synthetic-replicate",
                "synthetic cells",
            ]
        )
        + "\n"
    )
    guides, intended, _ = manifests()

    rows = normalize_experimental(
        native,
        "change-seq-be",
        "synthetic_HBG_guide",
        guides,
        intended,
        "SYNTHETIC",
        "CHANGE-seq-BE",
        "synthetic source",
    )

    assert rows[0]["site_role"] == "intended"
    assert rows[0]["score_name"] == "Nuclease_Read_Count"
    assert rows[0]["score_value"] == 42.0
    assert '"control_read_count":"2"' in rows[0]["source_metadata_json"]
    assert '"replicate_file":"synthetic-replicate"' in rows[0]["source_metadata_json"]


def test_bed6_adapter_uses_bed_coordinates_without_guessing(tmp_path):
    native = tmp_path / "circle-sites.bed"
    native.write_text("chr3\t900\t923\texperimental-site\t17\t-\n")
    guides, intended, _ = manifests()

    rows = normalize_experimental(
        native,
        "bed6",
        "synthetic_HBG_guide",
        guides,
        intended,
        "SYNTHETIC",
        "CIRCLE-seq",
        "synthetic source",
    )

    assert rows[0]["start_0based"] == 900
    assert rows[0]["end_0based_exclusive"] == 923
    assert rows[0]["site_role"] == "candidate"
    assert rows[0]["score_name"] == "BED_score"


def test_adapter_provenance_is_verified_by_candidate_ledger(tmp_path):
    native = tmp_path / "cas-offinder.txt"
    native.write_text(
        "synthetic_HBG_guide\tX\tACGTACGTACGTACGTACGTNNN\t"
        "ACGTACGTACGTACGTACGTAGG\tchr3\t900\t+\t0\t0\n"
    )
    guides, intended, profile = manifests()
    rows = normalize_cas_offinder(native, guides, intended, profile, "Cas-OFFinder test")
    normalized = tmp_path / "normalized.tsv"
    provenance = tmp_path / "source.json"
    reference = tmp_path / "synthetic.fa"
    target = "ACGTACGTACGTACGTACGTAGG"
    reference.write_text(">chr3\n" + "C" * 900 + target + "C" * 77 + "\n")
    Path(str(reference) + ".fai").write_text("chr3\t1000\t6\t1000\t1001\n")
    write_normalized(rows, normalized)
    write_provenance(
        provenance,
        normalized,
        native,
        "cas-offinder",
        "Cas-OFFinder",
        "Cas-OFFinder test",
        "cas-offinder input.txt G output.txt",
        "SYNTHETIC",
        reference,
        GUIDES_PATH,
        INTENDED_PATH,
        PROFILE_PATH,
        len(rows),
    )

    ledger = build_ledger(
        PROFILE_PATH,
        GUIDES_PATH,
        INTENDED_PATH,
        [normalized],
        source_provenance_paths=[provenance],
        reference_fasta_path=reference,
    )

    site = next(item for item in ledger["candidates"] if item["chrom"] == "chr3")
    assert site["site_role"] == "candidate"
    assert site["nominations"][0]["mismatch_count"] == 0
    source = ledger["nomination_source_provenance"][0]
    assert source["native_input"]["sha256"]
    assert source["source_run"]["command"] == "cas-offinder input.txt G output.txt"


def test_reference_validator_checks_sequence_and_reverse_strand(tmp_path):
    reference = tmp_path / "tiny.fa"
    reference.write_text(">1\nAACCGGTT\n")
    Path(str(reference) + ".fai").write_text("1\t8\t3\t8\t9\n")
    base = {
        "guide_id": "g",
        "genome_build": "GRCh38",
        "chrom": "chr1",
        "start_0based": 1,
        "end_0based_exclusive": 5,
        "source_record_id": "row",
    }
    plus = {**base, "strand": "+", "site_sequence": "ACCG"}
    minus = {**base, "strand": "-", "site_sequence": "CGGT", "source_record_id": "reverse"}

    validate_reference_sequences([plus, minus], reference)

    assert plus["reference_validation"]["status"] == "verified"
    assert minus["reference_validation"]["status"] == "verified"

    with pytest.raises(CandidatePanelError, match="does not match the declared reference"):
        validate_reference_sequences([{**plus, "site_sequence": "AAAA"}], reference)


def test_cas_offinder_cli_writes_normalized_rows_after_reference_validation(tmp_path, capsys):
    native = tmp_path / "native.tsv"
    native.write_text(
        "synthetic_HBG_guide\tX\tACGTACGTACGTACGTACGTNNN\t"
        "ACGTACGTACGTACGTACGTAGG\tchr3\t0\t+\t0\t0\n"
    )
    reference = tmp_path / "tiny.fa"
    reference.write_text(">chr3\nACGTACGTACGTACGTACGTAGG\n")
    Path(str(reference) + ".fai").write_text("chr3\t23\t6\t23\t24\n")
    output = tmp_path / "normalized.tsv"
    provenance = tmp_path / "source.json"

    status = main(
        [
            "cas-offinder",
            "--input",
            str(native),
            "--profile",
            str(PROFILE_PATH),
            "--guides",
            str(GUIDES_PATH),
            "--intended-loci",
            str(INTENDED_PATH),
            "--source-version",
            "synthetic test",
            "--source-command",
            "cas-offinder input.txt G output.tsv",
            "--reference-fasta",
            str(reference),
            "--output",
            str(output),
            "--provenance-output",
            str(provenance),
        ]
    )

    assert status == 0
    assert len(output.read_text().splitlines()) == 2
    assert provenance.exists()
    assert "wrote 1 normalized nominations" in capsys.readouterr().out


# Mock six- and seven-column Cas-OFFinder 2.4.1 inputs.
FIXTURE_241 = ROOT / "tests" / "fixtures" / "cas-offinder-2.4.1"


def fixture_guides(name="guides.tsv"):
    return load_guides(FIXTURE_241 / name)


def fixture_loci():
    return load_intended_loci(FIXTURE_241 / "intended-loci.tsv")


def test_cas_offinder_2_4_1_ids_normalizes_both_strands_with_provenance():
    _, _, profile = manifests()

    rows, accounting = normalize_cas_offinder_2_4_1(
        FIXTURE_241 / "valid-seven-field.tsv",
        fixture_guides(),
        fixture_loci(),
        profile,
        "Cas-OFFinder 2.4.1 (9816b94)",
        with_ids=True,
    )

    assert [row["guide_id"] for row in rows] == ["alpha_guide", "beta_guide", "alpha_guide"]
    assert [row["site_role"] for row in rows] == ["intended", "candidate", "candidate"]
    assert [row["strand"] for row in rows] == ["+", "-", "+"]
    assert [row["mismatch_count"] for row in rows] == [0, 1, 2]
    minus = rows[1]
    assert (minus["chrom"], minus["start_0based"], minus["end_0based_exclusive"]) == (
        "chrSYN1",
        200,
        223,
    )
    assert minus["site_sequence"] == "TTGCAAGTCCGTAGGAATCATGG"
    assert minus["pam"] == "TGG"
    assert rows[0]["end_0based_exclusive"] == 123
    assert rows[2]["mismatch_count"] == 2
    assert all(row["source_record_id"].startswith("casoff241:") for row in rows)
    assert all(row["dna_bulge"] == 0 and row["rna_bulge"] == 0 for row in rows)
    metadata = rows[1]["source_metadata_json"]
    assert '"native_format":"cas-offinder-2.4.1-ids"' in metadata
    assert '"native_query_id":"beta_guide"' in metadata
    assert accounting == {
        "declared_format": "cas-offinder-2.4.1-ids",
        "lines_read": 3,
        "blank_or_comment_lines": 0,
        "rows_normalized": 3,
        "query_binding": "native_query_id",
    }


def test_cas_offinder_2_4_1_six_field_binds_through_declared_query_map():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    rows, accounting = normalize_cas_offinder_2_4_1(
        FIXTURE_241 / "valid-six-field.tsv",
        guides,
        fixture_loci(),
        profile,
        "Cas-OFFinder 2.4.1 (9816b94)",
        with_ids=False,
        query_map=query_map,
    )

    assert [row["guide_id"] for row in rows] == ["alpha_guide", "beta_guide"]
    assert [row["site_role"] for row in rows] == ["intended", "candidate"]
    assert '"native_format":"cas-offinder-2.4.1"' in rows[0]["source_metadata_json"]
    assert "native_query_id" not in rows[0]["source_metadata_json"]
    assert accounting["query_binding"] == "query_map"
    assert accounting["rows_normalized"] == 2


def test_cas_offinder_2_4_1_six_field_accepts_single_declared_guide_override():
    _, _, profile = manifests()

    rows, accounting = normalize_cas_offinder_2_4_1(
        FIXTURE_241 / "valid-six-field-alpha.tsv",
        fixture_guides(),
        fixture_loci(),
        profile,
        "Cas-OFFinder 2.4.1 (9816b94)",
        with_ids=False,
        guide_override="alpha_guide",
    )

    assert [row["guide_id"] for row in rows] == ["alpha_guide"]
    assert accounting["query_binding"] == "guide_override"


def test_cas_offinder_2_4_1_six_field_requires_declared_binding():
    _, _, profile = manifests()

    with pytest.raises(CandidatePanelError, match="requires a declared"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "valid-six-field.tsv",
            fixture_guides(),
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
        )


def test_cas_offinder_2_4_1_shared_query_binds_only_by_declaration():
    _, _, profile = manifests()
    guides = fixture_guides("shared-guides.tsv")

    with pytest.raises(CandidatePanelError, match="requires a declared"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "valid-six-field-beta.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
        )

    query_map = load_query_map(FIXTURE_241 / "query-map-shared.json", guides)
    rows, _ = normalize_cas_offinder_2_4_1(
        FIXTURE_241 / "valid-six-field-beta.tsv",
        guides,
        fixture_loci(),
        profile,
        "Cas-OFFinder 2.4.1 (9816b94)",
        with_ids=False,
        query_map=query_map,
    )
    assert [row["guide_id"] for row in rows] == ["beta_guide_two"]


def test_cas_offinder_2_4_1_query_map_rejects_wrong_guide_declaration():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map-wrong-guide.json", guides)

    with pytest.raises(CandidatePanelError, match="query does not match guide 'beta_guide'"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "valid-six-field-alpha.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_query_map_rejects_unknown_guide():
    with pytest.raises(CandidatePanelError, match="query map guide 'gamma_guide' is absent"):
        load_query_map(FIXTURE_241 / "query-map-unknown-guide.json", fixture_guides())


def test_cas_offinder_2_4_1_rejects_query_absent_from_map():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map-alpha-only.json", guides)

    with pytest.raises(CandidatePanelError, match="absent from the declared query map"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "valid-six-field-beta.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_rejects_conflicting_bindings():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match="conflicting query bindings"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "valid-six-field.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
            guide_override="alpha_guide",
        )


def test_cas_offinder_2_4_1_rejects_mixed_width_rows():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match="declared 6-column.*row has 7 fields"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "mixed-width.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_rejects_nine_column_content_under_declared_version():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match="declared 6-column.*row has 9 fields"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "nine-column-row.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_nine_column_path_still_rejects_six_field_content():
    guides, intended, profile = manifests()

    with pytest.raises(CandidatePanelError, match="requires exactly 9 columns"):
        normalize_cas_offinder(
            FIXTURE_241 / "valid-six-field.tsv",
            guides,
            intended,
            profile,
            "Cas-OFFinder test",
        )


@pytest.mark.parametrize(
    ("fixture", "message"),
    [
        ("malformed-width.tsv", "row has 5 fields"),
        ("empty-field.tsv", "empty field"),
        ("bad-strand.tsv", "strand must be"),
        ("bad-coordinate.tsv", "non-negative"),
        ("bad-mismatch-text.tsv", "must be integers"),
        ("bad-mismatch-count.tsv", "declared mismatch count 3 does not match the 2"),
        ("bad-geometry.tsv", "declared spacer-plus-PAM geometry"),
        ("bad-pam.tsv", "observed PAM 'TTT' is outside the search profile"),
        ("duplicate-rows.tsv", "exact duplicate Cas-OFFinder row"),
    ],
)
def test_cas_offinder_2_4_1_rejects_invalid_rows(fixture, message):
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match=message):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / fixture,
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_rejects_unknown_query_id():
    _, _, profile = manifests()

    with pytest.raises(CandidatePanelError, match="ID 'gamma_guide' is absent from guides"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "unknown-id-seven-field.tsv",
            fixture_guides(),
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=True,
        )


@pytest.mark.parametrize("fixture", ["empty.tsv", "comment-only.tsv"])
def test_cas_offinder_2_4_1_empty_output_cannot_establish_completeness(fixture):
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match="cannot establish search completeness"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / fixture,
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_cli_writes_rows_and_row_accounting(tmp_path, capsys):
    output = tmp_path / "normalized.tsv"
    provenance = tmp_path / "source.json"

    status = main(
        [
            "cas-offinder",
            "--format",
            "2.4.1-ids",
            "--input",
            str(FIXTURE_241 / "valid-seven-field.tsv"),
            "--profile",
            str(PROFILE_PATH),
            "--guides",
            str(FIXTURE_241 / "guides.tsv"),
            "--intended-loci",
            str(FIXTURE_241 / "intended-loci.tsv"),
            "--source-version",
            "Cas-OFFinder 2.4.1 (9816b94c20c4cba2e79b039e1e2a6dee684b7b66)",
            "--source-command",
            "not-applicable: mock parser input",
            "--reference-fasta",
            str(FIXTURE_241 / "reference.fa"),
            "--output",
            str(output),
            "--provenance-output",
            str(provenance),
        ]
    )

    assert status == 0
    assert len(output.read_text().splitlines()) == 4
    payload = json.loads(provenance.read_text())
    assert payload["adapter"]["mode"] == "cas-offinder-2.4.1-ids"
    assert payload["row_accounting"] == {
        "declared_format": "cas-offinder-2.4.1-ids",
        "lines_read": 3,
        "blank_or_comment_lines": 0,
        "rows_normalized": 3,
        "query_binding": "native_query_id",
    }
    assert payload["native_input"]["sha256"]
    assert payload["normalized_output"]["row_count"] == 3
    assert "wrote 3 normalized nominations" in capsys.readouterr().out


def test_cas_offinder_2_4_1_rejects_uppercase_mismatch_claimed_zero():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match=r"lowercase markers \[\] do not match"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "uppercase_mismatch_claimed_zero.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_rejects_matching_base_marked_as_mismatch():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match=r"lowercase markers \[0\] do not match"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "matching_base_marked_mismatch.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_rejects_pam_region_mismatch():
    _, _, profile = manifests()
    guides = fixture_guides()

    with pytest.raises(CandidatePanelError, match="PAM-region mismatch cannot be normalized"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "pam-mismatch.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            guide_override="alpha_guide",
        )


def test_cas_offinder_2_4_1_rejects_extra_base_geometry():
    _, _, profile = manifests()
    guides = fixture_guides()

    with pytest.raises(CandidatePanelError, match=r"spacer-plus-PAM geometry \(20\+3 = 23"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "extra_base_geometry.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            guide_override="alpha_guide",
        )


def test_cas_offinder_2_4_1_rejects_conflicting_source_version():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match="conflicts with the explicit Cas-OFFinder"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "wrong_source_version.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 3.0",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_query_map_rejects_duplicate_serialized_keys():
    with pytest.raises(CandidatePanelError, match="duplicate key 'TTGCAAGTCCGTAGGAATCGNNN'"):
        load_query_map(FIXTURE_241 / "duplicate-query-map.json", fixture_guides("shared-guides.tsv"))


def test_cas_offinder_2_4_1_rejects_unresolved_target_ambiguity():
    _, _, profile = manifests()
    guides = fixture_guides()
    query_map = load_query_map(FIXTURE_241 / "query-map.json", guides)

    with pytest.raises(CandidatePanelError, match="unresolved target ambiguity"):
        normalize_cas_offinder_2_4_1(
            FIXTURE_241 / "ambiguous-target.tsv",
            guides,
            fixture_loci(),
            profile,
            "Cas-OFFinder 2.4.1 (9816b94)",
            with_ids=False,
            query_map=query_map,
        )


def test_cas_offinder_2_4_1_accepts_compatible_query_pam_mask():
    _, _, profile = manifests()

    rows, _ = normalize_cas_offinder_2_4_1(
        FIXTURE_241 / "valid-seven-field-ngg-mask.tsv",
        fixture_guides(),
        fixture_loci(),
        profile,
        "Cas-OFFinder 2.4.1 (9816b94)",
        with_ids=True,
    )

    assert [row["mismatch_count"] for row in rows] == [0, 1]
    assert [row["pam"] for row in rows] == ["AGG", "AGG"]
    assert [row["site_role"] for row in rows] == ["intended", "candidate"]


def test_cas_offinder_2_4_1_retains_line_locator_and_full_row_checksum():
    import hashlib

    _, _, profile = manifests()
    native = FIXTURE_241 / "valid-seven-field.tsv"
    rows, _ = normalize_cas_offinder_2_4_1(
        native,
        fixture_guides(),
        fixture_loci(),
        profile,
        "Cas-OFFinder 2.4.1 (9816b94)",
        with_ids=True,
    )
    expected = [
        hashlib.sha256(line.encode("utf-8")).hexdigest()
        for line in native.read_text().splitlines()
    ]

    for index, row in enumerate(rows):
        metadata = json.loads(row["source_metadata_json"])
        assert metadata["native_line"] == index + 1
        assert metadata["native_row_sha256"] == expected[index]


def test_cas_offinder_2_4_1_cli_pins_query_map_in_provenance(tmp_path):
    import hashlib

    output = tmp_path / "normalized.tsv"
    provenance = tmp_path / "source.json"

    status = main(
        [
            "cas-offinder",
            "--format",
            "2.4.1",
            "--query-map",
            str(FIXTURE_241 / "query-map.json"),
            "--input",
            str(FIXTURE_241 / "valid-six-field.tsv"),
            "--profile",
            str(PROFILE_PATH),
            "--guides",
            str(FIXTURE_241 / "guides.tsv"),
            "--intended-loci",
            str(FIXTURE_241 / "intended-loci.tsv"),
            "--source-version",
            "Cas-OFFinder 2.4.1 (9816b94c20c4cba2e79b039e1e2a6dee684b7b66)",
            "--source-command",
            "not-applicable: mock parser input",
            "--reference-fasta",
            str(FIXTURE_241 / "reference.fa"),
            "--output",
            str(output),
            "--provenance-output",
            str(provenance),
        ]
    )

    assert status == 0
    payload = json.loads(provenance.read_text())
    assert payload["query_map"] == {
        "file": "query-map.json",
        "sha256": hashlib.sha256((FIXTURE_241 / "query-map.json").read_bytes()).hexdigest(),
    }
    assert payload["row_accounting"]["query_binding"] == "query_map"
    assert payload["source_run"]["command"] == (
        "not-applicable: mock parser input"
    )
