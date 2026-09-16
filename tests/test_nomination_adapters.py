"""Native nomination adapters preserve rows, identity, and source-run provenance."""
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
    main,
    normalize_cas_offinder,
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
