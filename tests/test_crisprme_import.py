"""Check CRISPRme+ parsing with constructed rows and an optional local workbook."""
from __future__ import annotations

import csv
from hashlib import sha256
from html import escape
import json
import os
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import pytest

from base_edit_safety import __version__
from base_edit_safety.ingestion.crisprme import (
    CRISPRME_SCORING_MODEL,
    CRISPRME_TSV_FORMAT_SCHEMA,
    CRISPRME_VERSION,
    CRISPRME_XLSX_FORMAT_SCHEMA,
    REQUIRED_COLUMNS,
    CrisprmeImportResult,
    CrisprmeImportError,
    import_crisprme,
)
from base_edit_safety.ingestion.xlsx import iter_rows as iter_xlsx_rows
from base_edit_safety.ingestion.xlsx import sheet_names
from base_edit_safety.pipeline.candidate_panel import (
    CandidatePanelError,
    SearchProfile,
    _assess_allele_window,
    _comparison_is_armed,
    _why_crispor_missed,
    build_ledger,
    merge_evidence,
    validate_reference_sequences,
)
from base_edit_safety.pipeline.normalize_nominations import (
    NOMINATION_COLUMNS,
    main as normalize_main,
    normalize_crisprme,
    write_normalized,
    write_provenance,
)
from base_edit_safety.pipeline.datatypes import reverse_complement


# Constructed guide and variant row; no paper sequence, coordinates, scores or frequencies.
GUIDE_ID = "g1"
GUIDE = "GTAACAGTTGCTTTTATCAT"

FROZEN_REQUIRED_COLUMNS = frozenset(
    {
        "Spacer+PAM",
        "Chromosome",
        "Start_coordinate_(highest_CFD)",
        "Strand_(highest_CFD)",
        "Aligned_spacer+PAM_(highest_CFD)",
        "Aligned_protospacer+PAM_REF_(highest_CFD)",
        "Aligned_protospacer+PAM_ALT_(highest_CFD)",
        "PAM_(highest_CFD)",
        "Mismatches_(highest_CFD)",
        "Bulges_(highest_CFD)",
        "Bulge_type_(highest_CFD)",
        "REF/ALT_origin_(highest_CFD)",
        "PAM_creation_(highest_CFD)",
        "CFD_score_(highest_CFD)",
        "CFD_score_REF_(highest_CFD)",
        "CFD_score_ALT_(highest_CFD)",
        "Variant_info_spacer+PAM_(highest_CFD)",
        "Variant_info_genome_(highest_CFD)",
        "Variant_MAF_(highest_CFD)",
        "Variant_rsID_(highest_CFD)",
        "Not_found_in_REF",
    }
)


def test_crisprme_required_header_contract_is_frozen():
    assert REQUIRED_COLUMNS == FROZEN_REQUIRED_COLUMNS


def test_crisprme_plus_runtime_contract_is_pinned():
    assert CRISPRME_VERSION == "2.3.2"
    assert CRISPRME_SCORING_MODEL == "CFD"
    assert CRISPRME_TSV_FORMAT_SCHEMA == "crisprme-2.3.2-integrated-results-v1"


def base_row() -> dict[str, str]:
    return {
        "Spacer+PAM": GUIDE + "NNN",
        "Chromosome": "chr2",
        "Start_coordinate_(highest_CFD)": "1000",
        "Strand_(highest_CFD)": "-",
        "Aligned_spacer+PAM_(highest_CFD)": GUIDE + "NNN",
        "Aligned_protospacer+PAM_REF_(highest_CFD)": "aTAACAGcTGCcTTTATCATTGC",
        "Aligned_protospacer+PAM_ALT_(highest_CFD)": "aTAACAGcTGCcTTTATCATTGG",
        "PAM_(highest_CFD)": "TGG",
        "Mismatches_(highest_CFD)": "3",
        "Bulges_(highest_CFD)": "0",
        "Bulge_type_(highest_CFD)": "X",
        "REF/ALT_origin_(highest_CFD)": "alt",
        "PAM_creation_(highest_CFD)": "NA",
        "CFD_score_(highest_CFD)": "0.81250000000000000",
        "CFD_score_REF_(highest_CFD)": "3.1250000000000000E-2",
        "CFD_score_ALT_(highest_CFD)": "0.81250000000000000",
        "Variant_info_spacer+PAM_(highest_CFD)": "C23G",
        "Variant_info_genome_(highest_CFD)": "chr2_1001_G_C",
        "Variant_MAF_(highest_CFD)": "0.25,0.125",
        "Variant_rsID_(highest_CFD)": "NA,synthetic_variant_1",
        "Not_found_in_REF": "NA",
        "Variant_samples_(highest_CFD)": "SYNTHETIC_SAMPLE_MUST_DISAPPEAR",
    }


def write_tsv(path: Path, rows: list[dict[str, str]]) -> Path:
    header = sorted(set(REQUIRED_COLUMNS).union(*(row.keys() for row in rows)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def declaration_payload(
    native: Path,
    *,
    format_schema: str = CRISPRME_TSV_FORMAT_SCHEMA,
    coordinate_system: str = "zero_based_half_open",
    guide_id: str = GUIDE_ID,
    guide: str = GUIDE,
) -> dict:
    is_xlsx = format_schema == CRISPRME_XLSX_FORMAT_SCHEMA
    return {
        "record_type": "crisprme_import_declaration",
        "schema_version": "1",
        "format_schema": format_schema,
        "native_input": {
            "file": native.name,
            "sha256": sha256(native.read_bytes()).hexdigest(),
        },
        "selected_sheet": f"{guide_id}_integrated_results" if is_xlsx else None,
        "guide_id": guide_id,
        "guide_protospacer": guide,
        "genome_build": "GRCh38",
        "coordinate_system": coordinate_system,
        "crisprme_version": CRISPRME_VERSION,
        "scoring_model": CRISPRME_SCORING_MODEL,
        "variant_databases": [{"name": "synthetic variants", "release": "v1"}],
        "pam_patterns": ["NGG"],
        "max_mismatches": 6,
        "max_dna_bulge": 2,
        "max_rna_bulge": 2,
        "retention": {
            "completeness": "complete",
            "top_n": None,
            "score_cutoff": None,
            "merge_window_bp": 0,
            "all_candidate_coordinates_retained": True,
        },
        "source_command": "synthetic test fixture; CRISPRme was not executed",
    }


def crispor_declaration(
    *,
    pam_patterns: list[str] | None = None,
    max_mismatches: int = 4,
    completeness: str = "complete",
) -> dict:
    truncated = completeness == "truncated_top_n"
    return {
        "record_type": "crispor_search_declaration",
        "schema_version": "1",
        "guide_id": GUIDE_ID,
        "guide_protospacer": GUIDE,
        "genome_build": "GRCh38",
        "pam_patterns": pam_patterns or ["NGG"],
        "max_mismatches": max_mismatches,
        "max_dna_bulge": 0,
        "max_rna_bulge": 0,
        "retention": {
            "completeness": completeness,
            "top_n": 1 if truncated else None,
            "score_cutoff": None,
            "all_candidate_coordinates_retained": not truncated,
        },
    }


def write_declaration(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _column_name(index: int) -> str:
    result = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        result = chr(ord("A") + remainder) + result
    return result


def _sheet_xml(rows: list[list[str]]) -> str:
    rendered = []
    for row_index, row in enumerate(rows, start=1):
        cells = []
        for column_index, value in enumerate(row):
            reference = f"{_column_name(column_index)}{row_index}"
            cells.append(
                f'<c r="{reference}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'
            )
        rendered.append(f'<row r="{row_index}">{"".join(cells)}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(rendered)}</sheetData></worksheet>'
    )


def write_xlsx(path: Path, sheets: list[tuple[str, list[list[str]]]]) -> Path:
    workbook_sheets = "".join(
        f'<sheet name="{escape(name)}" sheetId="{index}" r:id="rId{index}"/>'
        for index, (name, _rows) in enumerate(sheets, start=1)
    )
    relationships = "".join(
        '<Relationship '
        f'Id="rId{index}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        f'Target="worksheets/sheet{index}.xml"/>'
        for index in range(1, len(sheets) + 1)
    )
    content_overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{index}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for index in range(1, len(sheets) + 1)
    )
    with ZipFile(path, "w", ZIP_DEFLATED) as book:
        book.writestr(
            "[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            f"{content_overrides}</Types>",
        )
        book.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
            'Target="xl/workbook.xml"/></Relationships>',
        )
        book.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f"<sheets>{workbook_sheets}</sheets></workbook>",
        )
        book.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f"{relationships}</Relationships>",
        )
        for index, (_name, rows) in enumerate(sheets, start=1):
            book.writestr(f"xl/worksheets/sheet{index}.xml", _sheet_xml(rows))
    return path


def test_native_tsv_normalizes_ref_and_alt_to_one_site_without_samples(tmp_path):
    native = write_tsv(tmp_path / "run.integrated_results.tsv", [base_row()])
    declaration = write_declaration(
        tmp_path / "run.crisprme-import.json", declaration_payload(native)
    )

    result = import_crisprme(native, declaration)

    assert result.import_summary == {
        "input_data_rows": 1,
        "emitted_nomination_rows": 1,
        "skipped_rows": 0,
        "skipped_by_reason": {},
    }
    assert len(result.rows) == 1
    row = result.rows[0]
    assert (row["chrom"], row["start_0based"], row["end_0based_exclusive"]) == (
        "chr2",
        1000,
        1023,
    )
    assert row["mismatch_count"] == 3
    assert [state["allele"] for state in row["allele_states"]] == ["ref", "alt"]
    assert [state["CFD"] for state in row["allele_states"]] == [0.03125, 0.8125]
    assert row["variant_maf_raw"] == "0.25,0.125"
    assert row["variant_rsid_raw"] == "NA,synthetic_variant_1"
    rendered = json.dumps(
        {
            "rows": result.rows,
            "declaration": result.declaration,
            "summary": result.import_summary,
        }
    )
    assert "SYNTHETIC_SAMPLE_MUST_DISAPPEAR" not in rendered
    assert "Variant_samples" not in rendered


def test_crisprme_plus_additive_columns_do_not_change_the_logical_contract(tmp_path):
    row = base_row()
    row.update(
        {
            "Seed_mismatches+bulges_REF_(highest_CFD)": "2",
            "Non_seed_mismatches+bulges_REF_(highest_CFD)": "1",
            "Susceptible_to_ABE": "y",
        }
    )
    native = write_tsv(
        tmp_path / "additive.integrated_results.tsv", [row]
    )
    declaration = write_declaration(
        tmp_path / "additive.crisprme-import.json", declaration_payload(native)
    )

    result = import_crisprme(native, declaration)

    assert len(result.rows) == 1
    rendered = json.dumps(result.rows)
    assert "Susceptible_to_ABE" not in rendered
    assert "Seed_mismatches" not in rendered


def test_native_tsv_requires_the_pinned_version_and_scoring_model(tmp_path):
    native = write_tsv(
        tmp_path / "pin.integrated_results.tsv", [base_row()]
    )
    wrong_version = declaration_payload(native)
    wrong_version["crisprme_version"] = "2.1.14"
    with pytest.raises(CrisprmeImportError, match="requires version 2.3.2"):
        import_crisprme(
            native,
            write_declaration(tmp_path / "wrong-version.json", wrong_version),
        )

    wrong_model = declaration_payload(native)
    wrong_model["scoring_model"] = "CFD++"
    with pytest.raises(CrisprmeImportError, match="scoring_model must be CFD"):
        import_crisprme(
            native,
            write_declaration(tmp_path / "wrong-model.json", wrong_model),
        )


def test_native_tsv_refuses_an_alternative_alignment_filename(tmp_path):
    native = write_tsv(tmp_path / "all_results_with_alternative_alignments.tsv", [base_row()])
    declaration = write_declaration(
        tmp_path / "alternative.json", declaration_payload(native)
    )

    with pytest.raises(CrisprmeImportError, match="integrated_results.tsv"):
        import_crisprme(native, declaration)


def test_declared_coordinate_systems_normalize_to_same_interval(tmp_path):
    zero_row = base_row()
    zero_native = write_tsv(tmp_path / "zero.integrated_results.tsv", [zero_row])
    zero_declaration = write_declaration(
        tmp_path / "zero.json", declaration_payload(zero_native)
    )

    one_row = base_row()
    one_row["Start_coordinate_(highest_CFD)"] = "1001"
    one_native = write_tsv(tmp_path / "one.integrated_results.tsv", [one_row])
    one_declaration = write_declaration(
        tmp_path / "one.json",
        declaration_payload(one_native, coordinate_system="one_based_inclusive"),
    )

    zero = import_crisprme(zero_native, zero_declaration).rows[0]
    one = import_crisprme(one_native, one_declaration).rows[0]
    assert (zero["start_0based"], zero["end_0based_exclusive"]) == (
        one["start_0based"],
        one["end_0based_exclusive"],
    )


@pytest.mark.parametrize("mismatch_count", [5, 6])
def test_deep_mismatch_rows_are_retained(tmp_path, mismatch_count):
    row = base_row()
    replacements = {"A": "C", "C": "A", "G": "T", "T": "G"}
    target_spacer = "".join(
        replacements[base] if index < mismatch_count else base
        for index, base in enumerate(GUIDE)
    )
    row.update(
        {
            "Aligned_protospacer+PAM_REF_(highest_CFD)": target_spacer + "TGG",
            "Aligned_protospacer+PAM_ALT_(highest_CFD)": "NA",
            "Mismatches_(highest_CFD)": str(mismatch_count),
            "REF/ALT_origin_(highest_CFD)": "ref",
            "CFD_score_(highest_CFD)": "0.5",
            "CFD_score_REF_(highest_CFD)": "0.5",
        }
    )
    native = write_tsv(
        tmp_path / f"mm{mismatch_count}.integrated_results.tsv", [row]
    )
    declaration = write_declaration(
        tmp_path / f"mm{mismatch_count}.json", declaration_payload(native)
    )

    imported = import_crisprme(native, declaration).rows[0]

    assert imported["mismatch_count"] == mismatch_count


@pytest.mark.parametrize(
    ("bulge_type", "guide_alignment", "target_alignment", "expected_span"),
    [
        (
            "DNA",
            GUIDE[:2] + "-" + GUIDE[2:] + "NNN",
            GUIDE[:2] + "C" + GUIDE[2:] + "TGG",
            24,
        ),
        (
            "RNA",
            GUIDE + "NNN",
            GUIDE[:2] + "-" + GUIDE[3:] + "TGG",
            22,
        ),
    ],
)
def test_bulge_rows_validate_gap_topology_and_interval(
    tmp_path, bulge_type, guide_alignment, target_alignment, expected_span
):
    row = base_row()
    row.update(
        {
            "Start_coordinate_(highest_CFD)": "100",
            "Aligned_spacer+PAM_(highest_CFD)": guide_alignment,
            "Aligned_protospacer+PAM_REF_(highest_CFD)": target_alignment,
            "Aligned_protospacer+PAM_ALT_(highest_CFD)": "NA",
            "Mismatches_(highest_CFD)": "0",
            "Bulges_(highest_CFD)": "1",
            "Bulge_type_(highest_CFD)": bulge_type,
            "REF/ALT_origin_(highest_CFD)": "ref",
            "CFD_score_(highest_CFD)": "0.5",
            "CFD_score_REF_(highest_CFD)": "0.5",
        }
    )
    native = write_tsv(
        tmp_path / f"{bulge_type}.integrated_results.tsv", [row]
    )
    declaration = write_declaration(
        tmp_path / f"{bulge_type}.json", declaration_payload(native)
    )

    imported = import_crisprme(native, declaration).rows[0]

    assert imported["end_0based_exclusive"] - imported["start_0based"] == expected_span
    assert imported["dna_bulge"] == (1 if bulge_type == "DNA" else 0)
    assert imported["rna_bulge"] == (1 if bulge_type == "RNA" else 0)


def test_inconsistent_bulge_gap_topology_is_skipped(tmp_path):
    bad = base_row()
    bad["Bulges_(highest_CFD)"] = "1"
    bad["Bulge_type_(highest_CFD)"] = "DNA"
    native = write_tsv(
        tmp_path / "bad-gap.integrated_results.tsv", [bad, base_row()]
    )
    declaration = write_declaration(
        tmp_path / "bad-gap.json", declaration_payload(native)
    )

    result = import_crisprme(native, declaration)

    assert result.import_summary["skipped_by_reason"] == {"inconsistent_span": 1}


def test_recoverable_bad_row_is_counted_without_row_values(tmp_path):
    bad = base_row()
    bad["Start_coordinate_(highest_CFD)"] = "SECRET_BAD_COORDINATE"
    native = write_tsv(
        tmp_path / "mixed.integrated_results.tsv", [bad, base_row()]
    )
    declaration = write_declaration(
        tmp_path / "mixed.json", declaration_payload(native)
    )

    result = import_crisprme(native, declaration)

    assert len(result.rows) == 1
    assert result.import_summary["skipped_rows"] == 1
    assert result.import_summary["skipped_by_reason"] == {"invalid_coordinate": 1}
    assert "SECRET_BAD_COORDINATE" not in json.dumps(result.import_summary)


def test_all_bad_rows_refuse_without_echoing_values(tmp_path):
    bad = base_row()
    bad["Start_coordinate_(highest_CFD)"] = "SECRET_BAD_COORDINATE"
    native = write_tsv(tmp_path / "bad.integrated_results.tsv", [bad])
    declaration = write_declaration(tmp_path / "bad.json", declaration_payload(native))

    with pytest.raises(CrisprmeImportError) as captured:
        import_crisprme(native, declaration)

    assert "no valid CRISPRme rows" in str(captured.value)
    assert "SECRET_BAD_COORDINATE" not in str(captured.value)


def test_mismatched_sidecar_checksum_is_refused(tmp_path):
    native = write_tsv(tmp_path / "run.integrated_results.tsv", [base_row()])
    payload = declaration_payload(native)
    payload["native_input"]["sha256"] = "0" * 64
    declaration = write_declaration(tmp_path / "run.json", payload)

    with pytest.raises(CrisprmeImportError, match="checksum does not match"):
        import_crisprme(native, declaration)


def test_missing_sidecar_is_refused(tmp_path):
    native = write_tsv(tmp_path / "run.integrated_results.tsv", [base_row()])

    with pytest.raises(CrisprmeImportError, match="cannot read CRISPRme import declaration"):
        import_crisprme(native, tmp_path / "missing-sidecar.json")


def test_sidecar_without_coordinate_system_is_refused(tmp_path):
    native = write_tsv(tmp_path / "run.integrated_results.tsv", [base_row()])
    payload = declaration_payload(native)
    del payload["coordinate_system"]
    declaration = write_declaration(tmp_path / "undeclared-coordinate-system.json", payload)

    with pytest.raises(CrisprmeImportError, match="must contain exactly"):
        import_crisprme(native, declaration)


def test_row_inconsistent_with_declared_coordinate_system_is_skipped_and_counted(tmp_path):
    inconsistent = base_row()
    inconsistent["Start_coordinate_(highest_CFD)"] = "0"
    valid = base_row()
    valid["Start_coordinate_(highest_CFD)"] = "1"
    native = write_tsv(
        tmp_path / "one-based.integrated_results.tsv", [inconsistent, valid]
    )
    declaration = write_declaration(
        tmp_path / "one-based.json",
        declaration_payload(native, coordinate_system="one_based_inclusive"),
    )

    result = import_crisprme(native, declaration)

    assert len(result.rows) == 1
    assert result.rows[0]["start_0based"] == 0
    assert result.import_summary["skipped_rows"] == 1
    assert result.import_summary["skipped_by_reason"] == {"invalid_coordinate": 1}


def test_validation_xlsx_reads_only_allowlisted_sheet(tmp_path):
    row = base_row()
    header = sorted(row)
    native = write_xlsx(
        tmp_path / "publication.xlsx",
        [
            (f"{GUIDE_ID}_integrated_results", [header, [row[field] for field in header]]),
            (
                f"{GUIDE_ID}_SYNTHETIC_PRIVATE_private_targets",
                [["Private_data"], ["PRIVATE_SHEET_VALUE_MUST_DISAPPEAR"]],
            ),
        ],
    )
    payload = declaration_payload(native, format_schema=CRISPRME_XLSX_FORMAT_SCHEMA)
    payload["retention"] = {
        "completeness": "truncated_top_n",
        "top_n": 1,
        "score_cutoff": None,
        "merge_window_bp": 3,
        "all_candidate_coordinates_retained": False,
    }
    declaration = write_declaration(tmp_path / "publication.json", payload)

    result = import_crisprme(native, declaration)

    assert len(result.rows) == 1
    assert result.sheets_ignored == {
        "count": 1,
        "categories": ["private_targets"],
    }
    assert result.privacy_filtering_applied is True
    rendered = json.dumps(
        {
            "rows": result.rows,
            "declaration": result.declaration,
            "summary": result.import_summary,
        }
    )
    assert "SYNTHETIC_SAMPLE_MUST_DISAPPEAR" not in rendered
    assert "PRIVATE_SHEET_VALUE_MUST_DISAPPEAR" not in rendered
    assert "private_targets" not in rendered


def test_validation_xlsx_refuses_private_sheet_selection(tmp_path):
    row = base_row()
    header = sorted(row)
    native = write_xlsx(
        tmp_path / "publication.xlsx",
        [
            (f"{GUIDE_ID}_integrated_results", [header, [row[field] for field in header]]),
            (f"{GUIDE_ID}_PRIVATE_private_targets", [["Private_data"], ["SECRET"]]),
        ],
    )
    payload = declaration_payload(native, format_schema=CRISPRME_XLSX_FORMAT_SCHEMA)
    payload["selected_sheet"] = f"{GUIDE_ID}_PRIVATE_private_targets"
    declaration = write_declaration(tmp_path / "publication.json", payload)

    with pytest.raises(CrisprmeImportError, match="integrated-results sheet") as captured:
        import_crisprme(native, declaration)

    assert "PRIVATE" not in str(captured.value)
    assert "SECRET" not in str(captured.value)


def write_manifest_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    guide_path = tmp_path / "guides.tsv"
    guide_path.write_text(
        "\t".join(
            [
                "guide_id",
                "spacer_sequence",
                "editor",
                "editor_class",
                "pam_class",
                "editing_window_start_1based",
                "editing_window_end_1based_inclusive",
                "intended_target_description",
                "locus_system",
            ]
        )
        + "\n"
        + "\t".join(
            [
                GUIDE_ID,
                GUIDE,
                "ABE8e",
                "ABE",
                "NGG",
                "4",
                "8",
                "synthetic HBG targets",
                "HBG1_HBG2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    intended_path = tmp_path / "intended.tsv"
    intended_path.write_text(
        "guide_id\tgenome_build\tchrom\tstart_0based\tend_0based_exclusive\tstrand\tlocus_name\n"
        f"{GUIDE_ID}\tGRCh38\tchr2\t100\t123\t+\tHBG1\n"
        f"{GUIDE_ID}\tGRCh38\tchr2\t5024\t5047\t+\tHBG2\n",
        encoding="utf-8",
    )
    profile_path = tmp_path / "profile.json"
    profile_path.write_text(
        json.dumps(
            {
                "name": "crisprme_synthetic_v1",
                "max_mismatches": 6,
                "max_dna_bulge": 2,
                "max_rna_bulge": 2,
                "pam_patterns": ["NGG"],
                "variant_aware": True,
                "intended_use": "Synthetic CRISPRme reconciliation test.",
                "reference_id": "GRCh38",
                "reference_provenance": "Synthetic chr2 sequence generated by this test.",
                "variant_source": "Synthetic variants v1.",
                "retention_rule": "Retain every row within the declared search limits.",
                "retention_score_name": None,
                "minimum_score": None,
            }
        ),
        encoding="utf-8",
    )
    return guide_path, intended_path, profile_path


def write_reference(tmp_path: Path) -> Path:
    constructed_ref = base_row()["Aligned_protospacer+PAM_REF_(highest_CFD)"].upper()
    ordinary = GUIDE + "TGG"
    sequence = list("A" * 5100)
    sequence[10:33] = reverse_complement(constructed_ref)
    sequence[34:57] = ordinary
    sequence[70:93] = ordinary
    sequence[100:123] = ordinary
    sequence[5024:5047] = ordinary
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr2\n" + "".join(sequence) + "\n", encoding="ascii")
    Path(str(reference) + ".fai").write_text(
        "chr2\t5100\t6\t5100\t5101\n", encoding="ascii"
    )
    return reference


def test_crisprme_manifest_records_structural_privacy_filtering_without_identifiers(tmp_path):
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    reference = write_reference(tmp_path)
    canary = "DONORXYZ123"
    private_sheet = f"{GUIDE_ID}_SYNTHETIC_SAMPLE_private_targets"
    row = base_row()
    row["Start_coordinate_(highest_CFD)"] = "10"
    row["Variant_samples_(highest_CFD)"] = canary
    header = sorted(row)
    native = write_xlsx(
        tmp_path / "privacy-canary.xlsx",
        [
            (
                f"{GUIDE_ID}_integrated_results",
                [header, [row[field] for field in header]],
            ),
            (f"{GUIDE_ID}_altMerge", [["private_value"], [canary]]),
            (private_sheet, [["private_value"], [canary]]),
        ],
    )
    payload = declaration_payload(native, format_schema=CRISPRME_XLSX_FORMAT_SCHEMA)
    payload["retention"] = {
        "completeness": "truncated_top_n",
        "top_n": 1,
        "score_cutoff": None,
        "merge_window_bp": 3,
        "all_candidate_coordinates_retained": False,
    }
    declaration = write_declaration(tmp_path / "privacy-canary.json", payload)
    normalized = tmp_path / "privacy-canary.normalized.tsv"
    provenance = tmp_path / "privacy-canary.provenance.json"

    assert normalize_main(
        [
            "crisprme",
            "--input",
            str(native),
            "--declaration",
            str(declaration),
            "--profile",
            str(profile_path),
            "--guides",
            str(guide_path),
            "--intended-loci",
            str(intended_path),
            "--reference-fasta",
            str(reference),
            "--output",
            str(normalized),
            "--provenance-output",
            str(provenance),
        ]
    ) == 0
    manifest = json.loads(provenance.read_text())
    ledger = build_ledger(
        profile_path,
        guide_path,
        intended_path,
        [normalized],
        source_provenance_paths=[provenance],
        reference_fasta_path=reference,
    )

    assert manifest["privacy_filtering_applied"] is True
    assert manifest["sheets_ignored"] == {
        "count": 2,
        "categories": ["altMerge", "private_targets"],
    }
    rendered = normalized.read_text() + provenance.read_text() + json.dumps(ledger)
    assert canary not in rendered
    assert private_sheet not in rendered


def ordinary_crisprme_row(start: int) -> dict[str, str]:
    row = base_row()
    row.update(
        {
            "Start_coordinate_(highest_CFD)": str(start),
            "Strand_(highest_CFD)": "+",
            "Aligned_protospacer+PAM_REF_(highest_CFD)": GUIDE + "TGG",
            "Aligned_protospacer+PAM_ALT_(highest_CFD)": "NA",
            "PAM_(highest_CFD)": "TGG",
            "Mismatches_(highest_CFD)": "0",
            "REF/ALT_origin_(highest_CFD)": "ref",
            "CFD_score_(highest_CFD)": "1",
            "CFD_score_REF_(highest_CFD)": "1",
            "CFD_score_ALT_(highest_CFD)": "1",
            "Variant_info_spacer+PAM_(highest_CFD)": "NA",
            "Variant_info_genome_(highest_CFD)": "NA",
            "Variant_MAF_(highest_CFD)": "NA",
            "Variant_rsID_(highest_CFD)": "NA",
        }
    )
    return row


def crispor_nomination(start: int, record_id: str) -> dict[str, object]:
    row = {column: "" for column in NOMINATION_COLUMNS}
    row.update(
        {
            "guide_id": GUIDE_ID,
            "nomination_kind": "computational",
            "source_method": "CRISPOR",
            "source_version_or_citation": "CRISPOR synthetic contract",
            "genome_build": "GRCh38",
            "chrom": "chr2",
            "start_0based": start,
            "end_0based_exclusive": start + 23,
            "strand": "+",
            "site_role": "candidate",
            "site_sequence": GUIDE + "TGG",
            "pam": "TGG",
            "mismatch_count": 0,
            "dna_bulge": 0,
            "rna_bulge": 0,
            "score_name": "CFD",
            "score_value": 1.0,
            "editable_base_in_window": True,
            "source_record_id": record_id,
            "search_profile": "crisprme_synthetic_v1",
        }
    )
    return row


def build_superset_test_ledger(tmp_path: Path, source_declaration: dict) -> dict:
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    reference = write_reference(tmp_path)
    native = write_tsv(
        tmp_path / "comparison.integrated_results.tsv",
        [ordinary_crisprme_row(34)],
    )
    declaration_path = write_declaration(
        tmp_path / "comparison.crisprme-import.json", declaration_payload(native)
    )
    crisprme_tsv = tmp_path / "comparison-crisprme.tsv"
    crisprme_provenance = tmp_path / "comparison-crisprme.json"
    assert normalize_main(
        [
            "crisprme",
            "--input",
            str(native),
            "--declaration",
            str(declaration_path),
            "--profile",
            str(profile_path),
            "--guides",
            str(guide_path),
            "--intended-loci",
            str(intended_path),
            "--reference-fasta",
            str(reference),
            "--output",
            str(crisprme_tsv),
            "--provenance-output",
            str(crisprme_provenance),
        ]
    ) == 0

    crispor_rows = [crispor_nomination(70, "crispor-missing-from-crisprme")]
    validate_reference_sequences(crispor_rows, reference)
    crispor_tsv = tmp_path / "comparison-crispor.tsv"
    write_normalized(crispor_rows, crispor_tsv)
    crispor_native = tmp_path / "comparison-crispor-native.tsv"
    crispor_native.write_text("synthetic native placeholder\n", encoding="utf-8")
    crispor_provenance = tmp_path / "comparison-crispor.json"
    write_provenance(
        crispor_provenance,
        crispor_tsv,
        crispor_native,
        "synthetic-crispor",
        "CRISPOR",
        "CRISPOR synthetic contract",
        "synthetic comparison fixture; CRISPOR was not executed",
        "GRCh38",
        reference,
        guide_path,
        intended_path,
        profile_path,
        len(crispor_rows),
        crispor_declaration=source_declaration,
    )
    return build_ledger(
        profile_path,
        guide_path,
        intended_path,
        [crisprme_tsv, crispor_tsv],
        source_provenance_paths=[crisprme_provenance, crispor_provenance],
        reference_fasta_path=reference,
    )


def test_truncated_crispor_provenance_does_not_arm_superset_alarm(tmp_path):
    ledger = build_superset_test_ledger(
        tmp_path, crispor_declaration(completeness="truncated_top_n")
    )

    site = next(item for item in ledger["candidates"] if item["start_0based"] == 70)
    assert site["source_bucket"] == "crispor_only"
    assert site["reconciliation"] == {
        "status": "comparison_not_armed",
        "warnings": [],
    }


def test_complete_compatible_crispor_provenance_arms_superset_alarm(tmp_path):
    ledger = build_superset_test_ledger(tmp_path, crispor_declaration())

    site = next(item for item in ledger["candidates"] if item["start_0based"] == 70)
    assert site["reconciliation"]["status"] == "qc_warning_missing_from_crisprme"
    assert len(site["reconciliation"]["warnings"]) == 1


def test_complete_crispor_provenance_with_incompatible_pam_does_not_arm_alarm(tmp_path):
    ledger = build_superset_test_ledger(
        tmp_path, crispor_declaration(pam_patterns=["NNG"])
    )

    site = next(item for item in ledger["candidates"] if item["start_0based"] == 70)
    assert site["reconciliation"] == {
        "status": "comparison_not_armed",
        "warnings": [],
    }


def test_complete_crispor_provenance_with_incompatible_mismatch_limit_does_not_arm_alarm(
    tmp_path,
):
    ledger = build_superset_test_ledger(
        tmp_path, crispor_declaration(max_mismatches=7)
    )

    site = next(item for item in ledger["candidates"] if item["start_0based"] == 70)
    assert site["reconciliation"] == {
        "status": "comparison_not_armed",
        "warnings": [],
    }


def test_variant_created_pam_positive_is_explained():
    metadata = {
        "canonical_interval_basis": "ref",
        "guide_alignment": GUIDE + "NNN",
        "selected_origin": "alt",
        "pam_creation_raw": "yes",
        "not_found_in_ref_raw": None,
        "allele_states": [
            {"allele": "ref", "alignment": GUIDE + "TGC"},
            {"allele": "alt", "alignment": GUIDE + "TGG"},
        ],
    }
    site = {
        "source_bucket": "crisprme_only",
        "nominations": [
            {
                "nomination_kind": "computational",
                "source_method": "CRISPRme",
                "mismatch_count": 0,
                "dna_bulge": 0,
                "rna_bulge": 0,
                "pam": "TGG",
                "score_name": "CFD",
                "score_value": 1.0,
                "source_metadata": {"crisprme": metadata},
            }
        ],
    }

    assert _why_crispor_missed(site, {"pam_class": "NGG"}) == [
        "variant_created_pam"
    ]


def test_minus_strand_row_normalization(tmp_path):
    native = write_tsv(tmp_path / "minus.integrated_results.tsv", [base_row()])
    declaration = write_declaration(
        tmp_path / "minus.json", declaration_payload(native)
    )

    row = import_crisprme(native, declaration).rows[0]

    assert row["strand"] == "-"
    assert (row["start_0based"], row["end_0based_exclusive"]) == (
        1000,
        1023,
    )
    assert {
        state["effective_interval"]["strand"] for state in row["allele_states"]
    } == {"-"}


def test_4924_bp_paralog_decoy_sites_do_not_merge():
    sites = merge_evidence(
        [
            crispor_nomination(100, "decoy-left"),
            crispor_nomination(5024, "decoy-right"),
        ],
        [],
        {},
        guides=[
            {
                "guide_id": GUIDE_ID,
                "spacer_sequence": GUIDE,
                "pam_class": "NGG",
                "editor_class": "ABE",
            }
        ],
    )

    assert [site["start_0based"] for site in sites] == [100, 5024]
    assert sites[1]["start_0based"] - sites[0]["start_0based"] == 4924
    assert sites[0]["site_id"] != sites[1]["site_id"]


@pytest.mark.parametrize(
    "format_schema",
    [CRISPRME_TSV_FORMAT_SCHEMA, CRISPRME_XLSX_FORMAT_SCHEMA],
    ids=["native", "legacy"],
)
def test_crisprme_normalization_builds_v3_reconciliation_and_provenance(
    tmp_path, format_schema
):
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    reference = write_reference(tmp_path)
    constructed_variant_row = base_row()
    constructed_variant_row["Start_coordinate_(highest_CFD)"] = "10"
    source_rows = [constructed_variant_row, ordinary_crisprme_row(34)]
    if format_schema == CRISPRME_XLSX_FORMAT_SCHEMA:
        header = sorted(constructed_variant_row)
        sheet_rows = [header] + [[row[field] for field in header] for row in source_rows]
        native = write_xlsx(
            tmp_path / "publication.xlsx",
            [(f"{GUIDE_ID}_integrated_results", sheet_rows)],
        )
        version = "publication workbook; version not encoded"
        expected_source = "CRISPRme publication workbook; version not encoded"
    else:
        native = write_tsv(tmp_path / "run.integrated_results.tsv", source_rows)
        version = "2.3.2"
        expected_source = "CRISPRme+ 2.3.2"
    declaration = declaration_payload(native, format_schema=format_schema)
    declaration["crisprme_version"] = version
    declaration_path = write_declaration(
        tmp_path / "run.crisprme-import.json", declaration
    )
    crisprme_tsv = tmp_path / "crisprme-nominations.tsv"
    crisprme_provenance = tmp_path / "crisprme-provenance.json"
    assert normalize_main(
        [
            "crisprme",
            "--input",
            str(native),
            "--declaration",
            str(declaration_path),
            "--profile",
            str(profile_path),
            "--guides",
            str(guide_path),
            "--intended-loci",
            str(intended_path),
            "--reference-fasta",
            str(reference),
            "--output",
            str(crisprme_tsv),
            "--provenance-output",
            str(crisprme_provenance),
        ]
    ) == 0

    crispor_rows = [crispor_nomination(34, "crispor-both"), crispor_nomination(70, "crispor-only")]
    validate_reference_sequences(crispor_rows, reference)
    crispor_tsv = tmp_path / "crispor-nominations.tsv"
    write_normalized(crispor_rows, crispor_tsv)
    crispor_native = tmp_path / "crispor-native.tsv"
    crispor_native.write_text("synthetic native placeholder\n", encoding="utf-8")
    crispor_provenance = tmp_path / "crispor-provenance.json"
    write_provenance(
        crispor_provenance,
        crispor_tsv,
        crispor_native,
        "synthetic-crispor",
        "CRISPOR",
        "CRISPOR synthetic contract",
        "synthetic test fixture; CRISPOR was not executed",
        "GRCh38",
        reference,
        guide_path,
        intended_path,
        profile_path,
        len(crispor_rows),
        crispor_declaration=crispor_declaration(),
    )

    ledger = build_ledger(
        profile_path,
        guide_path,
        intended_path,
        [crisprme_tsv, crispor_tsv],
        source_provenance_paths=[crisprme_provenance, crispor_provenance],
        reference_fasta_path=reference,
    )

    assert ledger["schema_version"] == "3"
    assert ledger["summary"]["by_source_bucket"] == {
        "both": 1,
        "crispor_only": 1,
        "crisprme_only": 1,
    }
    by_start = {site["start_0based"]: site for site in ledger["candidates"]}
    crisprme_only = by_start[10]
    assert crisprme_only["source_bucket"] == "crisprme_only"
    assert crisprme_only["why_crispor_missed"] == ["variant_created_pam"]
    assert [state["allele"] for state in crisprme_only["allele_states"]] == ["ref", "alt"]
    assert crisprme_only["editability"]["site_summary"] != "cannot_assess"
    assert by_start[34]["source_bucket"] == "both"
    assert by_start[34]["reconciliation"]["status"] == "matched"
    assert by_start[70]["source_bucket"] == "crispor_only"
    assert by_start[70]["reconciliation"]["status"] == "qc_warning_missing_from_crisprme"
    assert {
        site["start_0based"]
        for site in ledger["candidates"]
        if site["site_role"] == "intended"
    } == {100, 5024}
    assert ledger["locus_requirements"][0]["paralog_spacing_bp"] == 4924
    assert ledger["import_summaries"] == [
        {
            "normalized_output": "crisprme-nominations.tsv",
            "input_data_rows": 2,
            "emitted_nomination_rows": 2,
            "skipped_rows": 0,
            "skipped_by_reason": {},
        }
    ]
    provenance = json.loads(crisprme_provenance.read_text())
    assert provenance["schema_version"] == "2"
    assert provenance["source_run"]["version_or_citation"] == expected_source
    with crisprme_tsv.open(encoding="utf-8", newline="") as handle:
        assert {
            row["source_version_or_citation"]
            for row in csv.DictReader(handle, delimiter="\t")
        } == {expected_source}
    assert "SYNTHETIC_SAMPLE_MUST_DISAPPEAR" not in json.dumps(ledger)

    tampered_declaration = json.loads(crisprme_provenance.read_text())
    tampered_declaration["source_declarations"]["crisprme"]["unexpected"] = "value"
    tampered_declaration_path = tmp_path / "tampered-declaration.json"
    tampered_declaration_path.write_text(
        json.dumps(tampered_declaration), encoding="utf-8"
    )
    with pytest.raises(CandidatePanelError, match="invalid CRISPRme source declaration"):
        build_ledger(
            profile_path,
            guide_path,
            intended_path,
            [crisprme_tsv, crispor_tsv],
            source_provenance_paths=[tampered_declaration_path, crispor_provenance],
            reference_fasta_path=reference,
        )

    tampered_summary = json.loads(crisprme_provenance.read_text())
    tampered_summary["import_summary"] = {
        "input_data_rows": 3,
        "emitted_nomination_rows": 2,
        "skipped_rows": 1,
        "skipped_by_reason": {"not_a_supported_reason": 1},
    }
    tampered_summary_path = tmp_path / "tampered-summary.json"
    tampered_summary_path.write_text(json.dumps(tampered_summary), encoding="utf-8")
    with pytest.raises(CandidatePanelError, match="inconsistent CRISPRme import_summary"):
        build_ledger(
            profile_path,
            guide_path,
            intended_path,
            [crisprme_tsv, crispor_tsv],
            source_provenance_paths=[tampered_summary_path, crispor_provenance],
            reference_fasta_path=reference,
        )

    measurement_path = tmp_path / "measurements.tsv"
    measurement_path.write_text(
        "\t".join(
            [
                "guide_id",
                "genome_build",
                "chrom",
                "start_0based",
                "end_0based_exclusive",
                "assay",
                "source_version_or_citation",
                "status",
                "biological_context",
                "source_record_id",
            ]
        )
        + "\n"
        + "\t".join(
            [
                GUIDE_ID,
                "GRCh38",
                "chr2",
                "34",
                "57",
                "synthetic amplicon sequencing",
                "synthetic protocol v1",
                "edited",
                "synthetic cells",
                "measurement-both",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    measured_ledger = build_ledger(
        profile_path,
        guide_path,
        intended_path,
        [crisprme_tsv, crispor_tsv],
        measurement_path=measurement_path,
        source_provenance_paths=[crisprme_provenance, crispor_provenance],
        reference_fasta_path=reference,
    )
    before_states = {
        site["start_0based"]: site["evidence_state"] for site in ledger["candidates"]
    }
    after_states = {
        site["start_0based"]: site["evidence_state"]
        for site in measured_ledger["candidates"]
    }
    assert after_states[34] == "measurement_coverage_unresolved"
    measured_site = next(site for site in measured_ledger["candidates"] if site["start_0based"] == 34)
    assert measured_site["reported_measurement_statuses"] == ["edited"]
    assert {
        start: state for start, state in after_states.items() if start != 34
    } == {
        start: state for start, state in before_states.items() if start != 34
    }


def test_superset_warning_requires_complete_unmerged_retention(tmp_path):
    native = write_tsv(tmp_path / "run.integrated_results.tsv", [base_row()])
    declaration = declaration_payload(native)
    guide = {
        "guide_id": GUIDE_ID,
        "spacer_sequence": GUIDE,
        "pam_class": "NGG",
    }
    site = {"genome_build": "GRCh38"}

    assert _comparison_is_armed(
        site, guide, declaration, crispor_declaration()
    ) is True

    declaration["retention"] = {
        "completeness": "truncated_top_n",
        "top_n": 1000,
        "score_cutoff": None,
        "merge_window_bp": 3,
        "all_candidate_coordinates_retained": False,
    }
    assert _comparison_is_armed(
        site, guide, declaration, crispor_declaration()
    ) is False


def test_alt_selected_between_two_valid_pams_is_not_called_pam_creation():
    metadata = {
        "canonical_interval_basis": "ref",
        "guide_alignment": GUIDE + "NNN",
        "selected_origin": "alt",
        "pam_creation_raw": None,
        "not_found_in_ref_raw": None,
        "allele_states": [
            {"allele": "ref", "alignment": GUIDE + "AGG"},
            {"allele": "alt", "alignment": GUIDE + "TGG"},
        ],
    }
    site = {
        "source_bucket": "crisprme_only",
        "nominations": [
            {
                "nomination_kind": "computational",
                "source_method": "CRISPRme",
                "mismatch_count": 0,
                "dna_bulge": 0,
                "rna_bulge": 0,
                "pam": "TGG",
                "score_name": "CFD",
                "score_value": 1.0,
                "source_metadata": {"crisprme": metadata},
            }
        ],
    }

    reasons = _why_crispor_missed(site, {"pam_class": "NGG"})

    assert "variant_created_pam" not in reasons
    assert reasons == ["unexplained"]


def test_gap_mapping_consumes_guide_and_genome_axes_independently():
    dna_guide = "CT-AACAGTTGCTTTTATCACNNN"
    dna_target = "GTCAACAGTTGCTTTTATCATTGG"
    dna = _assess_allele_window(
        "ref",
        dna_guide,
        dna_target,
        {
            "start_0based": 100,
            "end_0based_exclusive": 124,
            "strand": "+",
        },
        "ABE", 3, 10,
    )
    assert next(base for base in dna["bases"] if base["position"] == 3)[
        "genomic_position_0based"
    ] == 103

    rna_guide = GUIDE + "NNN"
    rna_target = "CT-ACAGTTGCTTTTATCACTGG"
    rna = _assess_allele_window(
        "ref",
        rna_guide,
        rna_target,
        {
            "start_0based": 100,
            "end_0based_exclusive": 122,
            "strand": "+",
        },
        "ABE", 3, 10,
    )
    assert next(base for base in rna["bases"] if base["position"] == 4)[
        "genomic_position_0based"
    ] == 102

    minus = _assess_allele_window(
        "ref",
        GUIDE + "NNN",
        GUIDE + "TGG",
        {
            "start_0based": 100,
            "end_0based_exclusive": 123,
            "strand": "-",
        },
        "ABE", 3, 10,
    )
    assert next(base for base in minus["bases"] if base["position"] == 3)[
        "genomic_position_0based"
    ] == 120


GOLDEN_ZIP_SHA256 = "a48f261b3887af7285cdd03c345ce7cd2ddcb2e27c12d2d2c448cf7db12c76af"


@pytest.mark.skipif(
    not os.environ.get("BES_CRISPRME_GOLDEN_ZIP"),
    reason="set BES_CRISPRME_GOLDEN_ZIP for the opt-in deposited-workbook test",
)
def test_opt_in_legacy_deposited_sg1617_workbook(tmp_path):
    archive = Path(os.environ["BES_CRISPRME_GOLDEN_ZIP"])
    assert sha256(archive.read_bytes()).hexdigest() == GOLDEN_ZIP_SHA256
    with ZipFile(archive) as zipped:
        member = "File_S1_1617_search_results.xlsx"
        workbook = tmp_path / member
        workbook.write_bytes(zipped.read(member))

    names = sheet_names(workbook)
    assert "1617_integrated_results" in names
    assert "1617_altMerge" in names
    assert any(name.endswith("_private_targets") for name in names)
    source_rows = iter(iter_xlsx_rows(workbook, "1617_integrated_results"))
    header = next(source_rows)
    field_columns = {field: index for index, field in enumerate(header) if field in REQUIRED_COLUMNS}
    assert set(field_columns) == REQUIRED_COLUMNS
    native_rows = [
        {field: values[index] if index < len(values) else "" for field, index in field_columns.items()}
        for values in source_rows
    ]
    integrated_data_rows = len(native_rows)
    # Select a source example by sequence properties. Expected values stay in the local workbook.
    source_example = next(row for row in native_rows if (
        row["REF/ALT_origin_(highest_CFD)"].lower() == "alt"
        and row["Bulges_(highest_CFD)"] == "0"
        and len(row["Aligned_protospacer+PAM_REF_(highest_CFD)"]) == 23
        and len(row["Aligned_protospacer+PAM_ALT_(highest_CFD)"]) == 23
        and set(row["Aligned_protospacer+PAM_REF_(highest_CFD)"].upper()) <= set("ACGT")
        and set(row["Aligned_protospacer+PAM_ALT_(highest_CFD)"].upper()) <= set("ACGT")
        and row["Aligned_protospacer+PAM_ALT_(highest_CFD)"].upper().endswith("GG")
        and not row["Aligned_protospacer+PAM_REF_(highest_CFD)"].upper().endswith("GG")
    ))
    legacy_spacer = source_example["Spacer+PAM"][:-3].upper()

    payload = declaration_payload(
        workbook,
        format_schema=CRISPRME_XLSX_FORMAT_SCHEMA,
        guide_id="1617",
        guide=legacy_spacer,
    )
    payload["crisprme_version"] = "publication workbook; version not encoded"
    payload["retention"] = {
        "completeness": "truncated_top_n",
        "top_n": 1000,
        "score_cutoff": None,
        "merge_window_bp": 3,
        "all_candidate_coordinates_retained": False,
    }
    declaration = write_declaration(tmp_path / "sg1617.json", payload)

    result = import_crisprme(workbook, declaration)

    assert result.declaration["selected_sheet"] == "1617_integrated_results"
    assert result.declaration["retention"]["completeness"] == "truncated_top_n"
    assert result.sheets_ignored == {
        "count": 2,
        "categories": ["altMerge", "private_targets"],
    }
    assert result.privacy_filtering_applied is True
    assert result.import_summary["input_data_rows"] == integrated_data_rows
    assert result.import_summary["input_data_rows"] == 999
    assert result.import_summary["emitted_nomination_rows"] == 999
    source_start = int(source_example["Start_coordinate_(highest_CFD)"])
    source_ref = source_example["Aligned_protospacer+PAM_REF_(highest_CFD)"].upper()
    imported_examples = [row for row in result.rows if (
        row["chrom"], row["start_0based"], row["end_0based_exclusive"]
    ) == (source_example["Chromosome"], source_start, source_start + len(source_ref))]
    assert len(imported_examples) == 1
    imported_example = imported_examples[0]
    assert imported_example["strand"] == source_example["Strand_(highest_CFD)"]
    assert imported_example["mismatch_count"] == int(source_example["Mismatches_(highest_CFD)"])
    assert [state["allele"] for state in imported_example["allele_states"]] == ["ref", "alt"]
    assert [state["alignment"].upper() for state in imported_example["allele_states"]] == [
        source_example[f"Aligned_protospacer+PAM_{allele}_(highest_CFD)"].upper()
        for allele in ("REF", "ALT")
    ]
    assert [state["CFD"] for state in imported_example["allele_states"]] == [
        float(source_example[f"CFD_score_{allele}_(highest_CFD)"]) for allele in ("REF", "ALT")
    ]
    assert imported_example["variant_maf_raw"] == source_example["Variant_MAF_(highest_CFD)"]
    assert imported_example["variant_rsid_raw"] == source_example["Variant_rsID_(highest_CFD)"]
    rendered = json.dumps(result.rows)
    assert "Variant_samples" not in rendered
    assert "1617_altMerge" not in rendered
    assert not any(name in rendered for name in names if name.endswith("_private_targets"))

    # Exercise reconciliation with the locally read example and its actual truncated declaration.
    # The synthetic build label avoids requiring a local GRCh38 FASTA for this logic-only check.
    guide_id = "1617"
    guide_path = tmp_path / "golden-guides.tsv"
    guide_path.write_text(
        "guide_id\tspacer_sequence\teditor\teditor_class\tpam_class\t"
        "editing_window_start_1based\tediting_window_end_1based_inclusive\t"
        "intended_target_description\tlocus_system\n"
        f"{guide_id}\t{legacy_spacer}\tABE8e\tABE\tNGG\t4\t8\tgolden fixture\tother\n",
        encoding="utf-8",
    )
    intended_path = tmp_path / "golden-intended.tsv"
    intended_path.write_text(
        "guide_id\tgenome_build\tchrom\tstart_0based\tend_0based_exclusive\t"
        "strand\tlocus_name\n"
        f"{guide_id}\tSYNTHETIC\tchr2\t100\t123\t+\tcontrol\n",
        encoding="utf-8",
    )
    profile_path = tmp_path / "golden-profile.json"
    profile_path.write_text(
        json.dumps(
            {
                "name": "golden_truncated_fixture",
                "max_mismatches": 6,
                "max_dna_bulge": 2,
                "max_rna_bulge": 2,
                "pam_patterns": ["NGG"],
                "variant_aware": True,
                "intended_use": "Golden truncated-fixture reconciliation check.",
                "reference_id": "SYNTHETIC",
                "reference_provenance": "Logic-only coordinate fixture.",
                "variant_source": "Deposited workbook fields.",
                "retention_rule": "Retain normalized fixture rows.",
                "retention_score_name": None,
                "minimum_score": None,
            }
        ),
        encoding="utf-8",
    )
    guide = {
        "guide_id": guide_id,
        "spacer_sequence": legacy_spacer,
        "editor": "ABE8e",
        "editor_class": "ABE",
        "pam_class": "NGG",
        "editing_window_start_1based": 4,
        "editing_window_end_1based_inclusive": 8,
        "intended_target_description": "golden fixture",
        "locus_system": "other",
    }
    intended = {
        "guide_id": guide_id,
        "genome_build": "SYNTHETIC",
        "chrom": "chr2",
        "start_0based": 100,
        "end_0based_exclusive": 123,
        "strand": "+",
        "locus_name": "control",
    }
    fixture_declaration = json.loads(json.dumps(result.declaration))
    fixture_declaration["genome_build"] = "SYNTHETIC"
    fixture_row = {**imported_example, "genome_build": "SYNTHETIC"}
    fixture_import = CrisprmeImportResult(
        rows=(fixture_row,),
        declaration=fixture_declaration,
        import_summary={
            "input_data_rows": 1,
            "emitted_nomination_rows": 1,
            "skipped_rows": 0,
            "skipped_by_reason": {},
        },
        sheets_ignored=result.sheets_ignored,
        privacy_filtering_applied=True,
    )
    crisprme_rows = normalize_crisprme(
        fixture_import,
        [guide],
        [intended],
        SearchProfile.load(profile_path),
    )
    crisprme_tsv = tmp_path / "golden-crisprme.tsv"
    write_normalized(crisprme_rows, crisprme_tsv)

    crispor_row = crispor_nomination(70, "golden-crispor-only")
    crispor_row.update(
        {
            "guide_id": guide_id,
            "genome_build": "SYNTHETIC",
            "search_profile": "golden_truncated_fixture",
            "site_sequence": legacy_spacer + "TGG",
        }
    )
    crispor_tsv = tmp_path / "golden-crispor.tsv"
    write_normalized([crispor_row], crispor_tsv)
    crispor_native = tmp_path / "golden-crispor-native.tsv"
    crispor_native.write_text("synthetic comparator\n", encoding="utf-8")

    def digest(path: Path) -> str:
        return sha256(path.read_bytes()).hexdigest()

    common = {
        "record_type": "nomination_source_provenance",
        "guide_manifest": {"file": guide_path.name, "sha256": digest(guide_path)},
        "intended_locus_manifest": {
            "file": intended_path.name,
            "sha256": digest(intended_path),
        },
        "search_profile": {"file": profile_path.name, "sha256": digest(profile_path)},
        "reference_fasta": None,
        "generated_at_utc": "2026-08-31T00:00:00+00:00",
    }
    crisprme_provenance = tmp_path / "golden-crisprme-provenance.json"
    crisprme_provenance.write_text(
        json.dumps(
            {
                **common,
                "schema_version": "2",
                "normalized_output": {
                    "file": crisprme_tsv.name,
                    "sha256": digest(crisprme_tsv),
                    "row_count": 1,
                },
                "native_input": fixture_declaration["native_input"],
                "source_run": {
                    "method": "CRISPRme",
                    "version_or_citation": (
                        f"CRISPRme {fixture_declaration['crisprme_version']}"
                    ),
                    "command": fixture_declaration["source_command"],
                    "reference_id": "SYNTHETIC",
                    "reference_sha256": None,
                    "reference_checksum_status": "not_applicable_synthetic",
                },
                "adapter": {"mode": "crisprme", "base_edit_safety_version": __version__},
                "source_declarations": {"crisprme": fixture_declaration},
                "import_summary": fixture_import.import_summary,
                "sheets_ignored": fixture_import.sheets_ignored,
                "privacy_filtering_applied": True,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    complete_crispor_declaration = crispor_declaration()
    complete_crispor_declaration.update(
        {"guide_id": guide_id, "genome_build": "SYNTHETIC", "guide_protospacer": legacy_spacer}
    )
    crispor_provenance = tmp_path / "golden-crispor-provenance.json"
    crispor_provenance.write_text(
        json.dumps(
            {
                **common,
                "schema_version": "1",
                "normalized_output": {
                    "file": crispor_tsv.name,
                    "sha256": digest(crispor_tsv),
                    "row_count": 1,
                },
                "native_input": {
                    "file": crispor_native.name,
                    "sha256": digest(crispor_native),
                },
                "source_run": {
                    "method": "CRISPOR",
                    "version_or_citation": "synthetic golden comparator",
                    "command": "logic-only comparator",
                    "reference_id": "SYNTHETIC",
                    "reference_sha256": None,
                    "reference_checksum_status": "not_applicable_synthetic",
                },
                "adapter": {
                    "mode": "synthetic-crispor",
                    "base_edit_safety_version": __version__,
                },
                "source_declarations": {
                    "crispor": complete_crispor_declaration
                },
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    ledger = build_ledger(
        profile_path,
        guide_path,
        intended_path,
        [crisprme_tsv, crispor_tsv],
        source_provenance_paths=[crisprme_provenance, crispor_provenance],
    )
    crispor_only = next(
        site for site in ledger["candidates"] if site["source_bucket"] == "crispor_only"
    )
    assert crispor_only["reconciliation"] == {
        "status": "comparison_not_armed",
        "warnings": [],
    }
    assert not any(
        site["reconciliation"]["status"] == "qc_warning_missing_from_crisprme"
        for site in ledger["candidates"]
    )


def test_m2_editability_uses_the_declared_window():
    from base_edit_safety.pipeline.candidate_panel import _assess_site_windows

    spacer = "AGATATTTGCATTGAGATAG"
    site = {
        "chrom": "chrSYN", "start_0based": 100, "end_0based_exclusive": 123,
        "nominations": [{
            "nomination_kind": "computational", "source_method": "CRISPOR",
            "source_metadata": {}, "site_sequence": spacer + "TGT",
            "pam": "TGT", "strand": "+", "score_name": None,
        }],
    }
    guide = {
        "spacer_sequence": spacer, "editor_class": "ABE",
        "editing_window_start_1based": 4,
        "editing_window_end_1based_inclusive": 10,
    }
    result = _assess_site_windows(site, guide)
    assert [base["position"] for base in result["allele_states"][0]["bases"]] == [5]
    guide["editing_window_start_1based"] = 11
    guide["editing_window_end_1based_inclusive"] = 13
    result = _assess_site_windows(site, guide)
    assert [base["position"] for base in result["allele_states"][0]["bases"]] == [11]
    guide["editor_class"] = "CBE"
    result = _assess_site_windows(site, guide)
    assert result["site_summary"] == "no_base_in_declared_window"


def test_m2_unknown_target_base_keeps_editability_indeterminate():
    result = _assess_allele_window(
        "ref", GUIDE + "NNN", GUIDE[:4] + "N" + GUIDE[5:] + "TGG",
        {"start_0based": 100, "end_0based_exclusive": 123, "strand": "+"},
        "ABE", 4, 10,
    )
    assert result["summary"] == "cannot_assess"
    assert result["reason"] == "unknown_sequence"


@pytest.mark.parametrize("flag", ["False", "false", "0", "no", "n"])
def test_m2_false_native_flags_do_not_create_variant_explanations(flag):
    metadata = {
        "canonical_interval_basis": "ref", "guide_alignment": GUIDE + "NNN",
        "selected_origin": "ref", "pam_creation_raw": flag,
        "not_found_in_ref_raw": flag,
        "allele_states": [{"allele": "ref", "alignment": GUIDE + "TGG"}],
    }
    site = {"source_bucket": "crisprme_only", "nominations": [{
        "source_method": "CRISPRme", "nomination_kind": "computational",
        "mismatch_count": 0, "dna_bulge": 0, "rna_bulge": 0, "pam": "TGG",
        "source_metadata": {"crisprme": metadata},
    }]}
    assert _why_crispor_missed(site, {"pam_class": "NGG"}) == ["unexplained"]


def test_m2_skipped_rows_downgrade_completeness_without_rewriting_source(tmp_path):
    bad = base_row()
    bad["Start_coordinate_(highest_CFD)"] = "invalid"
    native = write_tsv(tmp_path / "incomplete.integrated_results.tsv", [base_row(), bad])
    original = declaration_payload(native)
    declaration = write_declaration(tmp_path / "incomplete.json", original)
    result = import_crisprme(native, declaration)
    assert result.declaration["retention"]["completeness"] == "filtered"
    assert result.declaration["retention"]["all_candidate_coordinates_retained"] is False
    assert result.original_declaration == original
    assert json.loads(declaration.read_text()) == original
    assert result.import_summary["skipped_by_reason"] == {"invalid_coordinate": 1}


def test_m2_distinct_alt_configurations_survive_import_and_ledger(tmp_path):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    first = ordinary_crisprme_row(34)
    second = ordinary_crisprme_row(34)
    first["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = GUIDE[:4] + "G" + GUIDE[5:] + "TGG"
    second["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = GUIDE[:4] + "T" + GUIDE[5:] + "TGG"
    native = write_tsv(tmp_path / "configurations.integrated_results.tsv", [first, second])
    declaration = write_declaration(tmp_path / "configurations.json", declaration_payload(native))
    imported = import_crisprme(native, declaration)
    assert imported.import_summary["emitted_nomination_rows"] == 2
    assert imported.import_summary["skipped_rows"] == 0
    assert len({row["source_record_id"] for row in imported.rows}) == 2
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    guides = load_guides(guide_path)
    rows = normalize_crisprme(imported, guides, load_intended_loci(intended_path), SearchProfile.load(profile_path))
    site = merge_evidence(rows, [], {}, guides=guides)[0]
    assert len(site["allele_configurations"]) == 2
    assert {state["alignment"] for state in site["allele_states"] if state["allele"] == "alt"} == {
        first["Aligned_protospacer+PAM_ALT_(highest_CFD)"],
        second["Aligned_protospacer+PAM_ALT_(highest_CFD)"],
    }
    assert len(site["editability"]["allele_states"]) == 4


@pytest.mark.parametrize("strand", ["+", "-"])
def test_m2_inserted_alt_bases_have_no_invented_reference_interval(tmp_path, strand):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    row = ordinary_crisprme_row(34)
    row.update({
        "Strand_(highest_CFD)": strand,
        "Aligned_spacer+PAM_(highest_CFD)": GUIDE[:2] + "-" + GUIDE[2:] + "NNN",
        "Aligned_protospacer+PAM_REF_(highest_CFD)": GUIDE[:2] + "-" + GUIDE[2:] + "TGG",
        "Aligned_protospacer+PAM_ALT_(highest_CFD)": GUIDE[:2] + "A" + GUIDE[2:] + "TGG",
        "REF/ALT_origin_(highest_CFD)": "alt",
        "Bulge_type_(highest_CFD)": "DNA", "Bulges_(highest_CFD)": "1",
    })
    native = write_tsv(tmp_path / "inserted.integrated_results.tsv", [row])
    declaration = write_declaration(tmp_path / "inserted.json", declaration_payload(native))
    imported = import_crisprme(native, declaration)
    ref, alt = imported.rows[0]["allele_states"]
    assert ref["effective_interval"]["end_0based_exclusive"] == 57
    assert alt["effective_interval"] is None
    assert alt["mapping_status"] == "unsupported_indel_alignment"
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    guides = load_guides(guide_path)
    rows = normalize_crisprme(imported, guides, load_intended_loci(intended_path), SearchProfile.load(profile_path))
    site = merge_evidence(rows, [], {}, guides=guides)[0]
    result = next(s for s in site["editability"]["allele_states"] if s["allele"] == "alt")
    assert result["summary"] == "cannot_assess"
    assert result["bases"] == []


@pytest.mark.parametrize(
    "strand, position, ref, alt",
    [("+", 1022, "C", "G"), ("-", 1000, "G", "C")],
)
def test_m2_allele_changes_are_genomic_forward(tmp_path, strand, position, ref, alt):
    row = base_row()
    row["Strand_(highest_CFD)"] = strand
    native = write_tsv(tmp_path / "forward.integrated_results.tsv", [row])
    declaration = write_declaration(tmp_path / "forward.json", declaration_payload(native))
    states = import_crisprme(native, declaration).rows[0]["allele_states"]
    assert states[0]["genomic_forward_changes"] == []
    assert states[1]["genomic_forward_changes"] == [{
        "chrom": "chr2", "position_0based": position, "ref": ref, "alt": alt,
    }]


def test_m2_alt_only_row_retains_only_a_reported_start_anchor(tmp_path):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci, load_nominations

    row = ordinary_crisprme_row(34)
    row["Aligned_protospacer+PAM_REF_(highest_CFD)"] = "NA"
    row["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = GUIDE + "TGG"
    row["REF/ALT_origin_(highest_CFD)"] = "alt"
    native = write_tsv(tmp_path / "alt-only.integrated_results.tsv", [row])
    declaration = write_declaration(tmp_path / "alt-only.json", declaration_payload(native))
    imported = import_crisprme(native, declaration)
    assert imported.rows[0]["end_0based_exclusive"] == 35
    assert imported.rows[0]["allele_states"][0]["effective_interval"] is None
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    profile = SearchProfile.load(profile_path)
    guides = load_guides(guide_path)
    rows = normalize_crisprme(imported, guides, load_intended_loci(intended_path), profile)
    assert rows[0]["site_sequence"] == ""
    validate_reference_sequences(rows, write_reference(tmp_path))
    assert rows[0]["reference_validation"]["status"] == "reported_start_anchor_only"
    normalized = tmp_path / "alt-only.normalized.tsv"
    write_normalized(rows, normalized)
    loaded = load_nominations([normalized], profile)
    site = merge_evidence(loaded, [], {}, guides=guides)[0]
    assert site["editability"]["site_summary"] == "cannot_assess"


def test_m2_unknown_reference_base_is_not_fully_verified(tmp_path):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    row = ordinary_crisprme_row(34)
    row["Aligned_protospacer+PAM_REF_(highest_CFD)"] = "N" + GUIDE[1:] + "TGG"
    native = write_tsv(tmp_path / "unknown.integrated_results.tsv", [row])
    declaration = write_declaration(tmp_path / "unknown.json", declaration_payload(native))
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    guides = load_guides(guide_path)
    rows = normalize_crisprme(import_crisprme(native, declaration), guides,
                             load_intended_loci(intended_path), SearchProfile.load(profile_path))
    validate_reference_sequences(rows, write_reference(tmp_path))
    assert rows[0]["reference_validation"]["status"] == "partial_unknown_sequence"
    site = merge_evidence(rows, [], {}, guides=guides)[0]
    assert site["editability"]["site_summary"] == "cannot_assess"


@pytest.mark.parametrize("raw_maf", ["0.8", "0.8,0.1", "NA"])
def test_m2_frequency_stays_opaque_without_an_allele_mapping(tmp_path, raw_maf):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    row = ordinary_crisprme_row(34)
    row["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = GUIDE[:4] + "A" + GUIDE[5:] + "TGG"
    row["Variant_MAF_(highest_CFD)"] = raw_maf
    row["Variant_info_genome_(highest_CFD)"] = "chr2_39_C_A,chr2_39_C_G"
    native = write_tsv(tmp_path / "frequency.integrated_results.tsv", [row])
    payload = declaration_payload(native)
    declaration = write_declaration(tmp_path / "frequency.json", payload)
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    guides = load_guides(guide_path)
    rows = normalize_crisprme(import_crisprme(native, declaration), guides,
                             load_intended_loci(intended_path), SearchProfile.load(profile_path))
    provenance = [{"schema_version": "2", "source_declarations": {"crisprme": payload}}]
    site = merge_evidence(rows, [], {}, guides=guides, source_provenance=provenance)[0]
    frequency = site["allele_configurations"][0]["frequency"]
    assert frequency == {
        "raw": None if raw_maf == "NA" else raw_maf,
        "semantics": "uninterpreted_upstream_MAF", "allele": None, "value": None,
        "population": None, "AC": None, "AN": None, "callability": None,
        "sources": payload["variant_databases"],
    }
    assert site["allele_states"][1]["variant_info_raw"]["genome"] == row["Variant_info_genome_(highest_CFD)"]


@pytest.mark.parametrize("calls", ["./.", "0/.", "0/1,0/1", "0|1,0|1", "0|1,1|0"])
def test_m2_missing_calls_and_cis_trans_are_not_inferred(tmp_path, calls):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    row = ordinary_crisprme_row(34)
    row["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = GUIDE[:4] + "A" + GUIDE[5:7] + "A" + GUIDE[8:] + "TGG"
    row["Variant_samples_(highest_CFD)"] = "SYNTHETIC_CALL_CANARY:" + calls
    row["Variant_MAF_(highest_CFD)"] = "0.5,0.5"
    native = write_tsv(tmp_path / "phase.integrated_results.tsv", [row])
    declaration = write_declaration(tmp_path / "phase.json", declaration_payload(native))
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    guides = load_guides(guide_path)
    imported = import_crisprme(native, declaration)
    rows = normalize_crisprme(imported, guides, load_intended_loci(intended_path),
                             SearchProfile.load(profile_path))
    site = merge_evidence(rows, [], {}, guides=guides)[0]
    block = site["allele_configurations"][0]
    assert block["phase_status"] == "not_provided_in_integrated_contract"
    assert block["genotype_status"] == "not_imported"
    assert block["configuration_coverage"] == "reported_highest_CFD_only"
    assert "SYNTHETIC_CALL_CANARY" not in json.dumps(site)
    alt = next(s for s in site["editability"]["allele_states"] if s["allele"] == "alt")
    assert alt["reason"] == "phase_not_provided_for_multiple_changes"
    assert alt["bases"] == []
    assert site["editability"]["site_summary"] == "cannot_assess"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("ref_pam,alt_pam,origin,expected", [
    ("TCG", "TGG", "alt", "created"),
    ("TGG", "TCG", "ref", "destroyed"),
    ("TGG", "AGG", "alt", "retained"),
    ("TNG", "TGG", "alt", "unknown_sequence"),
])
def test_m2_pam_changes_require_known_alleles(tmp_path, strand, ref_pam, alt_pam, origin, expected):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    row = ordinary_crisprme_row(34)
    row.update({
        "Strand_(highest_CFD)": strand,
        "Aligned_protospacer+PAM_REF_(highest_CFD)": GUIDE + ref_pam,
        "Aligned_protospacer+PAM_ALT_(highest_CFD)": GUIDE + alt_pam,
        "REF/ALT_origin_(highest_CFD)": origin,
        "PAM_(highest_CFD)": ref_pam if origin == "ref" else alt_pam,
    })
    native = write_tsv(tmp_path / "pam.integrated_results.tsv", [row])
    declaration = write_declaration(tmp_path / "pam.json", declaration_payload(native))
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    guides = load_guides(guide_path)
    rows = normalize_crisprme(import_crisprme(native, declaration), guides,
                             load_intended_loci(intended_path), SearchProfile.load(profile_path))
    site = merge_evidence(rows, [], {}, guides=guides)[0]
    assert site["allele_configurations"][0]["pam_change"] == expected
    assert ("variant_created_pam" in site["why_crispor_missed"]) == (expected == "created")
    if expected == "unknown_sequence":
        assert site["editability"]["site_summary"] == "cannot_assess"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_m2_deleted_alt_retains_sequence_but_no_reference_mapping(tmp_path, strand):
    row = ordinary_crisprme_row(34)
    row["Strand_(highest_CFD)"] = strand
    row["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = GUIDE[:4] + "-" + GUIDE[5:] + "TGG"
    native = write_tsv(tmp_path / "deletion.integrated_results.tsv", [row])
    declaration = write_declaration(tmp_path / "deletion.json", declaration_payload(native))
    states = import_crisprme(native, declaration).rows[0]["allele_states"]
    assert states[0]["effective_interval"]["end_0based_exclusive"] == 57
    assert states[1]["alignment"] == row["Aligned_protospacer+PAM_ALT_(highest_CFD)"]
    assert states[1]["effective_interval"] is None
    assert states[1]["genomic_forward_changes"] is None
    assert states[1]["mapping_status"] == "unsupported_indel_alignment"


def test_m2_unresolved_alternative_alignments_are_not_exact_edits(tmp_path):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    rows = []
    for gap, mismatches in [(2, 0), (3, 1)]:
        row = ordinary_crisprme_row(34)
        row.update({
            "Aligned_spacer+PAM_(highest_CFD)": GUIDE[:gap] + "-" + GUIDE[gap:] + "NNN",
            "Aligned_protospacer+PAM_REF_(highest_CFD)": "GTCAACAGTTGCTTTTATCATTGG",
            "Mismatches_(highest_CFD)": str(mismatches),
            "Bulge_type_(highest_CFD)": "DNA", "Bulges_(highest_CFD)": "1",
        })
        rows.append(row)
    native = write_tsv(tmp_path / "alignments.integrated_results.tsv", rows)
    declaration = write_declaration(tmp_path / "alignments.json", declaration_payload(native))
    imported = import_crisprme(native, declaration)
    assert imported.import_summary["emitted_nomination_rows"] == 2
    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    guides = load_guides(guide_path)
    normalized = normalize_crisprme(imported, guides, load_intended_loci(intended_path),
                                   SearchProfile.load(profile_path))
    site = merge_evidence(normalized, [], {}, guides=guides)[0]
    assert len(site["allele_configurations"]) == 2
    assert {s["reason"] for s in site["editability"]["allele_states"]} == {"unresolved_alternative_alignment"}
    assert all(s["bases"] == [] for s in site["editability"]["allele_states"])


def test_m2_extra_query_bases_are_rejected_and_counted(tmp_path):
    malformed = ordinary_crisprme_row(70)
    malformed["Aligned_spacer+PAM_(highest_CFD)"] = GUIDE + "NNNN"
    malformed["Aligned_protospacer+PAM_REF_(highest_CFD)"] = GUIDE + "ATGG"
    native = write_tsv(tmp_path / "extra.integrated_results.tsv", [ordinary_crisprme_row(34), malformed])
    declaration = write_declaration(tmp_path / "extra.json", declaration_payload(native))
    imported = import_crisprme(native, declaration)
    assert imported.import_summary["skipped_by_reason"] == {"invalid_alignment": 1}
    assert imported.import_summary["input_data_rows"] == 2
    assert imported.import_summary["emitted_nomination_rows"] == 1


def test_m2_report_shows_partial_import_and_unresolved_anchor(tmp_path):
    from base_edit_safety.pipeline.candidate_report import render_report

    anchor = ordinary_crisprme_row(34)
    anchor["Aligned_protospacer+PAM_REF_(highest_CFD)"] = "NA"
    anchor["Aligned_protospacer+PAM_ALT_(highest_CFD)"] = GUIDE + "TGG"
    anchor["REF/ALT_origin_(highest_CFD)"] = "alt"
    bad = ordinary_crisprme_row(70)
    bad["Start_coordinate_(highest_CFD)"] = "bad"
    native = write_tsv(tmp_path / "partial.integrated_results.tsv", [anchor, bad])
    original = declaration_payload(native)
    declaration = write_declaration(tmp_path / "partial.json", original)
    guides, loci, profile = write_manifest_inputs(tmp_path)
    reference = write_reference(tmp_path)
    normalized, provenance = tmp_path / "partial.tsv", tmp_path / "partial.source.json"
    assert normalize_main([
        "crisprme", "--input", str(native), "--declaration", str(declaration),
        "--guides", str(guides), "--intended-loci", str(loci), "--profile", str(profile),
        "--reference-fasta", str(reference), "--output", str(normalized),
        "--provenance-output", str(provenance),
    ]) == 0
    recorded = json.loads(provenance.read_text())
    assert recorded["original_source_declaration"] == original
    assert recorded["source_declarations"]["crisprme"]["retention"]["completeness"] == "filtered"
    ledger = build_ledger(profile, guides, loci, [normalized],
                          source_provenance_paths=[provenance], reference_fasta_path=reference)
    designs = {"record_type": "hbg_design_record_set", "schema_version": "7", "records": [{
        "record_type": "hbg_design_record", "schema_version": "7", "candidate": {"guide_id": GUIDE_ID},
    }]}
    report = render_report(designs, ledger)
    assert "Incomplete import" in report
    assert "rejected/skipped: 1" in report
    assert "reported start anchor; span unresolved" in report
    assert "highest-CFD" in report
    assert "not_provided_in_integrated_contract" in report
    assert "uninterpreted_upstream_MAF" in report
    assert "reported_sequence_only" in report


@pytest.mark.parametrize("fault", ["build", "strand", "coordinate", "guide"])
def test_m2_reference_identity_errors_cannot_be_normalized(tmp_path, fault):
    from base_edit_safety.pipeline.candidate_panel import load_guides, load_intended_loci

    row = base_row()
    row["Start_coordinate_(highest_CFD)"] = "10"
    if fault == "strand":
        row["Strand_(highest_CFD)"] = "+"
    elif fault == "coordinate":
        row["Start_coordinate_(highest_CFD)"] = "11"
    elif fault == "guide":
        row["Aligned_spacer+PAM_(highest_CFD)"] = "A" + GUIDE[1:] + "NNN"
    native = write_tsv(tmp_path / "identity.integrated_results.tsv", [row])
    payload = declaration_payload(native)
    if fault == "build":
        payload["genome_build"] = "GRCh37"
    declaration = write_declaration(tmp_path / "identity.json", payload)
    guides, loci, profile = write_manifest_inputs(tmp_path)
    with pytest.raises((CandidatePanelError, CrisprmeImportError)):
        normalized = normalize_crisprme(import_crisprme(native, declaration), load_guides(guides),
                                       load_intended_loci(loci), SearchProfile.load(profile))
        validate_reference_sequences(normalized, write_reference(tmp_path))


def test_m2_swapped_hbg_paralog_labels_are_rejected(tmp_path):
    from base_edit_safety.pipeline.candidate_panel import _validate_intended_loci, load_guides, load_intended_loci

    guide_path, intended_path, profile_path = write_manifest_inputs(tmp_path)
    loci = load_intended_loci(intended_path)
    loci[0]["locus_name"], loci[1]["locus_name"] = loci[1]["locus_name"], loci[0]["locus_name"]
    with pytest.raises(CandidatePanelError, match="4924"):
        _validate_intended_loci(load_guides(guide_path), loci, SearchProfile.load(profile_path))


@pytest.mark.parametrize("editor,window,genomic_base,coordinate", [
    ("ABE", 3, "T", 120), ("CBE", 1, "G", 122),
])
def test_m2_minus_strand_substrate_names_its_genomic_forward_base(editor, window, genomic_base, coordinate):
    guide = "C" + GUIDE[1:]  # Keep C at position 1 and A at position 3 for these strand checks.
    result = _assess_allele_window(
        "ref", guide + "NNN", guide + "TGG",
        {"start_0based": 100, "end_0based_exclusive": 123, "strand": "-"},
        editor, window, window,
    )
    assert len(result["bases"]) == 1
    assert result["bases"][0]["genomic_forward_base"] == genomic_base
    assert result["bases"][0]["genomic_position_0based"] == coordinate


def test_m2_one_explained_configuration_does_not_hide_an_unexplained_one():
    ref = "A" + GUIDE[1:] + "TGG"
    nominations = []
    for alt in (GUIDE + "TGG", "T" + GUIDE[1:] + "TGG"):
        nominations.append({
            "source_method": "CRISPRme", "mismatch_count": 1, "pam": "TGG",
            "source_metadata": {"crisprme": {
                "guide_alignment": GUIDE + "NNN", "canonical_interval_basis": "ref",
                "not_found_in_ref_raw": "False", "allele_states": [
                    {"allele": "ref", "alignment": ref}, {"allele": "alt", "alignment": alt},
                ],
            }},
        })
    reasons = _why_crispor_missed({"source_bucket": "crisprme_only", "nominations": nominations},
                                 {"pam_class": "NGG"})
    assert reasons == ["unexplained", "variant_reduced_mm"]


def test_m3_measurement_cli_keeps_hbg_copies_and_locus_checks_separate(tmp_path):
    """Synthetic measurements exercise the full CLI; no experimental result is asserted."""
    from base_edit_safety.pipeline.candidate_panel import main as panel_main
    from base_edit_safety.pipeline.candidate_report import main as report_main

    guides, loci, profile = write_manifest_inputs(tmp_path)
    reference = write_reference(tmp_path)
    native = write_tsv(tmp_path / "synthetic.integrated_results.tsv",
                       [ordinary_crisprme_row(100), ordinary_crisprme_row(5024)])
    declaration = write_declaration(tmp_path / "source.json", declaration_payload(native))
    normalized, provenance = tmp_path / "normalized.tsv", tmp_path / "provenance.json"
    shared = ["--guides", str(guides), "--intended-loci", str(loci), "--profile", str(profile),
              "--reference-fasta", str(reference)]
    assert normalize_main(["crisprme", "--input", str(native), "--declaration", str(declaration),
                           "--output", str(normalized), "--provenance-output", str(provenance), *shared]) == 0
    observed = {"guide_id": GUIDE_ID, "genome_build": "GRCh38", "chrom": "chr2", "start_0based": 100,
                "end_0based_exclusive": 123, "strand": "+", "assay": "synthetic amplicon assay",
                "source_version_or_citation": "synthetic fixture", "status": "not_detected",
                "limit_of_detection_value": 0.1, "limit_of_detection_unit": "percent",
                "biological_context": "synthetic cells", "source_record_id": "synthetic_reference_only",
                "editor": "ABE8e", "measurement_kind": "editing", "allele": "ref",
                "allele_sequence": GUIDE + "TGG", "allele_resolution": "resolved",
                "allele_evidence": "synthetic sequenced template", "qc_status": "passed",
                "treated_group": "treated", "control_group": "untreated"}
    measurements = tmp_path / "measurements.tsv"
    with measurements.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(observed), delimiter='\t')
        writer.writeheader()
        writer.writerow(observed)
    ledger_path, table = tmp_path / "ledger.json", tmp_path / "review.tsv"
    assert panel_main([*shared, "--nominations", str(normalized), "--source-provenance", str(provenance),
                       "--measurements", str(measurements), "--output", str(ledger_path),
                       "--panel-output", str(table)]) == 0
    ledger = json.loads(ledger_path.read_text())
    assert ledger["summary"]["measurement_rows"] == 1
    by_start = {site["start_0based"]: site for site in ledger["candidates"]}
    assert by_start[100]["measurement_summary"] == "not_detected"
    assert by_start[5024]["measurement_summary"] == "not_measured"
    for site in ledger["candidates"]:
        assert len(site["required_confirmation"]) == 2
        assert "locus-resolving" in " ".join(site["follow_up_questions"])
        assert site["panel_decision"]["status"] == "undecided"
    designs = tmp_path / "designs.json"
    designs.write_text(json.dumps({"record_type": "hbg_design_record_set", "schema_version": "7",
        "records": [{"record_type": "hbg_design_record", "schema_version": "7",
                     "candidate": {"guide_id": GUIDE_ID}, "evidence_tier": "predicted_only"}]}))
    report = tmp_path / "review.md"
    assert report_main(["--design-records", str(designs), "--candidate-ledger", str(ledger_path),
                        "--output", str(report)]) == 0
    assert "Untested for this DNA version" in report.read_text()
    assert "HBG locus checks" in report.read_text()
    assert "synthetic_reference_only" in report.read_text()
