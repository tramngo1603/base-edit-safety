"""Network-free fixture tests for the D1/D3/D4/D5 ingestion boundary."""
from __future__ import annotations

import csv
from dataclasses import replace
from datetime import date
import gzip
import hashlib
import json
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile

import pytest

from base_edit_safety.ingestion import d1, d3, d4, d5
from base_edit_safety.ingestion.common import (
    SourceFile,
    bytes_opener,
    download_source,
    verify_source_file,
)
from base_edit_safety.ingestion.xlsx import iter_rows, sheet_names


def _write_synthetic_d5_config(
    path: Path,
    *,
    report_date: str | None = "2030-01-02",
) -> d5.D5Config:
    payload = {
        "guide_sequence": "ACGTTGCAACGATCGTACGA",
        "loci": {
            "HBG1": {
                "outer_forward_primer": "TGCAGATCTACCG",
                "outer_reverse_primer": "AACCGTTAAGGCTAC",
                "inner_forward_primer": "GATCCGTAACGTT",
                "inner_reverse_primer": "CGTACCTAGGATC",
            },
            "HBG2": {
                "outer_forward_primer": "AGTCGATGCTCCA",
                "outer_reverse_primer": "AACCGTTAAGGCTAC",
                "inner_forward_primer": "GATCCGTAACGTT",
                "inner_reverse_primer": "CGTACCTAGGATC",
            },
        },
        "run_qc_notes": {"RUN-SENSITIVE": "synthetic technical note"},
        "sensitivity_runs": ["RUN-SENSITIVE"],
    }
    if report_date is not None:
        payload["report_date"] = report_date
    path.write_text(json.dumps(payload), encoding="utf-8")
    return d5.load_config(path)


def _column_name(index: int) -> str:
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(ord("A") + remainder) + result
    return result


def _write_xlsx(path: Path, sheets: dict[str, list[list[object]]]) -> None:
    """Write the minimal inline-string OOXML used by parser fixtures."""
    workbook_sheets = []
    relationships = []
    content_overrides = []
    with ZipFile(path, "w", ZIP_DEFLATED) as book:
        for sheet_index, (name, rows) in enumerate(sheets.items(), start=1):
            workbook_sheets.append(
                f'<sheet name="{escape(name)}" sheetId="{sheet_index}" r:id="rId{sheet_index}"/>'
            )
            relationships.append(
                f'<Relationship Id="rId{sheet_index}" '
                'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                f'Target="worksheets/sheet{sheet_index}.xml"/>'
            )
            content_overrides.append(
                f'<Override PartName="/xl/worksheets/sheet{sheet_index}.xml" '
                'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            )
            xml_rows = []
            for row_index, row in enumerate(rows, start=1):
                cells = []
                for column_index, value in enumerate(row):
                    reference = f"{_column_name(column_index)}{row_index}"
                    cells.append(
                        f'<c r="{reference}" t="inlineStr"><is><t>{escape(str(value))}</t></is></c>'
                    )
                xml_rows.append(f'<row r="{row_index}">{"".join(cells)}</row>')
            book.writestr(
                f"xl/worksheets/sheet{sheet_index}.xml",
                '<?xml version="1.0" encoding="UTF-8"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>',
            )
        book.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f'<sheets>{"".join(workbook_sheets)}</sheets></workbook>',
        )
        book.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{"".join(relationships)}</Relationships>',
        )
        book.writestr(
            "[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            f'{"".join(content_overrides)}</Types>',
        )


def _write_gzip_table(path: Path, header: list[str], rows: list[list[object]], delimiter: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def test_common_acquisition_uses_injected_bytes_and_checksum(tmp_path: Path):
    payload = b"fixture-only\n"
    source = SourceFile(
        "DX", "fixture", "raw.bin", "https://not-requested.invalid/raw.bin",
        hashlib.sha256(payload).hexdigest(), "2026-08-14",
    )
    downloaded = download_source(source, tmp_path, opener=bytes_opener(payload))
    assert downloaded.read_bytes() == payload
    assert verify_source_file(source, downloaded) == source.sha256
    assert download_source(source, tmp_path, opener=lambda _url: pytest.fail("network called")) == downloaded
    downloaded.write_bytes(b"tampered\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_source_file(source, downloaded)


def test_xlsx_inline_string_fixture(tmp_path: Path):
    workbook = tmp_path / "fixture.xlsx"
    _write_xlsx(workbook, {"First": [["a", "b"], ["1", "two"]]})
    assert sheet_names(workbook) == ["First"]
    assert list(iter_rows(workbook, "First")) == [["a", "b"], ["1", "two"]]


def test_d1_parsers_preserve_unresolved_crosswalks(tmp_path: Path):
    assert {row["accession"] for row in d1.REGISTRY_RECORDS} == {
        "GSE298535",
        "PRJNA1269627",
    }
    assert {row["accessed_on"] for row in d1.REGISTRY_RECORDS} == {"2026-08-14"}
    peak_header = [
        "PeakID (cmd=annotatePeaks.pl fixture hg38)", "Chr", "Start", "End", "Strand",
        "Peak Score", "Annotation",
    ]
    peak_rows = [["chrX:9-32", "chrX", "10", "32", "+", "12", "Intergenic"]]
    replicate_paths = []
    for replicate in (1, 2):
        path = tmp_path / f"rep{replicate}.tsv.gz"
        _write_gzip_table(path, peak_header, peak_rows, "\t")
        replicate_paths.append(path)
    cd40l = d1.parse_cd40l(replicate_paths)
    assert len(cd40l) == 2
    assert cd40l[0]["start_0based"] == 9
    assert cd40l[0]["start_1based"] == 10
    assert cd40l[0]["editor"] == "TBD"

    cbe_path = tmp_path / "cbe.csv.gz"
    cbe_header = [
        "sample", "chr", "start", "end", "Target_Sequence", "OffTarget_Sequence",
        "coordinate", "CHANGEseqBE_rep1", "CHANGEseqBE_rep2", "Annotation",
    ]
    _write_gzip_table(
        cbe_path,
        cbe_header,
        [["B2M", "chr1", "100", "123", "A" * 20 + "NGG", "A" * 23,
          "chr1:100-123", "7", "", "intron"]],
        ",",
    )
    cbe = d1.parse_cbe(cbe_path)
    assert cbe[0]["editor"] == "eA3A-BE3"
    assert cbe[0]["detected_by_changeseq_be"] == "true"

    table_1 = tmp_path / "table1.xlsx"
    _write_xlsx(
        table_1,
        {"MKSR": [["Genome_Editor", "Target_site", "Chromosome", "Start", "End",
                    "Genomic Coordinate", "MKSR_rep1", "MKSR_rep2", "Site_Sequence"],
                  ["ABE8e-NRCH", "MKSR", "chr2", "20", "43", "chr2:20-43", "1", "2", "C" * 23]]},
    )
    assert d1.parse_abe_table_1(table_1)[0]["editor"] == "ABE8e-NRCH"

    table_3 = tmp_path / "table3.xlsx"
    _write_xlsx(
        table_3,
        {"ABE": [["Genome_Editor", "Target", "Chromosome", "Start", "End",
                  "Genomic_Coordinate", "Strand", "CHANGEseq-BE_rep1_reads",
                  "CHANGEseq-BE_rep2_reads", "Target_Sequence", "Off-target_Sequence", "Annotation"],
                 ["ABE", "B2M", "chr3", "30", "53", "chr3:30-53", "-", "3", "4",
                  "G" * 20 + "NGG", "G" * 23, "exon"]]},
    )
    assert d1.parse_abe_table_3(table_3)[0]["editor"].startswith("TBD")


def test_d3_coordinate_normalization_and_product_flag(tmp_path: Path):
    workbook = tmp_path / "d3.xlsx"
    _write_xlsx(
        workbook,
        {"CHANGE-seq_Supp_Table_3": [
            ["chrom", "chromStart", "chromEnd", "name", "CHANGEseq_reads", "strand",
             "offtarget_sequence", "", "chromStart:chromEnd", "distance", "target"],
            ["chr4", "100", "123", "AAVS1_site_1", "9", "+", "A" * 23, "",
             "chr4:101-123:+", "2", "A" * 20 + "NGG"],
        ]},
    )
    sites = d3.parse_sites(workbook)
    assert sites[0]["start_0based"] == 100
    assert sites[0]["start_1based"] == 101
    assert sites[0]["product_identity_crosswalk"] == "unverified"


def test_d4_behive_allowlist_omits_prediction_and_comparison_fields(tmp_path: Path):
    path = tmp_path / "behive.csv"
    header = [
        "Name (unique)", "gRNA (20nt)", "Sequence context (56nt)",
        "Wildtype sequence context (56nt)", "Editing position",
        "Pathogenic nucleotide (gRNA orientation)", "Corrected nucleotide (gRNA orientation)",
        "hg19 chromosome", "hg19 coordinate", "hg19 strand", "ALLELEID", "gene",
        "pearsonr_obs_r1_obs_r2_HEK293T_BE4", "Pred total probability_HEK293T_BE4",
        "Total obs read count_HEK293T_BE4",
        "Obs gt correction precision among edited reads_HEK293T_BE4",
        "Obs aa correction precision among edited reads_HEK293T_BE4",
        "Obs aa correction precision among edited aas_HEK293T_BE4",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerow(
            ["target-1", "A" * 20, "A" * 56, "C" * 56, "6", "C", "T", "1", "99", "+",
             "7", "GENE", "0.99", "0.88", "1000", "0.5", "0.4", "0.3"]
        )
    rows = list(d4.iter_behive_measurements(path, source_table="Table S2"))
    assert len(rows) == 1
    assert set(rows[0]) == set(d4.BEHIVE_FIELDS)
    assert not any(key.lower().startswith(("pred", "pearson")) for key in rows[0])
    assert rows[0]["editor"] == "BE4"


def test_d4_crispron_outcome_accounting_fixture(tmp_path: Path):
    workbook = tmp_path / "crispron.xlsx"
    header = ["label", "gRNA", "outcome", "read", "total", "WTclass", "edited_read",
              "efficiency", "outcome_frequency"]
    rows = [
        header,
        ["1", "A" * 30, "A" * 30, "80", "100", "yes", "20", "20", "80"],
        ["1", "A" * 30, "C" + "A" * 29, "20", "100", "no", "20", "20", "20"],
    ]
    _write_xlsx(workbook, {"ABE7.10": rows, "BE4": rows})
    stats: dict[str, object] = {}
    outcomes = list(d4.iter_crispron_outcomes(workbook, stats=stats))
    assert len(outcomes) == 4
    assert stats["ABE7.10"]["outcome_read_accounting_errors"] == []
    assert stats["BE4"]["unique_targets"] == 1
    assert all(row["build"].startswith("not applicable") for row in outcomes)


def test_d4_bedict_split_manifest_reads_sequence_identity_only(tmp_path: Path):
    workbook = tmp_path / "bedict.xlsx"
    rows = [["ID", "Purpose", "Count", "Sequence"], ["x", "Train", "10", "A" * 20],
            ["y", "Validation/Test", "9", "C" * 20]]
    _write_xlsx(
        workbook,
        {sheet: rows for sheet in (
            "ABEmax_perbase", "CBE4max_perbase", "ABE8e_perbase", "Target-AID_perbase"
        )},
    )
    split = d4.bedict_development_sequences(workbook)
    assert split == {"Train": {"A" * 20}, "Validation/Test": {"C" * 20}}


def test_d5_reference_and_crispresso_normalizers(tmp_path: Path):
    config = _write_synthetic_d5_config(tmp_path / "d5_config.json")
    spec = config.loci["HBG1"]
    guide_sequence = config.guide_sequence
    amplicon = (
        spec.inner_reverse_primer
        + "AAAA"
        + d5.reverse_complement(guide_sequence)
        + "CCCC"
        + d5.reverse_complement(spec.inner_forward_primer)
    )
    sequence = (
        spec.outer_reverse_primer + "TTTT" + amplicon + "GGGG"
        + d5.reverse_complement(spec.outer_forward_primer)
    )
    fixture_filename = f"hg38_chr11_1001_{1000 + len(sequence)}.json"
    cache = tmp_path / fixture_filename
    cache.write_text(
        json.dumps({"genome_build": "hg38", "chrom": "chr11", "start_1based": 1001,
                    "end_1based": 1000 + len(sequence), "plus_sequence": sequence})
    )
    fixture_spec = replace(
        spec,
        cache_filename=fixture_filename,
        sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
    )
    reference = d5.derive_amplicon(
        fixture_spec, tmp_path, guide_sequence=guide_sequence
    )
    mapped_end_reference = d5.derive_mapped_end_reference(
        fixture_spec, tmp_path, guide_sequence=guide_sequence
    )
    with pytest.raises(ValueError, match="reference-sequence checksum mismatch"):
        d5.derive_amplicon(
            replace(fixture_spec, sequence_sha256="0" * 64),
            tmp_path,
            guide_sequence=guide_sequence,
        )
    assert reference.sequence_plus == amplicon
    assert mapped_end_reference.sequence_plus == sequence
    assert mapped_end_reference.start_0based == 1000
    assert mapped_end_reference.end_0based_exclusive == 1000 + len(sequence)
    assert mapped_end_reference.contig == "HBG1_mapped_end_gate"
    assert (
        mapped_end_reference.sequence_plus[
            mapped_end_reference.guide_start_0based:mapped_end_reference.guide_end_0based_exclusive
        ]
        == d5.reverse_complement(guide_sequence)
    )
    assert reference.guide_strand == "-"
    assert (
        reference.sequence_plus[
            reference.guide_start_0based:reference.guide_end_0based_exclusive
        ]
        == d5.reverse_complement(guide_sequence)
    )

    window_length = d5.ANALYSIS_WINDOW_LENGTH
    long_sequence = "A" * window_length + "C" + "T" * window_length
    small_reference = d5.AmpliconReference(
        locus="HBG1", contig="HBG1_amplicon", build=d5.REFERENCE_BUILD, chrom="chr11",
        start_0based=100, end_0based_exclusive=100 + len(long_sequence),
        sequence_plus=long_sequence,
        guide_start_0based=1, guide_end_0based_exclusive=3, guide_strand="-",
        source_cache="fixture", source_sha256="0" * 64,
        source_url="https://not-requested.invalid/reference",
        source_sequence_sha256="0" * 64,
    )
    window_dirs = {name: tmp_path / name for name in ("target", "distal")}
    for name, directory in window_dirs.items():
        directory.mkdir()
        (directory / "CRISPResso_mapping_statistics.txt").write_text(
            "READS IN INPUTS\tREADS AFTER PREPROCESSING\tREADS ALIGNED\n"
            "12\t11\t11\n"
        )
        base = "A" if name == "target" else "T"
        positions = "\t".join(str(index) for index in range(window_length))
        counts = {
            nucleotide: "\t".join(["10" if nucleotide == base else "0"] * window_length)
            for nucleotide in "ACGTN"
        }
        percentages = {
            nucleotide: "\t".join(["1" if nucleotide == base else "0"] * window_length)
            for nucleotide in "ACGTN"
        }
        counts["-"] = "\t".join(["1"] * window_length)
        percentages["-"] = "\t".join(["0.090909"] * window_length)
        (directory / "window.Nucleotide_frequency_table.txt").write_text(
            "Batch\t" + positions + "\n"
            + "".join(f"{nucleotide}\t{counts[nucleotide]}\n" for nucleotide in "ACGTN-")
        )
        (directory / "window.Nucleotide_percentage_table.txt").write_text(
            "Batch\t" + positions + "\n"
            + "".join(
                f"{nucleotide}\t{percentages[nucleotide]}\n" for nucleotide in "ACGTN-"
            )
        )
    for window_name, directory in window_dirs.items():
        allele_zip = directory / "Alleles_frequency_table.zip"
        with ZipFile(allele_zip, "w", ZIP_DEFLATED) as archive:
            archive.writestr(
                "alleles.txt",
                "Aligned_Sequence\tReference_Sequence\t#Reads\t%Reads\tn_deleted\tn_inserted\tn_mutated\n"
                f"AGCT\tAGCT\t10\t100\t0\t0\t0\n",
            )
    run = d5.RunSpec(
        "RUN-001", "STUDY-001", "SAMPLE-001", "EXPERIMENT-001", "editor", "HBG1", 1,
        "cell", "assay", "r1.fastq.gz", "r2.fastq.gz", "a", "b",
    )
    position_rows, minimum_aligned_depth, minimum_base_call_depth = d5.normalize_positions(
        run, small_reference, window_dirs, tmp_path / "positions.tsv"
    )
    allele_rows, allele_stats = d5.normalize_alleles(run, window_dirs, tmp_path / "alleles.tsv")
    assert (position_rows, minimum_aligned_depth, minimum_base_call_depth) == (
        2 * window_length + 1, 11, 10
    )
    assert allele_rows == 2
    assert allele_stats == {
        "target": {"allele_rows": 1, "accounted_reads": 10},
        "distal": {"allele_rows": 1, "accounted_reads": 10},
    }
    allele_output = list(csv.DictReader((tmp_path / "alleles.tsv").open(), delimiter="\t"))
    assert [row["analysis_window"] for row in allele_output] == ["target", "distal"]
    assert len({row["allele_key"] for row in allele_output}) == 2
    rows = list(csv.DictReader((tmp_path / "positions.tsv").open(), delimiter="\t"))
    assert rows[0]["deletion_count"] == "1"
    assert rows[0]["depth"] == "11"
    assert rows[0]["base_call_depth"] == "10"
    assert rows[window_length]["coverage_status"] == "unsequenced_between_mates"
    assert d5.crispresso_mapping_counts(window_dirs["target"]) == {
        "reads_input": 12,
        "reads_after_preprocessing": 11,
        "reads_aligned": 11,
        "alignment_fraction": 11 / 12,
    }
    assert d5._run_qc_failures(
        mapped_fraction=0.69,
        minimum_base_call_depth=999,
        position_count=len(long_sequence),
        expected_position_count=len(long_sequence),
        target_allele_reads=11,
        distal_allele_reads=10,
        target_aligned_reads=11,
        distal_aligned_reads=11,
    ) == [
        "mapped_end_fraction_below_minimum",
        "protospacer_depth_below_minimum",
        "distal_allele_read_accounting_mismatch",
    ]


def test_ingestion_modules_do_not_import_prediction_or_scoring_code():
    root = Path(__file__).parents[1] / "src" / "base_edit_safety" / "ingestion"
    text = "\n".join(path.read_text() for path in root.glob("*.py"))
    assert "base_edit_safety.adapters" not in text
    assert "base_edit_safety.scoring" not in text


def test_d5_qc_summary_is_technical_only(tmp_path: Path):
    config = _write_synthetic_d5_config(tmp_path / "d5_config.json")
    rendered = d5.render_qc_summary(
        [
            {
                "run_accession": "SRR1", "editor": "editor", "locus": "HBG1",
                "replicate": 1, "read_pairs": 2_000,
                "primary_mapped_end_fraction": 0.8,
                "target_crispresso_alignment_fraction": 0.75,
                "distal_crispresso_alignment_fraction": 0.7,
                "minimum_protospacer_base_call_depth": 1_500,
                "allele_rows": 20, "qc_result": "pass", "qc_failures": "",
                "mapped_end_gate_rerun_scope": "outcome_blind",
            }
        ],
        config,
    )
    assert "SRR1" in rendered
    assert "loaded no editing-frequency" in rendered
    assert "staged F1/R1 flank-inclusive geometry" in rendered
    assert "80.00%" in rendered
    assert "Generated: 2030-01-02" in rendered


def test_d5_config_uses_execution_date_without_override(tmp_path: Path):
    config = _write_synthetic_d5_config(
        tmp_path / "d5_config.json", report_date=None
    )
    assert config.effective_report_date == date.today().isoformat()


def test_d5_config_missing_file_has_clear_error(tmp_path: Path):
    with pytest.raises(d5.D5ConfigError, match="runtime configuration was not found"):
        d5.load_config(tmp_path / "missing.json")


def test_d5_outcome_blind_qc_recheck_changes_only_gate_fields(tmp_path: Path):
    config = _write_synthetic_d5_config(tmp_path / "d5_config.json")
    run = d5.RunSpec(
        "RUN-SENSITIVE", "STUDY-001", "SAMPLE-001", "EXPERIMENT-001", "editor", "HBG1", 1,
        "cell", "assay", "r1.fastq.gz", "r2.fastq.gz", "a", "b",
    )
    reference = d5.AmpliconReference(
        locus="HBG1", contig="HBG1_mapped_end_gate", build=d5.REFERENCE_BUILD,
        chrom="chr11", start_0based=100, end_0based_exclusive=740,
        sequence_plus="A" * 640, guide_start_0based=133,
        guide_end_0based_exclusive=153, guide_strand="-", source_cache="fixture",
        source_url="https://not-requested.invalid/reference", source_sha256="0" * 64,
        source_sequence_sha256="1" * 64,
    )
    run_dir = tmp_path / run.run_accession
    run_dir.mkdir()
    prior = {
        "run_accession": run.run_accession,
        "read_pairs": 500,
        "read_ends": 1_000,
        "minimum_protospacer_base_call_depth": 999,
        "position_rows": 307,
        "amplicon_length": 307,
        "target_allele_accounted_reads": 10,
        "distal_allele_accounted_reads": 11,
        "target_crispresso_reads_aligned": 10,
        "distal_crispresso_reads_aligned": 11,
        "untouched_technical_field": "preserved",
    }
    (run_dir / "run_qc.json").write_text(json.dumps(prior))
    updated = d5._recheck_existing_qc(
        run,
        reference,
        {"read_pairs": 500, "read_ends": 1_000, "md5_status": "pass"},
        {
            "mapped_read_ends": 800,
            "primary_mapped_read_ends": 800,
            "properly_paired_read_ends": 700,
        },
        run_dir,
        config,
    )
    assert updated["primary_mapped_end_fraction"] == 0.8
    assert updated["qc_failures"] == "protospacer_depth_below_minimum"
    assert updated["untouched_technical_field"] == "preserved"
    assert updated["mapped_end_reference_geometry"] == d5.MAPPED_END_REFERENCE_GEOMETRY
    assert updated["mapped_end_reference_length"] == 640
    assert updated["mapped_end_gate_rerun_scope"] == "outcome_blind"
    assert updated["qc_note"] == "synthetic technical note"
    assert updated["held_out_sensitivity_flag"] is True


def test_d5_mapped_end_rerun_refuses_partial_run_selection(tmp_path: Path):
    with pytest.raises(ValueError, match="must cover all 20 runs uniformly"):
        d5.rerun_mapped_end_qc(
            tmp_path,
            manifest_path=tmp_path / "not-read.tsv",
            fastq_dir=tmp_path,
            cache_dir=tmp_path,
            output_dir=tmp_path / "out",
            only_runs={"RUN-001"},
        )
