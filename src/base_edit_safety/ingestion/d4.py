"""Acquire and normalize measured D4 BE-Hive and CRISPRon-BE outcome tables.

The BE-Hive publisher tables contain measured values, predictions, and already-computed comparison
statistics side by side.  This parser uses an explicit measured-column allowlist and rejects any
attempt to emit the prediction or comparison columns.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
import json
from pathlib import Path
from typing import Iterable, Iterator

from .common import SourceFile, download_source, read_dict_rows, verify_source_file, write_tsv
from .xlsx import iter_rows


DATASET_ID = "D4"
ACCESSED_ON = "2026-08-14"
SOURCES = (
    SourceFile(
        DATASET_ID,
        "Arbab et al. 2020 Table S2; PMC7384975; PRJNA591007",
        "BEHive_Table_S2.csv",
        "https://pmc.ncbi.nlm.nih.gov/articles/instance/7384975/bin/"
        "NIHMS1599657-supplement-9.csv",
        "034e8e059b56c7b4e1670cd544f84210725f7692782f90b1893116c3f89b07cc",
        ACCESSED_ON,
    ),
    SourceFile(
        DATASET_ID,
        "Arbab et al. 2020 Table S3; PMC7384975; PRJNA591007",
        "BEHive_Table_S3.csv",
        "https://pmc.ncbi.nlm.nih.gov/articles/instance/7384975/bin/"
        "NIHMS1599657-supplement-10.csv",
        "88d456cf17745509efa39fa243cc58d30e9f3233e4f6232d113c27c43358e0aa",
        ACCESSED_ON,
    ),
    SourceFile(
        DATASET_ID,
        "Arbab et al. 2020 Table S8; PMC7384975; PRJNA591007",
        "BEHive_Table_S8.xlsx",
        "https://pmc.ncbi.nlm.nih.gov/articles/instance/7384975/bin/"
        "NIHMS1599657-supplement-15.xlsx",
        "99662791110b757e9f39090e7c8c9aed90e26b2dc5c76e67b564e7b222bec8e4",
        ACCESSED_ON,
    ),
    SourceFile(
        DATASET_ID,
        "Sun et al. 2025 Supplementary Data 1; PMC12595020; CNP0001031; PRJNA1321836",
        "CRISPRon_Supplementary_Data_1.xlsx",
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65200-5/"
        "MediaObjects/41467_2025_65200_MOESM3_ESM.xlsx",
        "a9413d9a01d6dc73e7fa783234ddc88f9dc4a08668de0597678ed4e7c2513cfa",
        ACCESSED_ON,
    ),
    SourceFile(
        DATASET_ID,
        "BE-DICT training-source audit; Marquart et al. 2021 Supplementary Dataset 1; PMC8387386",
        "BEDICT_Supplementary_Dataset_1.xlsx",
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-25375-z/"
        "MediaObjects/41467_2021_25375_MOESM2_ESM.xlsx",
        "bc8c1033d239305e256386cbe66f2cd19097ab84d9fd0ba35c183bf668e1ad92",
        ACCESSED_ON,
    ),
)

BEHIVE_FIELDS = [
    "dataset_id", "study", "source_table", "source_file", "source_row", "target_id",
    "guide_20nt", "sequence_context_56nt", "wildtype_sequence_context_56nt",
    "editing_position", "pathogenic_nucleotide_guide_orientation",
    "corrected_nucleotide_guide_orientation", "build", "chrom", "source_coordinate",
    "strand", "allele_id", "gene", "cell_type", "editor", "total_observed_read_count",
    "observed_genotype_correction_precision_among_edited_reads",
    "observed_amino_acid_correction_precision_among_edited_reads",
    "observed_amino_acid_correction_precision_among_edited_amino_acids",
]

CRISPRON_FIELDS = [
    "dataset_id", "study", "source_file", "source_sheet", "source_row", "editor",
    "cell_type", "library", "target_id", "target_sequence_30nt", "outcome_sequence_30nt",
    "outcome_read_count", "target_total_read_count", "wildtype_class", "edited_read_count",
    "measured_editing_efficiency_percent", "measured_outcome_frequency_percent", "build",
]


def acquire_sources(destination: str | Path) -> list[Path]:
    return [download_source(source, destination) for source in SOURCES]


def _profile_columns(header: Iterable[str]) -> list[str]:
    prefix = "Total obs read count_"
    profiles = [column[len(prefix):] for column in header if column.startswith(prefix)]
    if not profiles:
        raise ValueError("BE-Hive table has no measured observation profiles")
    return profiles


def iter_behive_measurements(
    path: str | Path,
    *,
    source_table: str,
    stats: dict[str, object] | None = None,
) -> Iterator[dict[str, object]]:
    source = Path(path)
    rows = read_dict_rows(source)
    try:
        first = next(rows)
    except StopIteration as exc:
        raise ValueError(f"empty BE-Hive table: {source}") from exc
    header = list(first)
    profiles = _profile_columns(header)
    required = {
        "Name (unique)", "gRNA (20nt)", "Sequence context (56nt)",
        "Wildtype sequence context (56nt)", "hg19 chromosome", "hg19 coordinate",
        "hg19 strand",
    }
    missing = required - set(header)
    if missing:
        raise ValueError(f"missing BE-Hive fields in {source}: {sorted(missing)}")

    profile_counts: Counter[str] = Counter()
    target_ids: set[str] = set()
    emitted = 0

    def source_rows():
        yield first
        yield from rows

    for source_row, row in enumerate(source_rows(), start=2):
        target_ids.add(row["Name (unique)"])
        for profile in profiles:
            total_column = f"Total obs read count_{profile}"
            gt_column = f"Obs gt correction precision among edited reads_{profile}"
            aa_reads_column = f"Obs aa correction precision among edited reads_{profile}"
            aa_edited_column = f"Obs aa correction precision among edited aas_{profile}"
            selected_columns = {total_column, gt_column, aa_reads_column, aa_edited_column}
            if any(column.startswith(("Pred ", "pearsonr_")) for column in selected_columns):
                raise AssertionError("non-measurement BE-Hive field entered the allowlist")
            values = [row.get(column, "") for column in selected_columns]
            if not any(value.strip() for value in values):
                continue
            cell_type, editor = profile.split("_", 1)
            emitted += 1
            profile_counts[profile] += 1
            yield {
                "dataset_id": DATASET_ID,
                "study": "BE-Hive / Arbab 2020",
                "source_table": source_table,
                "source_file": source.name,
                "source_row": source_row,
                "target_id": row["Name (unique)"],
                "guide_20nt": row["gRNA (20nt)"],
                "sequence_context_56nt": row["Sequence context (56nt)"],
                "wildtype_sequence_context_56nt": row["Wildtype sequence context (56nt)"],
                "editing_position": row.get("Editing position", ""),
                "pathogenic_nucleotide_guide_orientation": row.get(
                    "Pathogenic nucleotide (gRNA orientation)", ""
                ),
                "corrected_nucleotide_guide_orientation": row.get(
                    "Corrected nucleotide (gRNA orientation)", ""
                ),
                "build": "GRCh37/hg19",
                "chrom": row["hg19 chromosome"],
                "source_coordinate": row["hg19 coordinate"],
                "strand": row["hg19 strand"],
                "allele_id": row.get("ALLELEID", ""),
                "gene": row.get("gene", ""),
                "cell_type": cell_type,
                "editor": editor,
                "total_observed_read_count": row[total_column],
                "observed_genotype_correction_precision_among_edited_reads": row.get(
                    gt_column, ""
                ),
                "observed_amino_acid_correction_precision_among_edited_reads": row.get(
                    aa_reads_column, ""
                ),
                "observed_amino_acid_correction_precision_among_edited_amino_acids": row.get(
                    aa_edited_column, ""
                ),
            }
    if stats is not None:
        stats.update(
            {
                "source_rows": len(target_ids),
                "unique_target_ids": len(target_ids),
                "measurement_rows": emitted,
                "profile_rows": dict(sorted(profile_counts.items())),
                "profiles": sorted(profiles),
            }
        )


def _as_int(value: str, label: str) -> int:
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"non-integral {label}: {value}")
    return int(number)


def iter_crispron_outcomes(
    path: str | Path,
    *,
    stats: dict[str, object] | None = None,
) -> Iterator[dict[str, object]]:
    source = Path(path)
    sheet_stats: dict[str, object] = {}
    for sheet in ("ABE7.10", "BE4"):
        iterator = iter_rows(source, sheet)
        header = next(iterator)
        expected = {
            "label", "gRNA", "outcome", "read", "total", "WTclass", "edited_read",
            "efficiency", "outcome_frequency",
        }
        if set(header) != expected:
            raise ValueError(f"unexpected CRISPRon schema for {sheet}: {header}")
        target_totals: dict[str, int] = {}
        target_edited: dict[str, int] = {}
        outcome_sums: Counter[str] = Counter()
        outcome_rows = 0
        for source_row, values in enumerate(iterator, start=2):
            values += [""] * (len(header) - len(values))
            row = dict(zip(header, values))
            target = row["label"]
            total = _as_int(row["total"], "total")
            edited = _as_int(row["edited_read"], "edited_read")
            read_count = _as_int(row["read"], "read")
            if target in target_totals and target_totals[target] != total:
                raise ValueError(f"inconsistent total for {sheet} target {target}")
            if target in target_edited and target_edited[target] != edited:
                raise ValueError(f"inconsistent edited_read for {sheet} target {target}")
            target_totals[target] = total
            target_edited[target] = edited
            outcome_sums[target] += read_count
            outcome_rows += 1
            yield {
                "dataset_id": DATASET_ID,
                "study": "CRISPRon-BE / Sun 2025",
                "source_file": source.name,
                "source_sheet": sheet,
                "source_row": source_row,
                "editor": sheet,
                "cell_type": "HEK293T",
                "library": "SURRO-seq 12K gRNA-target pair library",
                "target_id": target,
                "target_sequence_30nt": row["gRNA"],
                "outcome_sequence_30nt": row["outcome"],
                "outcome_read_count": read_count,
                "target_total_read_count": total,
                "wildtype_class": row["WTclass"],
                "edited_read_count": edited,
                "measured_editing_efficiency_percent": row["efficiency"],
                "measured_outcome_frequency_percent": row["outcome_frequency"],
                "build": "not applicable: surrogate target library",
            }
        accounting_errors = sorted(
            target for target, total in target_totals.items() if outcome_sums[target] != total
        )
        sheet_stats[sheet] = {
            "outcome_rows": outcome_rows,
            "unique_targets": len(target_totals),
            "minimum_total_reads": min(target_totals.values()),
            "targets_below_100_total_reads": sum(total < 100 for total in target_totals.values()),
            "outcome_read_accounting_errors": accounting_errors,
        }
    if stats is not None:
        stats.update(sheet_stats)


def _read_oligo_manifest(path: Path) -> list[dict[str, object]]:
    iterator = iter_rows(path, "Sheet1")
    header = next(iterator)
    records = []
    for source_row, values in enumerate(iterator, start=2):
        values += [""] * (len(header) - len(values))
        row = dict(zip(header, values))
        if not any(row.values()):
            continue
        records.append({"source_row": source_row, **row})
    return records


def bedict_development_sequences(path: str | Path) -> dict[str, set[str]]:
    """Read only 20-nt sequence identities and publisher-assigned split purposes."""
    by_purpose: dict[str, set[str]] = defaultdict(set)
    for sheet in ("ABEmax_perbase", "CBE4max_perbase", "ABE8e_perbase", "Target-AID_perbase"):
        iterator = iter_rows(path, sheet)
        header = next(iterator)
        for values in iterator:
            values += [""] * (len(header) - len(values))
            row = dict(zip(header, values))
            sequence = row["Sequence"].upper()
            if len(sequence) != 20 or set(sequence) - set("ACGT"):
                raise ValueError(f"invalid BE-DICT development sequence in {sheet}")
            by_purpose[row["Purpose"]].add(sequence)
    return dict(by_purpose)


def _behive_guides(path: Path) -> set[str]:
    return {row["gRNA (20nt)"].upper() for row in read_dict_rows(path)}


def _crispron_protospacers(path: Path, sheet: str) -> set[str]:
    iterator = iter_rows(path, sheet)
    header = next(iterator)
    guides: set[str] = set()
    for values in iterator:
        values += [""] * (len(header) - len(values))
        row = dict(zip(header, values))
        target_30nt = row["gRNA"].upper()
        if len(target_30nt) != 30:
            raise ValueError(f"invalid CRISPRon 30-nt target in {sheet}")
        # Sun et al. define this field as 4-nt 5' flank + 20-nt protospacer + 3-nt PAM +
        # 3-nt 3' flank.  Only the explicit 20-nt sequence identity is used for this audit.
        guides.add(target_30nt[4:24])
    return guides


def ingest(source_dir: str | Path, output_dir: str | Path) -> dict[str, object]:
    source = Path(source_dir)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    ledger = [
        asdict(item) | {"local_sha256": verify_source_file(item, source / item.filename)}
        for item in SOURCES
    ]

    behive_stats: dict[str, dict[str, object]] = {}
    behive_generators = []
    for filename, table in (("BEHive_Table_S2.csv", "Table S2"), ("BEHive_Table_S3.csv", "Table S3")):
        table_stats: dict[str, object] = {}
        behive_stats[table] = table_stats
        behive_generators.append(
            iter_behive_measurements(source / filename, source_table=table, stats=table_stats)
        )

    def all_behive_rows():
        for generator in behive_generators:
            yield from generator

    behive_count = write_tsv(
        destination / "behive_measured_outcomes.tsv", BEHIVE_FIELDS, all_behive_rows()
    )

    crispron_stats: dict[str, object] = {}
    crispron_count = write_tsv(
        destination / "crispron_measured_outcomes.tsv",
        CRISPRON_FIELDS,
        iter_crispron_outcomes(
            source / "CRISPRon_Supplementary_Data_1.xlsx", stats=crispron_stats
        ),
    )
    oligos = _read_oligo_manifest(source / "BEHive_Table_S8.xlsx")
    write_tsv(
        destination / "behive_oligo_manifest.tsv",
        ["source_row", "Name", "Purpose", "Sequence"],
        oligos,
    )

    manifest_rows: list[dict[str, object]] = []
    behive_libraries = {
        "Table S2": "CBE-correctable disease-allele target library",
        "Table S3": "ABE-correctable disease-allele target library",
    }
    for table, table_stats in behive_stats.items():
        for profile, row_count in table_stats["profile_rows"].items():
            cell_type, editor = profile.split("_", 1)
            manifest_rows.append(
                {
                    "study": "BE-Hive / Arbab 2020",
                    "source_table_or_sheet": table,
                    "library": behive_libraries[table],
                    "cell_type": cell_type,
                    "editor": editor,
                    "normalized_profile_or_target_count": row_count,
                    "sample_mapping_status": (
                        "publisher table profile resolved; raw-run sample crosswalk not parsed"
                    ),
                }
            )
    for editor in ("ABE7.10", "BE4"):
        manifest_rows.append(
            {
                "study": "CRISPRon-BE / Sun 2025",
                "source_table_or_sheet": editor,
                "library": "SURRO-seq 12K gRNA-target pair library",
                "cell_type": "HEK293T",
                "editor": editor,
                "normalized_profile_or_target_count": crispron_stats[editor]["unique_targets"],
                "sample_mapping_status": (
                    "publisher processed-table context resolved; raw-run crosswalk not parsed"
                ),
            }
        )
    write_tsv(
        destination / "library_editor_sample_manifest.tsv",
        [
            "study", "source_table_or_sheet", "library", "cell_type", "editor",
            "normalized_profile_or_target_count", "sample_mapping_status",
        ],
        manifest_rows,
    )

    expected_targets = {"ABE7.10": 11_484, "BE4": 11_406}
    crispron_reconciliation = {
        editor: (
            "pass" if crispron_stats[editor]["unique_targets"] == expected else "fail"
        )
        for editor, expected in expected_targets.items()
    }
    summary: dict[str, object] = {
        "dataset_id": DATASET_ID,
        "behive_measurement_rows": behive_count,
        "behive_tables": behive_stats,
        "behive_oligo_rows": len(oligos),
        "library_editor_sample_manifest_rows": len(manifest_rows),
        "crispron_outcome_rows": crispron_count,
        "crispron_sheets": crispron_stats,
        "crispron_publication_reported_filtered_targets": expected_targets,
        "crispron_target_reconciliation": crispron_reconciliation,
        "build_handling": {
            "BE-Hive": "GRCh37/hg19 preserved; no liftover",
            "CRISPRon-BE": "not applicable to surrogate target sequences",
        },
        "column_scope": (
            "BE-Hive prediction and precomputed comparison columns excluded by allowlist; "
            "CRISPRon Supplementary Data 2 model-selection metrics not parsed"
        ),
        "scope_guard": "measured outcomes and source manifests only; no prediction records loaded",
    }
    (destination / "qc_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    development = bedict_development_sequences(source / "BEDICT_Supplementary_Dataset_1.xlsx")
    train_sequences = development.get("Train", set())
    validation_test_sequences = development.get("Validation/Test", set())
    audited_sets = {
        "BE-Hive Table S2": _behive_guides(source / "BEHive_Table_S2.csv"),
        "BE-Hive Table S3": _behive_guides(source / "BEHive_Table_S3.csv"),
        "CRISPRon ABE7.10": _crispron_protospacers(
            source / "CRISPRon_Supplementary_Data_1.xlsx", "ABE7.10"
        ),
        "CRISPRon BE4": _crispron_protospacers(
            source / "CRISPRon_Supplementary_Data_1.xlsx", "BE4"
        ),
    }
    sequence_audit = {
        name: {
            "unique_20nt_sequences": len(sequences),
            "exact_BE_DICT_train_matches": len(sequences & train_sequences),
            "exact_BE_DICT_validation_test_matches": len(sequences & validation_test_sequences),
        }
        for name, sequences in audited_sets.items()
    }
    overlap_review = {
        "repo_predictor": "BE-DICT",
        "pinned_revision": "d22e198d22dbcae2676c5b0493d45b3ec6d1156d",
        "development_manifest": "Marquart et al. 2021 Supplementary Dataset 1",
        "BE_DICT_unique_train_20nt_sequences": len(train_sequences),
        "BE_DICT_unique_validation_test_20nt_sequences": len(validation_test_sequences),
        "exact_sequence_audit": sequence_audit,
        "direct_training_overlap_result": (
            "no exact 20-nt matches to the 23,123 publisher-labeled BE-DICT Train sequences"
        ),
        "remaining_risks": [
            "BE-Hive Table S2 has 452 and Table S3 has 1,354 exact matches to BE-DICT's publisher-labeled Validation/Test sequences",
            "each CRISPRon editor table has one exact match to a BE-DICT Validation/Test sequence",
            "near-duplicate sequence similarity and undocumented pretraining are not resolved by an exact-identity audit",
        ],
        "scope": "sequence identities and publisher split labels only; no prediction output loaded",
    }
    (destination / "training_overlap_review.json").write_text(
        json.dumps(overlap_review, indent=2) + "\n"
    )
    (destination / "accession_ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", default="data/inputs/d4_outcomes")
    parser.add_argument("--output-dir", default="data/outputs/ingestion/d4")
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args(argv)
    if args.acquire:
        acquire_sources(args.source_dir)
    ingest(args.source_dir, args.output_dir)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
