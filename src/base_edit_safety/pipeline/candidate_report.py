"""Write a Markdown report from design-record v7 and candidate-ledger v3 files."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DISCLAIMER = "Organizes evidence; does not establish safety."


class CandidateReportError(ValueError):
    """Raised when report inputs use an unsupported format."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        handle = path.open("rb")
    except OSError as exc:
        raise CandidateReportError(f"cannot read report input {path}") from exc
    with handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path, label: str) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CandidateReportError(f"cannot read {label} {path}") from exc
    if not isinstance(payload, dict):
        raise CandidateReportError(f"{label} must be a JSON object")
    return payload


def _validate_inputs(design_records: dict, ledger: dict) -> None:
    if (
        design_records.get("record_type") != "hbg_design_record_set"
        or design_records.get("schema_version") != "7"
        or not isinstance(design_records.get("records"), list)
    ):
        raise CandidateReportError("design-record input must be an hbg_design_record_set v7")
    if (
        ledger.get("record_type") != "candidate_panel_ledger"
        or ledger.get("schema_version") != "3"
        or not isinstance(ledger.get("candidates"), list)
    ):
        raise CandidateReportError("candidate input must be a candidate_panel_ledger v3")
    design_guides = set()
    for record in design_records["records"]:
        if (
            not isinstance(record, dict)
            or record.get("record_type") != "hbg_design_record"
            or record.get("schema_version") != "7"
        ):
            raise CandidateReportError("design-record set contains an unsupported record")
        try:
            design_guides.add(record["candidate"]["guide_id"])
        except (KeyError, TypeError) as exc:
            raise CandidateReportError("design record is missing candidate.guide_id") from exc
    ledger_guides = {item.get("guide_id") for item in ledger.get("guides", [])}
    if None in ledger_guides or not ledger_guides:
        raise CandidateReportError("candidate ledger has no valid guide set")
    missing = sorted(ledger_guides - design_guides)
    if missing:
        raise CandidateReportError(
            "design-record input is missing ledger guides: " + ", ".join(missing)
        )


def _cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        text = "yes" if value else "no"
    elif isinstance(value, (list, dict)):
        text = json.dumps(value, sort_keys=True, separators=(",", ":"))
    else:
        text = str(value)
    return text.replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ")


def _joined(values) -> str:
    return "; ".join(sorted({str(value) for value in values if value not in {None, ""}}))


def _section(lines: list[str], title: str) -> None:
    lines.extend([f"## {title}", "", f"_{DISCLAIMER}_", ""])


def _design_rows(design_records: dict, ledger_guides: set[str]) -> list[list[object]]:
    rows = []
    for record in design_records["records"]:
        candidate = record["candidate"]
        if candidate["guide_id"] not in ledger_guides:
            continue
        intended = candidate.get("target", {}).get("intended_edits", [])
        offsets = _joined(item.get("promoter_offset") for item in intended)
        flags = record.get("hbg_paralog_flags", {})
        rows.append(
            [
                candidate["guide_id"],
                candidate.get("gene"),
                candidate.get("editor_name"),
                offsets,
                record.get("evidence_tier"),
                flags.get("flagged"),
            ]
        )
    return rows


def _candidate_rows(ledger: dict) -> list[list[object]]:
    rows = []
    for site in ledger["candidates"]:
        nominations = site.get("nominations", [])
        measurements = site.get("measurements", [])
        scores = _joined(
            f"{item.get('source_method')}:{item.get('score_name')}={item.get('score_value')}"
            for item in nominations
            if item.get("score_name") is not None and item.get("score_value") is not None
        )
        consequence = _joined(
            [
                *(
                    item.get("ch_driver_gene")
                    for item in nominations
                    if item.get("ch_driver_gene")
                ),
                *(item.get("site_role") for item in nominations),
            ]
        )
        rows.append(
            [
                site["site_id"],
                site["guide_id"],
                (f"{site['chrom']}:{site['start_0based']} (reported start anchor; span unresolved)"
                 if any(block.get("canonical_interval_basis") == "alt_only"
                        for block in site.get("allele_configurations", []))
                 else f"{site['chrom']}:{site['start_0based']}-{site['end_0based_exclusive']}"),
                site.get("source_bucket"),
                site.get("reconciliation", {}).get("status"),
                _joined(site.get("why_crispor_missed", [])),
                scores,
                (site.get("evidence_state") if "measurement_coverage" in site else "Legacy coverage not assessed"),
                (site.get("measurement_summary") if "measurement_coverage" in site else "Regenerate ledger"),
                site.get("editability", {}).get("site_summary"),
                consequence,
                site.get("panel_decision", {}).get("status"),
                site.get("panel_decision", {}).get("reason"),
                site.get("panel_decision", {}).get("note"),
                _joined(
                    item.get("assay")
                    for item in measurements
                    if item.get("assay")
                ),
            ]
        )
    return rows


def _table(lines: list[str], headers: list[str], rows: list[list[object]]) -> None:
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        lines.append("| " + " | ".join(_cell(value) for value in row) + " |")
    if not rows:
        lines.append("| " + " | ".join(["No records."] + [""] * (len(headers) - 1)) + " |")
    lines.append("")


def _measurement_tables(lines: list[str], ledger: dict) -> None:
    _section(lines, "Measurements and follow-up")
    lines.extend([
        "Coverage below applies only to the listed guide/editor, DNA sequence, assay and cell context. "
        "A reference-only negative does not test an alternate sequence. Failed or unreported quality "
        "checks, missing allele evidence and inconclusive results cannot supply a usable negative.",
        "",
        "Search scores and editable bases describe computational candidates. Editing measurements, "
        "consequence measurements and human panel choices remain separate. A design's evidence tier "
        "does not establish measurement coverage for its candidate sites.",
        "",
    ])
    measured = [site for site in ledger["candidates"] if site.get("measurements")]
    mapped = sum(bool(site.get("measurement_coverage")) and
                 all(s["mapping_reason"] is None for s in site["measurement_coverage"]) for site in measured)
    unsupported = sum(any(s.get("mapping_reason") == "experimental_candidate_model_unsupported"
                          for s in site.get("measurement_coverage", [])) for site in measured)
    if unsupported:
        lines.extend([
            f"Of {len(measured)} entries with source measurements, {mapped} have mapped candidate DNA; "
            f"{unsupported} lack a supported model for an experimental nomination; "
            f"{len(measured) - mapped - unsupported} have other unresolved models. "
            "These counts describe this workflow's representations, not experiment quality. "
            "Source reference identity and tested sample DNA require separate evidence.", "",
            "Missing-field reasons describe metadata supplied to this workflow. An unreported control "
            "or QC field does not establish that the experiment lacked controls or failed quality checks.", "",
        ])
    coverage_rows, measurement_rows = [], []
    for site in ledger["candidates"]:
        coverage = site.get("measurement_coverage", [])
        if not coverage:
            coverage_rows.append([site["guide_id"], site["site_id"], "Unknown", "Unknown", "Not assessed",
                                  "Regenerate this ledger to assess measurement coverage."])
        for scenario in coverage:
            allele = {"ref": "Reference", "alt": "Alternate"}.get(scenario["allele"], "Unknown DNA")
            dna = f"{allele}: {scenario.get('sequence_genomic_forward') or 'sequence unresolved'}"
            if scenario.get("mapping_reason") == "experimental_candidate_model_unsupported":
                dna = "candidate allele model unsupported; see nomination records"
            coverage_rows.append([
                site["guide_id"], site["site_id"], scenario["configuration_id"], dna,
                scenario["evidence_label"], "; ".join(scenario["follow_up_questions"]),
            ])
        for measurement in site.get("measurements", []):
            assessments = [assessment for scenario in coverage for assessment in scenario["measurements"]
                           if all(assessment.get(field) == measurement.get(field)
                                  for field in ("source_record_id", "source_version_or_citation", "assay"))]
            supported = sum(a["usable_for_editing"] for a in assessments)
            reasons = sorted({reason for a in assessments for reason in a["reasons"]})
            result = measurement["status"]
            if measurement.get("value") is not None:
                result += f"; {measurement['value']} {measurement.get('unit') or ''}"
            lod = measurement.get("limit_of_detection")
            measurement_rows.append([
                site["site_id"],
                f"{measurement.get('editor') or 'editor unknown'}; "
                f"{measurement.get('allele') or 'allele unknown'}; "
                f"{measurement.get('allele_sequence') or 'tested sequence unknown'}; "
                f"resolution: {measurement.get('allele_resolution') or 'unknown'}; "
                f"sequence evidence: {measurement.get('allele_evidence') or 'not supplied'}; "
                f"haplotype: {measurement.get('haplotype_id') or 'unknown'}; "
                f"phase: {measurement.get('phase_status') or 'unknown'}",
                f"{measurement['assay']}; {measurement.get('biological_context') or 'cells unknown'}; "
                f"measures: {measurement.get('measurement_kind') or 'not specified'}",
                result,
                f"{lod['value']} {lod['unit']}" if lod else "Not reported",
                f"{measurement.get('qc_status') or 'not reported'}; "
                f"{measurement.get('qc_note') or 'no quality note'}; "
                f"treated: {measurement.get('treated_group') or 'unknown'}; "
                f"control: {measurement.get('control_group') or 'unknown'}",
                f"Usable for {supported} listed DNA scenarios; " + "; ".join(reasons),
                f"{measurement['source_version_or_citation']}; {measurement['source_record_id']}",
            ])
    _table(lines, ["Guide", "Site", "Reported configuration", "DNA version (genomic forward)", "Editing evidence", "Follow-up question"],
           coverage_rows)
    _table(lines, ["Site", "Experiment identity", "Assay and cells", "Source result", "Detection limit",
                   "Quality and controls", "Use as editing evidence", "Source"], measurement_rows)
    requirements = sorted({
        f"{site['guide_id']}: {item['assay']}"
        for site in ledger["candidates"] for item in site.get("required_confirmation", [])
    })
    if requirements:
        lines.extend(["HBG locus checks still requiring review:", ""])
        lines.extend(f"- {_cell(item)}" for item in requirements)
        lines.append("")




def render_report(
    design_records: dict,
    ledger: dict,
    *,
    design_input: dict | None = None,
    ledger_input: dict | None = None,
) -> str:
    """Write the same Markdown for the same two input records without changing them."""
    _validate_inputs(design_records, ledger)
    guide_ids = {item["guide_id"] for item in ledger["guides"]}
    lines = [
        "# Base Edit Safety candidate review",
        "",
        f"_{DISCLAIMER}_",
        "",
        "Coordinates in candidate tables are zero-based and half-open.",
        "",
    ]

    _section(lines, "Design summary")
    _table(
        lines,
        ["Guide", "Gene", "Editor", "Intended offsets", "Evidence tier", "Paralog flag"],
        _design_rows(design_records, guide_ids),
    )

    _section(lines, "Nomination summary")
    summary = ledger.get("summary", {})
    lines.extend(
        [
            f"- Candidate sites: {_cell(summary.get('candidate_sites', 0))}",
            f"- Computational buckets: {_cell(summary.get('by_source_bucket', {}))}",
            f"- Reconciliation states: {_cell(summary.get('by_reconciliation_status', {}))}",
            (f"- Evidence states: {_cell(summary.get('by_evidence_state', {}))}"
             if all("measurement_coverage" in site for site in ledger["candidates"])
             else "- Legacy evidence summaries require regeneration for allele coverage."),
            "",
        ]
    )
    for item in ledger.get("import_summaries", []):
        status = "Incomplete import" if item["skipped_rows"] else "Import row accounting"
        lines.append(
            f"- {status}: {_cell(item['normalized_output'])}; input: {item['input_data_rows']}; "
            f"imported: {item['emitted_nomination_rows']}; rejected/skipped: {item['skipped_rows']}; "
            f"reasons: {_cell(item['skipped_by_reason'])}."
        )
    profile = ledger.get("search_profile", {})
    for field in ("reference_id", "intended_use", "variant_source", "retention_rule"):
        if profile.get(field):
            lines.append(f"- {_cell(field)}: {_cell(profile[field])}")
    for source in ledger.get("nomination_source_provenance", []):
        declaration = source.get("source_declarations", {}).get("crisprme")
        if declaration:
            lines.append(
                f"- CRISPRme+ retention for {_cell(declaration['guide_id'])}: "
                f"{_cell(declaration['retention'])}. Completeness applies only to the declared search."
            )
    lines.append("")
    warnings = [
        warning
        for site in ledger["candidates"]
        for warning in site.get("reconciliation", {}).get("warnings", [])
    ]
    if warnings:
        lines.append("Conditional QC warnings:")
        lines.append("")
        lines.extend(f"- {_cell(warning)}" for warning in sorted(set(warnings)))
        lines.append("")

    _section(lines, "Candidate evidence table")
    _table(
        lines,
        [
            "Site ID",
            "Guide",
            "Interval",
            "Bucket",
            "Reconciliation",
            "Possible source differences",
            "Source scores",
            "Evidence state",
            "Measurement",
            "Editability",
            "Imported gene/locus labels",
            "Decision",
            "Reason",
            "Note",
            "Assays",
        ],
        _candidate_rows(ledger),
    )

    configurations = [
        (site, block) for site in ledger["candidates"]
        for block in site.get("allele_configurations", [])
        if block.get("configuration_coverage") == "reported_highest_CFD_only"
    ]
    if configurations:
        lines.extend([
            "Reported highest-CFD configurations are not exhaustive alignments or haplotypes. "
            "Substrate positions describe only supported reported sequences in the declared editor window; "
            "they do not establish activity. No substrate in that window does not exclude editing elsewhere.",
            "",
            "Opaque upstream MAF is not converted to ALT AF, minor-allele frequency, pooled frequency, "
            "or patient prevalence. Missing genotype and phase information remains unknown. "
            "Source differences alone do not demonstrate a CRISPOR miss without a compatible comparison.",
            "",
            "[CRISPRme+ documents](https://github.com/pinellolab/crisprme-plus/blob/main/docs/INPUT_FORMATS.md#variant-driven-pam-creation-and-disruption) "
            "that alternate alleles destroying the required PAM are omitted. Their absence here does not "
            "establish that those variants are absent.",
            "",
        ])
        _table(lines, ["Guide", "Editor", "Class", "Declared window (1-based inclusive)"], [
            [guide["guide_id"], guide.get("editor"), guide.get("editor_class"),
             f"{guide.get('editing_window_start_1based')}-{guide.get('editing_window_end_1based_inclusive')}"]
            for guide in ledger["guides"]
        ])
        rows = []
        for site, block in configurations:
            assessments = [
                state for state in site.get("editability", {}).get("allele_states", [])
                if state.get("configuration_id") == block["configuration_id"]
            ]
            rows.append([
                site["site_id"], block["configuration_id"],
                [{"allele": state["allele"], "alignment": state["alignment"],
                  "mapping": state.get("mapping_status"),
                  "genomic_forward_changes": state.get("genomic_forward_changes")}
                 for state in block["allele_states"]],
                block.get("pam_change"), assessments, block.get("frequency"),
                block.get("phase_status"), block.get("interpretation"),
            ])
        _table(lines, ["Site ID", "Configuration", "Alleles and mapping", "PAM change",
                       "Window assessment", "Frequency and source", "Phase", "Interpretation"], rows)

    _measurement_tables(lines, ledger)

    _section(lines, "Panel decision summary")
    lines.extend(
        [
            f"- Decisions: {_cell(summary.get('by_panel_decision', {}))}",
            f"- Panel notes: {_cell(ledger.get('export_metadata', {}).get('panel_notes')) or 'None supplied.'}",
            "",
            "Panel exclusion records one proposed panel choice; it is not evidence that a site is safe.",
            "",
        ]
    )

    _section(lines, "Reproducibility and run diff")
    runtime = ledger.get("reproducibility", {})
    lines.extend(
        [
            "- Package: base-edit-safety "
            f"{_cell(runtime.get('base_edit_safety_version')) or 'version not recorded'}",
            f"- Python: {_cell(runtime.get('python_version')) or 'version not recorded'}",
        ]
    )
    if design_input is not None:
        lines.append(
            f"- Design records: `{_cell(design_input['file'])}` (`sha256:{_cell(design_input['sha256'])}`)"
        )
    if ledger_input is not None:
        lines.append(
            f"- Candidate ledger: `{_cell(ledger_input['file'])}` (`sha256:{_cell(ledger_input['sha256'])}`)"
        )
    for item in ledger.get("inputs", []):
        lines.append(
            f"- {_cell(item.get('kind'))}: `{_cell(item.get('file'))}` "
            f"(`sha256:{_cell(item.get('sha256'))}`)"
        )
    generated = ledger.get("reproducibility", {}).get("generated_at_utc")
    lines.extend(
        [
            f"- Ledger generated at: {_cell(generated)}",
            f"- Import summaries: {_cell(ledger.get('import_summaries', []))}",
            f"- Exact-identity diff: {_cell(ledger.get('diff_summary')) or 'Not requested.'}",
            "",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


def write_report(design_path: Path, ledger_path: Path, output_path: Path) -> None:
    design_records = _load_json(design_path, "design-record set")
    ledger = _load_json(ledger_path, "candidate ledger")
    rendered = render_report(
        design_records,
        ledger,
        design_input={"file": design_path.name, "sha256": _sha256(design_path)},
        ledger_input={"file": ledger_path.name, "sha256": _sha256(ledger_path)},
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(rendered, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design-records", type=Path, required=True)
    parser.add_argument("--candidate-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        write_report(args.design_records, args.candidate_ledger, args.output)
    except CandidateReportError as exc:
        parser.error(str(exc))
    print(f"wrote candidate review report to {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
