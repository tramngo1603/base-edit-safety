"""Render bounded sequence comparisons as Markdown and TSV files."""
from __future__ import annotations

import csv
import html
import io
import json
import re


COLUMNS = ("candidate_id", "sample_id", "contig", "start0", "end0", "strand",
           "target_sequence", "pam", "mismatches", "reference_target_sequence", "change",
           "comparison_basis", "model_carriage", "evidence_supported_carriage", "local_copies",
           "phase_status", "pairing_reason", "target_path", "variant_ids")
WARNING_COLUMNS = ("sample_id", "code", "detail", "count")


def _cell(value):
    if value is None:
        return "unknown"
    if isinstance(value, (dict, list, bool)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return str(value)


def _md(value):
    text = html.escape(_cell(value), quote=True)
    return re.sub(r"([\\`|*_{}\[\]()#!])", r"\\\1", text).replace("\r", " ").replace("\n", "<br>")


def _tsv(columns, rows):
    output = io.StringIO(newline="")
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    writer.writerows([_cell(row[column]) for column in columns] for row in rows)
    return output.getvalue()


def render(result):
    """Return report, candidate and warning text for a computed result, preserving row order."""
    scope, domain, summary = result["scope"], result["domain"], result["summary"]
    region, guide = scope["region"], scope["guide"]
    lines = ["# Bounded computational sequence comparison", "",
             f"Biological status: **{_md(domain['biological_status'])}**. Computational status: {_md(domain['computational_status'])}.", "",
             "Candidates are sequence matches within the declared region.", "",
             "Sample sequences use reference bases at positions without VCF records; genotypes at those positions remain unknown. "
             "Copy indexes label the modeled sequences. Links to physical chromosomes, donors, specimens and assays remain "
             "unknown. A single unphased heterozygote represents an unordered pair.", "",
             "Biological abstention reasons: " + ", ".join(_md(v) for v in domain["reasons"]) + ".", "",
             "## Inputs and scope", "",
             f"- Format: {_md(result['format_version'])}; purpose: {_md(result['purpose'])}; input kind: **{_md(result['input_kind'])}**.",
             f"- Run: {_md(result['run_id'])}.",
             f"- Assembly: {_md(scope['assembly'])}; region: {_md(region['contig'])} [{region['start']}, {region['end']}).",
             f"- Coordinates: {_md(scope['coordinate_convention'])}.",
             f"- Guide: {_md(guide['id'])}; spacer: {_md(guide['spacer'])}; PAM: {_md(guide['pam'])}; maximum mismatches: {guide['max_mismatches']}.",
             "- Both strands; full 23-base windows only; no alignment bulges; outside-region and boundary-spanning windows are not covered.",
             "- VCF sample IDs: " + ", ".join(_md(v) for v in scope["sample_ids"]) + ".",
             "- Genotype policy: " + _md(scope["genotype_policy"]) + ".", "",
             "| Input | SHA-256 | Bytes |", "| --- | --- | ---: |"]
    lines.extend(f"| {_md(row['role'])} | {_md(row['sha256'])} | {row['size_bytes']} |" for row in result["inputs"])
    lines.extend(["", "## Sequence results", "", "| Quantity | Value |", "| --- | --- |"])
    lines.extend(f"| {_md(key)} | {_md(value)} |" for key, value in summary.items())
    lines.extend(["", "Candidate counts refer to distinct sequence configurations. "
                  "Window totals sum the reference and both modeled copies of each sample. "
                  "The comparison basis is retained_callset_model; evidence-supported carriage stays unknown.", "",
                  "[Complete candidates TSV](candidates.tsv) preserves allele identities, reference counterparts, local copies, "
                  "phase and exact coordinate paths. [result.json](result.json) contains the same producer facts. "
                  "An indel path without an exact reference counterpart remains unpaired; a bounding interval does not replace its path.", ""])
    if result["candidates"]:
        lines.extend(["| Sample | Region / strand | Target + PAM | Mismatches | Change | Modeled carriage | Phase |",
                      "| --- | --- | --- | ---: | --- | --- | --- |"])
        for row in result["candidates"]:
            location = f"{row['contig']}:{row['start0']}-{row['end0']} {row['strand']}"
            lines.append("| " + " | ".join(_md(v) for v in (row["sample_id"], location,
                row["target_sequence"] + row["pam"], row["mismatches"], row["change"],
                row["model_carriage"], row["phase_status"])) + " |")
    else:
        lines.append("No supported sequence nomination was recovered. Unresolved windows remain excluded.")
    lines.extend(["", "## Warnings and unresolved inputs", "",
                  "Every warning is also in [warnings.tsv](warnings.tsv). Counts describe the unit stated in each detail; do not sum different units.", "",
                  "| Sample | Code | Count | Detail |", "| --- | --- | ---: | --- |"])
    for row in result["warnings"]:
        sample = "all / reference" if row["sample_id"] is None else row["sample_id"]
        lines.append("| " + " | ".join(_md(v) for v in (sample, row["code"], row["count"], row["detail"])) + " |")
    lines.extend(["", "## Execution limits and evidence", "", "| Limit | Value |", "| --- | ---: |"])
    lines.extend(f"| {_md(key)} | {_md(value)} |" for key, value in result["limits"].items())
    lines.extend(["", "Time and process peak-memory limits are checked cooperatively between bounded operations and during search. "
                  "[run-manifest.json](run-manifest.json) records the producing code, UTC completion time, actual host metrics, "
                  "input identities and every other output's checksum. Metadata timing may differ between repeated runs.", "",
                  "The evidence files retain normalized variants and original ALT indexes, genotype audit, sequence paths and gaps, "
                  "candidate configurations, search hits, sample-carriage models, nomination changes and search/population audits. "
                  "The reused raw record fields named donor_id refer only to the declared VCF sample IDs in this route. "
                  "No donor or specimen linkage was supplied.", ""])
    return {"report.md": "\n".join(lines), "candidates.tsv": _tsv(COLUMNS, result["candidates"]),
            "warnings.tsv": _tsv(WARNING_COLUMNS, result["warnings"])}
