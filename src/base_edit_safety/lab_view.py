"""Display labels and coverage summaries shared by the local interfaces."""
import json

STATUS_LABELS = {
    "retained": "Retained", "gained": "Gained", "lost": "Lost",
    "unresolved": "Unresolved", "unpaired_indel_path": "Unpaired indel path",
    "reference_only_run": "Reference only",
}

REASON_LABELS = {
    "unknown_phase": "Unknown phase",
    "unresolved_unphased": "Unresolved phase",
    "phase_scope": "Phase across sites unknown",
    "sample_contig_excluded": "Sample chromosome excluded",
    "sample_ploidy_undeclared": "Copy count not declared",
    "unknown_or_changing_ploidy": "Unknown or changing copy count",
    "unknown_ploidy": "Unknown copy count",
    "zero_ploidy": "Zero declared copies",
    "zero_declared_copies": "Zero declared copies",
    "overlapping_variants": "Overlapping variants",
    "overlapping_carried_edits": "Overlapping carried variants",
    "unsupported_carried_allele": "Unsupported carried allele",
    "unbounded_unsupported_allele": "Allele extent unresolved",
    "unresolved_spanning_deletion_extent": "Deletion extent unresolved",
    "local_reconstruction_resource_limit": "Local reconstruction limit",
    "sample_scope_or_ploidy_unknown": "Sample scope or copy count unknown",
    "sample_run_failed": "Sample analysis failed",
    "missing_genotype": "Missing genotype",
}


def reason_label(value):
    """Short display labels only; row details retain the original reason."""
    if not value:
        return "—"
    if value.startswith(("[", "{")):
        try:
            records = json.loads(value)
            records = records if isinstance(records, list) else [records]
            labels = list(dict.fromkeys(reason_label(row["reason"]) for row in records
                                        if isinstance(row, dict) and row.get("reason")))
            if labels:
                return "; ".join(labels)
        except (ValueError, TypeError):
            pass
    code = value.split(":", 1)[0]
    if code in REASON_LABELS:
        return REASON_LABELS[code]
    return value.replace("_", " ") if len(value) <= 70 else "See row details"


def coverage_brief(coverage, config):
    """Keep reference windows separate from conditional sample-modeled bases."""
    totals = coverage_totals(coverage)
    rows = coverage.get("chromosomes", [])
    reference = "Reference coverage not recorded"
    if rows and all("reference_evaluable_oriented_windows_per_guide" in row for row in rows):
        reference = "Reference windows: " + percent(totals["reference_evaluable_oriented_windows_per_guide"],
                                                    totals["reference_oriented_windows_per_guide"]) + " evaluable"
    if not config.get("sample"):
        return reference + " · No sample VCF"
    model = (percent(totals["modeled_reference_bases"], totals["reference_bases"]) + " of reference bases"
             if totals["reference_bases"] else "coverage not recorded")
    if totals["modeled_reference_bases"] < totals["reference_bases"] or coverage.get("exclusion_summary"):
        model += "; exclusions remain"
    return reference + " · Conditional sample model: " + model


def percent(part, whole):
    if not whole:
        return "—"
    value = 100 * part / whole
    if part < whole and round(value, 1) == 100:
        return "<100%"
    return f"{value:.1f}%"


def result_status(state):
    status = state.get("status", "")
    if status == "complete":
        return "Search complete"
    if status.startswith("complete_reference"):
        return "Reference complete; sample analysis incomplete"
    return "Run did not complete"


def coverage_totals(coverage):
    rows = coverage.get("chromosomes", [])
    keys = ("reference_bases", "modeled_reference_bases",
            "reference_evaluable_oriented_windows_per_guide",
            "reference_oriented_windows_per_guide", "comparable_oriented_windows_per_guide",
            "sample_evaluated_copy_oriented_windows_per_guide")
    return {key: sum(row.get(key, 0) or 0 for row in rows) for key in keys}


def donor_scope(coverage, config):
    if not config.get("sample"):
        return "Reference only · no sample VCF"
    totals = coverage_totals(coverage)
    scope = "Sample model"
    if totals["reference_bases"]:
        scope += " · " + percent(totals["modeled_reference_bases"], totals["reference_bases"]) + " of reference bases"
    unknown = list(dict.fromkeys(row["contig"] for row in
                                coverage.get("sample_undeclared_ploidy_intervals", [])))
    if unknown:
        # A ploidy interval can cover part of a chromosome; do not imply a whole-chromosome gap.
        lengths = {row["contig"]: row.get("reference_bases") for row in coverage.get("chromosomes", [])}
        whole = all(row.get("start") == 0 and row.get("end") == lengths.get(row["contig"])
                    for row in coverage.get("sample_undeclared_ploidy_intervals", []))
        scope += " · " + ", ".join(unknown) + (" unmodeled" if whole else " have undeclared-copy regions")
    if totals["modeled_reference_bases"] < totals["reference_bases"] or coverage.get("exclusion_summary"):
        scope += " · exclusions remain"
    return scope


def coverage_facts(coverage, config):
    totals = coverage_totals(coverage)
    bases = totals["reference_bases"]
    windows = totals["reference_oriented_windows_per_guide"]
    modeled = totals["modeled_reference_bases"]
    comparable = totals["comparable_oriented_windows_per_guide"]
    facts = [("Reference search scope", ", ".join(coverage.get("searched_contigs", [])) or "Not recorded")]
    if coverage.get("chromosomes"):
        if all("reference_evaluable_oriented_windows_per_guide" in row for row in coverage["chromosomes"]):
            facts.append(("Reference search windows", f"{totals['reference_evaluable_oriented_windows_per_guide']:,} evaluable / {windows:,} possible"))
        if config.get("sample"):
            facts.extend([
                ("Sample-modeled bases", f"{modeled:,} / {bases:,} reference bases ({percent(modeled, bases)})"),
                ("Comparable reference windows", f"{comparable:,} / {windows:,} ({percent(comparable, windows)})"),
                ("Evaluated sample copy windows", f"{totals['sample_evaluated_copy_oriented_windows_per_guide']:,} — a separate copy-specific denominator"),
            ])
    return facts
