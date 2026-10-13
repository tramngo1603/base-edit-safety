"""Read completed analyses and write separate candidate exports for follow-up."""
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil

from . import lab
from .population.common import checksum, save_json

COORDINATES = "0-based, half-open; end excluded"
SEARCH_SCOPE = {"spacer_length": 20, "pam": "NGG", "strands": "both", "gaps": False}
EXPORT_FIELDS = {
    "reference_build": "Assembly recorded in request.json.",
    "coordinate_convention": "0-based, half-open: start0 is included and end0 is excluded.",
    "spacer_sequence": "Input guide spacer, from the checksum-verified guide file; blank if unavailable.",
    "protospacer_sequence": "First 20 bases of the 23-base target, in guide orientation.",
    "pam_sequence": "Last 3 bases of that target, in guide orientation.",
    "locus_status": "Original locus comparison status; reference_only_run, retained, gained, lost, unresolved or unpaired_indel_path.",
    "sequence_changed": "Original comparison flag. A retained locus with a changed sequence is not a gained/lost locus.",
    "comparison_reason": "Original locus comparison reason, when recorded.",
    "sample_id": "Selected VCF sample; blank for reference-only analyses.",
    "max_mismatches": "Search ceiling from the original request, not the candidate mismatch count.",
    "pam_class": "NGG, the PAM supported by this search.",
    "searched_strands": "Both reference strands are searched.",
    "gaps_allowed": "False; this search does not nominate gapped guide-target alignments.",
}
RAW_FIELDS = {
    "configuration_id": "Stable identifier of this sequence configuration; unchanged from the source.",
    "locus_id": "Original locus identifier. Multiple sequences/copies can belong to a locus.",
    "guide_id": "Original guide identifier.",
    "contig": "Reference chromosome/contig name.",
    "start0": "Reference start, 0-based and included. Read coordinate_path for sample indel paths.",
    "end0": "Reference end, excluded. Read coordinate_path for sample indel paths.",
    "strand": "Target orientation relative to the reference (+ or -).",
    "target_sequence": "Original target sequence including PAM, in guide orientation.",
    "mismatches": "Existing engine mismatch count between the guide and target spacer.",
    "reference_present": "Whether this sequence configuration is present in the reference search.",
    "sample_copies": "JSON array of local sample copy indexes; not genome-wide haplotype identifiers.",
    "phase_status": "Original phase status records (JSON).",
    "coordinate_path": "Original coordinate/path records (JSON); retains indel and variant details where present.",
    "evidence_links": "Original linked nomination evidence (JSON). Empty means no linked evidence in this bundle, not experimentally negative.",
}


def row_matches(row, query="", filters=None):
    """Use identical search and filter rules for the table and full-file export."""
    for key, value in (filters or {}).items():
        if not value:
            continue
        if key == "difference":
            changed = row.get("sequence_changed") == "True"
            if value == "changed" and not (row.get("status") == "retained" and changed):
                return False
            if value == "differences" and row.get("status") == "retained" and not changed:
                return False
        elif row.get(key) != value:
            return False
    query = query.strip().casefold()
    return not query or query in "\t".join(str(value or "") for value in row.values()).casefold()


def read_analysis(folder):
    folder = Path(folder).expanduser().resolve()
    if not (folder / "run-manifest.json").is_file() and (folder / "results").is_dir():
        folder = folder / "results"
    try:
        manifest = json.loads((folder / "run-manifest.json").read_text())
        request = json.loads((folder / "request.json").read_text())
        coverage = json.loads((folder / "coverage.json").read_text())
    except (OSError, ValueError) as exc:
        raise ValueError("Choose a completed results folder (or its parent). This folder is missing readable run records.") from exc
    if manifest.get("status") not in {"complete", "complete_reference_sample_failed"}:
        raise ValueError("This run is incomplete. Its files remain available in the folder, but cannot be opened as completed results.")
    if not (folder / "candidates.tsv").is_file():
        raise ValueError("The candidate file is missing. Open the original complete result bundle.")
    return folder, manifest, request, coverage


def saved_status(folder):
    """Confirm completed output from the manifest written last by the runner."""
    try:
        folder, manifest, _, _ = read_analysis(folder)
        if manifest["status"] != "complete":
            return "Incomplete output", "Reference results are retained; sample analysis did not complete."
        if manifest.get("example_verification", {}).get("status") == "mismatched":
            return "Verification failed", "Output is retained, but the example did not match its accepted results."
        outputs = manifest.get("outputs", {})
        required = {"request.json", "coverage.json", "candidates.tsv", "locus-changes.tsv", "exclusions.tsv", "report.html"}
        if not required.issubset(outputs):
            return "Save not confirmed", "The completed output record is missing required file receipts."
        for name, receipt in outputs.items():
            path = folder / name
            if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size != receipt.get("bytes"):
                return "Output needs attention", f"Missing, unreadable or changed output: {name}."
        return "Saved", "Completion is recorded and the output files match their recorded sizes."
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return "Save not confirmed", str(exc)


def analysis_folder(folder):
    """Include prepared inputs only for the runner's established folder layout."""
    folder = Path(folder).resolve()
    return folder.parent if folder.name == "results" and (folder.parent / "inputs").is_dir() else folder


def export_copy(folder, destination, kind):
    """Copy an existing HTML report or analysis folder without changing the source."""
    folder, manifest, _, _ = read_analysis(folder)
    if kind not in {"report", "analysis"}:
        raise ValueError("Choose an HTML report or an analysis copy.")
    source = analysis_folder(folder)
    destination = Path(destination).expanduser().resolve()
    if destination == source or source in destination.parents:
        raise ValueError("Choose an export destination outside the original analysis folder.")
    if destination.exists() or not destination.parent.is_dir():
        raise ValueError("Choose an existing destination folder; existing exports are never overwritten.")
    report = folder / "report.html"
    if kind == "report" and not report.is_file():
        raise ValueError("This analysis has no HTML report. Its other results remain available.")
    names = ["report.html"] if kind == "report" else list(manifest.get("outputs", {}))
    for name in names:
        path = folder / name
        if not path.is_file():
            raise ValueError(f"The result bundle is missing {name}.")
        expected = manifest.get("outputs", {}).get(name, {}).get("sha256")
        if expected and checksum(path) != expected:
            raise ValueError(f"{name} no longer matches the run checksum. Use the original complete bundle.")
    temporary = destination.with_name(destination.name + ".partial")
    temporary.mkdir(mode=0o700)
    try:
        if kind == "report":
            shutil.copy2(report, temporary / "report.html")
        else:
            shutil.copytree(source, temporary, dirs_exist_ok=True, symlinks=True)
        copied_results = temporary / folder.relative_to(source) if kind == "analysis" else temporary
        for name in names:
            if checksum(folder / name) != checksum(copied_results / name):
                raise ValueError("A result file changed during copying. Reopen the source and try again.")
        if destination.exists():
            raise ValueError("The destination was created during export. Choose another folder.")
        temporary.rename(destination)
        return {"kind": kind, "source_results": str(folder)}
    except Exception:
        shutil.rmtree(temporary)
        raise


def pinned_input(folder, spec):
    """Find an original or relocated input only when its recorded checksum matches."""
    if not spec.get("path") or not spec.get("sha256"):
        return None
    original = Path(spec["path"])
    choices = [Path(folder) / original, Path(folder).parent / original, original,
               Path(folder).parent / "inputs" / original.name]
    for path in dict.fromkeys(choices):
        if path.is_file() and checksum(path) == spec["sha256"]:
            return path.resolve()
    return None


def guide_sequences(folder, request):
    spec = request.get("guides", {})
    path = pinned_input(folder, spec)
    if path is None:
        for name in ("examples/small/guides.tsv", "examples/research/guides.tsv"):
            candidate = lab.files() / name
            if candidate.is_file() and checksum(candidate) == spec.get("sha256"):
                path = candidate
                break
    return {row["guide_id"]: row["spacer_sequence"] for row in lab.read_guides(path.read_text())} if path else {}


def analysis_context(folder):
    path = Path(folder).parent / "analysis-context.json"
    if path.is_file():
        try:
            return json.loads(path.read_text())
        except (OSError, ValueError):
            return {"note": "The optional analysis context file could not be read."}
    return {}


def export_candidates(folder, destination, *, mode="all", query="", filters=None, selected_ids=()):
    """Stream all matching rows to a new folder without changing the result bundle."""
    folder, manifest, request, coverage = read_analysis(folder)
    if mode not in {"all", "filtered", "selected"}:
        raise ValueError("Choose all candidates, matching candidates, or selected rows.")
    selected = set(selected_ids)
    if mode == "selected" and not selected:
        raise ValueError("Select at least one candidate row first.")
    if mode != "filtered":
        query, filters = "", {}
    filters = {key: value for key, value in (filters or {}).items() if value}
    if set(filters) - {"guide_id", "contig"}:
        raise ValueError("Candidate exports support the guide and chromosome filters shown in Candidates.")
    destination = Path(destination).expanduser().resolve()
    if destination == folder or folder in destination.parents:
        raise ValueError("Save exports outside the original results folder to keep the complete bundle intact.")
    if destination.exists():
        raise ValueError("This export folder already exists. Choose a new name; existing files are preserved.")
    if not destination.parent.is_dir():
        raise ValueError("Choose an existing folder for the export.")
    source_names = ["run-manifest.json", "request.json", "candidates.tsv", "locus-changes.tsv", "coverage.json", "exclusions.tsv"]
    source_hashes = {}
    for name in source_names:
        path = folder / name
        if not path.is_file():
            raise ValueError(f"The result bundle is missing {name}. Use the original complete bundle.")
        digest = checksum(path)
        expected = manifest.get("outputs", {}).get(name, {}).get("sha256")
        if expected and digest != expected:
            raise ValueError(f"{name} no longer matches the run checksum. Use the original complete bundle.")
        source_hashes[name] = digest
    guides = guide_sequences(folder, request)
    statuses = {}
    with (folder / "locus-changes.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            statuses[row["locus_id"]] = (row["status"], row["sequence_changed"], row["reason"])
    # Keep a failed write distinct from a completed export and never overwrite a folder.
    temporary = destination.with_name(destination.name + ".partial")
    temporary.mkdir(mode=0o700)
    count, total, loci, found = 0, 0, set(), set()
    try:
        with (folder / "candidates.tsv").open(newline="") as source, (temporary / "candidates.tsv").open("w", newline="") as out:
            reader = csv.DictReader(source, delimiter="\t")
            raw_fields = reader.fieldnames or []
            if not {"configuration_id", "locus_id", "guide_id", "target_sequence"}.issubset(raw_fields):
                raise ValueError("The candidate table is missing required scientific identifiers or sequences.")
            fields = raw_fields + list(EXPORT_FIELDS)
            writer = csv.DictWriter(out, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                total += 1
                if mode == "filtered" and not row_matches(row, query, filters):
                    continue
                if mode == "selected" and row["configuration_id"] not in selected:
                    continue
                found.add(row["configuration_id"])
                loci.add(row["locus_id"])
                sequence = row["target_sequence"]
                valid_sequence = re.fullmatch(r"[ACGT]{23}", sequence)
                status, changed, reason = statuses.get(row["locus_id"], ("not_recorded", "", ""))
                row.update(reference_build=request.get("assembly", ""), coordinate_convention=COORDINATES,
                           spacer_sequence=guides.get(row["guide_id"], ""),
                           protospacer_sequence=sequence[:20] if valid_sequence else "",
                           pam_sequence=sequence[20:] if valid_sequence else "",
                           locus_status=status, sequence_changed=changed, comparison_reason=reason,
                           sample_id=request.get("sample", {}).get("id", ""),
                           max_mismatches=request.get("max_mismatches", ""),
                           pam_class="NGG", searched_strands="both", gaps_allowed="False")
                writer.writerow(row)
                count += 1
        if mode == "selected" and found != selected:
            raise ValueError("Some selected identifiers are missing from the source. Reopen the results and select again.")
        with (folder / "locus-changes.tsv").open(newline="") as source, (temporary / "locus-changes.tsv").open("w", newline="") as out:
            reader = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(out, reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row["locus_id"] in loci:
                    writer.writerow(row)
        for name in ("request.json", "run-manifest.json", "coverage.json", "exclusions.tsv"):
            shutil.copyfile(folder / name, temporary / name)
        # Detect a source changed while export was reading it.
        if any(checksum(folder / name) != digest for name, digest in source_hashes.items()):
            raise ValueError("A result file changed during export. Reopen the original complete result bundle and try again.")
        explanation = """Candidate export for experimental follow-up

candidates.tsv retains every original candidate field and adds the fields described below.
One row describes a sequence configuration; a locus can have several configurations.
Unresolved and unpaired loci remain distinct from lost loci. Empty evidence_links
means no linked evidence in this bundle; assay status is otherwise untested here.
At positions without VCF records, sample sequences use reference bases as a modeling
assumption. Genotype and callability at those positions remain unknown.

locus-changes.tsv contains the original comparison records for exported loci.
coverage.json and exclusions.tsv describe the WHOLE RUN, not the selected subset.
Reference evaluable windows / possible reference windows describe reference search
coverage. Sample modeled bases / reference bases and comparable reference windows
/ possible reference windows describe sample modeling. Evaluated sample-copy
windows are a separate copy-specific denominator. Callability remains unmeasured.
request.json and run-manifest.json are unchanged copies from the original run. Their paths and
output receipts refer to the ORIGINAL bundle, not to this derived export.
export.json records exact selection criteria, row counts, source hashes, search
parameters, sample source details and output hashes. Optional experiment labels are
context only and do not select an editor model or change search parameters.
Keep the original complete result bundle for reproduction and full accounting.

Field dictionary
"""
        explanation += "\n".join(f"{key}: {value}" for key, value in (RAW_FIELDS | EXPORT_FIELDS).items()) + "\n"
        (temporary / "READ-ME.txt").write_text(explanation)
        export = {"format_version": "1", "created_utc": datetime.now(timezone.utc).isoformat(),
                  "source_results": str(folder), "source_sha256": source_hashes,
                  "source_status": manifest.get("status"), "example_verification": manifest.get("example_verification"),
                  "selection": {"mode": mode, "query": query, "filters": filters,
                                "selected_configuration_ids": sorted(selected) if mode == "selected" else [],
                                "query_rule": "case-insensitive substring across original candidate fields; equality filters are combined with AND",
                                "display_limit_applied": False},
                  "candidate_rows": count, "source_candidate_rows": total, "loci": len(loci),
                  "guide_spacers": "checksum-verified input" if guides else "unavailable; spacer_sequence left blank",
                  "reference": request.get("reference"), "reference_build": request.get("assembly"),
                  "coordinate_convention": COORDINATES,
                  "search": {**SEARCH_SCOPE, "max_mismatches": request.get("max_mismatches"), "contigs": request.get("contigs")},
                  "sample": request.get("sample"), "reference_callability": coverage.get("reference_callability", "not_assessed"),
                  "sample_coverage_basis": coverage.get("sample_basis", manifest.get("coverage_basis")),
                  "analysis_context": analysis_context(folder)}
        export["files"] = {path.name: {"sha256": checksum(path), "bytes": path.stat().st_size} for path in temporary.iterdir()}
        save_json(temporary / "export.json", export)
        if destination.exists():
            raise ValueError("The destination was created during export. Choose a new folder name.")
        temporary.rename(destination)
        return export
    except Exception:
        shutil.rmtree(temporary)
        raise
