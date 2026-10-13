"""Versioned study and planner contracts with explicit unknown values."""
from __future__ import annotations

from importlib.resources import files as package_files
from functools import lru_cache
import math
import re

from .common import PopulationError, canonical, fields, identity, integer, nonempty, parse_json, unique

MODES = {"synthetic_test", "benchmark_simulation", "research_sample_planning"}
BASES = {"synthetic_truth", "retained_callset_model", "callable_supported"}
ALGORITHMS = {"greedy_cost", "greedy_coverage", "canonical_order", "random_order", "exact_tiny"}
LIMITS = {"max_donors": 32, "max_reference_bases": 100000,
          "max_reconstructed_bases": 5000000, "max_configurations": 100000,
          "max_hit_rows": 250000, "max_variant_rows": 100000}


def validate_schema(value, schema, location="request"):
    """Validate the closed JSON Schema vocabulary shipped with this package."""
    if "oneOf" in schema:
        matches = 0
        for option in schema["oneOf"]:
            try:
                validate_schema(value, option, location)
                matches += 1
            except PopulationError:
                pass
        if matches != 1:
            raise PopulationError(location + ": must match exactly one supported record contract")
        return
    if "const" in schema and (value != schema["const"] or type(value) is not type(schema["const"])):
        raise PopulationError(location + ": incorrect constant")
    if "enum" in schema and value not in schema["enum"]:
        raise PopulationError(location + ": unsupported value")
    types = schema.get("type", [])
    types = [types] if isinstance(types, str) else types
    actual = {dict: "object", list: "array", str: "string", int: "integer", float: "number",
              bool: "boolean", type(None): "null"}.get(type(value))
    if types and actual not in types and not (actual == "integer" and "number" in types):
        raise PopulationError(location + ": incorrect type")
    if type(value) is int:
        integer(value, location, schema.get("minimum", value), schema.get("maximum"))
    if type(value) is float and (not math.isfinite(value) or value < schema.get("minimum", value) or value > schema.get("maximum", value)):
        raise PopulationError(location + ": invalid number")
    if type(value) in (int, float) and "exclusiveMinimum" in schema and value <= schema["exclusiveMinimum"]:
        raise PopulationError(location + ": number must exceed its lower bound")
    if isinstance(value, str) and len(value) < schema.get("minLength", 0):
        raise PopulationError(location + ": empty string")
    if isinstance(value, str) and "pattern" in schema and re.search(schema["pattern"], value) is None:
        raise PopulationError(location + ": string does not match its required format")
    if isinstance(value, dict):
        props = schema.get("properties", {})
        fields(value, schema.get("required", []), props.keys() if schema.get("additionalProperties") is False
               else value.keys(), location)
        for key in value.keys() & props.keys():
            validate_schema(value[key], props[key], location + "." + key)
    if isinstance(value, list):
        if len(value) < schema.get("minItems", 0) or len(value) > schema.get("maxItems", len(value)):
            raise PopulationError(location + ": array length is outside its bounds")
        if schema.get("uniqueItems") and len({canonical(v) for v in value}) != len(value):
            raise PopulationError(location + ": duplicate values")
        for item in value:
            validate_schema(item, schema.get("items", {}), location + "[]")


@lru_cache(maxsize=None)
def planner_schema(kind):
    return parse_json(package_files(__package__).joinpath("schemas", f"donor-panel-{kind}.schema.json").read_text())


@lru_cache(maxsize=1)
def population_schema():
    return parse_json(package_files(__package__).joinpath("schemas", "population-record.schema.json").read_text())


@lru_cache(maxsize=1)
def study_schema():
    return parse_json(package_files(__package__).joinpath("schemas", "population-study.schema.json").read_text())


def validate_record(row):
    schema = next((s for s in population_schema()["oneOf"] if s["properties"]["record_type"]["const"] == row.get("record_type")), None)
    if schema is None:
        raise PopulationError("Unsupported population row type")
    validate_schema(row, schema, "population row")
    if row["record_type"] == "candidate_configuration":
        from .identity import bounding_interval, reverse_complement
        if identity("cfg", row["identity_payload"]) != row["configuration_id"] or row["identity_payload"]["target_path"] != row["target_path"]:
            raise PopulationError("Candidate identity does not match its sequence path")
        sequence = row["protospacer_guide_orientation"] + row["pam_guide_orientation"]
        forward = sequence if row["strand"] == "+" else reverse_complement(sequence)
        if forward != row["sequence_genomic_forward"] or bounding_interval(row["target_path"]) != row["bounding_interval"]:
            raise PopulationError("Candidate sequence orientation or bounding interval is inconsistent")
        for k, expected in (("reference_sha256", row["reference_id"]), ("contig", row["contig"]), ("strand", row["strand"]),
                            ("guide_definition_sha256", row["guide_definition_id"]), ("editor_definition_sha256", row["editor_definition_id"]),
                            ("protospacer_guide_orientation", row["protospacer_guide_orientation"]), ("pam_guide_orientation", row["pam_guide_orientation"])):
            if row["identity_payload"][k] != expected:
                raise PopulationError("Candidate identity payload differs from its record")
    if row["record_type"] == "sequence_configuration":
        if "".join(p["base"] for p in row["sequence_path"]) != row["sequence_genomic_forward"]:
            raise PopulationError("Materialized sequence and path disagree")
    return row


def validate_basis(mode, basis):
    if mode not in MODES or basis not in BASES:
        raise PopulationError("Unknown study mode or coverage basis")
    if basis == "retained_callset_model" and mode != "benchmark_simulation":
        raise PopulationError("Retained-callset modeling requires benchmark_simulation")
    if mode == "research_sample_planning" and basis != "callable_supported":
        raise PopulationError("Sample planning requires callable_supported")
    if mode == "synthetic_test" and basis != "synthetic_truth":
        raise PopulationError("Synthetic tests require synthetic_truth")


def validate_request(request):
    validate_schema(request, planner_schema("request"))
    validate_basis(request["mode"], request["coverage_basis"])
    units = unique(request["units"])
    donors = unique(request["donors"])
    unique(request["actions"])
    if request["coverage_basis"] == "retained_callset_model" and request["existing_usable_evidence"]:
        raise PopulationError("Conditional model planning cannot reduce demand using empirical evidence")
    for action in request["actions"]:
        if action["donor_id"] not in donors:
            raise PopulationError("Action donor does not resolve")
        covered, possible = set(action["covered_units"]), set(action["possible_only_units"])
        if covered & possible or (covered | possible) - units.keys():
            raise PopulationError("Invalid or overlapping action coverage links")
    for evidence in request["existing_usable_evidence"]:
        if evidence["unit_id"] not in units:
            raise PopulationError("Evidence unit does not resolve")
    return request


def _no_placeholders(value):
    if isinstance(value, str) and "REPLACE" in value:
        raise PopulationError("Unresolved REPLACE value; supply the named local source or policy")
    if isinstance(value, dict):
        for item in value.values():
            _no_placeholders(item)
    if isinstance(value, list):
        for item in value:
            _no_placeholders(item)


def intervals(rows, label):
    if not isinstance(rows, list) or not rows:
        raise PopulationError(label + " must list bounded intervals")
    for row in rows:
        fields(row, {"contig", "start", "end"}, label=label)
        nonempty(row["contig"], "contig")
        integer(row["start"], "interval start")
        integer(row["end"], "interval end", row["start"] + 1)
    ordered = sorted(rows, key=lambda r: (r["contig"], r["start"], r["end"]))
    for a, b in zip(ordered, ordered[1:]):
        if a["contig"] == b["contig"] and a["end"] > b["start"]:
            raise PopulationError(label + " intervals overlap")
    return ordered


def validate_config(c):
    _no_placeholders(c)
    validate_schema(c, study_schema(), "study config")
    fields(c, {"record_type", "schema_version", "study_id", "mode", "coverage_basis", "inputs",
               "input_sources", "reference_policy", "variant_policy", "genotype_policy", "scope",
               "search_profile", "assay_policy", "planning_policy", "benchmark_policy", "resource_limits"},
           label="study config")
    if c["record_type"] != "population_study_config" or c["schema_version"] != "1":
        raise PopulationError("Unsupported population study schema")
    nonempty(c["study_id"], "study_id")
    validate_basis(c["mode"], c["coverage_basis"])
    fields(c["inputs"], {"reference", "variants", "donors", "callable_intervals", "guides",
                         "intended_loci", "candidate_ledger", "measurement_linkage"}, label="inputs")
    for key, path in c["inputs"].items():
        if path is not None:
            nonempty(path, "input path")
        elif key in {"reference", "variants", "donors", "guides", "intended_loci"}:
            raise PopulationError("Missing required input " + key)
    if (c["inputs"]["candidate_ledger"] is None) != (c["inputs"]["measurement_linkage"] is None):
        raise PopulationError("Ledger and measurement linkage must be supplied together")
    fields(c["input_sources"], {k for k, v in c["inputs"].items() if v is not None}, label="input sources")
    for source in c["input_sources"].values():
        fields(source, {"source_id", "source_class", "release", "url", "use_notes", "sha256"},
               {"reference_call_contract"}, "source")
        for key in ("source_id", "source_class", "release", "use_notes"):
            nonempty(source[key], key)
        if source["url"] is not None:
            nonempty(source["url"], "source URL")
        if source["sha256"] is not None and not re.fullmatch(r"[0-9a-f]{64}", source["sha256"]):
            raise PopulationError("Expected SHA-256 must contain 64 lowercase hexadecimal characters")
    r = c["reference_policy"]
    fields(r, {"genome_build", "contig_alias_map", "permitted_contigs"}, label="reference policy")
    nonempty(r["genome_build"], "reference assembly")
    if c["coverage_basis"] != "synthetic_truth" and r["genome_build"] != "GRCh38":
        raise PopulationError("Real-data runs currently require GRCh38")
    if not isinstance(r["contig_alias_map"], dict) or not isinstance(r["permitted_contigs"], list):
        raise PopulationError("Explicit contig aliases and permitted contigs are required")
    if len(set(r["contig_alias_map"].values())) != len(r["contig_alias_map"]):
        raise PopulationError("Contig alias collision")
    if set(r["contig_alias_map"].values()) != set(r["permitted_contigs"]):
        raise PopulationError("Alias map must cover exactly the permitted contigs")
    if c["coverage_basis"] != "synthetic_truth" and any(
            not re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])", name) for name in r["permitted_contigs"]):
        raise PopulationError("Only diploid autosomes are supported")
    v = c["variant_policy"]
    fields(v, {"max_allele_length", "normalizer", "unsupported_footprint"}, label="variant policy")
    integer(v["max_allele_length"], "max allele length", 1, 50)
    nonempty(v["normalizer"], "normalizer and revision")
    if c["coverage_basis"] != "synthetic_truth" and not v["normalizer"].startswith("bcftools "):
        raise PopulationError("Public VCF inputs require a recorded bcftools normalizer version")
    if v["unsupported_footprint"] != "mark_affected_scope_unknown":
        raise PopulationError("Unsupported variant footprint policy")
    g = c["genotype_policy"]
    fields(g, {"ploidy", "accepted_site_filters", "accepted_sample_filters", "min_gq", "min_dp",
               "phase", "reference_evidence"}, label="genotype policy")
    if type(g["ploidy"]) is not int or g["ploidy"] != 2:
        raise PopulationError("Only diploid calls are supported")
    for key in ("accepted_site_filters", "accepted_sample_filters"):
        if not isinstance(g[key], list) or not g[key] or set(g[key]) - {"PASS", "."}:
            raise PopulationError("Declare accepted PASS and/or dot FILTER states")
    for key in ("min_gq", "min_dp"):
        if g[key] is not None:
            integer(g[key], key)
    if g["phase"] not in {"vcf_implicit_phase_set", "explicit_ps_only"}:
        raise PopulationError("Unsupported phase convention")
    if g["reference_evidence"] not in {"synthetic_complete", "callable_contract", "reference_imputed_scenario"}:
        raise PopulationError("Unknown reference evidence policy")
    if (g["reference_evidence"] == "synthetic_complete") != (c["coverage_basis"] == "synthetic_truth"):
        raise PopulationError("Synthetic reference evidence requires a synthetic truth study")
    if c["coverage_basis"] == "callable_supported" and (g["reference_evidence"] != "callable_contract"
                                                        or c["inputs"]["callable_intervals"] is None):
        raise PopulationError("Supported carriage requires a reference-callability contract")
    if c["inputs"]["callable_intervals"] is not None:
        nonempty(c["input_sources"]["callable_intervals"].get("reference_call_contract"),
                 "callability quality and nonreference representation contract")
    fields(c["scope"], {"requested_intervals", "materialized_intervals"}, label="scope")
    for key in c["scope"]:
        c["scope"][key] = intervals(c["scope"][key], key)
        if any(row["contig"] not in r["permitted_contigs"] for row in c["scope"][key]):
            raise PopulationError("Scope contig is not permitted")
    for row in c["scope"]["requested_intervals"]:
        if not any(m["contig"] == row["contig"] and m["start"] <= row["start"] and m["end"] >= row["end"]
                   for m in c["scope"]["materialized_intervals"]):
            raise PopulationError("Every requested interval must be materialized")
    s = c["search_profile"]
    fields(s, {"backend", "spacer_length", "pam_patterns", "max_mismatches", "max_dna_bulge",
               "max_rna_bulge", "retention"}, label="search profile")
    if s["backend"] != "regional_exhaustive_v1" or s["spacer_length"] != 20 \
            or s["max_dna_bulge"] != 0 or s["max_rna_bulge"] != 0 or s["retention"] != "all":
        raise PopulationError("Local scanner requires 20-base guides, zero alignment bulges and all rows")
    integer(s["max_mismatches"], "mismatch limit", 0, 20)
    integer(s["spacer_length"], "spacer length", 20, 20)
    integer(s["max_dna_bulge"], "DNA bulge limit", 0, 0)
    integer(s["max_rna_bulge"], "RNA bulge limit", 0, 0)
    if not isinstance(s["pam_patterns"], list) or not s["pam_patterns"] or any(
            not isinstance(p, str) or not re.fullmatch("[ACGTRYSWKMBDHVN]{3}", p) for p in s["pam_patterns"]):
        raise PopulationError("Declare three-base IUPAC PAM patterns")
    a = c["assay_policy"]
    fields(a, {"biological_context_id", "specimen_context_id", "measurement_class", "assay_id",
               "resolution_contract_id", "allele_resolved", "locus_resolved", "feasibility",
               "assay_cost", "control_costs_included", "extra_context_bases", "context", "guide_ids", "locus_scope"},
           {"max_detection_limit_fraction"}, label="assay policy")
    maximum = a.get("max_detection_limit_fraction")
    if maximum is not None and (type(maximum) not in (int, float) or not math.isfinite(maximum) or not 0 < maximum <= 1):
        raise PopulationError("Maximum detection limit must be a fraction greater than zero and at most one")
    fields(a["context"], {"cell_type", "delivery", "condition"}, label="biological context")
    for value in a["context"].values():
        nonempty(value, "biological context value")
        from ..pipeline.candidate_panel import _has_reported_text
        if not _has_reported_text(value):
            raise PopulationError("A declared assay context cannot use an unknown-value placeholder")
    if a["biological_context_id"] != identity("context", a["context"]):
        raise PopulationError("Biological context identity must hash its actual definition")
    if not isinstance(a["guide_ids"], list) or not a["guide_ids"] or len(set(a["guide_ids"])) != len(a["guide_ids"]):
        raise PopulationError("Assay requires unique declared guide IDs")
    for gid in a["guide_ids"]:
        nonempty(gid, "assay guide ID")
    fields(a["locus_scope"], {"mode", "loci"}, label="assay locus scope")
    scope = a["locus_scope"]
    if scope["mode"] == "all_admitted_loci_simulation":
        if c["mode"] == "research_sample_planning" or scope["loci"] != []:
            raise PopulationError("Physical planning requires explicit assay loci")
    elif scope["mode"] == "explicit_loci":
        if not isinstance(scope["loci"], list):
            raise PopulationError("Assay loci must be an array")
        for locus in scope["loci"]:
            fields(locus, {"contig", "start", "end", "strand"}, label="assay locus")
            integer(locus["start"], "assay locus start")
            integer(locus["end"], "assay locus end", locus["start"])
            if locus["strand"] not in {"+", "-"}:
                raise PopulationError("Assay locus requires a strand")
    else:
        raise PopulationError("Unknown assay locus scope")
    if a["resolution_contract_id"] != resolution_id(a):
        raise PopulationError("Assay resolution identity must hash its capabilities and scope")
    for key in ("biological_context_id", "specimen_context_id", "measurement_class", "assay_id", "resolution_contract_id"):
        nonempty(a[key], key)
    for key in ("allele_resolved", "locus_resolved", "control_costs_included"):
        if type(a[key]) is not bool:
            raise PopulationError("Assay capability must be boolean")
    if a["measurement_class"] != "editing" or a["feasibility"] not in {"declared_feasible", "unverified", "infeasible"}:
        raise PopulationError("Unsupported assay measurement or feasibility")
    integer(a["assay_cost"], "assay cost")
    integer(a["extra_context_bases"], "extra assay context", 0, 0)
    if a["extra_context_bases"] != 0:
        raise PopulationError("Additional assay sequence context is not yet supported")
    p = c["planning_policy"]
    fields(p, {"cost_unit", "max_cost", "max_donors", "max_assays", "fixed_overhead",
               "default_weight", "default_required_donors", "algorithm"}, label="planning policy")
    nonempty(p["cost_unit"], "cost unit")
    for key in ("max_cost", "fixed_overhead"):
        integer(p[key], key)
    if p["fixed_overhead"] > p["max_cost"]:
        raise PopulationError("Study overhead exceeds total budget")
    if not a["control_costs_included"] and p["fixed_overhead"] == 0:
        raise PopulationError("Declare control costs in assay packages or fixed overhead")
    for key in ("max_donors", "max_assays"):
        if p[key] is not None:
            integer(p[key], key)
    for key in ("default_weight", "default_required_donors"):
        integer(p[key], key, 1)
    if p["algorithm"] not in ALGORITHMS:
        raise PopulationError("Unsupported planning algorithm")
    validate_benchmark_policy(c["benchmark_policy"])
    fields(c["resource_limits"], LIMITS, label="resource limits")
    for key, value in c["resource_limits"].items():
        integer(value, key, 1)
    return c


def validate_benchmark_policy(b):
    fields(b, {"selectors", "seeds", "budget_grid", "split", "relatedness_policy", "endpoints"}, label="benchmark policy")
    if not isinstance(b["selectors"], list) or not b["selectors"] or any(not isinstance(s, str) or s not in ALGORITHMS for s in b["selectors"]) \
            or len(b["selectors"]) != len(set(b["selectors"])):
        raise PopulationError("Benchmark selectors must be a nonempty, unique list of supported algorithms")
    for key in ("seeds", "budget_grid"):
        if not isinstance(b[key], list) or not b[key] or len(b[key]) != len(set(b[key])):
            raise PopulationError("Benchmark values must be unique nonempty arrays")
        for value in b[key]:
            integer(value, key)
    if b["split"] != "transductive":
        raise PopulationError("External generalization is not implemented; freeze a separate family-aware protocol")
    nonempty(b["relatedness_policy"], "relatedness policy")
    if b["endpoints"] != ["nomination_changes", "planned_configuration_coverage"]:
        raise PopulationError("Computed nominations cannot supply empirical sensitivity endpoints")


def resolution_id(assay):
    payload = {k: assay[k] for k in ("assay_id", "biological_context_id", "specimen_context_id",
        "measurement_class", "allele_resolved", "locus_resolved", "feasibility", "guide_ids", "locus_scope", "extra_context_bases")}
    if assay.get("max_detection_limit_fraction") is not None:
        payload["max_detection_limit_fraction"] = assay["max_detection_limit_fraction"]
    return identity("resolution", payload)
