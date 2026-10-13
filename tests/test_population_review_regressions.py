"""Scientific review reproductions through reconstruction, evidence and assay planning."""
from copy import deepcopy
from pathlib import Path
import shutil

import pytest

from base_edit_safety.pipeline import candidate_panel as cp, normalize_nominations as nn, population as pipeline
from base_edit_safety.pipeline.donor_panel import execute
from base_edit_safety.population.common import checksum, json_text, read_json, read_jsonl
from base_edit_safety.population.provenance import resolve_study
from base_edit_safety.population.reporting import tsv
from base_edit_safety.population.variants import DONOR_COLUMNS, donors_from_tsv

ROOT = Path(__file__).resolve().parents[1]
GUIDE = "ACGTTGCAATGCTAGCTAAC"


def study_inputs(tmp_path):
    source = tmp_path/"inputs"
    shutil.copytree(ROOT/"examples/population/synthetic", source)
    return source, read_json(source/"study-config.json")


def workflow(source, config, root):
    (source/"study-config.json").write_text(json_text(config))
    study = resolve_study(source/"study-config.json")
    pipeline.build(study, root/"population")
    pipeline.search(study, root/"population", root/"search")
    pipeline.reconcile_stage(study, root/"search", root/"evidence")
    request = pipeline.actions_stage(study, root/"population", root/"search", root/"evidence", root/"panel-request.json")
    result = execute(root/"panel-request.json", root/"planning-manifest.json", "greedy_cost", root/"plan")
    return request, result


@pytest.mark.parametrize("alt_gt", ["1|1", "2|2", "1|2"])
def test_mixed_symbolic_row_cannot_create_an_assay_on_another_contig(tmp_path, alt_gt):
    source, config = study_inputs(tmp_path)
    contigs = ["toy_gain", "toy_mirror", "chrSYNcontrol"]
    (source/"reference.fa").write_text("".join(
        ">"+c+"\n"+GUIDE+("TGG" if c == "chrSYNcontrol" else "TGC")+"\n" for c in contigs))
    config["reference_policy"].update(permitted_contigs=contigs, contig_alias_map={c: c for c in contigs})
    config["scope"] = {k: [{"contig": c, "start": 0, "end": 23} for c in contigs]
                       for k in ("requested_intervals", "materialized_intervals")}
    (source/"variants.vcf").write_text(
        "##fileformat=VCFv4.3\n##reference=synthetic_v3_search_1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYN_REF\tSYN_ALT\tSYN_HET\tSYN_MISSING\n"
        "toy_gain\t23\tSYN_mixed\tC\tG,<DEL>\t.\tPASS\t.\tGT\t0|0\t"+alt_gt+"\t0/1\t./.\n")
    root = tmp_path/"run"
    request, result = workflow(source, config, root)
    candidates = read_jsonl(root/"search/candidate-configurations.jsonl")
    assert {c["contig"] for c in candidates if c["role"] == "candidate_offtarget"} == {"toy_gain"}
    assert len(request["units"]) == len(result["selected_actions"]) == result["objective_numerator"] == 1
    materials = read_jsonl(root/"population/sequence-configurations.jsonl")
    mirror = [m for m in materials if m["contig"] == "toy_mirror" and m["donor_id"] == "SYN_HET"]
    assert all(m["sequence_genomic_forward"] == GUIDE+"TGC" and not m["variant_intervals"] for m in mirror)
    if alt_gt != "1|1":
        # Carried unbounded alleles still mask every requested interval conservatively.
        assert all(m["excluded_intervals"] for m in materials if m["donor_id"] == "SYN_ALT")


def write_source_ledger(source, candidate, guide, config, measurement_changes=None):
    """Build a schema-3 ledger from mock BED6 and measurement inputs."""
    build = config["reference_policy"]["genome_build"]
    source_guide = deepcopy(guide)
    (source/"source-guides.tsv").write_text(tsv(sorted(cp.GUIDE_FIELDS), [source_guide]))
    locus = {"guide_id": guide["guide_id"], "genome_build": build, "chrom": "SYN_source_control",
             "start_0based": 0, "end_0based_exclusive": 23, "strand": "+", "locus_name": "source_control"}
    (source/"source-intended.tsv").write_text(tsv(sorted(cp.INTENDED_LOCUS_FIELDS), [locus]))
    profile = read_json(ROOT/"examples/candidate-panel/search-profile.json")
    profile.update(name="synthetic_source", max_mismatches=20, max_dna_bulge=0, max_rna_bulge=0,
                   variant_aware=False, reference_id=build, variant_source=None)
    (source/"source-profile.json").write_text(json_text(profile))
    measurement = {"editor": guide["editor"], "guide_id": guide["guide_id"], "genome_build": build,
        "chrom": candidate["contig"], "start_0based": 0, "end_0based_exclusive": 23, "strand": "+",
        "allele": "alt", "allele_sequence": candidate["sequence_genomic_forward"], "allele_resolution": "resolved",
        "allele_evidence": "synthetic resolved allele", "qc_status": "passed", "treated_group": "synthetic treated",
        "control_group": "synthetic control", "biological_context": "synthetic cells", "assay": "synthetic assay",
        "measurement_kind": "editing", "status": "not_detected", "limit_of_detection_value": 0.1,
        "limit_of_detection_unit": "percent", "source_record_id": "SYN_record", "source_version_or_citation": "synthetic"}
    measurement.update(measurement_changes or {})
    (source/"source-measurements.tsv").write_text(tsv(sorted(cp.MEASUREMENT_FIELDS), [measurement]))
    fasta, index, offset = "", "", 0
    for contig, sequence in [(candidate["contig"], candidate["sequence_genomic_forward"]),
                             ("SYN_source_control", guide["spacer_sequence"]+"TGG")]:
        header = ">"+contig+"\n"
        fasta += header+sequence+"\n"
        index += f"{contig}\t23\t{offset+len(header)}\t23\t24\n"
        offset = len(fasta)
    (source/"source-reference.fa").write_text(fasta)
    (source/"source-reference.fa.fai").write_text(index)
    (source/"source.bed").write_text(candidate["contig"]+"\t0\t23\tSYN_nomination\t0\t+\n")
    assert nn.main(["experimental", "--input", str(source/"source.bed"), "--guides", str(source/"source-guides.tsv"),
        "--intended-loci", str(source/"source-intended.tsv"), "--source-version", "synthetic",
        "--source-command", "not applicable: synthetic fixture", "--reference-fasta", str(source/"source-reference.fa"),
        "--output", str(source/"source-nominations.tsv"), "--provenance-output", str(source/"source-provenance.json"),
        "--format", "bed6", "--guide-id", guide["guide_id"], "--reference-id", build,
        "--source-method", "synthetic experimental assay"]) == 0
    ledger = cp.build_ledger(source/"source-profile.json", source/"source-guides.tsv", source/"source-intended.tsv",
        [source/"source-nominations.tsv"], source/"source-measurements.tsv", reference_fasta_path=source/"source-reference.fa",
        source_provenance_paths=[source/"source-provenance.json"])
    (source/"ledger.json").write_text(json_text(ledger))
    site = next(s for s in ledger["candidates"] if s["measurements"])
    link = {k: measurement[k] for k in ("source_version_or_citation", "assay", "source_record_id")}
    link.update(site_id=site["site_id"], source_person_id="SYN_historical", configuration_id=candidate["configuration_id"],
        target_path=candidate["target_path"], locus_resolved=True, source_assertion="synthetic truth",
        biological_context_id=config["assay_policy"]["biological_context_id"], context=config["assay_policy"]["context"])
    (source/"linkage.json").write_text(json_text({"record_type": "measurement_linkage", "schema_version": "1",
        "ledger_sha256": checksum(source/"ledger.json"), "people": [{"source_person_id": "SYN_historical",
        "canonical_donor_id": "SYN_historical", "identity_status": "resolved", "source_assertion": "synthetic person"}],
        "links": [link]}))


@pytest.mark.parametrize("wrong_spacer", [False, True])
def test_valid_source_ledger_cannot_transfer_credit_between_reused_guide_ids(tmp_path, wrong_spacer):
    source, config = study_inputs(tmp_path)
    workflow(source, config, tmp_path/"initial")
    candidate = next(c for c in read_jsonl(tmp_path/"initial/search/candidate-configurations.jsonl") if c["contig"] == "toy_gain")
    guide = cp.load_guides(source/"guides.tsv")[0]
    if wrong_spacer:
        guide["spacer_sequence"] = "C"*20
    write_source_ledger(source, candidate, guide, config)
    config["inputs"].update(candidate_ledger="ledger.json", measurement_linkage="linkage.json")
    for key in ("candidate_ledger", "measurement_linkage"):
        config["input_sources"][key] = {**config["input_sources"]["guides"], "source_id": "synthetic_"+key}
    root = tmp_path/"with_evidence"
    request, result = workflow(source, config, root)
    links = read_jsonl(root/"evidence/evidence-links.jsonl")
    assert len(links) == 1 and links[0]["usable"] is not wrong_spacer
    assert len(request["existing_usable_evidence"]) == (0 if wrong_spacer else 1)
    assert len(result["selected_actions"]) == (3 if wrong_spacer else 2)
    if wrong_spacer:
        assert "source_guide_definition_mismatch" in links[0]["reasons"]


@pytest.mark.parametrize("value,unit,maximum,usable", [
    (200, "percent", None, False), (0.1, "arbitrary dimension", None, False),
    (1.1, "fraction", None, False), (0.1, "percent", None, True),
    (0.001, "fraction", 0.001, True), (0.1, "percent", 0.0005, False),
])
def test_invalid_or_insensitive_non_detection_never_reduces_assay_demand(tmp_path, value, unit, maximum, usable):
    from base_edit_safety.population.contracts import resolution_id
    source, config = study_inputs(tmp_path)
    # Keep all three requirements affordable when sensitivity changes action tie order.
    config["planning_policy"]["max_cost"] = 8
    workflow(source, config, tmp_path/"initial")
    candidate = next(c for c in read_jsonl(tmp_path/"initial/search/candidate-configurations.jsonl") if c["contig"] == "toy_gain")
    guide = cp.load_guides(source/"guides.tsv")[0]
    if maximum is not None:
        config["assay_policy"]["max_detection_limit_fraction"] = maximum
        config["assay_policy"]["resolution_contract_id"] = resolution_id(config["assay_policy"])
    write_source_ledger(source, candidate, guide, config,
                        {"limit_of_detection_value": value, "limit_of_detection_unit": unit})
    config["inputs"].update(candidate_ledger="ledger.json", measurement_linkage="linkage.json")
    for key in ("candidate_ledger", "measurement_linkage"):
        config["input_sources"][key] = {**config["input_sources"]["guides"], "source_id": "synthetic_"+key}
    root = tmp_path/"with_evidence"
    request, result = workflow(source, config, root)
    links = read_jsonl(root/"evidence/evidence-links.jsonl")
    assert len(links) == 1 and links[0]["usable"] is usable
    assert links[0]["source_measurement"]["limit_of_detection"] == {"value": value, "unit": unit}
    assert len(request["existing_usable_evidence"]) == int(usable)
    assert len(result["selected_actions"]) == (2 if usable else 3)


def test_renamed_source_guide_has_an_explicit_unmatched_diagnostic(tmp_path):
    source, config = study_inputs(tmp_path)
    workflow(source, config, tmp_path/"initial")
    candidate = next(c for c in read_jsonl(tmp_path/"initial/search/candidate-configurations.jsonl") if c["contig"] == "toy_gain")
    guide = cp.load_guides(source/"guides.tsv")[0]
    write_source_ledger(source, candidate, {**guide, "guide_id": "SYN_historical_label"}, config)
    config["inputs"].update(candidate_ledger="ledger.json", measurement_linkage="linkage.json")
    for key in ("candidate_ledger", "measurement_linkage"):
        config["input_sources"][key] = {**config["input_sources"]["guides"], "source_id": "synthetic_"+key}
    root = tmp_path/"renamed"
    request, result = workflow(source, config, root)
    assert not read_jsonl(root/"evidence/evidence-links.jsonl")
    unmatched = read_json(root/"evidence/unmatched-measurements.json")["rows"]
    assert len(unmatched) == 1
    assert "source_guide_id_not_in_study" in unmatched[0]["reasons"]
    assert unmatched[0]["current_guide_ids_with_same_definitions"] == [guide["guide_id"]]
    assert unmatched[0]["source_guide"]["guide_id"] == "SYN_historical_label"
    assert not request["existing_usable_evidence"] and len(result["selected_actions"]) == 3


@pytest.mark.parametrize("percent,fraction,usable", [
    (0.07, 0.0007, True), (0.06999999999999999, 0.0006999999999999999, True),
    (0.07000000000000002, 0.0007000000000000002, False),
])
def test_equivalent_decimal_detection_limits_produce_the_same_assay_plan(tmp_path, percent, fraction, usable):
    from base_edit_safety.population.contracts import resolution_id
    plans = []
    for value, unit in [(percent, "percent"), (fraction, "fraction")]:
        folder = tmp_path/unit
        source, config = study_inputs(folder)
        config["planning_policy"]["max_cost"] = 8
        config["assay_policy"]["max_detection_limit_fraction"] = 0.0007
        config["assay_policy"]["resolution_contract_id"] = resolution_id(config["assay_policy"])
        workflow(source, config, folder/"initial")
        candidate = next(c for c in read_jsonl(folder/"initial/search/candidate-configurations.jsonl") if c["contig"] == "toy_gain")
        guide = cp.load_guides(source/"guides.tsv")[0]
        write_source_ledger(source, candidate, guide, config,
                            {"limit_of_detection_value": value, "limit_of_detection_unit": unit})
        config["inputs"].update(candidate_ledger="ledger.json", measurement_linkage="linkage.json")
        for key in ("candidate_ledger", "measurement_linkage"):
            config["input_sources"][key] = {**config["input_sources"]["guides"], "source_id": "synthetic_"+key}
        request, plan = workflow(source, config, folder/"with-evidence")
        assert len(request["existing_usable_evidence"]) == int(usable)
        assert len(plan["selected_actions"]) == (2 if usable else 3)
        plans.append({key: plan[key] for key in ("selected_actions", "selected_donors", "total_cost", "unit_coverage", "objective_numerator")})
    assert plans[0] == plans[1]


def test_multiple_specimens_do_not_multiply_genotypes_or_donor_credit(tmp_path):
    source, config = study_inputs(tmp_path)
    donors = donors_from_tsv(source/"donors.tsv", "synthetic_truth")
    donor = next(d for d in donors if d["donor_id"] == "SYN_ALT")
    specimens = [donor, {**donor, "specimen_id": "SYN_second"},
                 {**donor, "specimen_id": "SYN_wrong_context", "specimen_context_id": "other_cells"},
                 {**donor, "specimen_id": "SYN_no_use", "use_permitted": False}]
    rows = [{k: str(v).lower() if isinstance(v, bool) else v for k, v in d.items()} for d in specimens]
    (source/"donors.tsv").write_text(tsv(DONOR_COLUMNS, rows))
    config["resource_limits"]["max_donors"] = 1
    config["planning_policy"]["default_required_donors"] = 2
    root = tmp_path/"multi"
    request, result = workflow(source, config, root)
    assert len(request["donors"]) == 1
    assert len(read_jsonl(root/"population/donor-inventory.jsonl")) == 4
    assert all(v["frequency"]["AN"] == 2 for v in read_jsonl(root/"population/normalized-variants.jsonl"))
    assert len(read_jsonl(root/"population/sequence-configurations.jsonl")) == 12
    assert len(read_jsonl(root/"search/donor-carriage.jsonl")) == 4
    assert len(result["selected_actions"]) == 2 and result["total_cost"] == 3
    assert not result["fully_covered_units"]
    catalog = {a["id"]: a for a in read_jsonl(root/"assay-actions.jsonl")}
    assert all(catalog[a]["specimen_id"] not in {"SYN_wrong_context", "SYN_no_use"} for a in result["selected_actions"])
    assert all(len(u["planned_donors"]) <= 1 for u in result["unit_coverage"])
    import csv
    with (root/"plan/donor-site-plan.tsv").open() as handle:
        plan_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["specimen_id"] == catalog[row["action_id"]]["specimen_id"] for row in plan_rows)


@pytest.mark.parametrize("field,value", [("genotype_sample_id", "other_sample"), ("activation_cost", 2), ("family_id", "other_family")])
def test_conflicting_metadata_across_specimens_is_rejected(tmp_path, field, value):
    from base_edit_safety.population.common import PopulationError
    source, _ = study_inputs(tmp_path)
    donor = donors_from_tsv(source/"donors.tsv", "synthetic_truth")[0]
    rows = [donor, {**donor, "specimen_id": "SYN_second", field: value}]
    rows = [{k: str(v).lower() if isinstance(v, bool) else v for k, v in d.items()} for d in rows]
    (source/"donors.tsv").write_text(tsv(DONOR_COLUMNS, rows))
    with pytest.raises(PopulationError, match="disagree"):
        donors_from_tsv(source/"donors.tsv", "synthetic_truth")
