"""Independent sequence truths for phase, missing calls, strands and indel paths."""
from copy import deepcopy
import hashlib
from pathlib import Path

import pytest

from base_edit_safety.population.common import PopulationError, ResourceLimit, json_text, read_json
from base_edit_safety.population.haplotypes import build_materializations, local_status
from base_edit_safety.population.identity import canonical_path, candidate_identity, reverse_complement
from base_edit_safety.population.search import scan, nomination_changes
from base_edit_safety.population.variants import Reference, load_vcf, parse_call

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT/"tests/fixtures/population"
GUIDE = "ACGTTGCAATGCTAGCTAAC"


def experiment(tmp_path, sequence, variants, samples, *, start=0, guide=GUIDE, editor_class="ABE", config_change=None):
    """Write literal input alleles; expected sequences are supplied by each test."""
    config = read_json(ROOT/"examples/population/synthetic/study-config.json")
    interval = {"contig": "chr2", "start": start, "end": start+len(sequence)}
    config["reference_policy"] = {"genome_build": "GRCh38", "contig_alias_map": {"chr2": "chr2"}, "permitted_contigs": ["chr2"]}
    config["scope"] = {"requested_intervals": [interval], "materialized_intervals": [interval]}
    if config_change:
        config_change(config)
    reference_path = tmp_path/"reference.json"
    reference_path.write_text(json_text({"record_type": "reference_windows", "schema_version": "1", "genome_build": "GRCh38",
        "windows": [{**interval, "sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}]}))
    reference = Reference(reference_path, config["reference_policy"], True)
    donors = [{"donor_id": s, "genotype_sample_id": s} for s in samples]
    header = "##fileformat=VCFv4.3\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples)+"\n"
    text = header
    for i, v in enumerate(variants):
        text += "\t".join(["chr2", str(v["pos1"]), "synthetic_"+str(i), v["ref"], v["alt"], ".", v.get("filter", "PASS"),
                           v.get("info", "."), v.get("format", "GT"), *v["calls"]])+"\n"
    vcf = tmp_path/"variants.vcf"; vcf.write_text(text)
    records, audit, ingestion = load_vcf(vcf, config, reference, donors, "test:run")
    materials, counts = build_materializations(config, reference, records, donors, "test:run")
    g = {"guide_id": "synthetic_guide", "spacer_sequence": guide, "editor": "SpCas9" if editor_class == "nuclease" else "synthetic_ABE", "editor_class": editor_class, "pam_class": "NGG",
         "editing_window_start_1based": 4, "editing_window_end_1based_inclusive": 8}
    candidates, hits, carriage, search_audit = scan(config, reference, materials, [g], [], donors, [], "test:run")
    return locals()


def donor_hits(e, donor):
    return [h for h in e["hits"] if h["donor_id"] == donor]


@pytest.mark.parametrize("gt", ["./.", "0/.", "./1", ".", "0", "0/1/1"])
def test_missing_and_unsupported_ploidy_never_become_reference(tmp_path, gt):
    e = experiment(tmp_path, GUIDE+"TGG", [{"pos1": 23, "ref": "G", "alt": "C", "calls": [gt]}], ["SYN"])
    assert not donor_hits(e, "SYN")
    assert all(r["model_carriage"] == r["evidence_supported_carriage"] == "unknown" for r in e["carriage"])
    assert e["records"][0]["frequency"]["AN"] == 0
    assert e["records"][0]["frequency"]["AF"] == [None]


def test_new_pam_without_candidate_list_and_destroyed_pam(tmp_path):
    e = experiment(tmp_path, GUIDE+"TGC", [{"pos1": 23, "ref": "C", "alt": "G", "calls": ["0|0", "1|1", "0/1", "./."]}], ["REF", "ALT", "HET", "MISS"])
    assert not donor_hits(e, None)
    assert len(donor_hits(e, "ALT")) == 2 and len(donor_hits(e, "HET")) == 1
    changes = nomination_changes(e["config"], e["candidates"], e["materials"], e["donors"], [e["g"]], e["carriage"], "test:run")
    assert {r["change"] for r in changes if r["donor_id"] == "ALT"} == {"pam_gained"}
    assert {r["change"] for r in changes if r["donor_id"] == "MISS"} == {"unpaired_unevaluable"}
    assert e["records"][0]["frequency"] == {"AC": [3], "AN": 6, "AF": [0.5], "evaluable_donors": 3, "subset": "declared_donor_catalog"}
    e = experiment(tmp_path, GUIDE+"TGG", [{"pos1": 23, "ref": "G", "alt": "C", "calls": ["1|1"]}], ["ALT"])
    assert not donor_hits(e, "ALT") and donor_hits(e, None)
    changes = nomination_changes(e["config"], e["candidates"], e["materials"], e["donors"], [e["g"]], e["carriage"], "test:run")
    assert changes[0]["change"] == "pam_lost" and changes[0]["paired_evaluable"]


def test_multiallelic_original_indexes(tmp_path):
    e = experiment(tmp_path, GUIDE+"TGC", [{"pos1": 23, "ref": "C", "alt": "G,A", "calls": ["1|2"]}], ["SYN"])
    seqs = [m["sequence_genomic_forward"] for m in e["materials"] if m["donor_id"]]
    assert seqs == [GUIDE+"TGG", GUIDE+"TGA"]
    assert e["audit"][0]["selected_indexes"] == [1, 2]
    assert e["records"][0]["frequency"]["AC"] == [1, 1]


@pytest.mark.parametrize("gt,alleles", [("1|2", ["T", "TTTG"]), ("2|1", ["TTTG", "T"])])
def test_jointly_normalized_deletions_preserve_both_original_alleles(tmp_path, gt, alleles):
    # The shorter deletion can move within the repeat on its own. The other
    # ALT prevents that move when normalizing the complete unsplit record.
    e = experiment(tmp_path, "C"+"TTTGTTTG"+GUIDE+"TGG",
                   [{"pos1": 2, "ref": "TTTGTTTG", "alt": "T,TTTG", "calls": [gt]}], ["SYN"])
    copies = [m for m in e["materials"] if m["donor_id"] == "SYN"]
    assert [m["sequence_genomic_forward"] for m in copies] == ["C"+a+GUIDE+"TGG" for a in alleles]
    assert e["audit"][0]["selected_indexes"] == list(map(int, gt.split("|")))
    assert e["records"][0]["original_alts"] == ["T", "TTTG"]
    assert all(not m["excluded_intervals"] for m in copies)
    assert len(donor_hits(e, "SYN")) == 2


@pytest.mark.parametrize("ref,alt", [("TT", "T,TTT"), ("CAG", "CA,CAGG")])
def test_multiallelic_indel_still_requires_joint_normalization(tmp_path, ref, alt):
    with pytest.raises(PopulationError, match="jointly normalized"):
        experiment(tmp_path, "C"+ref+GUIDE+"TGG",
                   [{"pos1": 2, "ref": ref, "alt": alt, "calls": ["1|2"]}], ["SYN"])


def test_cis_trans_truth_and_unknown_phase(tmp_path):
    truth = read_json(FIXTURES/"phase-truth.json")
    variants = [{**v, "calls": [d["genotypes"][i] for d in truth["donors"]]} for i, v in enumerate(truth["variants"])]
    e = experiment(tmp_path, truth["reference_sequence"], variants, [d["id"] for d in truth["donors"]])
    for donor in truth["donors"][:2]:
        observed = [m["sequence_genomic_forward"] for m in e["materials"] if m["donor_id"] == donor["id"]]
        assert observed == donor["truth_haplotypes"]
        assert len(donor_hits(e, donor["id"])) == donor["exact_match_copies"]
    assert not donor_hits(e, "SYN_UNPHASED")
    assert all(r["model_carriage"] == "unknown" for r in e["carriage"] if r["donor_id"] == "SYN_UNPHASED")


@pytest.mark.parametrize("calls,format,expected", [(["0|1:1", "0|1:2"], "GT:PS", False),
    (["0/1:1", "0/1:1"], "GT:PS", False), (["0|1:1", "0|1:1"], "GT:PS", True),
    (["0|1", "0|1"], "GT", True), (["0|1:1", "0|1:."], "GT:PS", False)])
def test_phase_blocks_and_ps_semantics(tmp_path, calls, format, expected):
    truth = read_json(FIXTURES/"phase-truth.json")
    variants = [{**v, "format": format, "calls": [calls[i]]} for i, v in enumerate(truth["variants"])]
    e = experiment(tmp_path, truth["reference_sequence"], variants, ["SYN"])
    assert bool(donor_hits(e, "SYN")) == expected
    if not expected:
        material = next(m for m in e["materials"] if m["donor_id"])
        assert local_status(material, 10, 23)[0]


def test_phase_namespace_includes_source_donor_and_contig():
    policy = read_json(ROOT/"examples/population/synthetic/study-config.json")["genotype_policy"]
    keys = [parse_call("0|1:9", ["GT", "PS"], "PASS", policy, *args)["phase_key"]
            for args in [("s1", "d1", "chr1"), ("s2", "d1", "chr1"), ("s1", "d2", "chr1"), ("s1", "d1", "chr2")]]
    assert len({tuple(k) for k in keys}) == 4


@pytest.mark.parametrize("variant", [{"filter": "."}, {"format": "GT:GQ", "calls": ["0|0:."]}, {"format": "GT:FT", "calls": ["0|0:LowQual"]}])
def test_quality_never_invents_pass(tmp_path, variant):
    v = {"pos1": 23, "ref": "G", "alt": "C", "calls": ["0|0"], **variant}
    e = experiment(tmp_path, GUIDE+"TGG", [v], ["SYN"], config_change=lambda c: c["genotype_policy"].update(min_gq=10) if "GQ" in v.get("format", "") else None)
    assert e["audit"][0]["call_status"] == "filtered"
    assert not donor_hits(e, "SYN")


def test_unsupported_symbolic_allele_masks_only_affected_donors(tmp_path):
    e = experiment(tmp_path, GUIDE+"TGG", [{"pos1": 1, "ref": "A", "alt": "<DEL>", "info": "END=23", "calls": ["0|0", "0|1", "./."]}], ["REF", "ALT", "MISS"])
    assert len(donor_hits(e, "REF")) == 2
    assert not donor_hits(e, "ALT") and not donor_hits(e, "MISS")


def test_conditional_reference_fill_cannot_fill_retained_missing_calls(tmp_path):
    def conditional(c):
        c.update(mode="benchmark_simulation", coverage_basis="retained_callset_model")
        c["genotype_policy"]["reference_evidence"] = "reference_imputed_scenario"
    e = experiment(tmp_path, GUIDE+"TGG", [{"pos1": 23, "ref": "G", "alt": "C", "calls": ["0|0", "./."]}], ["REF", "MISS"], config_change=conditional)
    rows = {r["donor_id"]: r for r in e["carriage"]}
    assert rows["REF"]["model_carriage"] == "present" and rows["REF"]["evidence_supported_carriage"] == "unknown"
    assert rows["MISS"]["model_carriage"] == "unknown"


def test_paired_scope_counts_windows_even_when_neither_sequence_has_a_hit(tmp_path):
    e = experiment(tmp_path, "T"*24, [{"pos1": 1, "ref": "T", "alt": "A", "calls": ["0|1", "./."]}], ["OK", "MISSING"])
    assert not e["candidates"]
    rows = {r["donor_id"]: r for r in e["search_audit"]["paired_window_scope"]["rows"]}
    assert rows["OK"]["union_oriented_windows"] == rows["OK"]["paired_evaluable_oriented_windows"] == 4
    assert rows["OK"]["paired_reference_start_intervals"] == [{"start": 0, "end": 2}]
    assert rows["MISSING"]["paired_evaluable_oriented_windows"] == 2
    assert rows["MISSING"]["unpaired_reasons"] == {"unresolved_variant_footprint": 2}
    assert rows["MISSING"]["paired_reference_start_intervals"] == [{"start": 1, "end": 2}]


def test_paired_scope_distinguishes_callability_and_indel_paths(tmp_path):
    e = experiment(tmp_path, "T"*20, [{"pos1": 10, "ref": "T", "alt": "T"+GUIDE+"TGG", "calls": ["0/1"]}], ["SYN"])
    row = e["search_audit"]["paired_window_scope"]["rows"][0]
    assert row["reference_oriented_windows"] == row["paired_evaluable_oriented_windows"] == 0
    assert row["union_oriented_windows"] == row["unpaired_oriented_windows"] == 42
    assert row["unpaired_reasons"] == {"variant_only_indel_path": 42}
    changes = nomination_changes(e["config"], e["candidates"], e["materials"], e["donors"], [e["g"]], e["carriage"], "test:run")
    assert changes and all(not row["paired_evaluable"] and not row["pam_lost"] for row in changes)
    e = experiment(tmp_path, "T"*24, [], ["SYN"])
    e["config"].update(mode="research_sample_planning", coverage_basis="callable_supported")
    calls = [{"donor_id": "SYN", "contig": "chr2", "start_0based": 0, "end_0based_exclusive": 23}]
    audit = scan(e["config"], e["reference"], e["materials"], [e["g"]], [], e["donors"], calls, "test:run")[3]
    row = audit["paired_window_scope"]["rows"][0]
    assert row["paired_evaluable_oriented_windows"] == 2
    assert row["unpaired_reasons"] == {"reference_callability_unavailable": 2}


@pytest.mark.parametrize("genotype", ["0|1", "1|1"])
@pytest.mark.parametrize("kind", ["insertion", "deletion"])
def test_reference_nomination_and_window_audit_agree_across_indels(tmp_path, genotype, kind):
    reference = GUIDE+"TGGAC"
    ref, alt = ("T", "TC") if kind == "insertion" else (reference[9:11], "T")
    e = experiment(tmp_path, reference, [{"pos1": 10, "ref": ref, "alt": alt, "calls": [genotype]}], ["SYN"])
    expected_alt = reference[:10]+"C"+reference[10:] if kind == "insertion" else reference[:10]+reference[11:]
    copies = [m["sequence_genomic_forward"] for m in e["materials"] if m["donor_id"]]
    assert copies == ([reference, expected_alt] if genotype == "0|1" else [expected_alt, expected_alt])
    candidate = next(c for c in e["candidates"] if c["strand"] == "+" and c["target_path"] == [{"kind": "reference", "pos0": i} for i in range(23)])
    changes = nomination_changes(e["config"], e["candidates"], e["materials"], e["donors"], [e["g"]], e["carriage"], "test:run")
    change = next(c for c in changes if c["configuration_id"] == candidate["configuration_id"])
    assert not change["paired_evaluable"] and not change["pam_lost"]
    assert change["pairing_reason"] == "donor_indel_path_not_shared"
    assert e["search_audit"]["paired_window_scope"]["rows"][0]["paired_evaluable_oriented_windows"] == 0


@pytest.mark.parametrize("case_index", [0, 1])
def test_independent_indel_paths_and_insertion_only_locator(tmp_path, case_index):
    truth = read_json(FIXTURES/"indel-truth.json")
    case = truth["cases"][case_index]
    e = experiment(tmp_path, case["reference_sequence"], [{**case["variant"], "calls": ["1|1"]}], ["SYN"])
    material = next(m for m in e["materials"] if m["donor_id"])
    assert material["sequence_genomic_forward"] == case["truth_reconstructed_sequence"]
    c = next(c for c in e["candidates"] if c["configuration_id"] == donor_hits(e, "SYN")[0]["configuration_id"])
    assert c["target_path"] == case["target_path"]
    assert c["bounding_interval"] == case["bounding_interval"]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("at_boundary", [False, True])
def test_insertion_only_target_uses_its_variant_without_unrelated_flank(tmp_path, strand, at_boundary):
    target = GUIDE+"TGG"
    inserted = target if strand == "+" else reverse_complement(target)
    boundary = 20 if at_boundary else 10
    base = "T" if strand == "+" else "A"
    variants = [{"pos1": boundary, "ref": base, "alt": base+inserted, "calls": ["0/1"]}]
    if not at_boundary:
        variants.append({"pos1": boundary+1, "ref": base, "alt": "C", "calls": ["0/1"]})
    e = experiment(tmp_path, base*20, variants, ["SYN"])
    hits = donor_hits(e, "SYN")
    assert len(hits) == 1
    candidate = next(c for c in e["candidates"] if c["configuration_id"] == hits[0]["configuration_id"])
    assert candidate["strand"] == strand
    assert candidate["bounding_interval"] == [boundary, boundary]
    assert candidate["evaluation_interval"] == [boundary-1, boundary]
    assert e["carriage"][0]["supported_copy_count"] == 1
    assert len(candidate["variant_links"]) == 1
    # Only the insertion's establishing REF span is needed from a callable mask.
    e["config"].update(mode="benchmark_simulation", coverage_basis="retained_callset_model")
    rows = [{"donor_id": "SYN", "contig": "chr2", "start_0based": boundary-1, "end_0based_exclusive": boundary}]
    carriage = scan(e["config"], e["reference"], e["materials"], [e["g"]], [], e["donors"], rows, "test:run")[2]
    assert carriage[0]["evidence_supported_carriage"] == "confirmed"
    assert scan(e["config"], e["reference"], e["materials"], [e["g"]], [], e["donors"], [], "test:run")[2][0]["evidence_supported_carriage"] == "unknown"


def test_insertion_at_first_reference_boundary_needs_no_negative_coordinate(tmp_path):
    e = experiment(tmp_path, "T"*20, [{"pos1": 1, "ref": "T", "alt": GUIDE+"TGG"+"T", "calls": ["0/1"]}], ["SYN"])
    assert len(donor_hits(e, "SYN")) == 1
    assert e["candidates"][0]["evaluation_interval"] == [0, 1]
    assert e["carriage"][0]["supported_copy_count"] == 1


def test_minus_strand_new_pam_preserves_genomic_variant(tmp_path):
    e = experiment(tmp_path, reverse_complement(GUIDE+"TGC"), [{"pos1": 1, "ref": "G", "alt": "C", "calls": ["1/1"]}], ["SYN"])
    c = e["candidates"][0]
    assert c["strand"] == "-" and c["protospacer_guide_orientation"] == GUIDE and c["pam_guide_orientation"] == "TGG"
    assert c["target_path"][0] == {"kind": "reference", "pos0": 22}
    assert c["sequence_genomic_forward"] == reverse_complement(GUIDE+"TGG")
    assert e["records"][0]["ref"] == "G" and e["records"][0]["original_alts"] == ["C"]


def test_known_bcl11a_nuclease_regression_with_synthetic_genotypes(tmp_path):
    known = read_json(FIXTURES/"known-case-reference.json")
    v = known["variant"]
    e = experiment(tmp_path, known["sequence_plus_strand"], [{"pos1": v["pos_1based"], "ref": "G", "alt": "C", "calls": ["0|0", "0/1", "1|1", "./."]}],
        ["SYN_GG", "SYN_GC", "SYN_CC", "SYN_MISSING"], start=known["start"], guide=known["guide"]["spacer_5prime_to_3prime"],
        editor_class="nuclease", config_change=lambda c: c["search_profile"].update(max_mismatches=3))
    target = next(c for c in e["candidates"] if c["bounding_interval"] == [210530658, 210530681])
    assert target["strand"] == "-" and target["pam_guide_orientation"] == "TGG"
    assert target["sequence_genomic_forward"] == known["off_target"]["alt_sequence_plus_strand"]
    rows = {r["donor_id"]: r for r in e["carriage"] if r["configuration_id"] == target["configuration_id"]}
    assert [rows[d]["supported_copy_count"] for d in ["SYN_GG", "SYN_GC", "SYN_CC", "SYN_MISSING"]] == [0, 1, 2, None]
    hit = next(h for h in e["hits"] if h["configuration_id"] == target["configuration_id"])
    assert hit["mismatch_positions_1based"] == [1, 8, 12]


def test_bad_reference_and_discordant_duplicate_fail_before_search(tmp_path):
    for variants in [[{"pos1": 1, "ref": "T", "alt": "C", "calls": ["1|1"]}],
                     [{"pos1": 1, "ref": "A", "alt": "C", "calls": [gt]} for gt in ["0|0", "1|1"]]]:
        with pytest.raises(PopulationError):
            experiment(tmp_path, GUIDE+"TGG", variants, ["SYN"])


@pytest.mark.parametrize("info", ["END", "END=0", "END=23;END=2"])
def test_invalid_symbolic_footprint_cannot_shrink_unknown_scope(tmp_path, info):
    with pytest.raises(PopulationError):
        experiment(tmp_path, GUIDE+"TGG", [{"pos1": 1, "ref": "A", "alt": "<DEL>", "info": info, "calls": ["0|1"]}], ["SYN"])


def test_overlapping_records_are_explicit_unknown(tmp_path):
    e = experiment(tmp_path, GUIDE+"TGG", [{"pos1": 1, "ref": "AC", "alt": "GT", "calls": ["1|1"]},
        {"pos1": 2, "ref": "C", "alt": "A", "calls": ["1|1"]}], ["SYN"])
    assert not donor_hits(e, "SYN")
    assert e["carriage"][0]["model_carriage"] == "unknown"


def test_candidate_identity_ignores_padding_search_profile_and_donor(tmp_path):
    e = experiment(tmp_path, "TT"+GUIDE+"TGG"+"AA", [], ["SYN"])
    candidate = next(c for c in e["candidates"] if c["strand"] == "+")
    material = next(m for m in e["materials"] if m["donor_id"])
    padded = material["sequence_path"][2:25]
    assert candidate_identity(e["reference"].sha256, "chr2", "+", e["g"], GUIDE+"TGG", padded)[0] == candidate["configuration_id"]
    e["config"]["search_profile"]["max_mismatches"] = 2
    extra = scan(e["config"], e["reference"], e["materials"], [e["g"]], [], e["donors"], [], "different-run")[0]
    assert candidate["configuration_id"] in {c["configuration_id"] for c in extra}


def test_caps_fail_instead_of_truncating_with_complete_flag(tmp_path):
    with pytest.raises(ResourceLimit):
        experiment(tmp_path, GUIDE+"TGG", [], ["SYN"], config_change=lambda c: c["resource_limits"].update(max_hit_rows=1))
    with pytest.raises(ResourceLimit):
        experiment(tmp_path, GUIDE+"TGG", [], ["SYN"], config_change=lambda c: c["resource_limits"].update(max_reconstructed_bases=30))


def test_unsupported_unbounded_footprint_is_not_limited_to_its_anchor(tmp_path):
    e = experiment(tmp_path, GUIDE+"TGG"+"AAAA", [{"pos1": 27, "ref": "A", "alt": "<DEL>", "calls": ["0|0", "0|1"]}], ["REF", "ALT"])
    assert donor_hits(e, "REF") and not donor_hits(e, "ALT")
    assert all(r["model_carriage"] == "unknown" for r in e["carriage"] if r["donor_id"] == "ALT")


def test_uncarried_symbolic_allele_does_not_mask_another_partition(tmp_path):
    e = experiment(tmp_path, GUIDE+"TGC"+GUIDE+"TGG", [{"pos1": 23, "ref": "C", "alt": "G,<DEL>", "calls": ["1|1"]}], ["SYN"])
    e["config"]["scope"] = {key: [{"contig": "chr2", "start": start, "end": start+23} for start in (0, 23)]
                             for key in ("requested_intervals", "materialized_intervals")}
    materials, _ = build_materializations(e["config"], e["reference"], e["records"], e["donors"], "test:run")
    assert all(not m["excluded_intervals"] for m in materials)
    _, hits, carriage, _ = scan(e["config"], e["reference"], materials, [e["g"]], [], e["donors"], [], "test:run")
    assert len([h for h in hits if h["donor_id"]]) == 4
    assert all(r["supported_copy_count"] == 2 for r in carriage)


def test_callable_contract_requires_the_complete_target(tmp_path):
    def conditional(c):
        c.update(mode="benchmark_simulation", coverage_basis="retained_callset_model")
        c["genotype_policy"]["reference_evidence"] = "reference_imputed_scenario"
    e = experiment(tmp_path, GUIDE+"TGG", [], ["SYN"], config_change=conditional)
    for end, expected in [(22, "unknown"), (23, "confirmed")]:
        rows = [{"donor_id": "SYN", "contig": "chr2", "start_0based": 0, "end_0based_exclusive": end}]
        carriage = scan(e["config"], e["reference"], e["materials"], [e["g"]], [], e["donors"], rows, "test:run")[2]
        assert carriage[0]["evidence_supported_carriage"] == expected


def test_padding_and_input_record_order_preserve_candidate_ids(tmp_path):
    variants = [{"pos1": 23, "ref": "C", "alt": "G", "calls": ["1|1"]},
                {"pos1": 28, "ref": "A", "alt": "C", "calls": ["1|1"]}]
    e = experiment(tmp_path, GUIDE+"TGC"+"AAAAA", variants, ["SYN"])
    original = {c["configuration_id"] for c in e["candidates"]}
    shuffled = experiment(tmp_path, GUIDE+"TGC"+"AAAAA", variants[::-1], ["SYN"])
    assert {c["configuration_id"] for c in shuffled["candidates"]} == original
    # Use the same pinned reference with a smaller materialization interval.
    c = deepcopy(e["config"])
    c["scope"] = {"requested_intervals": [{"contig": "chr2", "start": 0, "end": 23}],
                  "materialized_intervals": [{"contig": "chr2", "start": 0, "end": 23}]}
    materials, _ = build_materializations(c, e["reference"], e["records"], e["donors"], "test:run")
    small = scan(c, e["reference"], materials, [e["g"]], [], e["donors"], [], "test:run")[0]
    assert {row["configuration_id"] for row in small} <= original
    assert small and len(small[0]["variant_links"]) == 1


def test_all_deletion_junction_windows_scan_and_boundary_cut_is_accounted(tmp_path):
    truth = read_json(FIXTURES/"indel-truth.json")["cases"][1]
    e = experiment(tmp_path, truth["reference_sequence"], [{**truth["variant"], "calls": ["1|1"]}], ["SYN"])
    assert len(donor_hits(e, "SYN")) == 2
    c = deepcopy(e["config"])
    c["scope"]["materialized_intervals"] = [{"contig": "chr2", "start": 11, "end": len(truth["reference_sequence"])}]
    materials, _ = build_materializations(c, e["reference"], e["records"], e["donors"], "test:run")
    assert all(m["excluded_intervals"][0]["reason"] == "variant_crosses_materialized_boundary" for m in materials if m["donor_id"])


def test_identical_hbg_copies_and_extra_perfect_match_keep_distinct_loci(tmp_path):
    target = GUIDE+"TGG"
    sequence = target + "N"*(4924-23) + target + "N"*(4924-23) + target
    e = experiment(tmp_path, sequence, [], ["SYN"])
    intended = [{"guide_id": "synthetic_guide", "chrom": "chr2", "start_0based": start,
                 "end_0based_exclusive": start+23, "strand": "+", "locus_name": gene}
                for start, gene in [(0, "HBG1"), (4924, "HBG2")]]
    candidates = scan(e["config"], e["reference"], e["materials"], [e["g"]], intended, e["donors"], [], "test:run")[0]
    assert len({c["configuration_id"] for c in candidates}) == 3
    assert {c["role"] for c in candidates} == {"intended_hbg1", "intended_hbg2", "candidate_offtarget"}
