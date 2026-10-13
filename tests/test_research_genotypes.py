"""Mock checks for streamed copy reconstruction and complete variant windows."""
import hashlib

import pytest

from base_edit_safety.population.common import PopulationError, ResourceLimit
from base_edit_safety.population.haplotypes import local_status
from base_edit_safety.population.research_genotypes import iter_variant_clusters, materialize_cluster
from base_edit_safety.population.variants import parse_call


class Reference:
    sequence = "ACGT" * 100
    build = "artificial"
    aliases = {"artificial_1": "artificial_1"}
    sha256 = hashlib.sha256(sequence.encode()).hexdigest()

    def fetch(self, contig, start, end):
        assert contig == "artificial_1"
        assert 0 <= start < end <= len(self.sequence)
        return self.sequence[start:end]

    def bounds(self, contig):
        assert contig == "artificial_1"
        return [[0, len(self.sequence)]]


POLICY = {"phase": "explicit_ps", "min_gq": None, "min_dp": None,
          "accepted_site_filters": ["PASS"], "accepted_sample_filters": ["."]}


def clusters(tmp_path, rows, copies=2, **kwargs):
    path = tmp_path / "artificial.vcf"
    path.write_text("##fileformat=VCFv4.2\n##reference=artificial\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n" +
                    "".join("artificial_1\t" + row + "\n" for row in rows))
    ploidy = kwargs.pop("ploidy", [{"contig": "artificial_1", "start": 0, "end": 400, "copies": copies}])
    return list(iter_variant_clusters(path, Reference(), "sample", POLICY, ploidy, **kwargs))


def test_haploid_is_one_actual_copy(tmp_path):
    cluster, = clusters(tmp_path, ["81\tv1\tA\tC\t.\tPASS\t.\tGT\t1"], copies=1)
    materials = materialize_cluster(Reference(), cluster, "sample", "run")
    assert len(materials) == 1
    assert materials[0]["local_copy"] == 0
    assert materials[0]["sequence_genomic_forward"][22] == "C"
    assert local_status(materials[0], 80, 103)[0]


def test_default_parse_call_stays_diploid():
    assert parse_call("1", ["GT"], "PASS", POLICY, "s", "d", "c")["call_status"] == "unsupported_ploidy"
    assert parse_call("1", ["GT"], "PASS", POLICY, "s", "d", "c", expected_ploidy=1)["selected_indexes"] == [1]


def test_multiallelic_indexes_and_diploid_copies(tmp_path):
    cluster, = clusters(tmp_path, ["81\tv1\tA\tC,G\t.\tPASS\t.\tGT:PS\t2|1:block"])
    variant, = cluster["variants"]
    assert variant["original_alts"] == ["C", "G"]
    assert variant["calls"]["sample"]["selected_indexes"] == [2, 1]
    mats = materialize_cluster(Reference(), cluster, "sample", "run")
    assert [m["sequence_genomic_forward"][22] for m in mats] == ["G", "C"]


@pytest.mark.parametrize("gt", ["./.", "0/."])
def test_missing_call_masks_footprint(tmp_path, gt):
    cluster, = clusters(tmp_path, [f"81\tv1\tA\tC\t.\tPASS\t.\tGT\t{gt}"])
    for mat in materialize_cluster(Reference(), cluster, "sample", "run"):
        assert mat["sequence_genomic_forward"][22] == "N"
        assert local_status(mat, 80, 90) == (False, "unresolved_variant_footprint")


@pytest.mark.parametrize("second,expected", [("0|1:block", True), ("0|1:other", False), ("0/1:.", False)])
def test_local_phase_requires_common_block(tmp_path, second, expected):
    cluster, = clusters(tmp_path, ["81\tv1\tA\tC\t.\tPASS\t.\tGT:PS\t0|1:block",
                                 f"86\tv2\tC\tG\t.\tPASS\t.\tGT:PS\t{second}"])
    mat = materialize_cluster(Reference(), cluster, "sample", "run")[0]
    assert local_status(mat, 80, 90)[0] is expected
    assert local_status(mat, 80, 85)[0] is True


def test_deletion_padding_uses_complete_ref_span(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tAC\tA\t.\tPASS\t.\tGT\t1/1",
                                 "105\tv2\tA\tC\t.\tPASS\t.\tGT\t1/1"])
    assert cluster["interval"] == {"contig": "artificial_1", "start": 58, "end": 127}
    mat = materialize_cluster(Reference(), cluster, "sample", "run")[0]
    assert len(mat["sequence_path"]) == 68
    assert 81 not in [p.get("pos0") for p in mat["sequence_path"]]
    assert mat["sequence_path"][21]["pos0"] == 79


def test_insertion_retains_exact_path(tmp_path):
    cluster, = clusters(tmp_path, ["81\tins\tA\tAGG\t.\tPASS\t.\tGT\t1"], copies=1)
    mat, = materialize_cluster(Reference(), cluster, "sample", "run")
    points = [p for p in mat["sequence_path"] if p["kind"] == "insertion"]
    assert [(p["boundary0"], p["offset0"], p["base"]) for p in points] == [(81, 0, "G"), (81, 1, "G")]


def test_unknown_and_zero_ploidy_are_excluded(tmp_path):
    for ploidy, reason in [([], "unknown_or_changing_ploidy"),
                           ([{"contig": "artificial_1", "start": 0, "end": 400, "copies": 0}], "zero_ploidy")]:
        cluster, = clusters(tmp_path, ["81\tv1\tA\tC\t.\tPASS\t.\tGT\t0/0"], ploidy=ploidy)
        assert cluster["copies"] == 0
        assert cluster["exclusions"][0]["reason"] == reason
        assert materialize_cluster(Reference(), cluster, "sample", "run") == []


def test_unbounded_unsupported_call_excludes_its_contig(tmp_path):
    cluster, = clusters(tmp_path, ["81\tv1\tA\t<DEL>\t.\tPASS\t.\tGT\t1/1"])
    exclusion, = cluster["unbounded_exclusions"]
    assert (exclusion["contig"], exclusion["start"], exclusion["end"]) == ("artificial_1", 0, 400)
    assert "affected interval" in exclusion["missing_information"]
    cluster, = clusters(tmp_path, ["81\tv1\tA\t<DEL>\t.\tPASS\tEND=95\tGT\t1/1"])
    mat = materialize_cluster(Reference(), cluster, "sample", "run")[0]
    assert local_status(mat, 90, 94) == (False, "unresolved_variant_footprint")


def test_sorted_stream_and_limits_fail_closed(tmp_path):
    rows = ["81\tv1\tA\tC\t.\tPASS\t.\tGT\t0/1", "85\tv2\tA\tC\t.\tPASS\t.\tGT\t0/1"]
    with pytest.raises(PopulationError, match="sorted"):
        clusters(tmp_path, rows[::-1])
    with pytest.raises(ResourceLimit, match="max_variant_rows"):
        clusters(tmp_path, rows, max_rows=1)
    cluster, = clusters(tmp_path, rows, max_cluster_bases=25)
    assert cluster["resource_excluded"]
    assert cluster["exclusions"][0]["reason"] == "local_reconstruction_resource_limit"
    assert materialize_cluster(Reference(), cluster, "sample", "run") == []


def test_identical_duplicate_retains_source_rows_once(tmp_path):
    row = "81\tv1\tA\tC\t.\tPASS\t.\tGT\t0/1"
    cluster, = clusters(tmp_path, [row, row])
    assert cluster["source_rows"] == [1, 2]
    assert len(cluster["variants"]) == 1


def test_filtered_call_and_ploidy_mismatch_remain_unknown(tmp_path):
    cluster, = clusters(tmp_path, ["81\tv1\tA\tC\t.\tLowQual\t.\tGT\t0/0"])
    for mat in materialize_cluster(Reference(), cluster, "sample", "run"):
        assert local_status(mat, 80, 90)[0] is False
    cluster, = clusters(tmp_path, ["81\tv1\tA\tC\t.\tPASS\t.\tGT\t0/1"], copies=1)
    mat, = materialize_cluster(Reference(), cluster, "sample", "run")
    assert mat["excluded_intervals"][0]["reason"] == "unsupported_ploidy"


def test_padding_crossing_ploidy_boundary_is_split(tmp_path):
    ploidy = [{"contig": "artificial_1", "start": 0, "end": 90, "copies": 1},
              {"contig": "artificial_1", "start": 90, "end": 400, "copies": 2}]
    left, right = clusters(tmp_path, ["81\tv1\tA\tC\t.\tPASS\t.\tGT\t1"], ploidy=ploidy)
    assert [(c["interval"]["start"], c["interval"]["end"], c["copies"]) for c in (left, right)] == [(58, 90, 1), (90, 103, 2)]
    assert left["source_rows"] == [1] and right["source_rows"] == []
    mat, = materialize_cluster(Reference(), left, "sample", "run")
    assert mat["sequence_genomic_forward"][22] == "C"
    assert local_status(mat, 58, 81)[0]
    assert not local_status(mat, 80, 103)[0]


def test_unbounded_missing_call_is_reported_without_reference_hits(tmp_path):
    cluster, = clusters(tmp_path, ["81\tv1\tA\t<DEL>\t.\tPASS\t.\tGT\t./."])
    assert cluster["unbounded_exclusions"][0]["reason"] == "unbounded_unsupported_allele"


def test_reference_mismatch_and_unnormalized_indel_are_rejected(tmp_path):
    with pytest.raises(PopulationError, match="REF disagrees"):
        clusters(tmp_path, ["81\tv1\tC\tG\t.\tPASS\t.\tGT\t0/1"])
    with pytest.raises(PopulationError, match="left normalized"):
        clusters(tmp_path, ["81\tins\tA\tAA\t.\tPASS\t.\tGT\t1/1"])


def test_budget_callback_runs_before_cluster_yield(tmp_path):
    checks = []
    def check():
        checks.append(True)
        if len(checks) == 7:
            raise ResourceLimit("streaming wall cap")
    rows = [f"{pos+1}\tv{pos}\t{Reference.sequence[pos]}\t{next(b for b in 'ACGT' if b != Reference.sequence[pos])}\t.\tPASS\t.\tGT\t0/1" for pos in range(80, 100)]
    with pytest.raises(ResourceLimit, match="streaming wall cap"):
        clusters(tmp_path, rows, check_limits=check)
    assert len(checks) == 7


def test_explicit_contig_mask_retains_star_and_ploidy_without_materialization(tmp_path):
    cluster, = clusters(tmp_path, ["81\tv1\tA\t*\t.\tPASS\t.\tGT\t0/1"],
                        excluded_contigs=["artificial_1"])
    assert cluster["copies"] == 2
    assert cluster["sample_contig_excluded"] is True
    assert cluster["source_rows"] == [1]
    variant, = cluster["variants"]
    assert variant["unbounded_footprint"] is True
    assert variant["calls"]["sample"]["selected_indexes"] == [0, 1]
    assert variant["calls"]["sample"]["call_status"] == "complete"
    assert materialize_cluster(Reference(), cluster, "sample", "run") == []


@pytest.mark.parametrize("alt", ["<DEL>", "A]artificial_1:100]", "*,<DEL>"])
def test_contig_mask_keeps_other_unbounded_exclusion_evidence(tmp_path, alt):
    cluster, = clusters(tmp_path, [f"81\tv1\tA\t{alt}\t.\tPASS\t.\tGT\t1/1"],
                        excluded_contigs=["artificial_1"])
    assert cluster["unbounded_exclusions"]
    assert materialize_cluster(Reference(), cluster, "sample", "run") == []


@pytest.mark.parametrize("deletion,marker,copies", [
    ("1|0:block", "2|1:block", 2),
    ("1/0:.", "1/2:.", 2),
    ("1|0:left", "1|2:right", 2),
    ("1/1:.", "2/2:.", 2),
    ("1:.", "2:.", 1),
])
def test_star_reuses_covering_deletion_without_changing_original_indexes(tmp_path, deletion, marker, copies):
    cluster, = clusters(tmp_path, [f"81\tdel\tACGT\tA\t.\tPASS\t.\tGT:PS\t{deletion}",
                                 f"82\tstar\tC\tT,*\t.\tPASS\t.\tGT:PS\t{marker}"], copies=copies)
    mats = materialize_cluster(Reference(), cluster, "sample", "run")
    lo, hi = cluster["interval"]["start"], cluster["interval"]["end"]
    deleted = Reference.sequence[lo:81] + Reference.sequence[84:hi]
    substituted = Reference.sequence[lo:81] + "T" + Reference.sequence[82:hi]
    assert sorted(m["sequence_genomic_forward"] for m in mats) == sorted(
        [deleted] * copies if deletion.startswith("1/1") or copies == 1 else [deleted, substituted])
    assert all(not m["excluded_intervals"] for m in mats)
    assert all(local_status(m, 80, 85)[0] for m in mats)
    assert cluster["variants"][1]["calls"]["sample"]["gt"] == marker.split(":")[0]
    assert cluster["variants"][1]["original_alts"] == ["T", "*"]
    assert cluster["variant_accounting"]["resolved_spanning_deletion_variant_records"] == 1
    assert not cluster["unbounded_exclusions"]


def test_star_requires_point_coverage_not_complete_ref_coverage(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tAC\tA\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "82\tstar\tCG\tTA,*\t.\tPASS\t.\tGT:PS\t2|1:block"])
    left, right = materialize_cluster(Reference(), cluster, "sample", "run")
    assert left["deletion_events"][0]["start"] == 81
    assert [p["base"] for p in right["sequence_path"] if p.get("pos0") in (81, 82)] == ["T", "A"]
    assert not left["excluded_intervals"] and not right["excluded_intervals"]


def test_star_constraint_does_not_merge_unrelated_phase_blocks(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tACGT\tA\t.\tPASS\t.\tGT:PS\t1|0:left",
                                 "82\tstar\tC\tT,*\t.\tPASS\t.\tGT:PS\t1|2:right",
                                 "90\tother\tC\tG\t.\tPASS\t.\tGT:PS\t0|1:unrelated"])
    for mat in materialize_cluster(Reference(), cluster, "sample", "run"):
        assert local_status(mat, 80, 85) == (True, "resolved_spanning_deletion_phase_constraint")
        assert local_status(mat, 80, 95) == (False, "unresolved_phase_blocks")
    assert [v["calls"]["sample"]["ps"] for v in cluster["variants"]] == ["left", "right", "unrelated"]


@pytest.mark.parametrize("marker", ["1|2:block", "0|2:block"])
def test_same_phase_contradiction_stays_masked_across_deletion_junction(tmp_path, marker):
    from base_edit_safety.population.search import evaluation_span
    cluster, = clusters(tmp_path, ["73\tleft\tA\tT\t.\tPASS\t.\tGT\t0/0",
                                 "81\tdel\tACGT\tA\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 f"82\tstar\tC\tT,*\t.\tPASS\t.\tGT:PS\t{marker}",
                                 "93\tright\tA\tT\t.\tPASS\t.\tGT\t0/0"])
    mat = materialize_cluster(Reference(), cluster, "sample", "run")[0]
    path = [p for p in mat["sequence_path"] if 70 <= p.get("pos0", -1) < 93]
    assert len(path) == 23
    assert not mat["deletion_events"]
    assert local_status(mat, *evaluation_span(path, mat)) == (False, "unresolved_variant_footprint")
    assert cluster["variant_accounting"]["unsupported_variant_records"] == 2


def test_wrong_copy_deletion_cannot_bound_an_unflanked_star(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tACGT\tA\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "82\tstar\tC\t*\t.\tPASS\t.\tGT:PS\t0|1:block"])
    assert [(g["start"], g["end"]) for g in cluster["unbounded_exclusions"]] == [(0, 400)]
    assert not cluster["variant_accounting"]["resolved_spanning_deletion_variant_records"]


def test_orphan_star_uses_present_bases_and_keeps_every_intervening_record(tmp_path):
    rows = ["41\tleft\tA\tC\t.\tPASS\t.\tGT\t1/1",
            "121\tmissing\tA\tC\t.\tPASS\t.\tGT\t./.",
            "201\tstar\tACG\tTCA,*\t.\tPASS\t.\tGT\t1/2",
            "202\tinside\tC\tG\t.\tPASS\t.\tGT\t0/0",
            "281\tright\tA\tC\t.\tPASS\t.\tGT\t0/0"]
    cluster, = clusters(tmp_path, rows)
    assert cluster["source_rows"] == [1, 2, 3, 4, 5]
    star = next(v for v in cluster["variants"] if v["original_id"] == "star")
    bounds = star["spanning_deletion_bounds"]
    assert (bounds["start"], bounds["end"]) == (41, 280)
    assert (bounds["left_anchor"]["source_record_index"], bounds["right_anchor"]["source_record_index"]) == (1, 5)
    assert not cluster["unbounded_exclusions"]
    known, unknown = materialize_cluster(Reference(), cluster, "sample", "run")
    assert not any(g["variant_id"] == star["variant_id"] for g in known["excluded_intervals"])
    assert [p["base"] for p in known["sequence_path"] if p.get("pos0") in (200, 201, 202)] == list("TCA")
    gap = next(g for g in unknown["excluded_intervals"] if g["variant_id"] == star["variant_id"])
    assert (gap["start"], gap["end"]) == (41, 280)
    assert local_status(known, 258, 281)[0]
    assert not local_status(unknown, 258, 281)[0]
    assert all(local_status(mat, 280, 303)[0] for mat in (known, unknown))


@pytest.mark.parametrize("gt", ["1|2:block", "1/2:.", "0|2:block", "0/2:."])
def test_bounded_orphan_masks_only_its_carried_copy(tmp_path, gt):
    cluster, = clusters(tmp_path, ["41\tleft\tA\tC\t.\tPASS\t.\tGT:PS\t0/0:.",
                                 f"81\tstar\tA\tG,*\t.\tPASS\t.\tGT:PS\t{gt}",
                                 "121\tright\tA\tC\t.\tPASS\t.\tGT:PS\t0/0:."])
    known, unknown = materialize_cluster(Reference(), cluster, "sample", "run")
    assert not cluster["unbounded_exclusions"]
    assert not known["excluded_intervals"]
    assert [p["base"] for p in known["sequence_path"] if p.get("pos0") == 80] == ["G" if gt[0] == "1" else "A"]
    assert [(g["start"], g["end"]) for g in unknown["excluded_intervals"]] == [(41, 120)]
    assert local_status(known, 70, 93)[0]
    assert not local_status(unknown, 70, 93)[0]
    marker = next(v for v in cluster["variants"] if v["original_id"] == "star")
    assert cluster["variant_dispositions"][marker["variant_id"]] == sorted([
        "applied" if gt[0] == "1" else "reference", "unsupported"])


@pytest.mark.parametrize("marker,other,expected", [("1|2:block", "0|1:block", True),
                                                  ("1|2:block", "0|1:other", False),
                                                  ("1/2:.", "0/1:.", False)])
def test_orphan_phase_dependency_covers_its_possible_deletion(tmp_path, marker, other, expected):
    cluster, = clusters(tmp_path, ["41\tleft\tA\tC\t.\tPASS\t.\tGT:PS\t0/0:.",
                                 f"81\tstar\tA\tG,*\t.\tPASS\t.\tGT:PS\t{marker}",
                                 f"101\tright\tA\tC\t.\tPASS\t.\tGT:PS\t{other}"])
    known, unknown = materialize_cluster(Reference(), cluster, "sample", "run")
    # This window is beyond the marker but crosses its unknown deletion bound.
    assert local_status(known, 90, 113)[0] is expected
    assert not local_status(unknown, 90, 113)[0]
    assert all(local_status(mat, 100, 123)[0] for mat in (known, unknown))


@pytest.mark.parametrize("left_gt,left_filter,right_gt", [("./.", "PASS", "0/0"),
                                                          ("0/0", "LowQual", "0/0"),
                                                          ("0/0", "PASS", "0/.")])
def test_missing_or_filtered_flanks_cannot_bound_an_orphan(tmp_path, left_gt, left_filter, right_gt):
    result = clusters(tmp_path, [f"41\tleft\tA\tC\t.\t{left_filter}\t.\tGT\t{left_gt}",
                                "81\tstar\tA\t*\t.\tPASS\t.\tGT\t0/1",
                                f"121\tright\tA\tC\t.\tPASS\t.\tGT\t{right_gt}"])
    assert any(c["unbounded_exclusions"] for c in result)


def test_unselected_deletion_does_not_bound_an_orphan_star(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tACGT\tA\t.\tPASS\t.\tGT\t0/0",
                                 "82\tstar\tC\t*\t.\tPASS\t.\tGT\t0/1"])
    assert cluster["unbounded_exclusions"]


def test_overlap_is_assessed_per_copy_and_changed_bases(tmp_path):
    cluster, = clusters(tmp_path, ["81\tone\tAC\tTC\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "81\ttwo\tAC\tAG\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "81\tthree\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:block"])
    mats = materialize_cluster(Reference(), cluster, "sample", "run")
    assert [[p["base"] for p in m["sequence_path"] if p.get("pos0") in (80, 81)] for m in mats] == [["T", "G"], ["G", "C"]]
    assert all(not m["excluded_intervals"] for m in mats)


def test_incompatible_carried_edits_exclude_only_the_affected_copy(tmp_path):
    cluster, = clusters(tmp_path, ["81\tone\tA\tC\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "81\ttwo\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:block"])
    a, b = materialize_cluster(Reference(), cluster, "sample", "run")
    assert {g["reason"] for g in a["excluded_intervals"]} == {"overlapping_carried_edits"}
    assert not b["excluded_intervals"]
    assert not local_status(a, 80, 103)[0] and local_status(b, 80, 103)[0]


def test_deletion_conflicts_with_same_copy_called_present_reference(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tACGT\tA\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "82\tsnp\tC\tT\t.\tPASS\t.\tGT:PS\t0|1:block"])
    deleted, present = materialize_cluster(Reference(), cluster, "sample", "run")
    assert not deleted["deletion_events"]
    assert {g["reason"] for g in deleted["excluded_intervals"]} == {"overlapping_carried_edits"}
    assert not local_status(deleted, 70, 93)[0]
    assert not present["excluded_intervals"]
    assert [p["base"] for p in present["sequence_path"] if p.get("pos0") == 81] == ["T"]
    assert cluster["variant_accounting"]["unsupported_variant_records"] == 2


def test_deletion_anchor_remains_present_for_an_overlapping_call(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tACGT\tA\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "81\tsnp\tA\tC\t.\tPASS\t.\tGT:PS\t0|1:block"])
    deleted, present = materialize_cluster(Reference(), cluster, "sample", "run")
    assert all(not mat["excluded_intervals"] for mat in (deleted, present))
    assert [(d["start"], d["end"]) for d in deleted["deletion_events"]] == [(81, 84)]
    assert [p["base"] for p in present["sequence_path"] if p.get("pos0") == 80] == ["C"]


def test_unsupported_literal_allele_with_unselected_star_stays_local(tmp_path):
    cluster, = clusters(tmp_path, ["81\tlong\tA\t" + "A" + "G"*60 + ",*\t.\tPASS\t.\tGT\t1/1"],
                        max_allele_length=50)
    assert not cluster["unbounded_exclusions"]
    for mat in materialize_cluster(Reference(), cluster, "sample", "run"):
        assert [(g["start"], g["end"]) for g in mat["excluded_intervals"]] == [(80, 81)]


def test_variant_spanning_ploidy_change_masks_both_pieces_once_in_accounting(tmp_path):
    ploidy = [{"contig": "artificial_1", "start": 0, "end": 82, "copies": 1},
              {"contig": "artificial_1", "start": 82, "end": 400, "copies": 2}]
    pieces = clusters(tmp_path, ["81\tdel\tACGT\tA\t.\tPASS\t.\tGT\t1"], ploidy=ploidy)
    assert sum(len(c["accounted_variant_ids"]) for c in pieces) == 1
    assert sum(len(c["source_rows"]) for c in pieces) == 1
    for cluster in pieces:
        for mat in materialize_cluster(Reference(), cluster, "sample", "run"):
            assert mat["excluded_intervals"]
            assert not mat["deletion_events"]


@pytest.mark.parametrize("gt", ["1/2", "./.", "0/."])
def test_star_bounds_do_not_bound_another_symbolic_allele(tmp_path, gt):
    result = clusters(tmp_path, ["41\tleft\tA\tC\t.\tPASS\t.\tGT\t0/0",
                                f"81\tstar\tA\t*,<DEL>\t.\tPASS\t.\tGT\t{gt}",
                                "121\tright\tA\tC\t.\tPASS\t.\tGT\t0/0"])
    assert any(c["unbounded_exclusions"] for c in result)


@pytest.mark.parametrize("gt", ["9/.", "./9", "0/9/.", "1//.", "0|1/."])
def test_missingness_does_not_hide_invalid_genotype_indexes_or_syntax(tmp_path, gt):
    with pytest.raises(PopulationError, match="genotype syntax or allele index"):
        clusters(tmp_path, [f"81\tv\tA\tC\t.\tPASS\t.\tGT\t{gt}"])


def test_unknown_literal_base_has_a_bounded_footprint(tmp_path):
    cluster, = clusters(tmp_path, ["81\tv\tA\tN\t.\tPASS\t.\tGT\t1/1"])
    assert not cluster["unbounded_exclusions"]
    for mat in materialize_cluster(Reference(), cluster, "sample", "run"):
        assert [(g["start"], g["end"]) for g in mat["excluded_intervals"]] == [(80, 81)]


def test_research_supports_exact_long_indel_and_its_star(tmp_path):
    ref = Reference.sequence[80:140]
    cluster, = clusters(tmp_path, [f"81\tdel\t{ref}\tA\t.\tPASS\t.\tGT\t1/1",
                                 "100\tstar\tT\t*\t.\tPASS\t.\tGT\t1/1"])
    assert not cluster["resource_excluded"]
    assert cluster["variant_accounting"]["resolved_spanning_deletion_variant_records"] == 1
    for mat in materialize_cluster(Reference(), cluster, "sample", "run"):
        assert not mat["excluded_intervals"]
        assert mat["sequence_genomic_forward"] == Reference.sequence[58:81] + Reference.sequence[140:162]


def test_supported_star_does_not_expand_to_distant_optional_anchors(tmp_path):
    result = clusters(tmp_path, ["1\tleft\tA\tC\t.\tPASS\t.\tGT\t0/0",
                                "161\tdel\tACGT\tA\t.\tPASS\t.\tGT\t1/1",
                                "162\tstar\tC\t*\t.\tPASS\t.\tGT\t1/1",
                                "321\tright\tA\tC\t.\tPASS\t.\tGT\t0/0"], max_cluster_bases=100)
    assert len(result) == 3
    assert all(not c["resource_excluded"] for c in result)
    assert result[1]["interval"] == {"contig": "artificial_1", "start": 138, "end": 186}


def test_large_local_resource_exclusion_preserves_source_accounting_and_later_region(tmp_path):
    ref = Reference.sequence[80:200]
    result = clusters(tmp_path, [f"81\tlarge\t{ref}\tA\t.\tPASS\t.\tGT\t1/1",
                                "321\tlater\tA\tC\t.\tPASS\t.\tGT\t1/1"], max_cluster_bases=100)
    assert result[0]["resource_excluded"]
    assert result[0]["copies"] == 2
    assert result[0]["source_rows"] == [1]
    assert result[0]["variant_accounting"]["unsupported_variant_records"] == 1
    assert not materialize_cluster(Reference(), result[0], "sample", "run")
    assert not result[1]["resource_excluded"]
    assert materialize_cluster(Reference(), result[1], "sample", "run")


@pytest.mark.parametrize("alt", ["A]artificial_1:100]", ".A", "A.", "<TRA>"])
def test_end_does_not_bound_an_unsupported_breakend(tmp_path, alt):
    cluster, = clusters(tmp_path, [f"81\tbnd\tA\t{alt}\t.\tPASS\tEND=95\tGT\t1/1"])
    assert cluster["unbounded_exclusions"]


@pytest.mark.parametrize("gt,excluded", [("1", True), ("0/1", True), ("./.", True),
                                        ("0/.", True), ("0", False), ("0/0", False)])
def test_unknown_regional_ploidy_does_not_contain_unbounded_alleles(tmp_path, gt, excluded):
    ploidy = [{"contig": "artificial_1", "start": 0, "end": 70, "copies": 2},
              {"contig": "artificial_1", "start": 100, "end": 400, "copies": 2}]
    result = clusters(tmp_path, [f"81\tunbounded\tA\t<DEL>\t.\tPASS\t.\tGT\t{gt}"], ploidy=ploidy)
    gaps = [gap for c in result for gap in c["unbounded_exclusions"]]
    assert bool(gaps) is excluded
    assert all((g["start"], g["end"]) == (0, 400) for g in gaps)
    unknown = next(c for c in result if c["interval"]["start"] == 70)
    assert unknown["copies"] == 0
    assert unknown["variants"][0]["calls"]["sample"]["selected_indexes"] is None


def test_rejected_call_cannot_prove_an_unbounded_allele_uncarried_without_ploidy(tmp_path):
    cluster, = clusters(tmp_path, ["81\tunbounded\tA\t<DEL>\t.\tLowQual\t.\tGT\t0/0"], ploidy=[])
    assert cluster["unbounded_exclusions"]


def test_explicit_zero_copies_remain_distinct_from_unknown_copies(tmp_path):
    cluster, = clusters(tmp_path, ["81\tunbounded\tA\t<DEL>\t.\tPASS\t.\tGT\t1/1"], copies=0)
    assert not cluster["unbounded_exclusions"]
    assert cluster["variants"][0]["calls"]["sample"]["call_status"] == "zero_ploidy"


@pytest.mark.parametrize("alt,gt,expected", [
    ("A]artificial_2:100]", "1/1", {"artificial_1", "artificial_2"}),
    ("A]artificial_2:100]", "0/0", set()),
    ("C,A]artificial_2:100]", "1/1", set()),
    ("A]outside_scope:100]", "1/1", {"artificial_1"}),
    ("A.", "1/1", {"artificial_1", "artificial_2", "artificial_3"}),
    ("<TRA>", "0/.", {"artificial_1", "artificial_2", "artificial_3"}),
])
def test_breakend_scope_includes_only_admitted_potentially_affected_chromosomes(tmp_path, alt, gt, expected):
    class ThreeChromosomes(Reference):
        aliases = {f"artificial_{i}": f"artificial_{i}" for i in (1, 2, 3)}

        def bounds(self, contig):
            assert contig in self.aliases
            return [[0, 400]]

    path = tmp_path / "breakends.vcf"
    path.write_text("##fileformat=VCFv4.3\n##reference=artificial\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                    f"artificial_1\t81\tbnd\tA\t{alt}\t.\tPASS\tEND=95\tGT\t{gt}\n")
    ploidy = [{"contig": c, "start": 0, "end": 400, "copies": 2} for c in ThreeChromosomes.aliases]
    cluster, = iter_variant_clusters(path, ThreeChromosomes(), "sample", POLICY, ploidy)
    assert {g["contig"] for g in cluster["unbounded_exclusions"]} == expected


def test_resolved_star_does_not_make_other_bounded_allele_unbounded(tmp_path):
    cluster, = clusters(tmp_path, ["81\tdel\tACGT\tA\t.\tPASS\t.\tGT:PS\t1|0:block",
                                 "82\tstar\tCGT\tAA,*\t.\tPASS\t.\tGT:PS\t2|1:block"])
    assert not cluster["unbounded_exclusions"]
    a, b = materialize_cluster(Reference(), cluster, "sample", "run")
    assert not a["excluded_intervals"]
    assert [(g["start"], g["end"]) for g in b["excluded_intervals"]] == [(81, 84)]
