"""Compare donor reconstruction with independently written mock sequences."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import pytest

from base_edit_safety.population.haplotypes import local_status
from base_edit_safety.population.research_coverage import Coverage
from base_edit_safety.population.research_genotypes import iter_variant_clusters, materialize_cluster
from base_edit_safety.population.search import evaluation_span


FIXTURES = Path(__file__).parent / "fixtures" / "donor-coverage"
TRUTH = json.loads((FIXTURES / "fixtures.json").read_text())
COMPLEMENT = str.maketrans("ACGT", "TGCA")


class FixtureReference:
    build = "artificial"

    def __init__(self, sequences):
        self.sequences = sequences
        self.aliases = {contig: contig for contig in sequences}
        self.sha256 = hashlib.sha256(json.dumps(sequences, sort_keys=True).encode()).hexdigest()

    def fetch(self, contig, start, end):
        assert 0 <= start < end <= len(self.sequences[contig])
        return self.sequences[contig][start:end]

    def bounds(self, contig):
        return [[0, len(self.sequences[contig])]]


def assemble_region(reference, region, copy, materials, exclusions):
    """Join padded neighborhoods with their unchanged reference context."""
    contig, start, end = (region[k] for k in ("contig", "start", "end"))
    result = {"interval": {"contig": contig, "start": start, "end": end},
              "sequence_path": [], "sequence_genomic_forward": "", "excluded_intervals": [],
              "heterozygous_sites": [], "phase_constraints": [], "deletion_events": [], "variant_intervals": []}
    sequence, cursor = [], start

    def append_reference(lo, hi):
        sequence.append(reference.sequences[contig][lo:hi])
        result["sequence_path"].extend({"kind": "reference", "pos0": pos, "contig": contig,
                                        "base": reference.sequences[contig][pos]}
                                       for pos in range(lo, hi))

    for material in materials:
        interval = material["interval"]
        if interval["contig"] != contig or material["local_copy"] != copy \
                or not start <= interval["start"] < interval["end"] <= end:
            continue
        assert cursor <= interval["start"]
        append_reference(cursor, interval["start"])
        sequence.append(material["sequence_genomic_forward"])
        for field in ("sequence_path", "excluded_intervals", "heterozygous_sites",
                      "phase_constraints", "deletion_events", "variant_intervals"):
            result[field].extend(material.get(field, []))
        cursor = interval["end"]
    append_reference(cursor, end)
    result["sequence_genomic_forward"] = "".join(sequence)
    result["excluded_intervals"].extend(g for g in exclusions if g["contig"] == contig)
    return result


def add_orientations(windows, key, sequence):
    windows[key + ("+", sequence)] += 1
    windows[key + ("-", sequence.translate(COMPLEMENT)[::-1])] += 1


@pytest.mark.parametrize("case", TRUTH["cases"], ids=lambda case: case["id"])
def test_independent_declared_windows(tmp_path, case):
    reference = FixtureReference(case["reference_sequences"])
    vcf = tmp_path / "independent.vcf"
    vcf.write_text(case["vcf"] + "\n")
    policy = {"phase": "vcf_implicit_phase_set" if case["id"] == "case4_multiallelic_shared_anchor" else "explicit_ps",
              "min_gq": None, "min_dp": None, "accepted_site_filters": ["."],
              "accepted_sample_filters": ["."]}
    clusters = list(iter_variant_clusters(vcf, reference, "S1", policy, case["ploidy_intervals"]))
    cluster_materials = [materialize_cluster(reference, cluster, "S1", "oracle") for cluster in clusters]
    materials = [material for group in cluster_materials for material in group]
    exclusions = [gap for cluster in clusters for gap in cluster["unbounded_exclusions"]]
    variants = {variant["source_record_index"]: variant for cluster in clusters for variant in cluster["variants"]}
    rows = [line.split("\t") for line in case["vcf"].splitlines() if not line.startswith("#")]
    assert len(variants) == len(rows)
    for index, row in enumerate(rows, 1):
        contig, pos, _, ref, alt = row[:5]
        raw = dict(zip(row[8].split(":"), row[9].split(":")))
        variant, call = variants[index], variants[index]["calls"]["S1"]
        assert reference.fetch(contig, int(pos)-1, int(pos)-1+len(ref)) == ref
        assert (variant["pos0"], variant["ref"], variant["original_alts"]) == (int(pos)-1, ref, alt.split(","))
        assert (call["gt"], call["ps"]) == (raw["GT"], raw.get("PS"))
        if call["call_status"] == "complete":
            assert call["selected_indexes"] == [int(value) for value in re.split(r"[/|]", raw["GT"])]

    # Unphased records describe unordered copies; their manual ordering is not phase evidence.
    unordered = any(call["phase_status"] == "unresolved_unphased"
                    for variant in variants.values() for call in variant["calls"].values())
    expected, actual = Counter(), Counter()
    declared = {(r["contig"], r["start"], r["end"], copy)
                for r in case["ploidy_intervals"] for copy in range(r["copies"])}
    provided = {(h["contig"], *h["ref_interval"], h["local_copy"]-1) for h in case["haplotypes"]}
    assert provided == declared
    assert case["ploidy_source"]
    for hap in case["haplotypes"]:
        start, end = hap["ref_interval"]
        copy = hap["local_copy"] - 1  # The independent fixtures use 1-based copy labels.
        domain = (hap["contig"], start, end, None if unordered else copy)
        for offset in range(len(hap["seq"])-22):
            if any(offset < hi and lo < offset+23 for lo, hi in hap.get("exclude", [])) \
                    or any(lo <= offset < hi for lo, hi in hap.get("excluded_window_starts", [])):
                continue
            window = hap["seq"][offset:offset+23]
            assert set(window) <= set("ACGT")
            add_orientations(expected, domain + (offset,), window)
        region = {"contig": hap["contig"], "start": start, "end": end}
        assembled = assemble_region(reference, region, copy, materials, exclusions)
        sequence, path = assembled["sequence_genomic_forward"], assembled["sequence_path"]
        for offset in range(len(sequence)-22):
            window = sequence[offset:offset+23]
            if set(window) <= set("ACGT") and local_status(assembled, *evaluation_span(path[offset:offset+23], assembled))[0]:
                add_orientations(actual, domain + (offset,), window)
    assert actual == expected
    coverage = Coverage([
        {"contig": contig, "bases": len(sequence), "unknown_bases": 0,
         "acgt_intervals": [(0, len(sequence))],
         "evaluable_oriented_windows_per_guide": 2*max(0, len(sequence)-22)}
        for contig, sequence in reference.sequences.items()
    ], {"ploidy_intervals": case["ploidy_intervals"]})
    for cluster, group in zip(clusters, cluster_materials):
        coverage.cluster(cluster, group)
    for gap in exclusions:
        assert (gap["start"], gap["end"]) == (0, len(reference.sequences[gap["contig"]]))
        coverage.exclude_contig(gap["contig"], gap["reason"])
    coverage_rows, _, _ = coverage.finish({})
    assert sum(row["sample_evaluated_copy_oriented_windows_per_guide"] for row in coverage_rows) == sum(expected.values())
    for row in coverage_rows:
        assert row["reference_oriented_windows_per_guide"] == 2*max(0, len(reference.sequences[row["contig"]])-22)
    candidates = sum(count for key, count in actual.items() if key[-1][21:] == "GG")
    hits = sum(count for key, count in actual.items() if key[-1][21:] == "GG"
               and sum(a != b for a, b in zip(key[-1][:20], TRUTH["guide"])) <= 6)
    assert (candidates, hits) == (case["expected_candidates"], case["expected_hits"])


def test_independent_oracle_counts():
    result = subprocess.run([sys.executable, "oracle.py"], cwd=FIXTURES, text=True,
                            capture_output=True, check=True)
    assert "TOTAL candidates=18 hits=18" in result.stdout
