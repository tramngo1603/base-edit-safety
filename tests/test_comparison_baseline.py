"""Preserve V3 allele order and reject malformed genotype inputs."""

import pytest

from base_edit_safety.population.common import PopulationError
from test_population_sequences import GUIDE, donor_hits, experiment


@pytest.mark.parametrize("alt", ["", "C,", ",C", "C,,A"])
def test_empty_alt_indexes_are_invalid_input(tmp_path, alt):
    with pytest.raises(PopulationError, match="Empty VCF ALT"):
        experiment(tmp_path, GUIDE + "TGG", [
            {"pos1": 23, "ref": "G", "alt": alt, "calls": ["1|1"]}
        ], ["SYN"])


@pytest.mark.parametrize("fmt,call", [
    ("", "1|2"),
    (":GT", ".:1|2"),
    ("GT:", "1|2:."),
    ("GT::PS", "1|2:.:7"),
    ("GT:GT", "1|2:1|2"),
])
def test_empty_and_duplicate_format_keys_are_invalid_input(tmp_path, fmt, call):
    with pytest.raises(PopulationError, match="VCF FORMAT key"):
        experiment(tmp_path, GUIDE + "TGG", [
            {"pos1": 23, "ref": "G", "alt": "C,A", "format": fmt, "calls": [call]}
        ], ["SYN"])


@pytest.mark.parametrize("call", ["1|2:", ":7", "1|2:7:extra"])
def test_empty_sample_values_and_excess_format_fields_are_invalid(tmp_path, call):
    with pytest.raises(PopulationError, match="VCF sample"):
        experiment(tmp_path, GUIDE + "TGG", [
            {"pos1": 23, "ref": "G", "alt": "C,A", "format": "GT:PS", "calls": [call]}
        ], ["SYN"])


@pytest.mark.parametrize("alts,gt,expected", [
    ("AA,AT", "1|2", ["AA", "AT"]),
    ("AA,AT", "2|1", ["AT", "AA"]),
    ("AT,AA", "2|1", ["AA", "AT"]),
])
def test_joint_insertions_keep_original_alt_index_and_phased_copy_order(tmp_path, alts, gt, expected):
    result = experiment(tmp_path, "CA" + GUIDE + "TGG", [
        {"pos1": 2, "ref": "A", "alt": alts, "format": "GT:PS", "calls": [gt + ":7"]}
    ], ["SYN"])
    copies = [row for row in result["materials"] if row["donor_id"] == "SYN"]
    assert [row["sequence_genomic_forward"] for row in copies] == [
        "C" + allele + GUIDE + "TGG" for allele in expected
    ]
    assert result["records"][0]["original_alts"] == alts.split(",")
    assert result["audit"][0]["selected_indexes"] == [int(value) for value in gt.split("|")]
    assert result["audit"][0]["phase_status"] == "phased_explicit_ps"
    assert result["audit"][0]["phase_key"][-2:] == ["explicit", "7"]
    assert all(not row["excluded_intervals"] for row in copies)


@pytest.mark.parametrize("call,status", [("./.", "missing"), ("0/.", "partial_missing")])
def test_missing_calls_remain_unknown_with_omitted_trailing_format(tmp_path, call, status):
    result = experiment(tmp_path, GUIDE + "TGG", [
        {"pos1": 23, "ref": "G", "alt": "C,A", "format": "GT:PS", "calls": [call]}
    ], ["SYN"])
    assert result["audit"][0]["call_status"] == status
    assert result["audit"][0]["selected_indexes"] is None
    assert result["records"][0]["frequency"]["AN"] == 0
    assert not donor_hits(result, "SYN")
    assert all(row["model_carriage"] == "unknown" for row in result["carriage"])


def test_unsupported_symbolic_alt_keeps_original_indexes_and_missingness(tmp_path):
    result = experiment(tmp_path, GUIDE + "TGG", [
        {"pos1": 23, "ref": "G", "alt": "C,<DEL>", "info": "END=23",
         "calls": ["1|1", "2|2", "./."]}
    ], ["LITERAL", "SYMBOLIC", "MISSING"])
    assert result["records"][0]["original_alts"] == ["C", "<DEL>"]
    assert result["records"][0]["edits"][2] is None
    assert result["records"][0]["calls"]["SYMBOLIC"]["selected_indexes"] == [2, 2]
    assert result["records"][0]["calls"]["MISSING"]["call_status"] == "missing"
    assert not donor_hits(result, "SYMBOLIC")
    assert not donor_hits(result, "MISSING")
    assert all(row["model_carriage"] == "unknown" for row in result["carriage"]
               if row["donor_id"] in {"SYMBOLIC", "MISSING"})
