"""The known-case recount preserves original indexes and excludes incomplete calls."""
import pytest

from base_edit_safety.population.common import PopulationError, checksum, read_json
from base_edit_safety.population.known_case import main, protocol, recount


def source_file(tmp_path, calls, extra_row=False, metadata="##reference=GRCh38\n"):
    path = tmp_path/"artificial.vcf"
    names = ["SYN_"+str(i) for i in range(len(calls))]
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n"
    row = "chr2\t210530659\tSYN_example\tG\tT,C\t.\t.\t.\tGT:FT\t" + "\t".join(calls) + "\n"
    path.write_text("##fileformat=VCFv4.3\n" + metadata + header + row * (2 if extra_row else 1))
    return path


def test_carrier_recount_uses_complete_calls_and_original_alt_indexes(tmp_path):
    path = source_file(tmp_path, ["0|2:.", "2/2:PASS", "0/0:.", "2|.:.", "0/2:low_quality", "0|1:."])
    summary, carriers = recount(path, protocol())
    assert (summary["AC"], summary["AN"], summary["carrier_donors"]) == (3, 8, 2)
    assert summary["site_filter_raw"] == "."
    assert summary["excluded_calls"] == {"missing_partial_or_non_diploid_gt": 1, "sample_filter": 1}
    assert [row["source_sample_id"] for row in carriers] == ["SYN_0", "SYN_1"]
    assert [row["original_gt"] for row in carriers] == ["0|2", "2/2"]


@pytest.mark.parametrize("call,duplicate", [("0|3:.", False), ("0|2:.", True)])
def test_carrier_recount_rejects_invalid_indexes_and_duplicate_alleles(tmp_path, call, duplicate):
    with pytest.raises(PopulationError):
        recount(source_file(tmp_path, [call], duplicate), protocol())


@pytest.mark.parametrize("metadata", [
    "##reference=GRCh37\n", "##reference=hg19\n", "##reference=unknown.fa\n",
    "##reference=GRCh38\n##reference=GRCh37\n",
    "##reference=GRCh38\n##contig=<ID=chr2,length=243199373>\n",
    '##reference=GRCh38\n##contig=<ID=chr2,length=242193529,assembly="GRCh37">\n',
])
def test_incompatible_or_unresolved_local_assembly_fails_even_with_a_checksum(tmp_path, metadata):
    source = source_file(tmp_path, ["0|2:."], metadata=metadata)
    output = tmp_path/"result"
    assert main(["--source-vcf", str(source), "--expected-sha256", checksum(source), "--output-dir", str(output)]) == 2
    assert not output.exists()


def test_missing_assembly_and_a_caller_checksum_do_not_establish_source(tmp_path):
    source = source_file(tmp_path, ["0|2:."], metadata="")
    with pytest.raises(PopulationError, match="no recorded acquisition contract"):
        recount(source, protocol())
    assert main(["--source-vcf", str(source), "--expected-sha256", checksum(source), "--output-dir", str(tmp_path/"result")]) == 2
    assert not (tmp_path/"result").exists()


@pytest.mark.parametrize("reference", ["GRCh38", "hg38", "file:///references/GRCh38.p14.fa"])
def test_compatible_local_assembly_does_not_invent_public_source_attribution(tmp_path, reference):
    source = source_file(tmp_path, ["0|2:.", "2/2:PASS", "0|2:.", "0|0:.", "2|.:.", "2:.", "2|2:low_quality"],
                         metadata="##reference="+reference+"\n")
    before = checksum(source)
    output = tmp_path/"result"
    assert main(["--source-vcf", str(source), "--expected-sha256", before, "--output-dir", str(output)]) == 0
    summary = read_json(output/"summary.json")
    assert (summary["AC"], summary["AN"], summary["carrier_donors"]) == (4, 8, 3)
    assert summary["assembly"] == "GRCh38" and summary["assembly_basis"] == "vcf_declaration"
    assert summary["source_verification"] == {"status": "unverified_local_source", "release": None,
        "source_vcf_url": None, "assembly_evidence_url": None}
    assert summary["expected_provenance"]["release"] == protocol()["release"]
    assert "SYN_" not in (output/"summary.json").read_text() and "/references/" not in (output/"summary.json").read_text()
    assert checksum(source) == checksum(output/"source.vcf") == before
    assert all(p.stat().st_mode & 0o777 == 0o600 for p in output.iterdir())


def test_recorded_acquisition_contract_binds_exact_bytes_when_metadata_is_absent(tmp_path):
    source = source_file(tmp_path, ["0|2:."], metadata="")
    # This invented contract tests binding without including any real donor genotypes.
    specification = protocol()
    specification["recorded_acquisition"]["source_subset_sha256"] = checksum(source)
    specification["release"], specification["source_vcf_url"] = "SYN_release", "https://example.invalid/synthetic.vcf"
    summary, _ = recount(source, specification)
    assert summary["assembly_basis"] == "source_documentation_and_acquisition"
    assert summary["source_verification"]["status"] == "recorded_public_acquisition"
    assert summary["source_verification"]["release"] == "SYN_release"
    source.write_text(source.read_text()+"\n")
    with pytest.raises(PopulationError, match="no recorded acquisition contract"):
        recount(source, specification)


@pytest.mark.parametrize("alt,call", [("T,C,C", "0|2:."), ("T,T,C", "0|3:."),
                                      ("T,C", "0|2:.:extra"), ("T,C", "0|2:low_quality:extra")])
def test_duplicate_alt_and_excess_format_values_are_rejected(tmp_path, alt, call):
    source = source_file(tmp_path, [call])
    source.write_text(source.read_text().replace("\tT,C\t", "\t"+alt+"\t"))
    with pytest.raises(PopulationError, match="duplicate or invalid ALT|excess sample FORMAT"):
        recount(source, protocol())
