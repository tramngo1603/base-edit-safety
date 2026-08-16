"""The helpers that record what a run depended on: checksums, cache state, tool paths."""
import datetime
import io
import json

from base_edit_safety.io import reference
from base_edit_safety.pipeline import environment


def test_sequence_checksum_stable_and_sensitive():
    c = reference.sequence_checksum("ACGTACGT")
    assert c == reference.sequence_checksum("ACGTACGT") and len(c) == 16
    assert reference.sequence_checksum("ACGTACGA") != c


def test_window_provenance_reports_uncached():
    prov = reference.window_provenance("chr11", 999_000, 999_130)  # not in the cache
    assert prov["cached"] is False
    assert prov["cache_key"] == "hg38_chr11_999000_999130.json"


def test_fetched_window_cache_and_provenance_carry_fetch_time(tmp_path, monkeypatch):
    monkeypatch.setattr(reference, "_CACHE_DIR", str(tmp_path))

    def fake_urlopen(url, timeout):
        assert "genome=hg38;chrom=chr11;start=99;end=104" in url
        assert timeout == 7
        return io.BytesIO(b'{"dna": "acgta"}')

    monkeypatch.setattr(reference.urllib.request, "urlopen", fake_urlopen)
    window = reference.fetch_window("chr11", 100, 104, timeout=7)

    cache_path = tmp_path / "hg38_chr11_100_104.json"
    payload = json.loads(cache_path.read_text())
    assert payload["fetched_utc"] == window.fetched_utc
    fetched = datetime.datetime.fromisoformat(payload["fetched_utc"])
    assert fetched.tzinfo is not None
    assert fetched.utcoffset() == datetime.timedelta(0)

    provenance = reference.window_provenance("chr11", 100, 104)
    assert provenance["fetched_utc"] == payload["fetched_utc"]
    assert provenance["checksum_sha256_16"] == reference.sequence_checksum("ACGTA")


def test_tool_environment_drops_missing_and_keeps_set_fields():
    env = environment.tool_environment("BE-DICT", repo_root=None, python_executable="/env/py")
    assert env["name"] == "BE-DICT" and env["python_executable"] == "/env/py"
    assert "repo_root" not in env and "genome_index" not in env  # empty fields are left out


def test_git_commit_never_raises():
    assert environment.git_commit(None) is None
    assert environment.git_commit("/no/such/path/xyz") is None


def test_failed_commit_lookup_is_written_down_not_dropped():
    env = environment.tool_environment("CRISPOR", repo_root="/no/such/repo/xyz", python_executable="/x/py")
    assert env["repo_root"] == "/no/such/repo/xyz"
    assert env["python_executable"] == "/x/py"
    assert env["commit"] is None and env["commit_lookup"] == "unavailable"


def test_reproducibility_block_shape():
    block = environment.reproducibility(
        tools=[{"name": "CRISPOR"}], editors={"ABE": "ABE8e"}, ensemble_runs=5,
        reference_windows={"HBG1": {"cached": False}},
    )
    assert block["genome_build"] == "GRCh38"
    assert block["ensemble_runs"] == 5
    assert block["tools"][0]["name"] == "CRISPOR"
    assert block["reference_windows"]["HBG1"]["cached"] is False
