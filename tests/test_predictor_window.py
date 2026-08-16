"""Share the BE-DICT edit window without changing the full-sequence model input.

BE-DICT's example configuration gives two `SeqProcessConfig` objects different jobs. The
`HaplotypeSeqProcessor` uses positions 3-10 to generate possible edited sequences, while the
`BEDICT_HaplotypeModel` uses the full guide-matched sequence as numeric input. The adapter passes the
3-10 window only to the first step.
"""
from __future__ import annotations

import ast
import pathlib

from base_edit_safety.adapters.bedict import BedictBystanderPredictor
from base_edit_safety.config import BedictConfig
from base_edit_safety.pipeline.datatypes import EDIT_WINDOW, Guide

RUNNER = (pathlib.Path(__file__).parent.parent
          / "src" / "base_edit_safety" / "adapters" / "_bedict_runner.py")


def _tree() -> ast.AST:
    return ast.parse(RUNNER.read_text())


def _calls(function_name: str) -> list[ast.Call]:
    return [node for node in ast.walk(_tree())
            if isinstance(node, ast.Call) and getattr(node.func, "id", None) == function_name]


def _assigned_value(name: str) -> ast.expr:
    for node in ast.walk(_tree()):
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
            return node.value
    raise AssertionError(f"no assignment to {name!r} found in {RUNNER.name}")


def _literal(node):
    """The literal value of an AST node, or None if it is not a literal."""
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError, SyntaxError):
        return None


GUIDE = Guide(
    guide_id="window_test",
    protospacer="AAAAAAAAAAAAAAAAAAAA",
    pam="AGG",
    strand="+",
    genome_build="GRCh38",
    chrom="chr11",
    gene="HBG1",
    gene_strand="-",
    tss_1based=100,
    editor_class="ABE",
    editor_name="ABE8e",
    edit_from="A",
    edit_to="G",
    proto_5p_genomic=100,
)


def test_adapter_request_carries_edit_window_bounds(monkeypatch):
    """Write the `EDIT_WINDOW` bounds into the external-tool request."""
    captured = []
    predictor = BedictBystanderPredictor(
        BedictConfig(repo_root="unused", python_executable="unused"), runs=2,
    )

    def fake_run(request):
        captured.append(request)
        return {
            "ensemble_runs": 2,
            "sequences": {
                GUIDE.guide_id: [{"output_seq": GUIDE.protospacer, "pred_score": 1.0}],
            },
        }

    monkeypatch.setattr(predictor, "_run", fake_run)
    predictor.predict(GUIDE)

    assert captured[0]["edit_window"] == [min(EDIT_WINDOW), max(EDIT_WINDOW)]


def test_runner_passes_the_requested_window_to_the_processor():
    """The external runner reads the request instead of repeating the 3-10 literal."""
    request_value = _assigned_value("edit_window")
    assert isinstance(request_value, ast.Call)
    assert getattr(request_value.func, "id", None) == "tuple"
    request_lookup = request_value.args[0]
    assert isinstance(request_lookup, ast.Subscript)
    assert getattr(request_lookup.value, "id", None) == "request"
    assert _literal(request_lookup.slice) == "edit_window"

    dataproc = _assigned_value("seqconfig_dataproc")
    assert isinstance(dataproc, ast.Call)
    assert getattr(dataproc.func, "id", None) == "SeqProcessConfig"
    assert getattr(dataproc.args[2], "id", None) == "edit_window"

    processor = _calls("HaplotypeSeqProcessor")
    assert len(processor) == 1
    assert getattr(processor[0].args[2], "id", None) == "seqconfig_dataproc"


def test_processor_window_and_tensor_configuration_are_intentionally_distinct():
    """Keep BE-DICT's full-sequence input separate from its edited-sequence window."""
    configs = _calls("SeqProcessConfig")
    assert len(configs) == 2

    datatensgen = _assigned_value("seqconfig_datatensgen")
    assert isinstance(datatensgen, ast.Call)
    assert getattr(datatensgen.func, "id", None) == "SeqProcessConfig"
    assert _literal(datatensgen.args[2]) == (1, 20)

    model = _calls("BEDICT_HaplotypeModel")
    assert len(model) == 1
    assert getattr(model[0].args[1], "id", None) == "seqconfig_datatensgen"


def test_edit_window_is_what_the_pipeline_thinks_it_is():
    """Do not widen predictions beyond BE-DICT's supported processor window."""
    assert list(EDIT_WINDOW) == [3, 4, 5, 6, 7, 8, 9, 10]
