"""Small launcher workflows using the existing preparation and search runner."""
from pathlib import Path
import shutil

from . import lab


def run_reference_example(root, backend, progress=print):
    """Run the bundled mock guides without adding the optional sample."""
    root, backend = Path(root).resolve(), Path(backend).resolve()
    lab.preflight(root.parent)
    if root.exists():
        raise ValueError("This run folder already exists. Choose a new folder; previous results are preserved.")
    root.mkdir(mode=0o700)
    inputs = root / "inputs"
    inputs.mkdir(mode=0o700)
    progress("Preparing the artificial reference example")
    for name in ("reference.fa", "guides.tsv"):
        shutil.copyfile(lab.files() / "examples/small" / name, inputs / name)
    request = lab.base_request(inputs / "reference.fa", inputs / "guides.tsv", backend, ["artificial_1"], small=True)
    config = inputs / "request.json"
    lab.save_json(config, request)
    result = lab.research.run(config, root / "results", progress=progress)
    return root / "results", result


def save_context(root, labels, original_vcf=None):
    """Keep optional experiment context outside the scientific request and results."""
    context = {key: value.strip() for key, value in labels.items() if value.strip()}
    if original_vcf:
        if lab.pin(original_vcf["path"]) == original_vcf:
            context["original_vcf"] = original_vcf
        else:
            context["original_vcf"] = {"path": original_vcf["path"],
                                       "status": "Source changed during the run; use the normalized VCF identity in request.json."}
    if context:
        context["role"] = "User-supplied context only; does not select an editor model or change search settings."
        lab.save_json(Path(root) / "analysis-context.json", context)
