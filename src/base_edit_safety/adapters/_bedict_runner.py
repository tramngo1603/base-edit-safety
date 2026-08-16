"""Runs BE-DICT. This file is executed by BE-DICT's own interpreter, not the pipeline's.

It therefore imports only BE-DICT and the standard library, not this package. It reads a
JSON request, runs the ensemble of models on CPU, and writes a JSON response holding the possible
output sequences and their probabilities.

Usage:  python _bedict_runner.py <bedict_repo_root> <request.json> <response.json>
"""
import json
import os
import sys

# BE-DICT uses aliases removed from the current NumPy release. Define them before importing BE-DICT;
# the change applies only to this process.
import numpy as _np

for _name, _builtin in {"bool": bool, "int": int, "float": float, "object": object, "str": str}.items():
    if not hasattr(_np, _name):
        setattr(_np, _name, _builtin)


def main() -> int:
    repo_root, request_path, response_path = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(request_path) as handle:
        request = json.load(handle)

    editor = request["editor"]
    conversion = (request["edit_from"], request["edit_to"])
    edit_window = tuple(request["edit_window"])
    sequences = request["sequences"]          # list of `seq_id` and `inp_seq` mappings
    runs = int(request.get("runs", 5))

    import pandas as pd
    from haplotype.data_preprocess import SeqProcessConfig, HaplotypeSeqProcessor
    from haplotype.utilities import get_device
    from haplotype.predict_model import BEDICT_HaplotypeModel

    device = get_device(False)  # CPU
    # BE-DICT uses this setting to generate output sequences within the supported edit window.
    seqconfig_dataproc = SeqProcessConfig(20, (1, 20), edit_window, 1)
    processor = HaplotypeSeqProcessor(editor, conversion, seqconfig_dataproc)
    # BE-DICT separately uses the full guide-matched sequence to create its numeric model input.
    seqconfig_datatensgen = SeqProcessConfig(20, (1, 20), (1, 20), 1)
    model = BEDICT_HaplotypeModel(processor, seqconfig_datatensgen, device)

    frame = pd.DataFrame(
        [{"seq_id": s["seq_id"], "Inp_seq": s["inp_seq"]} for s in sequences]
    )
    dloader = model.prepare_data(
        frame, ["seq_id", "Inp_seq"], outpseq_col=None, outcome_col=None,
        renormalize=False, batch_size=500,
    )

    run_frames = []
    for run_num in range(runs):
        model_dir = os.path.join(
            repo_root, "trained_models", "bystander", editor, "train_val", f"run_{run_num}"
        )
        pred = model.predict_from_dloader(dloader, model_dir, outcome_col=None)
        pred["run_num"] = run_num
        run_frames.append(pred)

    aggregated = model.compute_avg_predictions(pd.concat(run_frames, ignore_index=True))

    response = {"editor": editor, "ensemble_runs": runs, "sequences": {}}
    for seq_id, group in aggregated.groupby("seq_id"):
        response["sequences"][seq_id] = [
            {"output_seq": row["Outp_seq"], "pred_score": float(row["pred_score"])}
            for _, row in group.iterrows()
        ]
    with open(response_path, "w") as handle:
        json.dump(response, handle)
    return 0


if __name__ == "__main__":
    sys.exit(main())
