"""Check mock evaluation results with grouped splits and explicit assay denominators."""
from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import math
import random
import time

from .common import (PopulationError, checksum, digest, fields, integer, json_text,
                     nonempty, read_json, unique, write_bundle)
from .provenance import runtime
from .comparison_outcomes import OUTCOMES

MINIMUM_CLUSTERS = 3


def number(value, label, *, allow_float=False):
    if isinstance(value, float) and not allow_float:
        raise PopulationError(label + " must use an integer or decimal string; fractional JSON numbers are not exact")
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise PopulationError(label + " must be a finite number")
    try:
        result = Decimal(str(value))
    except InvalidOperation as exc:
        raise PopulationError(label + " must be a finite number") from exc
    if not result.is_finite():
        raise PopulationError(label + " must be a finite number")
    return result


def candidates_index(rows):
    if not isinstance(rows, list) or not rows:
        raise PopulationError("Candidate metadata must be a nonempty list")
    for row in rows:
        fields(row, ("candidate_id", "study_id", "guide_id", "family_id", "editor", "assay", "eligible"))
        for key in ("candidate_id", "study_id", "guide_id", "family_id", "editor", "assay"):
            nonempty(row[key], key)
            if row[key].strip().lower() in {"unknown", "none", "null", "na", "n/a", "missing"}:
                raise PopulationError("Missing identity/context cannot define evaluation groups")
        if type(row["eligible"]) is not bool:
            raise PopulationError("eligible must be boolean")
    return unique(rows, "candidate_id")


def dependency_groups(rows):
    """Connect shared studies, guides or families before any exclusion."""
    index = candidates_index(rows)
    parent = {key: key for key in index}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    seen = {}
    for key in sorted(index):
        for field in ("study_id", "guide_id", "family_id"):
            identity = (field, index[key][field])
            if identity in seen:
                parent[find(key)] = find(seen[identity])
            else:
                seen[identity] = key
    groups = {}
    for key in sorted(index):
        groups.setdefault(find(key), []).append(key)
    return {digest(members): members for members in sorted(groups.values())}


def make_splits(rows, seed=17, training_fraction="0.6", development_fraction="0.2"):
    integer(seed, "split seed")
    train = Fraction(number(training_fraction, "training fraction"))
    dev = Fraction(number(development_fraction, "development fraction"))
    if min(train, dev) < 0 or train + dev > 1:
        raise PopulationError("Split fractions must be nonnegative and sum to at most one")
    groups = dependency_groups(rows)
    order = sorted(groups, key=lambda group: digest([seed, group]))
    ntrain, ndev = int(len(order) * train), int(len(order) * dev)
    partitions = ["training"] * ntrain + ["development"] * ndev
    partitions += ["evaluation"] * (len(order) - len(partitions))
    return {"format_version": "1", "candidate_metadata_sha256": digest(sorted(rows, key=lambda r: r["candidate_id"])),
            "seed": seed, "training_fraction": str(number(training_fraction, "training fraction")),
            "development_fraction": str(number(development_fraction, "development fraction")),
            "groups": [{"group_id": group, "candidate_ids": groups[group], "partition": partition}
                       for group, partition in zip(order, partitions)]}


def validate_splits(rows, manifest):
    fields(manifest, ("format_version", "candidate_metadata_sha256", "seed", "training_fraction",
                      "development_fraction", "groups"))
    expected = make_splits(rows, manifest["seed"], manifest["training_fraction"], manifest["development_fraction"])
    if manifest != expected:
        raise PopulationError("Split manifest changed, leaks linked groups, or belongs to different metadata")
    return {key: (group["group_id"], group["partition"]) for group in manifest["groups"]
            for key in group["candidate_ids"]}


def validate_inputs(rows, labels, methods, protocol, splits):
    index = candidates_index(rows)
    assignments = validate_splits(rows, splits)
    fields(protocol, ("format_version", "data_kind", "protocol_id", "budget_per_editor_assay", "targets",
                      "bootstrap_seed", "bootstrap_replicates", "baseline_method"))
    if protocol["format_version"] != "1" or protocol["data_kind"] != "artificial":
        raise PopulationError("This draft harness accepts artificial evaluation only")
    nonempty(protocol["protocol_id"], "protocol_id")
    integer(protocol["budget_per_editor_assay"], "budget", maximum=1000000)
    integer(protocol["bootstrap_seed"], "bootstrap seed")
    integer(protocol["bootstrap_replicates"], "bootstrap replicates", minimum=100, maximum=10000)
    if not isinstance(protocol["targets"], list) or not protocol["targets"]:
        raise PopulationError("At least one matched-sensitivity target is required")
    targets = [number(x, "sensitivity target") for x in protocol["targets"]]
    if any(x < 0 or x > 1 for x in targets) or len(set(targets)) != len(targets):
        raise PopulationError("Sensitivity targets must be unique and between zero and one")
    if not isinstance(labels, list):
        raise PopulationError("Labels must be a list")
    for row in labels:
        fields(row, ("candidate_id", "adjudicated_outcome", "evidence_reference", "adjudication_reason"))
        if row["adjudicated_outcome"] not in OUTCOMES:
            raise PopulationError("Invalid adjudicated outcome")
        nonempty(row["evidence_reference"], "evidence reference")
        nonempty(row["adjudication_reason"], "adjudication reason")
    label_index = unique(labels, "candidate_id")
    if label_index.keys() != index.keys():
        raise PopulationError("Every candidate needs one explicit label, including untested candidates")
    if not isinstance(methods, list) or len(methods) < 2:
        raise PopulationError("At least two methods, including the baseline, are required")
    for method in methods:
        fields(method, ("method_id", "enumerated_ids", "scores", "runtime_seconds"))
        nonempty(method["method_id"], "method_id")
        elapsed = number(method["runtime_seconds"], "method runtime", allow_float=True)
        if elapsed < 0 or not math.isfinite(float(elapsed)):
            raise PopulationError("Method runtime must be a finite nonnegative number")
        enumerated = method["enumerated_ids"]
        if not isinstance(enumerated, list) or any(not isinstance(x, str) for x in enumerated):
            raise PopulationError("enumerated_ids must be a list of candidate identifiers")
        if len(set(enumerated)) != len(enumerated) or not set(enumerated) <= index.keys():
            raise PopulationError("Enumeration has duplicate or unknown candidates")
        if not isinstance(method["scores"], list):
            raise PopulationError("Scores must be a list")
        for score in method["scores"]:
            fields(score, ("candidate_id", "score"))
            number(score["score"], "score")
        scores = unique(method["scores"], "candidate_id")
        if not scores.keys() <= set(enumerated):
            raise PopulationError("Scored candidates must be enumerated by that method")
    method_index = unique(methods, "method_id")
    if protocol["baseline_method"] not in method_index:
        raise PopulationError("Baseline method is missing")
    return index, label_index, assignments


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def recall_difference(current, baseline):
    """Subtract selected positive counts before dividing by their shared denominator."""
    return ratio(current["at_budget"]["positive"] - baseline["at_budget"]["positive"],
                 current["eligible_positives"])


def metrics(ids, labels, method, budget, targets, weights=None, curve=False):
    """Rerank the resampled population; repeated copies consume assay units."""
    weights = weights if weights is not None else {key: 1 for key in ids}
    count = lambda keys: sum(weights.get(key, 0) for key in keys)
    positives = {key for key in ids if labels[key]["adjudicated_outcome"] == "observed_positive"}
    positive_count = count(positives)
    enumerated = set(method["enumerated_ids"]) & ids
    scored = {row["candidate_id"] for row in method["scores"]} & ids
    ranking = sorted((row for row in method["scores"] if row["candidate_id"] in ids),
                     key=lambda row: (number(row["score"], "score").copy_negate(), row["candidate_id"]))
    ranked = [row["candidate_id"] for row in ranking for _ in range(weights.get(row["candidate_id"], 0))]

    def at(k):
        selected = Counter(labels[key]["adjudicated_outcome"] for key in ranked[:k])
        positive = selected["observed_positive"]
        conclusive = positive + selected["tested_no_signal"]
        return {"assays_used": min(k, len(ranked)), "positive": positive,
                "tested_no_signal": selected["tested_no_signal"], "inconclusive": selected["inconclusive"],
                "untested": selected["untested"], "conclusive_count": conclusive,
                "recall": ratio(positive, positive_count), "precision_conclusive": ratio(positive, conclusive),
                "conclusive_outcome_coverage": ratio(conclusive, min(k, len(ranked)))}

    matched = {}
    positions = [i for i, key in enumerate(ranked, 1) if key in positives]
    for target in targets:
        numeric = number(target, "target")
        required = math.ceil(Fraction(numeric) * positive_count)
        if numeric == 0:
            status, cost = "reached", 0
        elif not positive_count:
            status, cost = "undefined_no_eligible_positives", None
        elif len(positions) < required:
            status, cost = "unreachable", None
        else:
            status, cost = "reached", positions[required - 1]
        matched[str(target)] = {"status": status, "required_positives": required, "assays": cost}
    result = {"eligible_candidates": count(ids), "eligible_positives": positive_count,
              "enumerated_candidates": count(enumerated), "enumerated_positives": count(enumerated & positives),
              "enumeration_sensitivity": ratio(count(enumerated & positives), positive_count),
              "missed_enumeration_positives": count(positives - enumerated),
              "enumerated_unscored_candidates": count(enumerated - scored),
              "enumerated_unscored_positives": count((enumerated - scored) & positives),
              "ranked_candidates": len(ranked), "budget": budget, "budget_unfilled": max(0, budget - len(ranked)),
              "at_budget": at(budget), "matched_sensitivity": matched}
    if curve:
        # Accumulate once rather than rescanning every prefix.
        cumulative = Counter()
        result["precision_recall"] = [at(0)]
        for i, key in enumerate(ranked, 1):
            cumulative[labels[key]["adjudicated_outcome"]] += 1
            positive = cumulative["observed_positive"]
            conclusive = positive + cumulative["tested_no_signal"]
            result["precision_recall"].append({"assays_used": i, "positive": positive,
                "conclusive_count": conclusive, "inconclusive": cumulative["inconclusive"],
                "untested": cumulative["untested"], "tested_no_signal": cumulative["tested_no_signal"],
                "recall": ratio(positive, positive_count), "precision_conclusive": ratio(positive, conclusive),
                "conclusive_outcome_coverage": ratio(conclusive, i)})
    return result


def interval(values, independent_clusters, *, original_defined=True):
    valid = sorted(x for x in values if x is not None)
    result = {"independent_clusters": independent_clusters, "replicates": len(values),
              "defined_replicates": len(valid), "undefined_replicates": len(values) - len(valid),
              "lower": None, "upper": None, "reason": None}
    if not original_defined:
        result["reason"] = "original_metric_undefined_or_unreachable"
    elif independent_clusters < MINIMUM_CLUSTERS:
        result["reason"] = "fewer_than_three_relevant_dependency_groups"
    elif not valid or len(valid) != len(values):
        result["reason"] = "undefined_replicates_present"
    else:
        # Nearest-rank percentile interval, with no interpolation.
        result["lower"] = valid[max(0, math.ceil(0.025 * len(valid)) - 1)]
        result["upper"] = valid[math.ceil(0.975 * len(valid)) - 1]
    return result


def evaluate(rows, labels, methods, protocol, splits):
    index, labels, assignments = validate_inputs(rows, labels, methods, protocol, splits)
    eligible = {key for key, row in index.items() if row["eligible"] and assignments[key][1] == "evaluation"}
    if not eligible:
        raise PopulationError("No eligible candidates in the evaluation partition")
    groups = sorted({assignments[key][0] for key in eligible})
    rng = random.Random(protocol["bootstrap_seed"])
    draws = [Counter(rng.choices(groups, k=len(groups))) for _ in range(protocol["bootstrap_replicates"])]
    budget, targets = protocol["budget_per_editor_assay"], protocol["targets"]
    strata = sorted({(index[key]["editor"], index[key]["assay"]) for key in eligible})
    report = {"format_version": "3", "data_kind": "artificial", "protocol_id": protocol["protocol_id"],
              "assay_count_basis": "simulated_testing_slots_including_currently_untested",
              "completed_assays": None, "completed_assays_reason": "not_recorded_in_evaluation_inputs",
              "scope": "software demonstration; no measured accuracy or assay savings claim",
              "candidate_metadata_sha256": splits["candidate_metadata_sha256"], "split_sha256": digest(splits),
              "protocol_sha256": digest(protocol), "label_sha256": digest(list(labels.values())),
              "methods_sha256": digest(methods), "evaluation_candidates": len(eligible),
              "excluded_candidates": len(index) - len(eligible), "evaluation_dependency_groups": len(groups),
              "bootstrap": {"seed": protocol["bootstrap_seed"], "replicates": len(draws),
                  "draw_sha256": digest([dict(sorted(draw.items())) for draw in draws]),
                  "unit": "transitive shared study/guide/family component", "minimum_relevant_groups": MINIMUM_CLUSTERS,
                  "method": "paired cluster bootstrap; rerank repeated copies at fixed budget; nearest-rank 95% percentile"},
              "strata": []}
    for editor, assay in strata:
        ids = {key for key in eligible if (index[key]["editor"], index[key]["assay"]) == (editor, assay)}
        positive_groups = {assignments[key][0] for key in ids if labels[key]["adjudicated_outcome"] == "observed_positive"}
        outputs, resamples = {}, {}
        for method in sorted(methods, key=lambda row: row["method_id"]):
            name = method["method_id"]
            outputs[name] = metrics(ids, labels, method, budget, targets, curve=True)
            outputs[name]["declared_method_runtime_seconds_total"] = float(number(method["runtime_seconds"], "runtime", allow_float=True))
            resamples[name] = [metrics(ids, labels, method, budget, targets,
                {key: draw[assignments[key][0]] for key in ids}) for draw in draws]
            selected_order = sorted((row for row in method["scores"] if row["candidate_id"] in ids),
                                    key=lambda row: (number(row["score"], "score").copy_negate(), row["candidate_id"]))[:budget]
            conclusive_groups = {assignments[row["candidate_id"]][0] for row in selected_order
                if labels[row["candidate_id"]]["adjudicated_outcome"] in {"observed_positive", "tested_no_signal"}}
            outputs[name]["uncertainty"] = {
                "enumeration_sensitivity": interval([r["enumeration_sensitivity"] for r in resamples[name]], len(positive_groups)),
                "recall_at_budget": interval([r["at_budget"]["recall"] for r in resamples[name]], len(positive_groups)),
                "precision_conclusive_at_budget": interval([r["at_budget"]["precision_conclusive"] for r in resamples[name]], len(conclusive_groups)),
                "assays_at_sensitivity": {str(target): interval([r["matched_sensitivity"][str(target)]["assays"]
                    for r in resamples[name]], len(positive_groups),
                    original_defined=outputs[name]["matched_sensitivity"][str(target)]["assays"] is not None) for target in targets}}
        baseline = protocol["baseline_method"]
        comparisons = {}
        recall_comparisons = {}
        for name in outputs:
            if name == baseline:
                continue
            current, reference = outputs[name], outputs[baseline]
            difference = recall_difference(current, reference)
            recall_comparisons[name] = {"baseline_method": baseline,
                "status": "comparable" if difference is not None else "not_comparable",
                "budget": budget, "eligible_positives": current["eligible_positives"],
                "method_recall": current["at_budget"]["recall"],
                "baseline_recall": reference["at_budget"]["recall"],
                "recall_difference_method_minus_baseline": difference,
                "uncertainty": interval([recall_difference(a, b)
                    for a, b in zip(resamples[name], resamples[baseline])],
                    len(positive_groups), original_defined=difference is not None)}
            comparisons[name] = {}
            for target in targets:
                key = str(target)
                left = outputs[name]["matched_sensitivity"][key]
                right = outputs[baseline]["matched_sensitivity"][key]
                comparable = left["assays"] is not None and right["assays"] is not None
                differences = []
                for current, reference in zip(resamples[name], resamples[baseline]):
                    a = current["matched_sensitivity"][key]["assays"]
                    b = reference["matched_sensitivity"][key]["assays"]
                    differences.append(a - b if a is not None and b is not None else None)
                comparisons[name][key] = {"baseline_method": baseline,
                    "status": "comparable" if comparable else "not_comparable",
                    "method_status": left["status"], "baseline_status": right["status"],
                    "assay_difference_method_minus_baseline": left["assays"] - right["assays"] if comparable else None,
                    "uncertainty": interval(differences, len(positive_groups), original_defined=comparable)}
        report["strata"].append({"editor": editor, "assay": assay,
            "dependency_groups": len({assignments[key][0] for key in ids}),
            "positive_dependency_groups": len(positive_groups), "methods": outputs,
            "recall_at_budget_comparisons": recall_comparisons,
            "matched_sensitivity_comparisons": comparisons})
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    split = commands.add_parser("split", help="Make label-free connected-group assignments")
    split.add_argument("--candidates", required=True)
    split.add_argument("--seed", type=int, default=17)
    split.add_argument("--training-fraction", default="0.6")
    split.add_argument("--development-fraction", default="0.2")
    run = commands.add_parser("evaluate", help="Evaluate artificial frozen labels and method outputs")
    for flag in ("candidates", "labels", "methods", "protocol", "splits"):
        run.add_argument("--" + flag, required=True)
    for command in (split, run):
        command.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    started = time.perf_counter()
    try:
        roles = ("candidates",) if args.command == "split" else ("candidates", "labels", "methods", "protocol", "splits")
        inputs = {role: getattr(args, role) for role in roles}
        values = {role: read_json(path) for role, path in inputs.items()}
        if args.command == "split":
            result = make_splits(values["candidates"], args.seed, args.training_fraction, args.development_fraction)
            name = "splits.json"
        else:
            result = evaluate(values["candidates"], values["labels"], values["methods"], values["protocol"], values["splits"])
            name = "evaluation.json"
        receipt = {"runtime": runtime(), "elapsed_seconds": time.perf_counter() - started,
                   "inputs": {role: {"path": path, "sha256": checksum(path)} for role, path in inputs.items()},
                   "semantic_result_sha256": digest(result)}
        write_bundle(args.output, {name: json_text(result), "receipt.json": json_text(receipt)}, inputs=inputs.values())
    except (PopulationError, OSError) as exc:
        parser.exit(2, str(exc) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
