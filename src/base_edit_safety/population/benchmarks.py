"""Compare selectors on one frozen universe without treating nominations as labels."""
from __future__ import annotations

from copy import deepcopy

from .common import PopulationError, record
from .contracts import validate_benchmark_policy
from .planning import eligible_actions, plan


STRATA = {
    "unknown_denominator": "AN = 0",
    "not_observed": "AC = 0 with AN > 0",
    "singleton": "AC = 1",
    "rare_below_1_percent": "AC >= 2 and AC/AN < 0.01",
    "frequency_1_to_5_percent": "0.01 <= AC/AN < 0.05, excluding singletons",
    "frequency_at_least_5_percent": "AC/AN >= 0.05, excluding singletons",
    "no_linked_retained_variant": "Configuration has no linked retained variant",
}


def variant_strata(variants, candidates, units):
    """Use fixed descriptive bins from called alleles; never interpret them as prevalence."""
    by_variant, counts = {}, {s: 0 for s in STRATA if s != "no_linked_retained_variant"}
    for variant in variants:
        f = variant["frequency"]
        strata = set()
        for ac in f["AC"]:
            an = f["AN"]
            label = "unknown_denominator" if not an else "not_observed" if not ac else "singleton" if ac == 1 else (
                "rare_below_1_percent" if 100*ac < an else "frequency_1_to_5_percent" if 20*ac < an else "frequency_at_least_5_percent")
            counts[label] += 1
            strata.add(label)
        by_variant[variant["variant_id"]] = strata
    by_candidate = {c["configuration_id"]: c for c in candidates}
    by_unit = {}
    for unit in units:
        c = by_candidate[unit["configuration_id"]]
        strata = set().union(*(by_variant[v] for v in c["variant_links"]))
        by_unit[unit["id"]] = strata or {"no_linked_retained_variant"}
    return {"status": "reported", "definitions": STRATA, "retained_alt_alleles": counts,
        "scope": "Accepted complete diploid calls in the declared donor subset. Linked-variant strata may overlap across units; "
                 "a link does not establish that each ALT creates the configuration. Source filtering limits rare-allele completeness."}, by_unit


def benchmark(config, request, run_id, *, variants=None, candidates=None):
    policy = config["benchmark_policy"]
    validate_benchmark_policy(policy)
    actionable = {u for a in eligible_actions(request) for u in a["covered_units"]}
    strata, unit_strata = ({"status": "source_artifacts_not_supplied"}, {}) if variants is None or candidates is None else \
        variant_strata(variants, candidates, request["units"])
    rows = []
    for total_budget in policy["budget_grid"]:
        if total_budget < config["planning_policy"]["fixed_overhead"]:
            raise PopulationError("Benchmark budget is below fixed study overhead")
        r = deepcopy(request)
        r["budget"]["max_cost"] = total_budget - config["planning_policy"]["fixed_overhead"]
        for algorithm in policy["selectors"]:
            for seed in policy["seeds"] if algorithm == "random_order" else [None]:
                result, order = plan(r, algorithm, seed)
                rows.append({"algorithm": algorithm, "seed": seed, "visit_order": order,
                    "study_budget": total_budget, "activation_and_assay_cost": result["total_cost"],
                    "fixed_overhead": config["planning_policy"]["fixed_overhead"],
                    "unused_budget": result["unused_budget"], "plan_status": result["plan_status"],
                    "objective_numerator": result["objective_numerator"], "objective_denominator": result["objective_denominator"],
                    "all_units": len(request["units"]), "actionable_units": len(actionable),
                    "fully_satisfied_all_units": len(result["fully_covered_units"]),
                    "fully_satisfied_actionable_units": len(set(result["fully_covered_units"]) & actionable),
                    "eligible_donors": len({a["donor_id"] for a in eligible_actions(r)}),
                    "excluded_donors": len(request["donors"]) - len({a["donor_id"] for a in eligible_actions(r)}),
                    "selected_actions": result["selected_actions"],
                    "variant_strata": [{"stratum": label,
                        "all_units": sum(label in bins for bins in unit_strata.values()),
                        "actionable_units": sum(label in unit_strata[u] for u in actionable),
                        "fully_satisfied_all_units": sum(label in unit_strata[u] for u in result["fully_covered_units"]),
                        "fully_satisfied_actionable_units": sum(label in unit_strata[u] for u in set(result["fully_covered_units"]) & actionable)}
                        for label in STRATA] if unit_strata or strata["status"] == "reported" else []})
    expected_rows = len(policy["budget_grid"]) * sum(len(policy["seeds"]) if a == "random_order" else 1 for a in policy["selectors"])
    if not rows or len(rows) != expected_rows:
        raise PopulationError("Benchmark did not produce every requested selector/budget/seed comparison")
    return record("benchmark_summary", run_id, mode=config["mode"], coverage_basis=config["coverage_basis"],
                  universe_id=request["universe_id"], split="transductive", rows=rows,
                  interpretation="Planned configuration coverage; no empirical sensitivity or editing outcome is inferred.",
                  rare_variant_strata=strata)
