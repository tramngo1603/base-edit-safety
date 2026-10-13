"""Coverage follows actual carrier sequences and stays within one search."""
from itertools import combinations
from types import SimpleNamespace

from base_edit_safety.population.experiment_planning import _Coverage


def catalog():
    configurations = [
        {"record_type": "comparison_design_configuration", "design_configuration_id": identity,
         "design_id": design, "sequence_configuration_id": sequence, "role": role}
        for identity, design, sequence, role in [
            ("off-a", "design-a", "carrier-sequence", "off_target"),
            ("off-b", "design-b", "carrier-sequence", "off_target"),
            ("second-a", "design-a", "second-sequence", "off_target"),
            ("intended-a", "design-a", "carrier-sequence", "intended_target"),
            ("unresolved-a", "design-a", "unresolved-sequence", "off_target"),
        ]]
    carriers = [
        {"record_type": "comparison_donor_configuration", "donor_id": donor,
         "design_configuration_id": identity, "carrier_status": status}
        for donor, identity, status in [
            ("carrier", "off-a", "supported_carrier"),
            ("carrier", "off-b", "modeled_carrier"),
            ("carrier", "second-a", "modeled_carrier"),
            ("carrier", "intended-a", "supported_carrier"),
            ("carrier", "unresolved-a", "unresolved"),
            ("noncarrier", "off-a", "supported_noncarrier"),
        ]]
    actions = {}
    for aid, donor, design, sequence, operation in [
        ("single-a", "carrier", "design-a", "carrier-sequence", "locus_assay"),
        ("single-b", "carrier", "design-b", "carrier-sequence", "locus_assay"),
        ("panel", "carrier", "design-a", "carrier-sequence", "measurement_panel"),
        ("noncarrier", "noncarrier", "design-a", "carrier-sequence", "locus_assay"),
        ("wrong-sequence", "carrier", "design-a", "noncarrier-sequence", "locus_assay"),
        ("unresolved", "carrier", "design-a", "unresolved-sequence", "locus_assay"),
        ("control", "carrier", "design-a", "carrier-sequence", "delivery_control"),
    ]:
        actions[aid] = {"donor_id": donor, "design_id": {"value": design},
                        "target_sequence_id": {"value": sequence}, "operation_type": operation}
    actions["panel"]["target_sequences"] = [{"sequence_configuration_ids": ["second-sequence"]}]
    return SimpleNamespace(records=configurations + carriers, actions=actions)


def test_coverage_cache_preserves_actual_sequence_and_design_membership():
    qset = catalog()
    coverage = _Coverage(qset)
    expected = {"single-a": {"off-a"}, "single-b": {"off-b"},
                "panel": {"off-a", "second-a"}, "noncarrier": set(),
                "wrong-sequence": set(), "unresolved": set(), "control": set()}
    for size in range(len(expected) + 1):
        for selected in combinations(expected, size):
            assert coverage.units(selected) == set().union(*(expected[aid] for aid in selected))
    assert len(coverage.by_action) == len(expected)
    returned = coverage.units(["panel"])
    returned.clear()
    assert coverage.units(["panel"]) == {"off-a", "second-a"}


def test_coverage_cache_does_not_scan_records_again_or_cross_snapshots():
    class CountedRows(list):
        reads = 0

        def __iter__(self):
            self.reads += 1
            return super().__iter__()

    first = catalog()
    first.records = CountedRows(first.records)
    coverage = _Coverage(first)
    for _ in range(3):
        assert coverage.units(["single-a", "panel"]) == {"off-a", "second-a"}
    assert first.records.reads == 1
    second = catalog()
    second.actions["panel"]["donor_id"] = "noncarrier"
    assert _Coverage(second).units(["panel"]) == set()
    assert coverage.units(["panel"]) == {"off-a", "second-a"}
