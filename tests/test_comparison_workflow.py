"""Exercise the offline commands, immutable history and protected failure paths."""
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import shutil
import stat

import pytest

from base_edit_safety.pipeline.comparison import main
from base_edit_safety.population.common import checksum
from base_edit_safety.population.comparison_provenance import resolve_study, stage_key

ROOT = Path(__file__).resolve().parents[1]


def invoke(*arguments, expected=0):
    stdout, stderr = io.StringIO(), io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        status = main([str(argument) for argument in arguments])
    assert status == expected, stderr.getvalue()
    text = stdout.getvalue() + stderr.getvalue()
    assert 'synthetic-donor-C' not in text and 'synthetic-sample-C' not in text
    assert str(ROOT) not in text
    return json.loads(stdout.getvalue() if status == 0 else stderr.getvalue())


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')
    return path


@pytest.fixture(scope='module')
def workflow(tmp_path_factory):
    root = tmp_path_factory.mktemp('comparison-workflow')
    inputs = root / 'inputs'
    shutil.copytree(ROOT / 'examples/comparison/workflow', inputs)
    config = inputs / 'study.json'
    invoke('validate', '--config', config)
    for command, parents in (
        ('search', ()), ('observations', ()),
        ('evidence', ('search', 'observations')), ('compare', ('search', 'evidence')),
        ('questions', ('comparison',)), ('plan', ('questions', 'evidence')),
        ('report', ('comparison', 'questions', 'plan')),
    ):
        arguments = [command, '--config', config, '--output-dir', root / command]
        for parent in parents:
            arguments += ['--' + parent, root / ('compare' if parent == 'comparison' else parent)]
        invoke(*arguments)
    return root


def test_complete_workflow_exports_exact_costs_and_17_unit_shortfall(workflow):
    plan = json.loads((workflow / 'plan/plan.json').read_text())
    assert plan['accounting']['total_cost'] == 18
    assert plan['unweighted_completion'] == 1 and plan['meets_request']
    assert len(plan['selected_action_ids']) == 10
    assert all(row['evidence_status'] == 'no_measurements' for row in plan['questions'].values())
    required = {'comparison-report.md', 'design-comparison.tsv', 'donor-design-matrix.tsv',
        'observations.jsonl', 'evidence-applicability.tsv', 'out-of-scope-evidence.tsv',
        'questions.tsv', 'experimental-plan.tsv', 'resources.tsv', 'unmet-requirements.tsv',
        'plan.json', 'search-audit.json', 'run-manifest.json', 'plan-records.jsonl',
        'report-format.json'}
    assert required <= {path.name for path in (workflow / 'report').iterdir()}
    manifest = json.loads((workflow / 'report/run-manifest.json').read_text())
    assert 'report-format.json' in manifest['outputs']
    report_format = json.loads((workflow / 'report/report-format.json').read_text())
    assert report_format['format_version'] == '2' and report_format['producing_code'] == manifest['code']
    assert {key: row['status'] for key, row in report_format['human_records'].items()} == {
        'reviewer_attestation': 'not_recorded', 'funding_disclosure': 'not_recorded',
        'customer_disclosure': 'not_recorded'}
    report = (workflow / 'report/comparison-report.md').read_text()
    assert 'synthetic-donor-C' not in report and 'synthetic-sample-C' not in report
    assert str(workflow) not in report
    assert '9.8000000000000007' in (workflow / 'report/observations.jsonl').read_text()
    config = workflow / 'inputs/study-budget17.json'
    invoke('plan', '--config', config, '--questions', workflow / 'questions',
           '--evidence', workflow / 'evidence', '--output', workflow / 'plan17')
    short = json.loads((workflow / 'plan17/plan.json').read_text())
    assert not short['meets_request'] and short['unweighted_completion'] == 0
    assert short['status'] == 'infeasible_proven'
    invoke('report', '--config', config, '--comparison', workflow / 'compare',
           '--questions', workflow / 'questions', '--plan', workflow / 'plan17',
           '--output', workflow / 'report17')


def test_reuse_reorder_and_changed_bytes_are_checked(workflow, tmp_path):
    config = workflow / 'inputs/study.json'
    invoke('search', '--config', config, '--output', workflow / 'search', '--reuse')
    original = resolve_study(config)
    reordered = deepcopy(original['config'])
    reordered['catalogs'].reverse()
    reordered['input_refs'].reverse()
    path = write(config.parent / 'reordered.json', reordered)
    assert stage_key(original, 'search') == stage_key(resolve_study(path), 'search')
    invoke('search', '--config', path, '--output', workflow / 'search', '--reuse')
    copied = tmp_path / 'protected-source'
    shutil.copytree(config.parent, copied)
    (copied / 'variants.vcf').write_text((copied / 'variants.vcf').read_text() + '\n')
    result = invoke('validate', '--config', copied / 'study.json', expected=2)
    assert 'checksum' in result['reason'] and str(copied) not in result['reason']


def test_protected_writes_prerequisites_and_no_overwrite(workflow, tmp_path):
    for path in (workflow / 'report').iterdir():
        assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert stat.S_IMODE((workflow / 'report').stat().st_mode) == 0o700
    config = workflow / 'inputs/study.json'
    result = invoke('search', '--config', config, '--output', workflow / 'search', expected=2)
    assert '--reuse' in result['reason']
    invoke('questions', '--config', config, '--comparison', tmp_path / 'missing-private-parent',
           '--output', tmp_path / 'questions', expected=2)
    assert not (tmp_path / 'questions').exists()
    invoke('observations', '--config', config, '--output', workflow / 'inputs', expected=2)
    invoke('search', '--config', config, '--backend', 'private-unrecognized-value',
           '--output', tmp_path / 'search', expected=2)


def test_plan_revision_binds_verified_history_without_rewriting_it(workflow):
    old_config = workflow / 'inputs/study.json'
    old_plan = workflow / 'plan/plan.json'
    old_hash = checksum(old_plan)
    config = json.loads(old_config.read_text())
    config['policies']['planning']['total_budget'] = 17
    for label, source in [('configuration', old_config), ('manifest', workflow / 'plan/run-manifest.json'),
                          ('audit', old_plan)]:
        config['input_refs'].append({'input_id': 'previous-' + label, 'role': 'plan_history',
            'path': str(source), 'sha256': checksum(source), 'format': 'json', 'stages': ['planning']})
    path = write(old_config.parent / 'revision.json', config)
    args = ['plan', '--config', path, '--questions', workflow / 'questions', '--evidence', workflow / 'evidence',
            '--previous-plan', workflow / 'plan', '--previous-config', old_config,
            '--output', workflow / 'revision']
    invoke(*args)
    revised = json.loads((workflow / 'revision/plan.json').read_text())
    original = json.loads(old_plan.read_text())
    assert revised['previous_plan_id'] == original['plan_id']
    assert revised['plan_id'] != original['plan_id']
    assert checksum(old_plan) == old_hash
    invoke(*args, '--reuse')


def test_rehashed_false_plan_is_rejected_before_reporting(workflow):
    folder = workflow / 'false-plan'
    shutil.copytree(workflow / 'plan', folder)
    path = folder / 'plan.json'
    plan = json.loads(path.read_text())
    plan['accounting']['total_cost'] -= 1
    write(path, plan)
    manifest_path = folder / 'run-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['outputs']['plan.json'] = checksum(path)
    write(manifest_path, manifest)
    result = invoke('report', '--config', workflow / 'inputs/study.json',
        '--comparison', workflow / 'compare', '--questions', workflow / 'questions',
        '--plan', folder, '--output', workflow / 'false-report', expected=2)
    assert 'recomputation' in result['reason']
    assert not (workflow / 'false-report').exists()


def test_resource_cap_retains_partial_status_without_completed_manifest(workflow):
    original = workflow / 'inputs/study.json'
    config = json.loads(original.read_text())
    config['policies']['search']['execution_limits'] = {'local': {'max_resident_bytes': 1}}
    path = write(original.parent / 'memory-cap.json', config)
    output = workflow / 'capped'
    result = invoke('search', '--config', path, '--output', output, expected=4)
    assert result['status'] == 'partial' and result['reason_code'] == 'resource_limit'
    audit = json.loads((output / 'stage-status.json').read_text())
    assert audit['execution_audit']['reason_codes'] == ['local_memory_cap']
    assert not (output / 'run-manifest.json').exists()
    invoke('search', '--config', path, '--output', output, '--reuse', expected=2)


def test_atomic_interruption_does_not_publish_complete_stage(workflow, monkeypatch):
    from base_edit_safety.population import common
    output = workflow / 'interrupted'
    original = common.os.replace
    calls = 0
    def interrupt_once(source, target):
        nonlocal calls
        if Path(target) == output and calls == 0:
            calls += 1
            raise KeyboardInterrupt()
        return original(source, target)
    monkeypatch.setattr(common.os, 'replace', interrupt_once)
    result = invoke('observations', '--config', workflow / 'inputs/study.json',
                    '--output', output, expected=4)
    assert result['status'] == 'interrupted'
    assert not (output / 'run-manifest.json').exists()
    assert not (output / 'observations.jsonl').exists()
    assert not list(workflow.glob('.bes-*'))


def _remap(value, mapping):
    if isinstance(value, str):
        return mapping.get(value, value)
    if isinstance(value, list):
        return [_remap(item, mapping) for item in value]
    if isinstance(value, dict):
        return {key: _remap(item, mapping) for key, item in value.items()}
    return value


def _reidentify(config):
    from base_edit_safety.population.comparison_contracts import record_identity
    for _ in range(32):
        mapping = {row['record_id']: record_identity(row) for row in config['catalogs']
                   if row['record_id'] != record_identity(row)}
        if not mapping:
            return config
        config = _remap(config, mapping)
        for row in config['catalogs']:
            if row['record_type'] == 'comparison_design':
                row['definition_hash'] = row['record_id'].rsplit(':', 1)[1]
    raise AssertionError('Changed catalog identities did not converge')


@pytest.mark.parametrize('change,reusable', [
    ('budget', {'search', 'observations', 'evidence', 'compare', 'questions'}),
    ('profile', {'observations'}),
    ('context', {'search'}),
    ('controls', {'search'}),
    ('annotations', {'search', 'observations', 'evidence'}),
    ('report', {'search', 'observations', 'evidence', 'compare', 'questions', 'plan'}),
    ('disabled_prediction', {'search', 'observations', 'evidence', 'compare', 'questions', 'plan', 'report'}),
])
def test_actual_stage_cache_invalidation_matrix(workflow, change, reusable):
    from base_edit_safety.population.comparison_provenance import load_stage
    from base_edit_safety.population.common import PopulationError
    config = json.loads((workflow / 'inputs/study.json').read_text())
    if change == 'budget':
        config['policies']['planning']['total_budget'] = 17
    elif change == 'profile':
        profile = next(row for row in config['catalogs'] if row['record_type'] == 'comparison_search_profile')
        profile['max_mismatches'] = 0
        config = _reidentify(config)
    elif change == 'context':
        context = next(row for row in config['catalogs'] if row['record_type'] == 'comparison_experimental_context')
        context['cell_type'] = 'Another explicitly artificial cell type'
        config = _reidentify(config)
    elif change == 'controls':
        config['policies']['evidence']['require_controls'] = False
    elif change == 'annotations':
        source = write(workflow / 'inputs/annotations.json', {'schema_version': '1', 'assembly': 'SYNTHETIC_TEST',
            'source': 'Artificial annotation', 'release': 'artificial-v1', 'features': []})
        config['input_refs'].append({'input_id': 'artificial-annotation', 'role': 'annotations',
            'path': source.name, 'sha256': checksum(source), 'format': 'json', 'stages': ['comparison']})
    elif change == 'report':
        config['policies']['report']['policy_version'] = 'revised-local-report-policy'
    else:
        config['policies']['predictions']['policy_version'] = 'different-unselected-model-policy'
    path = write(workflow / ('inputs/cache-' + change + '.json'), config)
    study = resolve_study(path)
    roles = {'search': 'search', 'observations': 'observation_normalization', 'evidence': 'evidence_applicability',
             'compare': 'comparison', 'questions': 'questions', 'plan': 'planning', 'report': 'report'}
    for name, role in roles.items():
        if name in reusable:
            load_stage(study, workflow / name, role)
        else:
            with pytest.raises(PopulationError, match='stale'):
                load_stage(study, workflow / name, role)


def test_report_only_policy_change_can_reuse_verified_scientific_plan(workflow):
    config = json.loads((workflow / 'inputs/study.json').read_text())
    config['policies']['report']['policy_version'] = 'changed-presentation-only'
    path = write(workflow / 'inputs/report-only.json', config)
    invoke('report', '--config', path, '--comparison', workflow / 'compare',
           '--questions', workflow / 'questions', '--plan', workflow / 'plan',
           '--output', workflow / 'report-new-presentation')
