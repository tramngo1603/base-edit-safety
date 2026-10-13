"""Whole questions require actual donors, complete panels and independent repeats."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from base_edit_safety.population.common import PopulationError, ResourceLimit, known, unknown
from base_edit_safety.population.comparison_contracts import make_record, record_identity, record_schema, validate_catalog, validate_record
from base_edit_safety.population.comparison_provenance import resolve_study, effective_inputs
from base_edit_safety.population.comparison_questions import (prepare_question_set,evaluate_requirements,eligible_donor_bundles,generate_questions,dependency_closure)

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('question_fixture',ROOT/'examples/comparison/questions/build_fixture.py')
fixture_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixture_module)
build_fixture=fixture_module.build_fixture


def replace(value, mapping):
    if isinstance(value,str):return mapping.get(value,value)
    if isinstance(value,list):return [replace(item,mapping) for item in value]
    if isinstance(value,dict):return {key:replace(item,mapping) for key,item in value.items()}
    return value


def mutated(bundle, change):
    config=deepcopy(bundle['study']['config']);config['catalogs']+=deepcopy(bundle['records'])
    change(config['catalogs'])
    labels=deepcopy(bundle['labels'])
    for _ in range(30):
        mapping={row['record_id']:record_identity(row) for row in config['catalogs'] if row['record_id']!=record_identity(row)}
        if not mapping:break
        config=replace(config,mapping);labels=replace(labels,mapping)
        for row in config['catalogs']:
            if row['record_type']=='comparison_design':row['definition_hash']=row['record_id'].rsplit(':',1)[1]
    else:raise AssertionError('Identity graph did not converge')
    path=bundle['study']['config_path'];path.write_text(json.dumps(config))
    study=resolve_study(path)
    return prepare_question_set(study,[]),labels


def policy_changed(bundle, change):
    path=bundle['study']['config_path'];config=deepcopy(bundle['study']['config'])
    policy_path=path.parent/'experiment-policy.json';policy=json.loads(policy_path.read_text());change(policy)
    policy_path.write_text(json.dumps(policy))
    next(row for row in config['input_refs'] if row['role']=='experiment_policy')['sha256']=hashlib.sha256(policy_path.read_bytes()).hexdigest()
    path.write_text(json.dumps(config));study=resolve_study(path)
    output=generate_questions(study,[])
    return prepare_question_set(study,output['question-records.jsonl']),output


@pytest.fixture
def fixture(tmp_path):return build_fixture(tmp_path/'fixture')


def state(qset, selected=None):
    ids=[identifier for identifier in qset.actions if qset.eligibility[identifier]['eligible']] if selected is None else selected
    return evaluate_requirements(qset,ids)


def test_full_fixture_has_separate_endpoints_actual_noncarrier_and_no_empirical_result(fixture):
    qset=prepare_question_set(fixture['study'],fixture['records']);result=state(qset)
    assert len(qset.actions)==10 and len(qset.slots)==20
    assert sum(row['cost'] for row in qset.resources.values())==18
    assert result['unweighted_completion']==1
    assert next(iter(result['questions'].values()))['planning_status']=='fully_planned'
    assert next(iter(result['questions'].values()))['evidence_status']=='no_measurements'
    labels=fixture['labels'];noncarrier=labels['donors']['N']
    for row in qset.actions.values():
        if row['donor_id']==noncarrier:
            assert next(t for t in row['target_sequences'] if t['locus_id']==labels['loci']['L'])['sequence_configuration_ids']==[labels['sequences']['L_ref']]
    assert fixture['study']['config']['coverage_basis']=='synthetic_truth'


@pytest.mark.parametrize('template',['candidate_validation','paired_design_comparison','variant_association','design_by_variant'])
def test_each_nonreplication_template_builds_complete_requirements(fixture,template):
    qset,output=policy_changed(fixture,lambda p:p['rules'][0].update(templates=[template]))
    assert {row['template'] for row in qset.questions.values()}=={template}
    assert state(qset)['unweighted_completion']>=1
    for row in qset.questions.values():
        assert all(arm['required_independent_donors']==1 for arm in row['arms'])
        assert row['technical_replicates']==1
        assert row['pairing']==('within_donor' if template in {'paired_design_comparison','design_by_variant'} else 'none')


@pytest.mark.parametrize('field',['required_donors_per_arm','technical_replicates'])
def test_missing_scientific_count_is_a_field_error(fixture,field):
    with pytest.raises(PopulationError,match=field):policy_changed(fixture,lambda p:p['rules'][0].pop(field))


def test_dropping_b_side_slots_cannot_create_a_complete_question(fixture):
    def change(rows):
        design=fixture['labels']['designs']['B'];removed={r['record_id'] for r in rows if r['record_type']=='comparison_requirement_slot' and r['design_id']['value']==design}
        rows[:]=[r for r in rows if r['record_id'] not in removed]
        for row in rows:
            if row['record_type']=='comparison_question':row['requirement_slot_ids']=[v for v in row['requirement_slot_ids'] if v not in removed]
            if row['record_type']=='comparison_action':row['memberships']=[v for v in row['memberships'] if v['slot_id'] not in removed]
    qset,_=mutated(fixture,change)
    assert state(qset)['unweighted_completion']==0
    assert all('required_paired_or_control_slot_omitted:intended_activity' in b['reason_codes'] for b in eligible_donor_bundles(qset,next(iter(qset.questions))))


def test_missing_required_intended_assessment_blocks_original_carrier(fixture):
    donor=fixture['labels']['donors']['C']
    qset,_=mutated(fixture,lambda rows:rows.__setitem__(slice(None),[r for r in rows if not(r['record_type']=='comparison_intended_assessment' and r['donor_id']==donor)]))
    assert state(qset)['unweighted_completion']==0
    assert any('required_intended_locus_unresolved' in r for b in eligible_donor_bundles(qset,next(iter(qset.questions))) for r in b['reason_codes'])


def test_alternate_carrier_restores_question_without_changing_historical_input(tmp_path):
    bundle=build_fixture(tmp_path/'alternate',alternate_carrier=True);before=deepcopy(bundle['records']);donor=bundle['labels']['donors']['C']
    qset,labels=mutated(bundle,lambda rows:rows.__setitem__(slice(None),[r for r in rows if not(r['record_type']=='comparison_intended_assessment' and r['donor_id']==donor)]))
    result=state(qset);assert result['unweighted_completion']==1
    assert labels['donors']['C2'] in next(iter(result['questions'].values()))['selected_donor_ids']
    assert bundle['records']==before


def test_unknown_noncarrier_copy_does_not_become_reference(fixture):
    donor=fixture['labels']['donors']['N'];sequence=fixture['labels']['sequences']['L_ref']
    def change(rows):
        removed={r['record_id'] for r in rows if r['record_type']=='comparison_sequence_occurrence' and r['donor_id']==donor and r['sequence_configuration_id']==sequence and r['copy_index']==1}
        rows[:]=[r for r in rows if r['record_id'] not in removed and not(r['record_type']=='comparison_donor_configuration' and r['sequence_occurrence_id'] in removed)]
    qset,_=mutated(fixture,change);assert state(qset)['unweighted_completion']==0
    assert any('both_actual_copies_are_not_resolved' in v['reason_codes'] for v in qset.donor_eligibility.values())


def test_related_carrier_and_noncarrier_are_not_independent(fixture):
    qset,_=mutated(fixture,lambda rows:[r.update(independence_group_id=known('same-family')) for r in rows if r['record_type']=='comparison_donor'])
    assert state(qset)['unweighted_completion']==0


def test_unknown_relationship_strict_blocks_and_explicit_simulation_labels_assumption(fixture):
    def change(rows):
        for row in rows:
            if row['record_type']=='comparison_donor':row.update(independence_group_id=unknown('not measured'),relatedness_policy='strict')
    qset,_=mutated(fixture,change);assert state(qset)['unweighted_completion']==0
    bundle={'study':qset.study,'records':[],'labels':fixture['labels']}
    qset,_=mutated(bundle,lambda rows:[r.update(relatedness_policy='declared_simulation_assumption') for r in rows if r['record_type']=='comparison_question'])
    assert state(qset)['unweighted_completion']==1
    assert all(v['assumptions'] for v in qset.donor_eligibility.values())


def test_unmatched_mock_or_qc_operations_cannot_fill_prerequisites(fixture):
    def change(rows):
        for row in rows:
            if row['record_type']=='comparison_action' and row['operation_type']=='measurement_panel':
                row['dependency_action_ids']=[i for i in row['dependency_action_ids'] if next(r for r in rows if r['record_id']==i)['operation_type']!='mock_control']
    qset,_=mutated(fixture,change);assert state(qset)['unweighted_completion']==0
    assert any('compatible_control_prerequisite_missing:mock' in v['reason_codes'] for v in qset.eligibility.values())


def test_cycle_rejected_before_optimization(fixture):
    def change(rows):
        action=next(r for r in rows if r['record_type']=='comparison_action');action['dependency_action_ids'].append(action['record_id'])
    with pytest.raises(PopulationError,match='[Cc]ycl'):mutated(fixture,change)


def test_actual_target_cannot_be_copied_from_another_donor(fixture):
    def change(rows):
        action=next(r for r in rows if r['record_type']=='comparison_action' and r['donor_id']==fixture['labels']['donors']['N'])
        next(t for t in action['target_sequences'] if t['locus_id']==fixture['labels']['loci']['L'])['sequence_configuration_ids']=[fixture['labels']['sequences']['L_alt']]
    with pytest.raises(PopulationError,match='actual donor'):mutated(fixture,change)


def test_shortlist_does_not_remove_scientific_questions(fixture):
    qset,output=policy_changed(fixture,lambda p:p.update(max_action_questions=0))
    assert len(qset.questions)==1 and not qset.actions
    assert output['question-audit.json']['question_exclusions']
    assert state(qset)['unweighted_completion']==0


def test_specimen_use_must_be_explicit(fixture):
    qset,_=mutated(fixture,lambda rows:[r.update(use_authorized=False) for r in rows if r['record_type']=='comparison_specimen'])
    assert not any(v['eligible'] for v in qset.eligibility.values())
    assert state(qset)['unweighted_completion']==0


def test_research_mode_refuses_simulation_specimen_creation(fixture):
    config=deepcopy(fixture['study']['config']);config.update(mode='research_sample_planning',coverage_basis='callable_supported');path=fixture['study']['config_path'];path.write_text(json.dumps(config))
    with pytest.raises(PopulationError,match='Research specimens'):generate_questions(resolve_study(path),[])


def test_alternative_resources_cannot_combine_distinct_package_types(fixture):
    def change(rows):
        panel=next(r for r in rows if r['record_type']=='comparison_action' and r['operation_type']=='measurement_panel')
        resources=[r for r in rows if r['record_type']=='comparison_resource']
        seq=next(r for r in resources if r['resource_type']=='sequencing_batch');mock=next(r for r in resources if r['resource_type']=='mock_control')
        panel['resource_demands'].append({'eligible_resource_ids':[seq['record_id'],mock['record_id']],'assignment_key':'mixed','dimension':'uses','quantity':1})
    with pytest.raises(PopulationError,match='same resource type'):mutated(fixture,change)


def test_dependency_selection_must_be_closed(fixture):
    qset=prepare_question_set(fixture['study'],fixture['records']);panel=next(row for row in qset.actions.values() if row['operation_type']=='measurement_panel')
    with pytest.raises(PopulationError,match='prerequisite'):evaluate_requirements(qset,[panel['record_id']])
    assert len(dependency_closure(qset,[panel['record_id']]))==3


def test_policy_bytes_and_dependency_projection_are_bound(fixture):
    before=effective_inputs(fixture['study'],'questions');assert any(row['role']=='experiment_policy' for row in before['sources'].values())
    path=fixture['study']['config_path'].parent/'experiment-policy.json';path.write_text(path.read_text()+' ')
    with pytest.raises(PopulationError,match='source changed'):prepare_question_set(fixture['study'],fixture['records'])


def observation_records(bundle, *, donor='C', design='A', locus='L', experiment='completed-experiment', replicate='1', export='one', credit=True, executed=None):
    templates=json.loads((ROOT/'tests/fixtures/comparison/canonical-records.json').read_text())
    labels=bundle['labels'];source_template=deepcopy(next(r for r in templates if r['record_type']=='comparison_source_observation'))
    def fields(row):return {k:v for k,v in row.items() if k not in {'record_type','schema_version','record_id',record_schema(row['record_type'])['x-id-field']}}
    artifact=next(r for r in bundle['study']['config']['catalogs'] if r['record_type']=='comparison_source_artifact')
    source_template.update(source_artifact_id=artifact['record_id'],measurement_key=donor+':'+str(design)+':'+locus+':'+export,source_experiment=experiment,source_replicate=replicate,locator={'kind':'file','member':'artificial-observation-'+export},raw_row_text='Artificial measurement: 0.2 fraction')
    source=make_record('source_observation',**fields(source_template))
    row=deepcopy(next(r for r in templates if r['record_type']=='comparison_observation_revision'))
    specimen=next((r for r in bundle['records'] if r['record_type']=='comparison_specimen' and r['donor_id']==labels['donors'][donor]),None) if executed else None
    row.update(source_record_id=source['record_id'],design_id=known(labels['designs'][design]) if design else unknown('Design-neutral completed mock'),unresolved_source_design=unknown('Resolved design') if design else known('mock'),context_id=known(labels['context_id']),sequence_configuration_id=known(labels['sequences']['T' if locus=='T' else 'L_ref' if donor=='N' else 'L_alt']),locus_ids=[labels['loci'][locus]],endpoint_id=labels['endpoints'][locus],executed_action_id=known(executed) if executed else unknown('Historical measurement without operation record'))
    row['subject'].update(donor_id=known(labels['donors'][donor]),source_subject_id=known('synthetic-'+donor),specimen_id=known(specimen['record_id']) if specimen else unknown('Historical specimen not supplied'))
    row['applicability']['empirical_donor_credit']=credit;row['applicability']['reason_codes']=[] if credit else ['no_empirical_credit']
    observation=make_record('observation_revision',**fields(row))
    state_record=make_record('observation_state',source_record_id=source['record_id'],active_revision_id=observation['record_id'],revision_ids=[observation['record_id']],supersession_reason='Initial artificial test normalization')
    return [source,observation,state_record]


def test_evidence_replication_is_generated_and_excludes_its_source(fixture):
    observations=observation_records(fixture,credit=False)
    fixture['study']['config']['catalogs'].extend(observations)
    fixture['study']['config_path'].write_text(json.dumps(fixture['study']['config']))
    fixture['study']=resolve_study(fixture['study']['config_path'])
    qset,_=policy_changed(fixture,lambda p:p['rules'][0].update(templates=['evidence_replication'],trigger='existing_observation',replication_dimensions=['experiment']))
    assert qset.questions and all(q['motivating_observation_ids'] for q in qset.questions.values())
    assert state(qset)['unweighted_completion']>=1
    for question in qset.questions.values():assert set(question['motivating_observation_ids']) <= set(question['evidence_reuse']['excluded_observation_ids'])


def test_one_actual_measurement_and_its_duplicate_export_cannot_fill_two_repeats(fixture):
    qset,output=policy_changed(fixture,lambda p:p['rules'][0].update(technical_replicates=2,allow_existing=True))
    bundle={'study':qset.study,'records':output['question-records.jsonl'],'labels':fixture['labels']}
    observations=observation_records(bundle,export='one')+observation_records(bundle,export='duplicate')
    observation_ids=[r['record_id'] for r in observations if r['record_type']=='comparison_observation_revision']
    def change(rows):
        rows.extend(observations)
        for row in rows:
            if row['record_type']=='comparison_requirement_slot' and row['donor_id']==fixture['labels']['donors']['C'] and row['design_id']['value']==fixture['labels']['designs']['A'] and row['control_role']=='measurement':row['allowed_observation_ids']=observation_ids
    qset,_=mutated(bundle,change)
    assignments=evaluate_requirements(qset)['slot_assignments']
    assert sum('observation_id' in value for value in assignments.values())==1
    assert len({qset.observation_tokens[identifier] for identifier in observation_ids})==1
    explicit={identifier:observation_ids[0] for identifier,slot in qset.slots.items() if slot['allowed_observation_ids']}
    with pytest.raises(PopulationError,match='reuse'):evaluate_requirements(qset,observation_assignments=explicit)


def test_distinct_completed_measurements_can_fill_distinct_repeat_slots(fixture):
    qset,output=policy_changed(fixture,lambda p:p['rules'][0].update(technical_replicates=2,allow_existing=True))
    bundle={'study':qset.study,'records':output['question-records.jsonl'],'labels':fixture['labels']}
    observations=observation_records(bundle,export='one',replicate='1')+observation_records(bundle,export='two',replicate='2')
    ids=[r['record_id'] for r in observations if r['record_type']=='comparison_observation_revision']
    def change(rows):
        rows.extend(observations)
        for row in rows:
            if row['record_type']=='comparison_requirement_slot' and row['donor_id']==fixture['labels']['donors']['C'] and row['design_id']['value']==fixture['labels']['designs']['A'] and row['control_role']=='measurement':row['allowed_observation_ids']=ids
    qset,_=mutated(bundle,change)
    assert sum('observation_id' in value for value in evaluate_requirements(qset)['slot_assignments'].values())==2


def test_planned_controls_never_supply_passed_observation_qc(fixture):
    observations=observation_records(fixture,credit=False);observation=next(r for r in observations if r['record_type']=='comparison_observation_revision');observation['qc_status']='failed'
    # A failed empirical record remains rejected even with a fully paid future plan.
    observation['record_id']=record_identity(observation);observation['observation_id']=observation['record_id']
    old=next(r for r in observations if r['record_type']=='comparison_observation_state');old['active_revision_id']=observation['record_id'];old['revision_ids']=[observation['record_id']];old['record_id']=record_identity(old);old['observation_state_id']=old['record_id']
    qset=prepare_question_set(fixture['study'],fixture['records']+observations)
    assert next(iter(state(qset)['questions'].values()))['evidence_status']=='no_measurements'
    assert next(iter(state(qset)['questions'].values()))['planning_status']=='fully_planned'


@pytest.mark.parametrize('dimension',['experiment','batch','specimen','donor'])
def test_replication_rejects_original_dimension(fixture,dimension):
    # Bind the motivating source to a known prior operation, then request a repeat.
    operation=next(r for r in fixture['records'] if r['record_type']=='comparison_action' and r['operation_type']=='measurement_panel' and r['donor_id']==fixture['labels']['donors']['C'] and r['design_id']['value']==fixture['labels']['designs']['A'])
    observations=observation_records(fixture,experiment=operation['experiment_id'],executed=operation['record_id'],credit=False)
    def change(rows):
        rows.extend(observations)
        for question in [r for r in rows if r['record_type']=='comparison_question']:
            # Keep only A and one all-resolved arm, preserving controls and donor alternatives.
            question.update(template='evidence_replication',design_ids=[fixture['labels']['designs']['A']],pairing='none',arms=[{'arm_id':'all_resolved','membership':'all_resolved','required_independent_donors':1}])
            question['genotype_predicate'].update(kind='any_resolved',variants=[],sequence_configuration_id=unknown('any resolved donor'))
            source=observations[0]['record_id'];question['motivating_observation_ids']=[source];question['evidence_reuse'].update(excluded_observation_ids=[source],excluded_experiment_ids=[operation['experiment_id']],require_new=[dimension])
            # A new experiment is required in all replication templates; hold only
            # the separately requested dimension equal to its motivating source.
            for action in [r for r in rows if r['record_type']=='comparison_action']:
                if dimension!='experiment':action['experiment_id']='new-replication-experiment'
            keep={r['record_id'] for r in rows if r['record_type']=='comparison_requirement_slot' and r['design_id']['value'] in {None,fixture['labels']['designs']['A']}}
            rows[:]=[r for r in rows if r['record_type']!='comparison_requirement_slot' or r['record_id'] in keep]
            question['requirement_slot_ids']=sorted(keep)
            for slot in [r for r in rows if r['record_type']=='comparison_requirement_slot']:slot['arm_id']='all_resolved'
            for action in [r for r in rows if r['record_type']=='comparison_action']:action['memberships']=[m for m in action['memberships'] if m['slot_id'] in keep]
    qset,_=mutated(fixture,change)
    carrier_slots=[s for s in qset.slots.values() if s['donor_id']==fixture['labels']['donors']['C']]
    assert all(not qset.options[slot['record_id']]['action_ids'] for slot in carrier_slots)


def test_joint_observation_assignment_prioritizes_hard_question_over_optional_weight(fixture):
    observations=observation_records(fixture)
    observation_id=next(r['record_id'] for r in observations if r['record_type']=='comparison_observation_revision')
    def change(rows):
        rows.extend(observations)
        hard=next(r for r in rows if r['record_type']=='comparison_question');hard['evidence_reuse'].update(allow_existing=True,allow_across_questions=False)
        optional=deepcopy(hard);optional.update(hypothesis='Optional competing measurement requirement',priority_weight=1000,hard_required=False)
        optional['record_id']=record_identity(optional);optional['question_id']=optional['record_id']
        original_actions=[r for r in rows if r['record_type']=='comparison_action'];new_actions=deepcopy(original_actions)
        mapping={}
        for row in new_actions:
            old=row['record_id'];row['experiment_id']='optional-prospective-experiment';row['record_id']=record_identity(row);row['action_id']=row['record_id'];row['memberships']=[];mapping[old]=row['record_id']
        for row in new_actions:row['dependency_action_ids']=[mapping[v] for v in row['dependency_action_ids']]
        original_slots=[r for r in rows if r['record_type']=='comparison_requirement_slot']
        for row in original_slots:
            if row['donor_id']==fixture['labels']['donors']['C'] and row['design_id']['value']==fixture['labels']['designs']['A'] and row['control_role']=='measurement':row.update(allowed_action_ids=[],allowed_observation_ids=[observation_id])
        new_slots=deepcopy(original_slots)
        for row in new_slots:
            row['question_id']=optional['record_id'];row['record_id']=record_identity(row);row['slot_id']=row['record_id'];row['allowed_action_ids']=[mapping[v] for v in row['allowed_action_ids']]
            for action in new_actions:
                if action['record_id'] in row['allowed_action_ids']:action['memberships'].append({'question_id':optional['record_id'],'slot_id':row['record_id']})
        optional['requirement_slot_ids']=[row['record_id'] for row in new_slots];rows.extend([optional,*new_slots,*new_actions])
    qset,_=mutated(fixture,change);result=state(qset)
    assert result['unweighted_completion']==1
    assert all(qset.questions[identifier]['hard_required'] for identifier in result['completed_question_ids'])
    assert result['weighted_completion']==1


def test_completed_panels_and_matched_completed_mocks_are_already_satisfied(fixture):
    observations=[]
    for action in fixture['records']:
        if action['record_type']!='comparison_action' or action['operation_type'] not in {'measurement_panel','mock_control'}:continue
        donor=next(name for name,identifier in fixture['labels']['donors'].items() if identifier==action['donor_id'])
        design=next((name for name,identifier in fixture['labels']['designs'].items() if identifier==action['design_id']['value']),None)
        for locus in ('T','L'):observations.extend(observation_records(fixture,donor=donor,design=design,locus=locus,experiment=action['experiment_id'],export=donor+str(design)+locus,executed=action['record_id']))
    def change(rows):
        rows.extend(observations)
        for row in rows:
            if row['record_type']=='comparison_question':row['evidence_reuse']['allow_existing']=True
            if row['record_type']=='comparison_requirement_slot':
                row['allowed_observation_ids']=[obs['record_id'] for obs in observations if obs['record_type']=='comparison_observation_revision' and obs['subject']['donor_id']['value']==row['donor_id'] and obs['design_id']['value']==row['design_id']['value'] and obs['endpoint_id']==row['endpoint_id']]
    qset,_=mutated(fixture,change);result=evaluate_requirements(qset)
    assert result['unweighted_completion']==1
    question=next(iter(result['questions'].values()))
    assert question['evidence_status']=='evaluable' and question['planning_status']=='already_satisfied'
    assert all('observation_id' in value for value in result['slot_assignments'].values())


def test_excluded_questions_have_no_empty_cartesian_slot_catalog(fixture):
    qset,output=policy_changed(fixture,lambda p:p.update(max_action_questions=0))
    question=next(iter(qset.questions.values()))
    assert question['construction_status']=='excluded_by_declared_limit'
    assert question['requirement_slot_ids']==[] and not qset.slots
    assert not qset.donor_eligibility
    assert next(iter(evaluate_requirements(qset)['questions'].values()))['remaining_gaps']==['requirement_construction_shortlist_exclusion']


def test_supported_lost_pam_intended_dna_remains_an_eligible_assay(fixture):
    qset,_=mutated(fixture,lambda rows:[r.update(pam_compatible=known(False),reason_codes=['supported_altered_intended_PAM']) for r in rows if r['record_type']=='comparison_intended_assessment'])
    assert state(qset)['unweighted_completion']==1


def test_declared_covariate_matching_is_checked_across_arms(fixture):
    def change(rows):
        for row in rows:
            if row['record_type']=='comparison_donor':row['covariates']=[{'name':'declared_block','value':known(row['canonical_label'])}]
            if row['record_type']=='comparison_question':row['matching_covariates']=['declared_block']
    qset,_=mutated(fixture,change)
    assert state(qset)['unweighted_completion']==0


def test_mock_sharing_policy_adds_real_control_packages(fixture):
    qset,_=policy_changed(fixture,lambda p:p['assay_policy'].update(mock_shared_across_designs=False))
    assert len([r for r in qset.actions.values() if r['operation_type']=='mock_control'])==4
    assert sum(r['cost'] for r in qset.resources.values())==20
    assert state(qset)['unweighted_completion']==1


def test_a_treatment_label_cannot_replace_priced_assay_qc(fixture):
    def change(rows):
        for row in rows:
            if row['record_type']=='comparison_action' and row['operation_type']=='measurement_panel':row['control_requirements']=[]
    qset,_=mutated(fixture,change)
    assert all(not value['eligible'] for identifier,value in qset.eligibility.items() if qset.actions[identifier]['operation_type']=='measurement_panel')
    assert state(qset)['unweighted_completion']==0
