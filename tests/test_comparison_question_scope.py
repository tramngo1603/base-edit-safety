"""A common readout remains common when designs have different intended loci."""
from copy import deepcopy
import hashlib
import json

from base_edit_safety.population.comparison_contracts import make_record,record_identity,record_registry
from base_edit_safety.population.comparison_provenance import resolve_study
from base_edit_safety.population.comparison_questions import generate_questions,prepare_question_set,evaluate_requirements
from base_edit_safety.population.experiment_planning import plan_experiments
from test_comparison_questions import build_fixture,replace


def distinct_intended_fixture(directory):
    fixture=build_fixture(directory)
    config=deepcopy(fixture['study']['config']);labels=deepcopy(fixture['labels'])
    policy_path=directory/'experiment-policy.json';policy=json.loads(policy_path.read_text())
    rows=config['catalogs'];index={r['record_id']:r for r in rows}
    intended=deepcopy(index[labels['loci']['T']]);old_candidate=labels['loci']['L']
    for key in ('record_type','schema_version','record_id',record_registry()[intended['record_type']]):intended.pop(key)
    intended.update(start=200,end=223,locus_name='Distinct artificial intended locus L',path=[{**p,'pos0':p['pos0']+100} for p in intended['path']])
    new_locus=make_record('intended_locus',**intended)
    rows=replace([r for r in rows if r['record_id']!=old_candidate]+[new_locus],{old_candidate:new_locus['record_id']});labels=replace(labels,{old_candidate:new_locus['record_id']})
    for row in rows:
        if row['record_id']==labels['designs']['A']:
            row['intended_locus_ids']=[new_locus['record_id']]
            for edit in row['desired_edit_definitions']:edit.update(locus_id=new_locus['record_id'],pos0=200)
        if row['record_type']=='comparison_design_configuration' and row['design_id']==labels['designs']['A']:
            row['role']='intended' if row['sequence_configuration_id'] in (labels['sequences']['L_ref'],labels['sequences']['L_alt']) else 'off_target'
        if row['record_type']=='comparison_intended_assessment' and row['design_id']==labels['designs']['A']:
            row['locus_id']=new_locus['record_id'];row['sequence_configuration_id']['value']=labels['sequences']['L_ref' if row['donor_id']==labels['donors']['N'] else 'L_alt']
    for _ in range(30):
        mapping={r['record_id']:record_identity(r) for r in rows if r['record_id']!=record_identity(r)}
        if not mapping:break
        rows,labels,policy=replace(rows,mapping),replace(labels,mapping),replace(policy,mapping)
        config['policies']=replace(config['policies'],mapping)
        for row in rows:
            if row['record_type']=='comparison_design':row['definition_hash']=row['record_id'].rsplit(':',1)[1]
    else:raise AssertionError('Source identity graph did not converge')
    policy['rules'][0]['trigger']='variant_intended';policy['max_action_questions']=10
    policy_path.write_text(json.dumps(policy,indent=2)+'\n')
    next(r for r in config['input_refs'] if r['role']=='experiment_policy')['sha256']=hashlib.sha256(policy_path.read_bytes()).hexdigest()
    config['catalogs']=rows;(directory/'study.json').write_text(json.dumps(config,indent=2)+'\n')
    study=resolve_study(directory/'study.json');output=generate_questions(study,[])
    return prepare_question_set(study,output['question-records.jsonl']),labels,output


def test_every_design_measures_the_declared_common_endpoint_at_other_designs_intended_locus(tmp_path):
    qset,labels,_=distinct_intended_fixture(tmp_path/'source-case')
    expected={(labels['endpoints']['L'],labels['loci']['L']),
              (labels['endpoints']['T'],labels['loci']['L']),
              (labels['endpoints']['T'],labels['loci']['T'])}
    for donor in labels['donors'].values():
        for design in labels['designs'].values():
            panel=next(r for r in qset.actions.values() if r['operation_type']=='measurement_panel' and r['donor_id']==donor and r['design_id']['value']==design)
            assert {(p['endpoint_id'],p['locus_id']) for p in panel['endpoint_locus_pairs']}==expected
            assert any(s['donor_id']==donor and s['design_id']['value']==design and s['endpoint_id']==labels['endpoints']['L'] and s['locus_id']['value']==labels['loci']['L'] and s['control_role']=='measurement' for s in qset.slots.values())
    for name,design in labels['designs'].items():
        intended_locus=labels['loci']['L' if name=='A' else 'T']
        intended_slots=[s for s in qset.slots.values() if s['design_id']['value']==design and s['control_role']=='intended_activity']
        assert intended_slots and {(s['endpoint_id'],s['locus_id']['value']) for s in intended_slots}=={(labels['endpoints']['T'],intended_locus)}
    result=evaluate_requirements(qset,qset.actions)
    assert result['unweighted_completion']==1
    assert next(iter(result['questions'].values()))['evidence_status']=='no_measurements'


def test_removing_common_readout_slots_cannot_redefine_the_scientific_contrast(tmp_path):
    qset,labels,_=distinct_intended_fixture(tmp_path/'missing-readout')
    removed={s['record_id'] for s in qset.slots.values() if s['design_id']['value']==labels['designs']['B'] and s['endpoint_id']==labels['endpoints']['L'] and s['locus_id']['value']==labels['loci']['L']}
    rows=[deepcopy(r) for r in qset.records if r['record_id'] not in removed]
    for row in rows:
        if row['record_type']=='comparison_question':row['requirement_slot_ids']=[v for v in row['requirement_slot_ids'] if v not in removed]
        if row['record_type']=='comparison_action':row['memberships']=[v for v in row['memberships'] if v['slot_id'] not in removed]
    checked=prepare_question_set(qset.study,rows)
    result=evaluate_requirements(checked,checked.actions)
    assert result['unweighted_completion']==0
    assert any('required_paired_or_control_slot_omitted:measurement' in donor['reason_codes'] for q in result['questions'].values() for donor in q['donors'])


def test_original_section9_complete_cost_and_budget_boundary_remain_unchanged(tmp_path):
    fixture=build_fixture(tmp_path/'section9');qset=prepare_question_set(fixture['study'],fixture['records'])
    assert len(qset.actions)==10 and len(qset.slots)==20
    assert sum(r['cost'] for r in qset.resources.values())==18
    full=plan_experiments(qset,budget=18,selector='exact_tiny')
    short=plan_experiments(qset,budget=17,selector='exact_tiny')
    assert full['meets_request'] and full['weighted_completion']==1 and full['accounting']['total_cost']==18
    assert not short['meets_request'] and short['weighted_completion']==0
