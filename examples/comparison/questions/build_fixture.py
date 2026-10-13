"""Export a mock carrier/noncarrier question with an 18-unit full panel."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from base_edit_safety.population.comparison_contracts import make_record, validate_study, validate_catalog
from base_edit_safety.population.common import known, unknown
from base_edit_safety.population.comparison_provenance import resolve_study
from base_edit_safety.population.comparison_questions import generate_questions

ROOT=Path(__file__).resolve().parents[3]


def build_fixture(directory, *, alternate_carrier=False):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    base=json.loads((ROOT/'examples/comparison/contracts/study.json').read_text())
    catalog=base['catalogs']; bykind={row['record_type'].removeprefix('comparison_'):row for row in catalog}
    def add(kind, **fields):
        row=make_record(kind,**fields);catalog.append(row);return row
    def clone(row, **changes):
        fields={key:value for key,value in row.items() if key not in {'record_type','schema_version','record_id',row['record_type'].removeprefix('comparison_')+'_id','definition_hash'}}
        # Some named primary IDs differ from the record kind.
        for key in ('context_id','endpoint_id','intended_locus_id'):fields.pop(key,None)
        fields.update(changes);return add(row['record_type'].removeprefix('comparison_'),**fields)
    reference=(ROOT/'examples/comparison/contracts/reference.fa').read_bytes();(directory/'reference.fa').write_bytes(reference)
    designs={'A':bykind['design']}
    editor=clone(bykind['editor_definition'],display_name='Synthetic second editor',variant_name='synthetic-editor-B',editing_window={'value':None,'reason':'Unreported biological window; search and assay remain evaluable'})
    designs['B']=clone(bykind['design'],display_name='Synthetic design B',editor_definition_id=editor['record_id'])
    profile=clone(bykind['search_profile'],design_id=designs['B']['record_id'])
    base['policies']['search']['design_ids'].append(designs['B']['record_id']);base['policies']['search']['profile_ids'].append(profile['record_id'])
    context=bykind['experimental_context']; intended=bykind['intended_locus']; intended_endpoint=bykind['endpoint_definition']
    candidate=add('locus',reference_bundle_id=bykind['reference_bundle']['record_id'],contig='synthetic_chr11',start=200,end=223,strand='+',locus_name='L',gene=unknown('Artificial candidate locus'),path=[])
    endpoint=clone(intended_endpoint,name='unintended_edit_fraction',biological_meaning='Measured editing at the separate candidate locus L',favorable_direction='lower')
    donors={'C':bykind['donor']}
    for name in ['N']+(['C2'] if alternate_carrier else []):
        donors[name]=clone(bykind['donor'],canonical_label='synthetic-donor-'+name,sample_ids=['synthetic-sample-'+name],independence_group_id=known('synthetic-family-'+name))
    raw={'kind':'artificial_diploid_truth','intended_sequence':'ACGTACGTACGTACGTACGTAGG','candidate_reference':'ACGTACGTACGTACGTACGTAGG','candidate_alternate':'ACGTGCGTACGTACGTACGTAGG','donors':{name:{'intended_copies':['reference','reference'],'candidate_copies':['alternate','reference'] if name!='N' else ['reference','reference']} for name in donors}}
    (directory/'synthetic-genomes.json').write_text(json.dumps(raw,indent=2)+'\n')
    sha=hashlib.sha256((directory/'synthetic-genomes.json').read_bytes()).hexdigest()
    artifact=add('source_artifact',sha256=sha,source_id='artificial-diploid-truth',release='section9-v1',uri='urn:bes:synthetic:section9',media_type='application/json',scope='complete_file',parent_artifact_id=unknown('Original artificial inputs'),locator={'kind':'file','member':'synthetic-genomes.json'},access='synthetic',use_terms='Artificial software fixture')
    refs=[{'source_artifact_id':artifact['record_id'],'locator':{'kind':'file','member':'synthetic-genomes.json'},'note':'Artificial copy sequences are literal source truth; no human specimens or observations'}]
    variant={'contig':'synthetic_chr11','pos0':204,'reference':'A','alternates':['G'],'selected_alt_index':1,'source_artifact_id':artifact['record_id'],'source_row':1}
    sequences={}
    for key,locus,seq,variants in [('T',intended,raw['intended_sequence'],[]),('L_ref',candidate,raw['candidate_reference'],[]),('L_alt',candidate,raw['candidate_alternate'],[variant])]:
        sequences[key]=add('sequence_configuration',reference_bundle_id=bykind['reference_bundle']['record_id'],locus_id=locus['record_id'],contig='synthetic_chr11',strand='+',sequence=seq,path=[{'kind':'reference','pos0':locus['start']+i,'base':base} for i,base in enumerate(seq)],variants=variants,support_status='supported',support_reason=None)
    configs={}
    for label,design in designs.items():
        for key,sequence in sequences.items():configs[(label,key)]=add('design_configuration',design_id=design['record_id'],sequence_configuration_id=sequence['record_id'],role='intended' if key=='T' else 'off_target')
    for name,donor in donors.items():
        for copy in (0,1):
            for locus_key in ('T','L_alt' if copy==0 and name!='N' else 'L_ref'):
                sequence=sequences[locus_key]
                occurrence=add('sequence_occurrence',sequence_configuration_id=sequence['record_id'],donor_id=donor['record_id'],sample_id=donor['sample_ids'][0],copy_index=copy,phase_block=known('artificial-phase'),phase_status='confirmed',source_version='section9-v1',source_refs=refs)
                for design in designs:
                    add('donor_configuration',donor_id=donor['record_id'],sample_id=donor['sample_ids'][0],sequence_occurrence_id=occurrence['record_id'],design_configuration_id=configs[(design,locus_key)]['record_id'],copy_index=copy,carrier_status='supported_carrier',carrier_of_sequence_id=sequence['record_id'],genotype_support='called',phase_support='confirmed',reason_codes=[],source_refs=refs)
        for design in designs.values():add('intended_assessment',donor_id=donor['record_id'],design_id=design['record_id'],locus_id=intended['record_id'],sequence_configuration_id=known(sequences['T']['record_id']),support_status='supported',pam_compatible=known(True),alignment_ids=[],reference_differences=[],reason_codes=[])
    policy={'record_type':'comparison_experiment_policy','schema_version':'1','policy_version':'section9-declared-v1','context_id':context['record_id'],'experiment_id':'section9-new-experiment','batch_id':'section9-new-batch','donor_ids':[row['record_id'] for row in donors.values()],'rules':[{'rule_id':'complete-design-by-variant','trigger':'variant_candidate','templates':['design_by_variant'],'design_ids':[row['record_id'] for row in designs.values()],'candidate_endpoint_id':endpoint['record_id'],'intended_endpoint_id':intended_endpoint['record_id'],'required_donors_per_arm':1,'technical_replicates':1,'priority_weight':1,'hard_required':True,'allow_existing':False,'allow_across_questions':True,'replication_dimensions':[],'matching_covariates':[]}],'max_action_questions':1,'max_donor_alternatives_per_arm':2,'relatedness_policy':'strict','specimen_policy':{'create_simulation_specimens':True,'permitted_use':'Artificial software fixture only','use_authorized':True,'capacity':100,'assumption':'Artificial specimens exist by fixture construction'},'assay_policy':{'feasibility':'simulation_assumption','resolution':'copy_specific','allele_resolved':True,'qc_requirements':['assay_qc'],'mock_shared_across_designs':True,'package_basis':'Each 2-unit treated panel measures T and L with its assay controls/QC; each 1-unit donor mock measures T and L with its assay controls/QC'},'costs':{'cost_unit':'simulation_units','cost_scale':1,'study_setup':2,'donor_activation':1,'specimen_preparation':0,'treatment':1,'treated_panel':2,'mock_panel':1,'sequencing_batch':0,'additional_replicate':0},'resources':{'study_setup_activation':'nonempty_plan','specimen_capacity':100,'panel_capacity':1,'mock_capacity':1,'sequencing_capacity':6 if not alternate_carrier else 9,'max_sequencing_instances':1}}
    (directory/'experiment-policy.json').write_text(json.dumps(policy,indent=2)+'\n')
    base['study_id']='artificial-section9-complete-panel';base['input_refs'] += [{'input_id':'artificial-genomes','role':'variants','path':'synthetic-genomes.json','sha256':sha,'format':'json','stages':['search']},{'input_id':'prospective-policy','role':'experiment_policy','path':'experiment-policy.json','sha256':hashlib.sha256((directory/'experiment-policy.json').read_bytes()).hexdigest(),'format':'json','stages':['questions']}]
    base['policies']['questions'].update(templates=['design_by_variant'],relatedness_policy='strict',control_roles=['intended_activity','mock','assay_qc'])
    base['policies']['planning'].update(max_expansions=300000,beam_width=64,max_alternatives_per_question=16,exact_action_cap=16)
    validate_study(base);(directory/'study.json').write_text(json.dumps(base,indent=2)+'\n')
    study=resolve_study(directory/'study.json');output=generate_questions(study,[])
    records=output['question-records.jsonl'];validate_catalog(base['catalogs']+records)
    labels={'donors':{k:v['record_id'] for k,v in donors.items()},'designs':{k:v['record_id'] for k,v in designs.items()},'loci':{'T':intended['record_id'],'L':candidate['record_id']},'endpoints':{'T':intended_endpoint['record_id'],'L':endpoint['record_id']},'sequences':{k:v['record_id'] for k,v in sequences.items()},'context_id':context['record_id'],'expected':{'complete_cost':18,'insufficient_budget':17,'treated_panels':4,'mock_panels':2,'treatments':4,'independent_donors_per_arm':1,'technical_replicates':1,'evidence_status':'no_measurements','planning_status':'fully_planned'}}
    labels['questions']=[row['record_id'] for row in records if row['record_type']=='comparison_question']
    labels['actions']={name+':'+('mock' if row['design_id']['value'] is None else next(label for label,design in designs.items() if design['record_id']==row['design_id']['value']))+':'+row['operation_type']:row['record_id'] for name,donor in donors.items() for row in records if row['record_type']=='comparison_action' and row['donor_id']==donor['record_id']}
    (directory/'records.json').write_text(json.dumps(records,indent=2)+'\n');(directory/'labels.json').write_text(json.dumps(labels,indent=2)+'\n')
    (directory/'question-audit.json').write_text(json.dumps(output['question-audit.json'],indent=2)+'\n')
    return {'study':study,'records':records,'labels':labels}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path(__file__).parent);args=parser.parse_args()
    value=build_fixture(args.output);print(json.dumps(value['labels']['expected'],indent=2))
