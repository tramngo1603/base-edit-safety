"""Generate prospective questions and priced operations from declared packages."""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from itertools import combinations
from pathlib import Path

from .common import PopulationError, known, read_json, unknown
from .comparison_contracts import active_observations, make_record, record_identity, validate_catalog, validate_schema
from .comparison_identity import canonical_bytes, content_id
from .comparison_outcomes import outcome_reasons
from .comparison_provenance import verify_sources
from .comparison_questions import (_QuestionIndex, _actual_copies, _kind, _pairs,
    _observation_slot_reasons, donor_membership, evaluate_requirements, prepare_question_set,
    required_design_measurements)


def read_experiment_policy(study, index):
    sources = [row for row in study["sources"].values() if row["role"] == "experiment_policy"]
    if len(sources) != 1:
        raise PopulationError("Question generation requires exactly one experiment_policy input")
    source = sources[0]
    policy = read_json(source["path"])
    schema = read_json(Path(__file__).with_name("schemas") / "comparison-experiment-policy.schema.json")
    validate_schema(policy, schema, "experiment_policy")
    if len({rule["rule_id"] for rule in policy["rules"]}) != len(policy["rules"]):
        raise PopulationError("experiment_policy.rules: rule_id must be distinct")
    wanted = [(policy["context_id"], "experimental_context"), *[(v, "donor") for v in policy["donor_ids"]]]
    for rule in policy["rules"]:
        wanted.extend((v, "design") for v in rule["design_ids"])
        wanted.extend((rule[key], "endpoint_definition") for key in ("candidate_endpoint_id", "intended_endpoint_id"))
        if "evidence_replication" in rule["templates"] and not rule["replication_dimensions"]:
            raise PopulationError("experiment_policy.rules.replication_dimensions: replication requires an explicit new dimension")
    for identifier, kind in wanted:
        if identifier not in index or index[identifier]["record_type"] != "comparison_" + kind:
            raise PopulationError("experiment_policy has an unresolved " + kind + " reference")
    if policy["costs"]["cost_unit"] != study["config"]["policies"]["planning"]["cost_unit"]:
        raise PopulationError("experiment_policy costs use another planning cost unit")
    if study["config"]["mode"] == "research_sample_planning" and policy["specimen_policy"]["create_simulation_specimens"]:
        raise PopulationError("Research specimens cannot be created from public donor genomes")
    return policy, source


def _events(index, rule, *, synthetic=False):
    """Suggest comparison questions from recorded source facts."""
    by_locus = defaultdict(list)
    for configuration in _kind(index, "design_configuration").values():
        if configuration["design_id"] in rule["design_ids"]:
            sequence = index[configuration["sequence_configuration_id"]]
            by_locus[sequence["locus_id"]].append(configuration)
    result = []
    summaries = list(_kind(index, "donor_design_summary").values())
    nominations = {identifier for row in summaries for identifier in row["nominated_configuration_ids"]}
    active_ids = {row["record_id"] for row in active_observations(list(index.values()))}
    for locus_id, configurations in sorted(by_locus.items()):
        configurations.sort(key=lambda row: row["record_id"])
        off_target = [row for row in configurations if row["role"] == "off_target" and (row["record_id"] in nominations or synthetic and not summaries)]
        variant_rows = [row for row in configurations if index[row["sequence_configuration_id"]]["variants"]]
        intended_rows = [row for row in variant_rows if row["role"] == "intended"]
        candidate_rows = [row for row in variant_rows if row in off_target]
        observations = [row for row in _kind(index, "observation_revision").values() if row["record_id"] in active_ids and locus_id in row["locus_ids"] and row["design_id"]["value"] in rule["design_ids"]]
        predictions = [row for row in _kind(index, "prediction").values() if row["pre_edit_sequence_id"] in {v["sequence_configuration_id"] for v in configurations}]
        trigger = rule["trigger"]
        chosen = configurations
        if trigger == "candidate_nomination": chosen = off_target
        elif trigger == "variant_candidate": chosen = candidate_rows
        elif trigger == "variant_intended": chosen = intended_rows
        elif trigger == "search_profile_difference":
            chosen = off_target if {row["design_id"] for row in off_target} != set(rule["design_ids"]) else []
        elif trigger == "existing_observation": chosen = configurations if observations else []
        elif trigger == "prediction_followup": chosen = configurations if predictions else []
        elif trigger == "missing_functional_evidence":
            chosen = off_target if not any(row["endpoint_id"] == rule["candidate_endpoint_id"] and row["normalized_value"]["value"] is not None and not outcome_reasons(row) for row in observations) else []
        elif trigger == "observed_conflict":
            # Differing replicate values alone do not establish a contradiction.
            conflicts = {identifier for row in _kind(index, "decision").values()
                         if "conflicting_observations" in row["reason_codes"]
                         for identifier in row["observation_ids"]}
            chosen = configurations if any(row["record_id"] in conflicts or
                row["source_reported_outcome"]["value"] == "conflicting" for row in observations) else []
        if not chosen:
            continue
        # Separate different allele predicates at one locus instead of pooling them.
        by_sequence = defaultdict(list)
        for row in chosen:
            by_sequence[row["sequence_configuration_id"]].append(row)
        for sequence_id, rows in sorted(by_sequence.items()):
            result.append({"locus_id":locus_id,"sequence_id":sequence_id,
                "configuration_ids":[row["record_id"] for row in rows],
                "observation_ids":sorted({row["source_record_id"] for row in observations}),
                "prediction_ids":sorted(row["record_id"] for row in predictions)})
    return result


class _Builder:
    def __init__(self, study, records, policy, source):
        self.study, self.policy = study, policy
        self.index = _QuestionIndex(records)
        self.initial = set(self.index)
        self.audit = {"policy_sha256":source["sha256"],"suggestions":[],"question_exclusions":[],"donor_exclusions":[],"assumptions":[],"sequence_preparation":{}}
        artifact = self.add("source_artifact", sha256=source["sha256"], source_id="prospective-experiment-policy", release=policy["policy_version"], uri="urn:bes:experiment-policy:"+source["sha256"], media_type="application/json", scope="complete_file",parent_artifact_id=unknown("original declared policy"),locator={"kind":"file","member":Path(source["path"]).name},access="synthetic" if study["config"]["mode"] == "synthetic_test" else "restricted",use_terms="Declared prospective software planning assumptions")
        self.refs = [{"source_artifact_id":artifact["record_id"],"locator":{"kind":"file","member":Path(source["path"]).name},"note":"Declared prospective capability, cost and specimen assumptions; no observed outcome"}]
        self.assays, self.specimens, self.resource_cache = {}, {}, {}
        self.questions, self.slots, self.actions = {}, {}, {}
        self.setup = self.resource("study_setup", "study", policy["costs"]["study_setup"], unbounded=True, activation=policy["resources"]["study_setup_activation"])
        self.sequencing = [self.resource("sequencing_batch", "sequencing", policy["costs"]["sequencing_batch"], capacity=policy["resources"]["sequencing_capacity"],pool="declared-sequencing:"+policy["policy_version"],instance=i,maximum=policy["resources"]["max_sequencing_instances"]) for i in range(policy["resources"]["max_sequencing_instances"])]

    def add(self, kind, **fields):
        row = make_record(kind, **fields)
        old = self.index.get(row["record_id"])
        if old is not None and old != row:
            raise PopulationError("Generated canonical record conflicts with an existing definition")
        self.index[row["record_id"]] = row
        return row

    def resource(self, kind, key, cost, *, donor=None, design=None, assay=None, capacity=1, unbounded=False, activation="on_use", pool=None, instance=0, maximum=1):
        compatibility = {name:[] for name in ("donor_ids","design_ids","context_ids","batch_ids","assay_capability_ids","timepoints","endpoint_ids","locus_ids")}
        compatibility.update(donor_ids=[donor] if donor else [],design_ids=[design] if design else [], context_ids=[self.policy["context_id"]],batch_ids=[self.policy["batch_id"]],assay_capability_ids=[assay] if assay else [])
        return self.add("resource", resource_type=kind,cost=cost,cost_unit=self.policy["costs"]["cost_unit"],cost_scale=self.policy["costs"]["cost_scale"],sharing_key=self.policy["policy_version"]+":"+kind+":"+key,capacities=[{"dimension":"uses","capacity":{"kind":"unbounded","justification":"Explicit nonconsumable activation package; no extra assay capacity"} if unbounded else {"kind":"finite","quantity":capacity}}],activation_rule=activation,compatibility=compatibility,pool_id=known(pool) if pool else unknown("one named instance"),instance_index=instance,max_pool_instances=maximum)

    def assay(self, pairs):
        key = canonical_bytes(sorted(pairs, key=canonical_bytes))
        if key not in self.assays:
            policy = self.policy["assay_policy"]
            self.assays[key] = self.add("assay_capability",name=self.policy["policy_version"]+" declared panel "+content_id("comparison_panel",pairs).rsplit(":",1)[1][:16],endpoint_ids=sorted({p["endpoint_id"] for p in pairs}),locus_ids=sorted({p["locus_id"] for p in pairs}),endpoint_locus_pairs=pairs,context_ids=[self.policy["context_id"]],resolution=policy["resolution"],allele_resolved=policy["allele_resolved"],feasibility=policy["feasibility"],required_control_roles=["intended_activity","mock","assay_qc"],qc_requirements=policy["qc_requirements"],max_targets=len({p["locus_id"] for p in pairs}),source_refs=self.refs)
        return self.assays[key]

    def specimen(self, donor):
        if donor in self.specimens:
            return self.specimens[donor]
        candidates = [row for row in self.index.values() if row["record_type"] == "comparison_specimen" and row["donor_id"] == donor and self.policy["context_id"] in row["context_ids"]]
        if candidates:
            permitted = [row for row in candidates if row["availability"] in ({"verified"} if self.study["config"]["mode"] == "research_sample_planning" else {"verified","simulation_assumption"}) and row.get("use_authorized") is True]
            specimen = min(permitted or candidates,key=lambda row:row["record_id"])
        elif self.policy["specimen_policy"]["create_simulation_specimens"]:
            policy = self.policy["specimen_policy"]
            specimen = self.add("specimen",donor_id=donor,source_specimen_id="declared-simulation:"+self.policy["policy_version"]+":"+donor,context_ids=[self.policy["context_id"]],availability="simulation_assumption",permitted_use=known(policy["permitted_use"]),use_authorized=policy["use_authorized"],capacity=policy["capacity"],source_refs=self.refs)
            self.audit["assumptions"].append({"donor_id":donor,"specimen_id":specimen["record_id"],"reason":policy["assumption"]})
        else:
            return None
        self.specimens[donor] = specimen
        return specimen

    def controls(self, pairs, assay, intended_pairs=()):
        result = []
        intended = {canonical_bytes(pair) for pair in intended_pairs}
        for pair in pairs:
            for role in (["intended_activity"] if canonical_bytes(pair) in intended else []) + ["mock","assay_qc"]:
                result.append({"role":role,"endpoint_ids":[pair["endpoint_id"]],"locus_ids":[pair["locus_id"]],"assay_capability_id":assay["record_id"],"required_qc":self.policy["assay_policy"]["qc_requirements"]})
        return result

    def question(self, rule, event, template, designs):
        intended = {locus for design in designs for locus in self.index[design]["intended_locus_ids"]}
        pairs = [{"endpoint_id":rule["candidate_endpoint_id"],"locus_id":event["locus_id"]}]
        intended_pairs = [{"endpoint_id":rule["intended_endpoint_id"],"locus_id":locus} for locus in sorted(intended)]
        pairs += intended_pairs
        pairs = list({canonical_bytes(pair):pair for pair in pairs}.values())
        assay = self.assay(pairs)
        sequence = self.index[event["sequence_id"]]
        contrasted = template in {"variant_association","design_by_variant"}
        predicate = {"kind":"allele" if len(sequence["variants"]) == 1 else "haplotype" if sequence["variants"] else "configuration", "locus_id":event["locus_id"],"sequence_configuration_id":known(sequence["record_id"]) if not sequence["variants"] else unknown("literal allele predicate"),"variants":sequence["variants"],"dosage":unknown("one or two matching copies establish carriage"),"phase_requirement":"confirmed_cis" if len(sequence["variants"]) > 1 else "not_required","coverage_basis":self.study["config"]["coverage_basis"]}
        if not contrasted:
            predicate.update(kind="any_resolved",sequence_configuration_id=unknown("all donors with resolved actual DNA"),variants=[],phase_requirement="not_required")
        motivations = event["observation_ids"] if template == "evidence_replication" else []
        if template == "evidence_replication" and not motivations:
            return None
        fields = dict(template=template,hypothesis="Evaluate "+template+" at "+event["locus_id"]+" under declared rule "+rule["rule_id"],decision="Measure the declared endpoint; genotype contrasts establish association, not variant causality",design_ids=designs,locus_ids=sorted({p["locus_id"] for p in pairs}),endpoint_ids=sorted({p["endpoint_id"] for p in pairs}),endpoint_locus_pairs=pairs,context_ids=[self.policy["context_id"]],genotype_predicate=predicate,comparability={"require_same_context":True,"require_same_endpoint":True,"require_same_timepoint":True},arms=[{"arm_id":name,"membership":name,"required_independent_donors":rule["required_donors_per_arm"]} for name in (["carrier","noncarrier"] if contrasted else ["all_resolved"])],pairing="within_donor" if len(designs)==2 else "none",relatedness_policy=self.policy["relatedness_policy"],technical_replicates=rule["technical_replicates"],control_requirements=self.controls(pairs,assay,intended_pairs),requirement_slot_ids=[],evidence_reuse={"allow_existing":rule["allow_existing"],"allow_across_questions":rule["allow_across_questions"],"excluded_observation_ids":motivations,"excluded_experiment_ids":sorted({self.index[identifier]["source_experiment"] for identifier in motivations}),"require_new":rule["replication_dimensions"] if template=="evidence_replication" else []},matching_covariates=rule["matching_covariates"],priority_weight=rule["priority_weight"],hard_required=rule["hard_required"],rationale="Source-bound suggestion from "+rule["trigger"]+"; "+self.policy["assay_policy"]["package_basis"],policy_version=self.policy["policy_version"]+":"+rule["rule_id"],motivating_observation_ids=motivations,motivating_configuration_ids=event["configuration_ids"],motivating_prediction_ids=event["prediction_ids"],planning_status="unplanned",evidence_status="no_measurements",remaining_gaps=["prospective_measurements_required"])
        row = {"record_type":"comparison_question","schema_version":"1",**fields}
        identifier = record_identity(row)
        row.update(record_id=identifier,question_id=identifier)
        if identifier in self.questions:
            previous = self.questions[identifier]
            for field in ("motivating_observation_ids", "motivating_configuration_ids", "motivating_prediction_ids"):
                row[field] = sorted(set(row[field]) | set(previous[field]))
        self.questions[identifier] = row
        self.audit["suggestions"].append({"question_id":identifier,"trigger":rule["trigger"],"rule_id":rule["rule_id"],"configuration_ids":event["configuration_ids"],"observation_ids":event["observation_ids"],"prediction_ids":event["prediction_ids"]})
        return row

    def action(self, operation, donor, specimen, design, pairs, assay, repeat, dependencies=(), *, intended_pairs=()):
        costs, limits = self.policy["costs"], self.policy["resources"]
        targets = []
        for locus in sorted({pair["locus_id"] for pair in pairs}):
            copies, reasons = _actual_copies(self.index, donor, locus, self.study["config"]["coverage_basis"])
            if reasons:
                return None
            targets.append({"locus_id":locus,"sequence_configuration_ids":sorted({copy["sequence"]["record_id"] for copy in copies.values()})})
        key = donor+":"+str(design)+":"+assay["record_id"]+":"+str(repeat)
        activation = self.resource("donor_activation",donor,costs["donor_activation"],donor=donor,unbounded=True)
        preparation = self.resource("specimen_preparation",specimen["record_id"],costs["specimen_preparation"],donor=donor,capacity=min(limits["specimen_capacity"],specimen["capacity"]))
        demands = [{"resource_id":row["record_id"],"dimension":"uses","quantity":1} for row in (self.setup,activation,preparation)]
        if operation == "treatment":
            resource = self.resource("treatment",key,costs["treatment"],donor=donor,design=design,capacity=1)
        elif operation == "mock_control":
            resource = self.resource("mock_control",key,costs["mock_panel"],donor=donor,design=design,assay=assay["record_id"],capacity=limits["mock_capacity"])
        else:
            resource = self.resource("assay_library",key,costs["treated_panel"],donor=donor,design=design,assay=assay["record_id"],capacity=limits["panel_capacity"])
        demands.append({"resource_id":resource["record_id"],"dimension":"uses","quantity":1})
        if operation != "treatment":
            # Each library consumes one slot in one explicit paid batch instance.
            if len(self.sequencing)==1:
                demands.append({"resource_id":self.sequencing[0]["record_id"],"dimension":"uses","quantity":1})
            else:
                demands.append({"eligible_resource_ids":[row["record_id"] for row in self.sequencing],"assignment_key":"sequencing","dimension":"uses","quantity":1})
            if repeat > 1:
                extra=self.resource("additional_replicate",key,costs["additional_replicate"],donor=donor,capacity=1)
                demands.append({"resource_id":extra["record_id"],"dimension":"uses","quantity":1})
        controls = self.controls(pairs,assay,intended_pairs) if operation == "measurement_panel" else []
        roles = [] if operation=="treatment" else ["mock","assay_qc"] if operation=="mock_control" else ["intended_activity","assay_qc"]
        row = make_record("action",operation_type=operation,donor_id=donor,specimen_id=specimen["record_id"],design_id=known(design) if design else unknown("shared mock is design neutral"),context_id=self.policy["context_id"],assay_capability_id=assay["record_id"],endpoint_ids=sorted({p["endpoint_id"] for p in pairs}),locus_ids=sorted({p["locus_id"] for p in pairs}),endpoint_locus_pairs=pairs,target_sequence_id=unknown("exact per-locus diploid sequences are listed"),target_sequences=targets,replicate_position=repeat,batch_id=self.policy["batch_id"],experiment_id=self.policy["experiment_id"],dependency_action_ids=sorted(set(dependencies)),control_requirements=controls,resource_demands=demands,eligibility_status="eligible",reason_codes=[],memberships=[],performed_control_roles=roles,performed_qc=[] if operation=="treatment" else self.policy["assay_policy"]["qc_requirements"])
        old=self.actions.get(row["record_id"])
        if old and {k:v for k,v in old.items() if k!='memberships'} != {k:v for k,v in row.items() if k!='memberships'}:
            raise PopulationError("Same biological action has contradictory prospective packages")
        self.actions.setdefault(row["record_id"],row)
        return self.actions[row["record_id"]]

    def fill(self, question, allow_actions):
        if not allow_actions:
            question.update(construction_status="excluded_by_declared_limit",planning_status="blocked",remaining_gaps=["requirement_construction_shortlist_exclusion"])
            return
        question["construction_status"] = "complete"
        membership = {donor:donor_membership(self.index,question,donor,self.study["config"]["mode"],self.study["config"]["coverage_basis"]) for donor in self.policy["donor_ids"]}
        pairs = _pairs(question); assay = self.assay(pairs)
        selected=set()
        for arm in question["arms"]:
            candidates=sorted(donor for donor,state in membership.items() if state["eligible"] and (arm["membership"]=="all_resolved" or state["membership"]==arm["membership"]))
            selected.update(candidates[:self.policy["max_donor_alternatives_per_arm"]] if allow_actions else [])
        for donor,state in sorted(membership.items()):
            arm_id=state["membership"] if question["template"] in {"variant_association","design_by_variant"} else "all_resolved"
            arm_id=arm_id if arm_id in {arm["arm_id"] for arm in question["arms"]} else question["arms"][0]["arm_id"]
            specimen = self.specimen(donor) if donor in selected else None
            if donor not in selected or specimen is None:
                self.audit["donor_exclusions"].append({"question_id":question["record_id"],"donor_id":donor,"reason_codes":state["reason_codes"] or ["canonical_action_shortlist_exclusion" if donor not in selected else "no_declared_specimen"]})
            if donor not in selected:
                continue
            for repeat in range(1,question["technical_replicates"]+1):
                mock_shared=None
                if specimen and self.policy["assay_policy"]["mock_shared_across_designs"]:
                    mock_shared=self.action("mock_control",donor,specimen,None,pairs,assay,repeat)
                for design in question["design_ids"]:
                    measurements = required_design_measurements(self.index, question, design)
                    covered = [{key:pair[key] for key in ("endpoint_id","locus_id")} for pair in measurements]
                    intended_pairs = [pair for pair, requirement in zip(covered, measurements) if requirement["control_role"] == "intended_activity"]
                    mock=mock_shared
                    if specimen and not self.policy["assay_policy"]["mock_shared_across_designs"]:
                        mock=self.action("mock_control",donor,specimen,design,covered,assay,repeat)
                    treatment=self.action("treatment",donor,specimen,design,covered,assay,repeat) if specimen else None
                    panel=self.action("measurement_panel",donor,specimen,design,covered,assay,repeat,[treatment["record_id"],mock["record_id"]],intended_pairs=intended_pairs) if treatment and mock else None
                    for pair in measurements:
                        main_role=pair["control_role"]
                        for role,action,slot_design in ((main_role,panel,design),("assay_qc",panel,design),("mock",mock,None if self.policy["assay_policy"]["mock_shared_across_designs"] else design)):
                            row=make_record("requirement_slot",question_id=question["record_id"],donor_id=donor,arm_id=arm_id,design_id=known(slot_design) if slot_design else unknown("shared design-neutral mock"),endpoint_id=pair["endpoint_id"],locus_id=known(pair["locus_id"]),context_id=self.policy["context_id"],control_role=role,replicate_position=repeat,allowed_observation_ids=[],allowed_action_ids=[action["record_id"]] if action else [])
                            self.slots[row["record_id"]]=row
                            question["requirement_slot_ids"].append(row["record_id"])
                            if action:
                                entry={"question_id":question["record_id"],"slot_id":row["record_id"]}
                                if entry not in action["memberships"]:action["memberships"].append(entry)
        question["requirement_slot_ids"]=sorted(set(question["requirement_slot_ids"]))
        if not question["requirement_slot_ids"]:
            question.update(planning_status="blocked",remaining_gaps=["no_eligible_declared_donor"])


def generate_questions(study, records):
    """Return immutable canonical inputs plus a source-bound eligibility audit."""
    verify_sources(study)
    merged={row["record_id"]:deepcopy(row) for row in study["config"]["catalogs"]}
    for row in records:
        if row["record_id"] in merged and merged[row["record_id"]]!=row:
            raise PopulationError("Conflicting parent record during question generation")
        merged[row["record_id"]]=deepcopy(row)
    parent_ids = set(merged)
    for source in study["sources"].values():
        if source["role"]=="question_inventory":
            for row in read_json(source["path"]):
                if row["record_id"] in merged and merged[row["record_id"]]!=row:raise PopulationError("Conflicting question inventory record")
                merged[row["record_id"]]=row
    index=validate_catalog(list(merged.values()))
    policy,source=read_experiment_policy(study,index)
    builder=_Builder(study,list(index.values()),policy,source)
    builder.initial = parent_ids
    for rule in policy["rules"]:
        for event in _events(index,rule,synthetic=study["config"]["mode"]=="synthetic_test"):
            for template in rule["templates"]:
                sets=list(combinations(sorted(rule["design_ids"]),2)) if template in {"paired_design_comparison","design_by_variant"} else [(design,) for design in sorted(rule["design_ids"])]
                for designs in sets:builder.question(rule,event,template,list(designs))
    if not builder.questions:
        return {"question-records.jsonl":[],"question-audit.json":{**builder.audit,"reason_codes":["no_source_facts_match_declared_rules"]},"question-states.json":{}}
    chosen=set(sorted(builder.questions)[:policy["max_action_questions"]])
    for identifier,question in sorted(builder.questions.items()):
        if identifier not in chosen:builder.audit["question_exclusions"].append({"question_id":identifier,"reason":"canonical_question_shortlist_exclusion"})
        builder.fill(question,identifier in chosen)
    rows=list(builder.index.values())+list(builder.questions.values())+list(builder.slots.values())+list(builder.actions.values())
    qset=prepare_question_set(study,rows)
    # Offer only source-compatible active measurements; explicit reuse is still checked jointly.
    active = {row["record_id"] for row in active_observations(rows)}
    for slot in builder.slots.values():
        slot["allowed_observation_ids"]=[identifier for identifier,row in _kind(qset.index,"observation_revision").items() if identifier in active and not _observation_slot_reasons(qset,row,slot)]
    rows=list(builder.index.values())+list(builder.questions.values())+list(builder.slots.values())+list(builder.actions.values())
    qset=prepare_question_set(study,rows)
    states=evaluate_requirements(qset)
    builder.audit["action_eligibility"]=qset.eligibility
    builder.audit["donor_eligibility"]=[{"question_id":key[0],"donor_id":key[1],**value} for key,value in sorted(qset.donor_eligibility.items())]
    builder.audit["input_snapshot_run_id"]=study["snapshot_run_id"]
    builder.audit["code"]=study["code"]
    verify_sources(study)
    generated=[row for row in rows if row["record_id"] not in builder.initial]
    return {"question-records.jsonl":sorted(generated,key=lambda row:(row["record_type"],row["record_id"])),"question-audit.json":builder.audit,"question-states.json":states}
