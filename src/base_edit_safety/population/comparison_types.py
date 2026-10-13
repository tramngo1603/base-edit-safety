"""Static input types for the closed comparison schemas.

Validate JSON with comparison_contracts before treating it as these records.
Reason-bearing values retain an explicit reason whenever their value is null.
"""
from __future__ import annotations

from typing import Literal, NotRequired, TypedDict, Union

class Interval(TypedDict):
    contig: str
    start: int
    end: int


class PathPointKnown(TypedDict):
    kind: Literal['reference']
    pos0: int
    base: str


class PathPointUnknown(TypedDict):
    kind: Literal['insertion']
    boundary0: int
    offset0: int
    base: str


PathPoint = Union[PathPointKnown, PathPointUnknown]


class LocatorKnown(TypedDict):
    kind: Literal['table']
    sheet: str
    table: str
    row: int
    column: str
    cell: str


class LocatorUnknown(TypedDict):
    kind: Literal['text']
    page: int
    line_start: int
    line_end: int


class Locator2(TypedDict):
    kind: Literal['file']
    member: str


Locator = Union[LocatorKnown, LocatorUnknown, Locator2]


class SourceRef(TypedDict):
    source_artifact_id: str
    locator: Locator
    note: str


SourceRefs = list[SourceRef]


class ResolvedTextKnown(TypedDict):
    value: str
    reason: None


class ResolvedTextUnknown(TypedDict):
    value: None
    reason: str


ResolvedText = Union[ResolvedTextKnown, ResolvedTextUnknown]


class QuantityValueKnown(TypedDict):
    value: Union[int, float]
    reason: None


class QuantityValueUnknown(TypedDict):
    value: None
    reason: str


class Quantity(TypedDict):
    value: Union[QuantityValueKnown, QuantityValueUnknown]
    unit: str


class Edit(TypedDict):
    locus_id: str
    pos0: int
    reference_base: str
    alternate_base: str
    strand: Literal['+', '-']
    purpose: Literal['intended', 'bystander', 'scenario']


class Variant(TypedDict):
    contig: str
    pos0: int
    reference: str
    alternates: list[str]
    selected_alt_index: int
    source_artifact_id: str
    source_row: int


class Assessment(TypedDict):
    status: Literal['applicable', 'related', 'unresolved', 'incompatible', 'not_evaluated']
    reason_codes: list[str]


class ReusePolicy(TypedDict):
    allow_existing: bool
    allow_across_questions: bool
    require_new: list[Literal['experiment', 'batch', 'specimen', 'donor']]
    excluded_observation_ids: list[str]
    excluded_experiment_ids: list[str]


class GenotypePredicateSequenceConfigurationIdKnown(TypedDict):
    value: str
    reason: None


class GenotypePredicateSequenceConfigurationIdUnknown(TypedDict):
    value: None
    reason: str


class GenotypePredicateDosageKnown(TypedDict):
    value: int
    reason: None


class GenotypePredicateDosageUnknown(TypedDict):
    value: None
    reason: str


class GenotypePredicate(TypedDict):
    kind: Literal['any_resolved', 'allele', 'haplotype', 'configuration']
    locus_id: str
    sequence_configuration_id: Union[GenotypePredicateSequenceConfigurationIdKnown, GenotypePredicateSequenceConfigurationIdUnknown]
    variants: list[Variant]
    dosage: Union[GenotypePredicateDosageKnown, GenotypePredicateDosageUnknown]
    phase_requirement: Literal['not_required', 'confirmed_cis', 'confirmed_trans']
    coverage_basis: Literal['synthetic_truth', 'callable_supported', 'retained_callset_model']


class Compatibility(TypedDict):
    donor_ids: list[str]
    design_ids: list[str]
    context_ids: list[str]
    batch_ids: list[str]
    assay_capability_ids: list[str]
    timepoints: list[str]
    endpoint_ids: list[str]
    locus_ids: list[str]


class ResourceDemand(TypedDict):
    resource_id: str
    dimension: str
    quantity: int


class ControlRequirement(TypedDict):
    role: str
    endpoint_ids: list[str]
    locus_ids: list[str]
    assay_capability_id: str
    required_qc: list[str]


class Membership(TypedDict):
    question_id: str
    slot_id: str


class Applicability(TypedDict):
    identity: Assessment
    context: Assessment
    endpoint: Assessment
    quality: Assessment
    classification: Literal['exact_donor_configuration', 'exact_sequence_context', 'related_design_context', 'unresolved', 'incompatible']
    empirical_donor_credit: bool
    reason_codes: list[str]


class SourceArtifactParentArtifactIdKnown(TypedDict):
    value: str
    reason: None


class SourceArtifactParentArtifactIdUnknown(TypedDict):
    value: None
    reason: str


class SourceArtifact(TypedDict):
    record_type: Literal['comparison_source_artifact']
    schema_version: Literal['1']
    record_id: str
    source_artifact_id: str
    sha256: str
    source_id: str
    release: str
    uri: str
    media_type: str
    scope: Literal['complete_file', 'indexed_window', 'extracted_member']
    parent_artifact_id: Union[SourceArtifactParentArtifactIdKnown, SourceArtifactParentArtifactIdUnknown]
    locator: Locator
    access: Literal['public', 'restricted', 'synthetic']
    use_terms: str


class ReferenceBundleWindowsItem(TypedDict):
    interval: Interval
    source_sha256: str


class ReferenceBundleCompleteReferenceSha256Known(TypedDict):
    value: str
    reason: None


class ReferenceBundleCompleteReferenceSha256Unknown(TypedDict):
    value: None
    reason: str


class ReferenceBundle(TypedDict):
    record_type: Literal['comparison_reference_bundle']
    schema_version: Literal['1']
    record_id: str
    reference_bundle_id: str
    assembly: str
    identity_kind: Literal['complete_reference_files', 'indexed_windows']
    source_artifact_ids: list[str]
    source_checksums: list[str]
    windows: list[ReferenceBundleWindowsItem]
    complete_reference_sha256: Union[ReferenceBundleCompleteReferenceSha256Known, ReferenceBundleCompleteReferenceSha256Unknown]
    transformation_version: str


class EditorDefinitionGuideGeometry(TypedDict):
    spacer_length: int
    pam_side: Literal['3prime', '5prime']
    pam_length: int


class EditorDefinitionConversion(TypedDict):
    substrate_base: str
    product_base: str
    strand_convention: Literal['guide_oriented', 'genomic_forward']


class EditorDefinitionEditingWindowKnown(TypedDict):
    start_1based: int
    end_1based: int
    basis: str


class EditorDefinitionEditingWindowUnknown(TypedDict):
    value: None
    reason: str


class EditorDefinition(TypedDict):
    record_type: Literal['comparison_editor_definition']
    schema_version: Literal['1']
    record_id: str
    editor_definition_id: str
    display_name: str
    variant_name: str
    editor_class: Literal['ABE', 'CBE', 'dual', 'other_registered']
    cas_component: str
    pam_patterns: list[str]
    guide_geometry: EditorDefinitionGuideGeometry
    conversion: EditorDefinitionConversion
    editing_window: Union[EditorDefinitionEditingWindowKnown, EditorDefinitionEditingWindowUnknown]
    source_refs: SourceRefs


class LocusGeneKnown(TypedDict):
    value: str
    reason: None


class LocusGeneUnknown(TypedDict):
    value: None
    reason: str


class Locus(TypedDict):
    record_type: Literal['comparison_locus']
    schema_version: Literal['1']
    record_id: str
    locus_id: str
    reference_bundle_id: str
    contig: str
    start: int
    end: int
    strand: Literal['+', '-']
    locus_name: str
    gene: Union[LocusGeneKnown, LocusGeneUnknown]
    path: list[PathPoint]


class IntendedLocus(TypedDict):
    record_type: Literal['comparison_intended_locus']
    schema_version: Literal['1']
    record_id: str
    intended_locus_id: str
    reference_bundle_id: str
    contig: str
    start: int
    end: int
    strand: Literal['+', '-']
    gene: str
    locus_name: str
    path: list[PathPoint]
    intended_endpoint: str
    copy_resolution: Literal['copy_specific', 'combined_copies']
    source_refs: SourceRefs


class Design(TypedDict):
    record_type: Literal['comparison_design']
    schema_version: Literal['1']
    record_id: str
    design_id: str
    definition_hash: str
    display_name: str
    guide_sequence: str
    editor_definition_id: str
    intended_locus_ids: list[str]
    desired_edit_definitions: list[Edit]
    source_status: Literal['published_verified', 'proposed', 'synthetic']
    source_refs: SourceRefs


class ExperimentalContextCellStateKnown(TypedDict):
    value: str
    reason: None


class ExperimentalContextCellStateUnknown(TypedDict):
    value: None
    reason: str


class ExperimentalContextCellProvenanceKnown(TypedDict):
    value: str
    reason: None


class ExperimentalContextCellProvenanceUnknown(TypedDict):
    value: None
    reason: str


class ExperimentalContext(TypedDict):
    record_type: Literal['comparison_experimental_context']
    schema_version: Literal['1']
    record_id: str
    context_id: str
    organism: str
    cell_type: str
    cell_state: Union[ExperimentalContextCellStateKnown, ExperimentalContextCellStateUnknown]
    cell_provenance: Union[ExperimentalContextCellProvenanceKnown, ExperimentalContextCellProvenanceUnknown]
    delivery_method: str
    exposure: Quantity
    dose: Quantity
    condition: str
    measurement_timepoint: Quantity
    source_refs: SourceRefs


class EndpointDefinition(TypedDict):
    record_type: Literal['comparison_endpoint_definition']
    schema_version: Literal['1']
    record_id: str
    endpoint_id: str
    name: Literal['intended_edit_fraction', 'unintended_edit_fraction', 'bystander_edit_fraction', 'indel_fraction', 'HBG_expression', 'HbF_fraction', 'F_cell_fraction', 'globin_chain_ratio', 'cell_viability', 'other_registered_endpoint']
    biological_meaning: str
    unit: str
    denominator: str
    measurement_level: Literal['reads', 'cells', 'molecules', 'donor', 'pooled_group', 'study_summary']
    favorable_direction: Literal['higher', 'lower', 'not_defined']
    locus_resolution: Literal['copy_specific', 'locus_specific', 'combined_copies', 'not_applicable']
    allele_resolution: Literal['required', 'not_required']
    requires_sequence: bool
    required_context_ids: list[str]


class DesignAlias(TypedDict):
    record_type: Literal['comparison_design_alias']
    schema_version: Literal['1']
    record_id: str
    alias_id: str
    source_name: str
    design_id: str
    source_location: SourceRef
    justification: str


class Domain(TypedDict):
    record_type: Literal['comparison_domain']
    schema_version: Literal['1']
    record_id: str
    domain_id: str
    reference_bundle_id: str
    intervals: list[Interval]
    selection_basis: Literal['autosome_wide', 'prespecified_bounded', 'conditioned_on_nominations', 'synthetic']
    description: str


class SearchProfileIndexIdKnown(TypedDict):
    value: str
    reason: None


class SearchProfileIndexIdUnknown(TypedDict):
    value: None
    reason: str


class SearchProfile(TypedDict):
    record_type: Literal['comparison_search_profile']
    schema_version: Literal['1']
    record_id: str
    search_profile_id: str
    design_id: str
    domain_id: str
    pam_patterns: list[str]
    max_mismatches: int
    max_dna_bulges: int
    max_rna_bulges: int
    max_total_edits: int
    total_edit_accounting: Literal['mismatch_plus_bulged_bases', 'mismatch_plus_bulge_events']
    strands: list[Literal['+', '-']]
    variant_domain: Literal['literal_snv_indel', 'synthetic_literal']
    backend_name: str
    backend_version: str
    backend_digest: str
    index_id: Union[SearchProfileIndexIdKnown, SearchProfileIndexIdUnknown]
    source_panel_id: str
    source_equivalence: Literal['verified', 'candidate_only', 'not_applicable_synthetic']
    max_rows: int
    max_window_variants: int
    retention: Literal['all']


class SequenceConfiguration(TypedDict):
    record_type: Literal['comparison_sequence_configuration']
    schema_version: Literal['1']
    record_id: str
    sequence_configuration_id: str
    reference_bundle_id: str
    locus_id: str
    contig: str
    strand: Literal['+', '-']
    sequence: str
    path: list[PathPoint]
    variants: list[Variant]
    support_status: Literal['supported', 'modeled', 'putative', 'unresolved']
    support_reason: Union[str, None]


class SequenceOccurrencePhaseBlockKnown(TypedDict):
    value: str
    reason: None


class SequenceOccurrencePhaseBlockUnknown(TypedDict):
    value: None
    reason: str


class SequenceOccurrence(TypedDict):
    record_type: Literal['comparison_sequence_occurrence']
    schema_version: Literal['1']
    record_id: str
    sequence_occurrence_id: str
    sequence_configuration_id: str
    donor_id: str
    sample_id: str
    copy_index: int
    phase_block: Union[SequenceOccurrencePhaseBlockKnown, SequenceOccurrencePhaseBlockUnknown]
    phase_status: Literal['confirmed', 'putative', 'unphased', 'unresolved']
    source_version: str
    source_refs: SourceRefs


class DesignConfiguration(TypedDict):
    record_type: Literal['comparison_design_configuration']
    schema_version: Literal['1']
    record_id: str
    design_configuration_id: str
    design_id: str
    sequence_configuration_id: str
    role: Literal['intended', 'off_target', 'contextual']


class AlignmentColumnsItem(TypedDict):
    guide: Union[str, None]
    target: Union[str, None]
    guide_position_1based: Union[int, None]
    target_path_index: Union[int, None]


class Alignment(TypedDict):
    record_type: Literal['comparison_alignment']
    schema_version: Literal['1']
    record_id: str
    alignment_id: str
    design_configuration_id: str
    columns: list[AlignmentColumnsItem]
    mismatch_positions: list[int]
    dna_bulge_positions: list[int]
    rna_bulge_positions: list[int]
    pam_sequence: str
    pam_path: list[PathPoint]
    native_row_refs: SourceRefs


class DonorIndependenceGroupIdKnown(TypedDict):
    value: str
    reason: None


class DonorIndependenceGroupIdUnknown(TypedDict):
    value: None
    reason: str


class DonorCovariatesItemValueKnown(TypedDict):
    value: Union[str, Union[int, float]]
    reason: None


class DonorCovariatesItemValueUnknown(TypedDict):
    value: None
    reason: str


class DonorCovariatesItem(TypedDict):
    name: str
    value: Union[DonorCovariatesItemValueKnown, DonorCovariatesItemValueUnknown]


class Donor(TypedDict):
    record_type: Literal['comparison_donor']
    schema_version: Literal['1']
    record_id: str
    donor_id: str
    canonical_label: str
    source_panel_id: str
    sample_ids: list[str]
    independence_group_id: Union[DonorIndependenceGroupIdKnown, DonorIndependenceGroupIdUnknown]
    relatedness_policy: Literal['strict', 'declared_simulation_assumption']
    relatedness_basis: str
    source_refs: SourceRefs
    covariates: NotRequired[list[DonorCovariatesItem]]


class DonorConfiguration(TypedDict):
    record_type: Literal['comparison_donor_configuration']
    schema_version: Literal['1']
    record_id: str
    donor_configuration_id: str
    donor_id: str
    sample_id: str
    sequence_occurrence_id: str
    design_configuration_id: str
    copy_index: int
    carrier_status: Literal['supported_carrier', 'supported_noncarrier', 'modeled_carrier', 'modeled_noncarrier', 'unresolved', 'not_evaluated']
    carrier_of_sequence_id: str
    genotype_support: Literal['called', 'callable_reference', 'retained_callset_model', 'missing', 'filtered', 'not_evaluated']
    phase_support: Literal['confirmed', 'putative', 'unphased', 'unresolved']
    reason_codes: list[str]
    source_refs: SourceRefs


class SearchScopeDonorIdKnown(TypedDict):
    value: str
    reason: None


class SearchScopeDonorIdUnknown(TypedDict):
    value: None
    reason: str


class SearchScopeSkippedIntervalsItem(TypedDict):
    interval: Interval
    reason: str


class SearchScopeBackendLimits(TypedDict):
    max_rows: int
    max_window_variants: int
    timeout_seconds: int


class SearchScopeExitCodeKnown(TypedDict):
    value: int
    reason: None


class SearchScopeExitCodeUnknown(TypedDict):
    value: None
    reason: str


class SearchScope(TypedDict):
    record_type: Literal['comparison_search_scope']
    schema_version: Literal['1']
    record_id: str
    search_scope_id: str
    design_id: str
    donor_id: Union[SearchScopeDonorIdKnown, SearchScopeDonorIdUnknown]
    source_panel_id: str
    search_profile_id: str
    requested_intervals: list[Interval]
    evaluated_intervals: list[Interval]
    skipped_intervals: list[SearchScopeSkippedIntervalsItem]
    execution_status: Literal['complete', 'partial', 'failed', 'not_evaluated']
    retention_status: Literal['complete', 'truncated', 'unknown']
    backend_limits: SearchScopeBackendLimits
    exit_code: Union[SearchScopeExitCodeKnown, SearchScopeExitCodeUnknown]
    reason_codes: list[str]


class VariantEffectDonorConfigurationIdKnown(TypedDict):
    value: str
    reason: None


class VariantEffectDonorConfigurationIdUnknown(TypedDict):
    value: None
    reason: str


class VariantEffectReferenceConfigurationIdKnown(TypedDict):
    value: str
    reason: None


class VariantEffectReferenceConfigurationIdUnknown(TypedDict):
    value: None
    reason: str


class VariantEffectDonorAlignmentIdKnown(TypedDict):
    value: str
    reason: None


class VariantEffectDonorAlignmentIdUnknown(TypedDict):
    value: None
    reason: str


class VariantEffectReferenceAlignmentIdKnown(TypedDict):
    value: str
    reason: None


class VariantEffectReferenceAlignmentIdUnknown(TypedDict):
    value: None
    reason: str


class VariantEffect(TypedDict):
    record_type: Literal['comparison_variant_effect']
    schema_version: Literal['1']
    record_id: str
    variant_effect_id: str
    donor_id: str
    donor_configuration_id: Union[VariantEffectDonorConfigurationIdKnown, VariantEffectDonorConfigurationIdUnknown]
    reference_configuration_id: Union[VariantEffectReferenceConfigurationIdKnown, VariantEffectReferenceConfigurationIdUnknown]
    donor_alignment_id: Union[VariantEffectDonorAlignmentIdKnown, VariantEffectDonorAlignmentIdUnknown]
    reference_alignment_id: Union[VariantEffectReferenceAlignmentIdKnown, VariantEffectReferenceAlignmentIdUnknown]
    status: Literal['created', 'lost', 'changed', 'unpaired', 'unknown', 'unchanged']
    reason: str
    source_variants: list[Variant]


class IntendedAssessmentSequenceConfigurationIdKnown(TypedDict):
    value: str
    reason: None


class IntendedAssessmentSequenceConfigurationIdUnknown(TypedDict):
    value: None
    reason: str


class IntendedAssessmentPamCompatibleKnown(TypedDict):
    value: bool
    reason: None


class IntendedAssessmentPamCompatibleUnknown(TypedDict):
    value: None
    reason: str


class IntendedAssessment(TypedDict):
    record_type: Literal['comparison_intended_assessment']
    schema_version: Literal['1']
    record_id: str
    intended_assessment_id: str
    donor_id: str
    design_id: str
    locus_id: str
    sequence_configuration_id: Union[IntendedAssessmentSequenceConfigurationIdKnown, IntendedAssessmentSequenceConfigurationIdUnknown]
    support_status: Literal['supported', 'modeled', 'unresolved', 'not_evaluated']
    pam_compatible: Union[IntendedAssessmentPamCompatibleKnown, IntendedAssessmentPamCompatibleUnknown]
    alignment_ids: list[str]
    reference_differences: list[Variant]
    reason_codes: list[str]


class SourceObservationRawValueKnown(TypedDict):
    value: Union[Union[int, float], str]
    reason: None


class SourceObservationRawValueUnknown(TypedDict):
    value: None
    reason: str


class SourceObservationSourceReportedLabelKnown(TypedDict):
    value: str
    reason: None


class SourceObservationSourceReportedLabelUnknown(TypedDict):
    value: None
    reason: str


class SourceObservation(TypedDict):
    record_type: Literal['comparison_source_observation']
    schema_version: Literal['1']
    record_id: str
    source_record_id: str
    source_artifact_id: str
    locator: Locator
    measurement_key: str
    source_study: str
    source_experiment: str
    source_replicate: str
    raw_row_text: str
    raw_value: Union[SourceObservationRawValueKnown, SourceObservationRawValueUnknown]
    raw_unit: str
    source_reported_label: Union[SourceObservationSourceReportedLabelKnown, SourceObservationSourceReportedLabelUnknown]


class ObservationRevisionSubjectSourceSubjectIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionSubjectSourceSubjectIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionSubjectDonorIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionSubjectDonorIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionSubjectSpecimenIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionSubjectSpecimenIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionSubject(TypedDict):
    level: Literal['donor', 'specimen', 'cell_line', 'pooled_group', 'study_summary']
    source_subject_id: Union[ObservationRevisionSubjectSourceSubjectIdKnown, ObservationRevisionSubjectSourceSubjectIdUnknown]
    donor_id: Union[ObservationRevisionSubjectDonorIdKnown, ObservationRevisionSubjectDonorIdUnknown]
    specimen_id: Union[ObservationRevisionSubjectSpecimenIdKnown, ObservationRevisionSubjectSpecimenIdUnknown]
    resolution_status: Literal['resolved', 'unresolved', 'pooled', 'summary']
    reason: Union[str, None]


class ObservationRevisionDesignIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionDesignIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionUnresolvedSourceDesignKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionUnresolvedSourceDesignUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionContextIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionContextIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionSequenceConfigurationIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionSequenceConfigurationIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionNormalizedValueKnown(TypedDict):
    value: Union[Union[int, float], str]
    reason: None


class ObservationRevisionNormalizedValueUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionConversionKnownValue(TypedDict):
    method: str
    factor: Union[int, float]
    offset: Union[int, float]


class ObservationRevisionConversionKnown(TypedDict):
    value: ObservationRevisionConversionKnownValue
    reason: None


class ObservationRevisionConversionUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionIndependentUnitIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionIndependentUnitIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionRelationshipGroupKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionRelationshipGroupUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionPairingKeyKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionPairingKeyUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionReadCountKnown(TypedDict):
    value: int
    reason: None


class ObservationRevisionReadCountUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionDenominatorKnown(TypedDict):
    value: int
    reason: None


class ObservationRevisionDenominatorUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionSampleSizeKnown(TypedDict):
    value: int
    reason: None


class ObservationRevisionSampleSizeUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionUncertaintyKnownValueLevelKnown(TypedDict):
    value: Union[int, float]
    reason: None


class ObservationRevisionUncertaintyKnownValueLevelUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionUncertaintyKnownValue(TypedDict):
    kind: Literal['sd', 'sem', 'confidence_interval', 'range', 'other_reported']
    lower: Union[int, float]
    upper: Union[int, float]
    level: Union[ObservationRevisionUncertaintyKnownValueLevelKnown, ObservationRevisionUncertaintyKnownValueLevelUnknown]


class ObservationRevisionUncertaintyKnown(TypedDict):
    value: ObservationRevisionUncertaintyKnownValue
    reason: None


class ObservationRevisionUncertaintyUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionAssayIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionAssayIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionSourceReportedOutcomeKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionSourceReportedOutcomeUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionExecutedActionIdKnown(TypedDict):
    value: str
    reason: None


class ObservationRevisionExecutedActionIdUnknown(TypedDict):
    value: None
    reason: str


class ObservationRevisionFields(TypedDict):
    record_type: Literal['comparison_observation_revision']
    record_id: str
    observation_id: str
    source_record_id: str
    extraction_version: str
    mapping_version: str
    extraction_map_sha256: str
    subject: ObservationRevisionSubject
    design_id: Union[ObservationRevisionDesignIdKnown, ObservationRevisionDesignIdUnknown]
    unresolved_source_design: Union[ObservationRevisionUnresolvedSourceDesignKnown, ObservationRevisionUnresolvedSourceDesignUnknown]
    context_id: Union[ObservationRevisionContextIdKnown, ObservationRevisionContextIdUnknown]
    sequence_configuration_id: Union[ObservationRevisionSequenceConfigurationIdKnown, ObservationRevisionSequenceConfigurationIdUnknown]
    locus_ids: list[str]
    endpoint_id: str
    normalized_value: Union[ObservationRevisionNormalizedValueKnown, ObservationRevisionNormalizedValueUnknown]
    normalized_unit: str
    conversion: Union[ObservationRevisionConversionKnown, ObservationRevisionConversionUnknown]
    measurement_kind: Literal['individual', 'technical_replicate', 'biological_replicate', 'aggregate', 'censored']
    independent_unit_id: Union[ObservationRevisionIndependentUnitIdKnown, ObservationRevisionIndependentUnitIdUnknown]
    relationship_group: Union[ObservationRevisionRelationshipGroupKnown, ObservationRevisionRelationshipGroupUnknown]
    pairing_key: Union[ObservationRevisionPairingKeyKnown, ObservationRevisionPairingKeyUnknown]
    read_count: Union[ObservationRevisionReadCountKnown, ObservationRevisionReadCountUnknown]
    denominator: Union[ObservationRevisionDenominatorKnown, ObservationRevisionDenominatorUnknown]
    sample_size: Union[ObservationRevisionSampleSizeKnown, ObservationRevisionSampleSizeUnknown]
    uncertainty: Union[ObservationRevisionUncertaintyKnown, ObservationRevisionUncertaintyUnknown]
    assay_id: Union[ObservationRevisionAssayIdKnown, ObservationRevisionAssayIdUnknown]
    control_status: Literal['passed', 'failed', 'unknown', 'not_applicable']
    qc_status: Literal['passed', 'failed', 'unknown']
    detection_limit: Quantity
    quantification_limit: Quantity
    source_reported_outcome: Union[ObservationRevisionSourceReportedOutcomeKnown, ObservationRevisionSourceReportedOutcomeUnknown]
    applicability: Applicability
    executed_action_id: Union[ObservationRevisionExecutedActionIdKnown, ObservationRevisionExecutedActionIdUnknown]


class LegacyObservationRevision(ObservationRevisionFields):
    schema_version: Literal["1"]
    normalization_schema_version: Literal["1"]


class ObservationRevisionV2(ObservationRevisionFields):
    schema_version: Literal["2"]
    normalization_schema_version: Literal["2"]
    adjudicated_outcome: Literal["observed_positive", "tested_no_signal", "inconclusive", "untested"]
    adjudication_reason: str


ObservationRevision = Union[LegacyObservationRevision, ObservationRevisionV2]

class ObservationState(TypedDict):
    record_type: Literal['comparison_observation_state']
    schema_version: Literal['1']
    record_id: str
    observation_state_id: str
    source_record_id: str
    active_revision_id: str
    revision_ids: list[str]
    supersession_reason: str


class ObservationEquivalence(TypedDict):
    record_type: Literal['comparison_observation_equivalence']
    schema_version: Literal['1']
    record_id: str
    observation_equivalence_id: str
    source_record_ids: list[str]
    canonical_source_record_id: str
    basis: str
    source_refs: list[SourceRef]


class EvidenceApplicabilityDonorIdKnown(TypedDict):
    value: str
    reason: None


class EvidenceApplicabilityDonorIdUnknown(TypedDict):
    value: None
    reason: str


class EvidenceApplicabilitySequenceConfigurationIdKnown(TypedDict):
    value: str
    reason: None


class EvidenceApplicabilitySequenceConfigurationIdUnknown(TypedDict):
    value: None
    reason: str


class EvidenceApplicability(TypedDict):
    record_type: Literal['comparison_evidence_applicability']
    schema_version: Literal['1']
    record_id: str
    evidence_applicability_id: str
    observation_id: str
    design_id: str
    context_id: str
    endpoint_id: str
    donor_id: Union[EvidenceApplicabilityDonorIdKnown, EvidenceApplicabilityDonorIdUnknown]
    sequence_configuration_id: Union[EvidenceApplicabilitySequenceConfigurationIdKnown, EvidenceApplicabilitySequenceConfigurationIdUnknown]
    assessment: Applicability
    policy_version: str
    coverage_basis: Literal['synthetic_truth', 'callable_supported', 'retained_callset_model']


class PredictionPostEditScenario(TypedDict):
    sequence: str
    edits: list[Edit]
    scenario_id: str


class PredictionWindow(TypedDict):
    reference_bundle_id: str
    interval: Interval
    sequence: str
    path: list[PathPoint]
    construction_basis: Literal['complete_donor_reconstruction', 'partial_donor_reconstruction', 'reference_modeled']
    padding_left: int
    padding_right: int


class PredictionCellTypeKnown(TypedDict):
    value: str
    reason: None


class PredictionCellTypeUnknown(TypedDict):
    value: None
    reason: str


class PredictionTissueKnown(TypedDict):
    value: str
    reason: None


class PredictionTissueUnknown(TypedDict):
    value: None
    reason: str


class PredictionScoreKnown(TypedDict):
    value: Union[int, float]
    reason: None


class PredictionScoreUnknown(TypedDict):
    value: None
    reason: str


class PredictionUncertaintyKnown(TypedDict):
    value: Union[int, float]
    reason: None


class PredictionUncertaintyUnknown(TypedDict):
    value: None
    reason: str


class Prediction(TypedDict):
    record_type: Literal['comparison_prediction']
    schema_version: Literal['1']
    record_id: str
    prediction_id: str
    provider: str
    model: str
    model_version: str
    request_sha256: str
    response_sha256: str
    pre_edit_sequence_id: str
    post_edit_scenario: PredictionPostEditScenario
    window: PredictionWindow
    track: str
    cell_type: Union[PredictionCellTypeKnown, PredictionCellTypeUnknown]
    tissue: Union[PredictionTissueKnown, PredictionTissueUnknown]
    endpoint_id: str
    score: Union[PredictionScoreKnown, PredictionScoreUnknown]
    unit: str
    applicability: Assessment
    uncertainty: Union[PredictionUncertaintyKnown, PredictionUncertaintyUnknown]
    terms_uri: str
    terms_version: str
    redistribution: Literal['permitted', 'restricted', 'unknown']


class AssayCapabilityEndpointLocusPairsItem(TypedDict):
    endpoint_id: str
    locus_id: str


class AssayCapability(TypedDict):
    record_type: Literal['comparison_assay_capability']
    schema_version: Literal['1']
    record_id: str
    assay_capability_id: str
    name: str
    endpoint_ids: list[str]
    locus_ids: list[str]
    context_ids: list[str]
    resolution: Literal['copy_specific', 'locus_specific', 'combined_copies', 'not_applicable']
    allele_resolved: bool
    feasibility: Literal['verified', 'simulation_assumption', 'unverified', 'infeasible']
    required_control_roles: list[str]
    qc_requirements: list[str]
    max_targets: int
    source_refs: SourceRefs
    endpoint_locus_pairs: NotRequired[list[AssayCapabilityEndpointLocusPairsItem]]


class SpecimenPermittedUseKnown(TypedDict):
    value: str
    reason: None


class SpecimenPermittedUseUnknown(TypedDict):
    value: None
    reason: str


class Specimen(TypedDict):
    record_type: Literal['comparison_specimen']
    schema_version: Literal['1']
    record_id: str
    specimen_id: str
    donor_id: str
    source_specimen_id: str
    context_ids: list[str]
    availability: Literal['verified', 'simulation_assumption', 'unavailable', 'unresolved']
    permitted_use: Union[SpecimenPermittedUseKnown, SpecimenPermittedUseUnknown]
    capacity: int
    source_refs: SourceRefs
    use_authorized: NotRequired[bool]


class ResourceCapacitiesItemCapacityKnown(TypedDict):
    kind: Literal['finite']
    quantity: int


class ResourceCapacitiesItemCapacityUnknown(TypedDict):
    kind: Literal['unbounded']
    justification: str


class ResourceCapacitiesItem(TypedDict):
    dimension: str
    capacity: Union[ResourceCapacitiesItemCapacityKnown, ResourceCapacitiesItemCapacityUnknown]


class ResourcePoolIdKnown(TypedDict):
    value: str
    reason: None


class ResourcePoolIdUnknown(TypedDict):
    value: None
    reason: str


class Resource(TypedDict):
    record_type: Literal['comparison_resource']
    schema_version: Literal['1']
    record_id: str
    resource_id: str
    resource_type: Literal['study_setup', 'donor_activation', 'specimen_preparation', 'treatment', 'mock_control', 'assay_library', 'sequencing_batch', 'additional_target', 'additional_replicate']
    cost: int
    cost_unit: str
    cost_scale: int
    sharing_key: str
    capacities: list[ResourceCapacitiesItem]
    activation_rule: Literal['always', 'nonempty_plan', 'on_use']
    compatibility: Compatibility
    pool_id: Union[ResourcePoolIdKnown, ResourcePoolIdUnknown]
    instance_index: int
    max_pool_instances: int


class QuestionComparability(TypedDict):
    require_same_context: bool
    require_same_endpoint: bool
    require_same_timepoint: bool


class QuestionArmsItem(TypedDict):
    arm_id: str
    membership: Literal['all_resolved', 'carrier', 'noncarrier']
    required_independent_donors: int


class QuestionEndpointLocusPairsItem(TypedDict):
    endpoint_id: str
    locus_id: str


class Question(TypedDict):
    record_type: Literal['comparison_question']
    schema_version: Literal['1']
    record_id: str
    question_id: str
    template: Literal['candidate_validation', 'paired_design_comparison', 'variant_association', 'design_by_variant', 'evidence_replication']
    hypothesis: str
    decision: str
    design_ids: list[str]
    locus_ids: list[str]
    endpoint_ids: list[str]
    context_ids: list[str]
    genotype_predicate: GenotypePredicate
    comparability: QuestionComparability
    arms: list[QuestionArmsItem]
    pairing: Literal['within_donor', 'none']
    relatedness_policy: Literal['strict', 'declared_simulation_assumption']
    technical_replicates: int
    control_requirements: list[ControlRequirement]
    requirement_slot_ids: list[str]
    evidence_reuse: ReusePolicy
    matching_covariates: list[str]
    priority_weight: int
    hard_required: bool
    rationale: str
    policy_version: str
    motivating_observation_ids: list[str]
    motivating_configuration_ids: list[str]
    motivating_prediction_ids: list[str]
    planning_status: Literal['unplanned', 'partially_planned', 'fully_planned', 'already_satisfied', 'infeasible_proven', 'no_feasible_bundle_found', 'blocked']
    evidence_status: Literal['no_measurements', 'partial_measurements', 'measurements_available', 'evaluable', 'conflicting', 'insufficient_quality']
    remaining_gaps: list[str]
    endpoint_locus_pairs: NotRequired[list[QuestionEndpointLocusPairsItem]]
    construction_status: NotRequired[Literal['complete', 'excluded_by_declared_limit']]


class RequirementSlotDesignIdKnown(TypedDict):
    value: str
    reason: None


class RequirementSlotDesignIdUnknown(TypedDict):
    value: None
    reason: str


class RequirementSlotLocusIdKnown(TypedDict):
    value: str
    reason: None


class RequirementSlotLocusIdUnknown(TypedDict):
    value: None
    reason: str


class RequirementSlot(TypedDict):
    record_type: Literal['comparison_requirement_slot']
    schema_version: Literal['1']
    record_id: str
    slot_id: str
    question_id: str
    donor_id: str
    arm_id: str
    design_id: Union[RequirementSlotDesignIdKnown, RequirementSlotDesignIdUnknown]
    endpoint_id: str
    context_id: str
    locus_id: Union[RequirementSlotLocusIdKnown, RequirementSlotLocusIdUnknown]
    control_role: str
    replicate_position: int
    allowed_observation_ids: list[str]
    allowed_action_ids: list[str]


class ActionDesignIdKnown(TypedDict):
    value: str
    reason: None


class ActionDesignIdUnknown(TypedDict):
    value: None
    reason: str


class ActionTargetSequenceIdKnown(TypedDict):
    value: str
    reason: None


class ActionTargetSequenceIdUnknown(TypedDict):
    value: None
    reason: str


class ActionResourceDemandsItemUnknown(TypedDict):
    eligible_resource_ids: list[str]
    assignment_key: str
    dimension: str
    quantity: int


class ActionEndpointLocusPairsItem(TypedDict):
    endpoint_id: str
    locus_id: str


class ActionTargetSequencesItem(TypedDict):
    locus_id: str
    sequence_configuration_ids: list[str]


class Action(TypedDict):
    record_type: Literal['comparison_action']
    schema_version: Literal['1']
    record_id: str
    action_id: str
    operation_type: Literal['treatment', 'locus_assay', 'intended_assay', 'functional_assay', 'mock_control', 'assay_control', 'preparation', 'measurement_panel']
    donor_id: str
    specimen_id: str
    design_id: Union[ActionDesignIdKnown, ActionDesignIdUnknown]
    context_id: str
    assay_capability_id: str
    endpoint_ids: list[str]
    target_sequence_id: Union[ActionTargetSequenceIdKnown, ActionTargetSequenceIdUnknown]
    locus_ids: list[str]
    replicate_position: int
    batch_id: str
    experiment_id: str
    dependency_action_ids: list[str]
    control_requirements: list[ControlRequirement]
    resource_demands: list[Union[ResourceDemand, ActionResourceDemandsItemUnknown]]
    eligibility_status: Literal['eligible', 'ineligible', 'unresolved']
    reason_codes: list[str]
    memberships: list[Membership]
    endpoint_locus_pairs: NotRequired[list[ActionEndpointLocusPairsItem]]
    target_sequences: NotRequired[list[ActionTargetSequencesItem]]
    performed_control_roles: NotRequired[list[str]]
    performed_qc: NotRequired[list[str]]


class IdentityBridge(TypedDict):
    record_type: Literal['comparison_identity_bridge']
    schema_version: Literal['1']
    record_id: str
    bridge_id: str
    legacy_family: Literal['population', 'candidate_ledger', 'design']
    legacy_schema_version: Literal['1', '3', '7']
    legacy_record_id: str
    legacy_payload_json: str
    legacy_payload_sha256: str
    new_record_id: str
    relationship: Literal['semantic_bridge', 'incompatible_identity_semantics']
    mapping_version: str
    reason: str


class DesignSummaryCompleteLocusCountKnown(TypedDict):
    value: int
    reason: None


class DesignSummaryCompleteLocusCountUnknown(TypedDict):
    value: None
    reason: str


class DesignSummaryCompleteSequenceConfigurationCountKnown(TypedDict):
    value: int
    reason: None


class DesignSummaryCompleteSequenceConfigurationCountUnknown(TypedDict):
    value: None
    reason: str


class DesignSummaryCompleteAlignmentCountKnown(TypedDict):
    value: int
    reason: None


class DesignSummaryCompleteAlignmentCountUnknown(TypedDict):
    value: None
    reason: str


class DesignSummary(TypedDict):
    record_type: Literal['comparison_design_summary']
    schema_version: Literal['1']
    record_id: str
    design_summary_id: str
    design_id: str
    domain_id: str
    search_scope_ids: list[str]
    donor_count: int
    intended_evaluable_donors: int
    intended_unresolved_donors: int
    search_evaluable_donors: int
    search_unresolved_donors: int
    observed_locus_count: int
    observed_sequence_configuration_count: int
    observed_alignment_count: int
    complete_locus_count: Union[DesignSummaryCompleteLocusCountKnown, DesignSummaryCompleteLocusCountUnknown]
    complete_sequence_configuration_count: Union[DesignSummaryCompleteSequenceConfigurationCountKnown, DesignSummaryCompleteSequenceConfigurationCountUnknown]
    complete_alignment_count: Union[DesignSummaryCompleteAlignmentCountKnown, DesignSummaryCompleteAlignmentCountUnknown]
    supported_carrier_donors: int
    modeled_carrier_donors: int
    variant_effect_ids: list[str]
    evidence_observation_ids: list[str]
    reason_codes: list[str]


class DonorDesignSummaryCompleteConfigurationCountKnown(TypedDict):
    value: int
    reason: None


class DonorDesignSummaryCompleteConfigurationCountUnknown(TypedDict):
    value: None
    reason: str


class DonorDesignSummary(TypedDict):
    record_type: Literal['comparison_donor_design_summary']
    schema_version: Literal['1']
    record_id: str
    donor_design_summary_id: str
    donor_id: str
    design_id: str
    domain_id: str
    intended_assessment_ids: list[str]
    search_scope_ids: list[str]
    intended_evaluability: Literal['supported', 'modeled', 'unresolved']
    search_status: Literal['complete', 'partial', 'not_evaluated']
    nominated_configuration_ids: list[str]
    nominated_alignment_ids: list[str]
    observed_configuration_count: int
    observed_alignment_count: int
    complete_configuration_count: Union[DonorDesignSummaryCompleteConfigurationCountKnown, DonorDesignSummaryCompleteConfigurationCountUnknown]
    reason_codes: list[str]
    evidence_observation_ids: list[str]
    evidence_availability: Literal['not_assessed', 'none_in_current_sources', 'contextual_only', 'exact_donor_evidence']


class Decision(TypedDict):
    record_type: Literal['comparison_decision']
    schema_version: Literal['1']
    record_id: str
    decision_id: str
    design_ids: list[str]
    domain_id: str
    state: Literal['observed_difference_in_matched_context', 'sequence_based_hypothesis', 'tradeoff_requires_experiment', 'insufficient_comparable_evidence', 'no_difference_detected_within_scope']
    observation_ids: list[str]
    prediction_ids: list[str]
    reason_codes: list[str]
    question: str


class ExperimentPlanPreviousPlanIdKnown(TypedDict):
    value: str
    reason: None


class ExperimentPlanPreviousPlanIdUnknown(TypedDict):
    value: None
    reason: str


class ExperimentPlanLimits(TypedDict):
    beam_width: int
    max_expansions: int
    max_alternatives_per_question: int
    exact_action_cap: int
    exact_resource_cap: int
    exact_state_cap: int


class ExperimentPlanResourceAssignmentsItem(TypedDict):
    action_id: str
    assignment_key: str
    resource_id: str


class ExperimentPlanMandatorySeedCostKnown(TypedDict):
    value: int
    reason: None


class ExperimentPlanMandatorySeedCostUnknown(TypedDict):
    value: None
    reason: str


class ExperimentPlanQuestionStatesItem(TypedDict):
    question_id: str
    planning_status: Literal['unplanned', 'partially_planned', 'fully_planned', 'already_satisfied', 'infeasible_proven', 'no_feasible_bundle_found', 'blocked']
    evidence_status: Literal['no_measurements', 'partial_measurements', 'measurements_available', 'evaluable', 'conflicting', 'insufficient_quality']


class ExperimentPlan(TypedDict):
    record_type: Literal['comparison_experiment_plan']
    schema_version: Literal['1']
    record_id: str
    plan_id: str
    planner_version: str
    source_sha256: str
    implementation_sha256: str
    previous_plan_id: Union[ExperimentPlanPreviousPlanIdKnown, ExperimentPlanPreviousPlanIdUnknown]
    selector: Literal['complete_bundle_greedy', 'donor_first', 'adapted_coverage', 'canonical_bundle', 'random_bundle', 'exact_tiny']
    random_seed: int
    limits: ExperimentPlanLimits
    budget: int
    fill_policy: Literal['none']
    status: Literal['optimal', 'heuristic_complete', 'feasible_search_incomplete', 'infeasible_proven', 'no_feasible_bundle_found', 'invalid_plan']
    optimality_established: bool
    meets_request: bool
    selected_action_ids: list[str]
    resource_assignments: list[ExperimentPlanResourceAssignmentsItem]
    total_cost: int
    cost_unit: str
    cost_scale: int
    active_resource_ids: list[str]
    mandatory_seed_cost: Union[ExperimentPlanMandatorySeedCostKnown, ExperimentPlanMandatorySeedCostUnknown]
    weighted_completion: int
    unweighted_completion: int
    existing_evidence_weight: int
    newly_planned_weight: int
    question_states: list[ExperimentPlanQuestionStatesItem]
    audit_sha256: str


ComparisonRecord = Union[
    SourceArtifact,
    ReferenceBundle,
    EditorDefinition,
    Locus,
    IntendedLocus,
    Design,
    ExperimentalContext,
    EndpointDefinition,
    DesignAlias,
    Domain,
    SearchProfile,
    SequenceConfiguration,
    SequenceOccurrence,
    DesignConfiguration,
    Alignment,
    Donor,
    DonorConfiguration,
    SearchScope,
    VariantEffect,
    IntendedAssessment,
    SourceObservation,
    ObservationRevision,
    ObservationState,
    ObservationEquivalence,
    EvidenceApplicability,
    Prediction,
    AssayCapability,
    Specimen,
    Resource,
    Question,
    RequirementSlot,
    Action,
    IdentityBridge,
    DesignSummary,
    DonorDesignSummary,
    Decision,
    ExperimentPlan,
]
