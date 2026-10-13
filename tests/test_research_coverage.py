"""Independent arithmetic checks of full denominators and indel window counts."""
from collections import defaultdict, Counter

import pytest

from base_edit_safety.population.research_coverage import Coverage, merge, subtract
from base_edit_safety.population.research_genotypes import iter_variant_clusters, materialize_cluster
from test_research_genotypes import POLICY, Reference


def counted(tmp_path, rows, ploidy=None):
    reference=Reference()
    reference.sequence='ACGT'*25
    vcf=tmp_path/'calls.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##reference=artificial\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'+''.join('artificial_1\t'+row+'\n' for row in rows))
    ploidy=ploidy or [{'contig':'artificial_1','start':0,'end':100,'copies':2}]
    counter=Coverage([{'contig':'artificial_1','bases':100,'unknown_bases':0,'acgt_intervals':[(0,100)],'evaluable_oriented_windows_per_guide':156}], {'ploidy_intervals':ploidy})
    counts=defaultdict(Counter)
    for cluster in iter_variant_clusters(vcf,reference,'sample',POLICY,ploidy):
        counter.cluster(cluster,materialize_cluster(reference,cluster,'sample','invented'))
        for exclusion in cluster.get('unbounded_exclusions',[]):
            counter.exclude_contig('artificial_1',exclusion['reason'])
    rows,regions,windows=counter.finish(counts)
    assert sum(r['end']-r['start'] for r in regions)==100
    assert sum(r['end']-r['start'] for r in windows)==78
    assert rows[0]['comparable_oriented_windows_per_guide']+rows[0]['excluded_comparison_oriented_windows_per_guide']==156
    return rows[0]


def test_missing_base_keeps_full_denominator(tmp_path):
    row=counted(tmp_path,['41\tmissing\tA\tC\t.\tPASS\t.\tGT\t./.'])
    assert row['modeled_reference_bases']==99
    assert row['comparable_oriented_windows_per_guide']==2*(78-23)
    assert row['sample_evaluated_copy_oriented_windows_per_guide']==2*2*(78-23)
    assert row['nominal_reference_copy_oriented_windows_per_guide']==312
    assert row['confidently_genotyped_bases'] is None


def test_phase_excludes_only_windows_touching_both_sites(tmp_path):
    row=counted(tmp_path,['41\ta\tA\tC\t.\tPASS\t.\tGT\t0/1','46\tb\tC\tG\t.\tPASS\t.\tGT\t0/1'])
    assert row['modeled_reference_bases']==100
    # Starts 23 through 40 inclusive span both positions 40 and 45.
    assert row['comparable_oriented_windows_per_guide']==2*(78-18)
    assert row['sample_evaluated_copy_oriented_windows_per_guide']==4*(78-18)


@pytest.mark.parametrize('record,physical,paired,indel',[
    ('41\tdel\tAC\tA\t.\tPASS\t.\tGT\t1/1',77,55,22),
    ('41\tins\tA\tAGG\t.\tPASS\t.\tGT\t1/1',80,56,24),
])
def test_indels_keep_reference_grid_and_actual_copy_windows_separate(tmp_path,record,physical,paired,indel):
    row=counted(tmp_path,[record])
    assert row['reference_oriented_windows_per_guide']==156
    assert row['nominal_reference_copy_oriented_windows_per_guide']==312
    assert row['conditional_model_oriented_windows_per_guide']==4*physical
    assert row['sample_evaluated_copy_oriented_windows_per_guide']==4*physical
    assert row['comparable_oriented_windows_per_guide']==2*paired
    assert row['evaluated_indel_copy_oriented_windows_per_guide']==4*indel
    assert row['modeled_reference_bases']==100


def test_ploidy_transition_excludes_crossing_windows_without_variants(tmp_path):
    p=[{'contig':'artificial_1','start':0,'end':50,'copies':1},{'contig':'artificial_1','start':50,'end':100,'copies':2}]
    row=counted(tmp_path,[],p)
    assert row['comparable_oriented_windows_per_guide']==2*(28+28)
    assert row['sample_evaluated_copy_oriented_windows_per_guide']==2*(28+2*28)
    assert row['reference_oriented_windows_per_guide']==156


def test_interval_subtraction_preserves_overlapping_denominator():
    assert merge([(2,6),(3,8),(8,9)])==[(2,9)]
    assert subtract([(0,10)],[(2,6),(3,8)])==[(0,2),(8,10)]


def test_haploid_missing_call_and_no_right_star_anchor_do_not_shrink_scope(tmp_path):
    p=[{'contig':'artificial_1','start':0,'end':100,'copies':1}]
    row=counted(tmp_path,['41\tmissing\tA\tC\t.\tPASS\t.\tGT\t.'],p)
    assert row['nominal_reference_copy_oriented_windows_per_guide']==156
    assert row['sample_evaluated_copy_oriented_windows_per_guide']==110
    row=counted(tmp_path,['41\tstar\tA\t*\t.\tPASS\t.\tGT\t1'],p)
    assert row['reference_bases']==100
    assert row['reference_oriented_windows_per_guide']==156
    assert row['modeled_reference_bases']==0
    assert row['comparable_oriented_windows_per_guide']==0
    assert row['sample_evaluated_copy_oriented_windows_per_guide']==0


def test_symbolic_local_gap_preserves_flanking_regions(tmp_path):
    row=counted(tmp_path,['41\tdel\tA\t<DEL>\t.\tPASS\tEND=50\tGT\t1/1'])
    assert row['modeled_reference_bases']==90
    # [40,50) masks starts 18..49 inclusive.
    assert row['comparable_oriented_windows_per_guide']==2*(78-32)
    assert row['sample_evaluated_copy_oriented_windows_per_guide']==4*(78-32)


def test_local_resource_gap_removes_boundary_windows_from_sample_counts():
    sample={'ploidy_intervals':[{'contig':'artificial_1','start':0,'end':100,'copies':2}]}
    counter=Coverage([{'contig':'artificial_1','bases':100,'unknown_bases':0,'acgt_intervals':[(0,100)],'evaluable_oriented_windows_per_guide':156}],sample)
    counter.cluster({'interval':{'contig':'artificial_1','start':30,'end':70},'copies':2,'resource_excluded':True},[])
    rows,_,_=counter.finish({})
    assert rows[0]['modeled_reference_bases']==60
    assert rows[0]['comparable_oriented_windows_per_guide']==2*(78-62)
    assert rows[0]['sample_evaluated_copy_oriented_windows_per_guide']==4*(78-62)


def test_neighboring_resource_gaps_count_shared_context_once():
    sample={'ploidy_intervals':[{'contig':'artificial_1','start':0,'end':300,'copies':2}]}
    counter=Coverage([{'contig':'artificial_1','bases':300,'unknown_bases':0,'acgt_intervals':[(0,300)],'evaluable_oriented_windows_per_guide':556}],sample)
    for start,end in [(30,100),(110,180)]:
        counter.cluster({'interval':{'contig':'artificial_1','start':start,'end':end},'copies':2,'resource_excluded':True},[])
    rows,_,_=counter.finish({})
    # Boundary masks [8,100) and [88,180) have a 172-start union.
    assert rows[0]['comparable_oriented_windows_per_guide']==2*(278-172)
    assert rows[0]['sample_evaluated_copy_oriented_windows_per_guide']==4*(278-172)
