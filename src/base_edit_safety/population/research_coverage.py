"""Count conditional sequence coverage on fixed reference and copy denominators."""
from bisect import bisect_right
from collections import Counter, defaultdict
import re

from .haplotypes import local_status
from .search import evaluation_span


def merge(intervals):
    result = []
    for lo, hi in sorted(intervals):
        if hi <= lo:
            continue
        if result and lo <= result[-1][1]:
            result[-1] = (result[-1][0], max(hi, result[-1][1]))
        else:
            result.append((lo, hi))
    return result


def subtract(intervals, masks):
    masks = merge(masks)
    result, index = [], 0
    for lo, hi in merge(intervals):
        while index < len(masks) and masks[index][1] <= lo:
            index += 1
        cursor, j = lo, index
        while j < len(masks) and masks[j][0] < hi:
            a, b = masks[j]
            if cursor < a:
                result.append((cursor, min(a, hi)))
            cursor = max(cursor, b)
            j += 1
        if cursor < hi:
            result.append((cursor, hi))
    return result


def size(intervals):
    return sum(b-a for a, b in intervals)


def intersect(a, b):
    a, b = merge(a), merge(b)
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return out


def _sequence_windows(sequence):
    return sum(max(0, len(m[0])-22) for m in re.finditer('[ACGT]+', sequence))


def _status_ranges(material, starts):
    """Local status changes only where a 23-base window meets a site boundary."""
    sites = material['excluded_intervals'] + material['heterozygous_sites']
    points = sorted({p for h in sites for p in (h['start']-22, h['end'])})
    valid, reasons = [], defaultdict(list)
    for lo, hi in starts:
        breaks = [lo] + [p for p in points if lo < p < hi] + [hi]
        for a, b in zip(breaks, breaks[1:]):
            ok, reason = local_status(material, a, a+23)
            (valid if ok else reasons[reason]).append((a, b))
    return merge(valid), reasons


def material_windows(material):
    """Return exact reference paths and counts of all reconstructed copy windows."""
    path, sequence = material['sequence_path'], material['sequence_genomic_forward']
    runs, lo, last = [], None, None
    for point, base in zip(path, sequence):
        pos = point.get('pos0') if point['kind'] == 'reference' and base in 'ACGT' else None
        if pos is None or last is not None and pos != last+1:
            if lo is not None and last-lo+1 >= 23:
                runs.append((lo, last-21))
            lo = None
        if pos is not None and lo is None:
            lo = pos
        last = pos
    if lo is not None and last-lo+1 >= 23:
        runs.append((lo, last-21))
    exact, reasons = _status_ranges(material, runs)
    total = max(0, len(sequence)-22)
    ok, _ = local_status(material, material['interval']['start'], material['interval']['end'])
    if ok and not material.get("phase_constraints"):
        evaluated = _sequence_windows(sequence)
    else:
        # This branch is limited to neighborhoods with gaps or incompatible phase.
        evaluated = 0
        for offset in range(total):
            if 'N' in sequence[offset:offset+23]:
                continue
            lo, hi = evaluation_span(path[offset:offset+23], material)
            evaluated += int(local_status(material, lo, hi)[0])
    return exact, reasons, total, evaluated


class Coverage:
    """Accumulate masks without scanning billions of reference windows in Python."""
    def __init__(self, reference_rows, sample, confidence=None):
        self.reference = {r['contig']: r for r in reference_rows}
        self.sample = sample
        self.confidence = confidence
        self.base_masks = defaultdict(lambda: defaultdict(list))
        self.window_masks = defaultdict(lambda: defaultdict(list))
        self.sample_counts = defaultdict(Counter)
        self.valid_starts = {}
        self.valid_ends = {}
        for c, row in self.reference.items():
            runs = row['acgt_intervals']
            starts = [(a, b-22) for a, b in runs if b-a >= 23]
            self.valid_starts[c] = starts
            self.valid_ends[c] = [b for a, b in starts]
            unknown = subtract([(0, row['bases'])], runs)
            self.base_masks[c]['reference_sequence_unknown'].extend(unknown)
            self.window_masks[c]['reference_sequence_unknown'].extend(
                subtract([(0, max(0, row['bases']-22))], starts))

    def _reference_windows(self, contig, lo, hi):
        runs = self.valid_starts[contig]
        i = bisect_right(self.valid_ends[contig], lo)
        total = 0
        while i < len(runs) and runs[i][0] < hi:
            total += max(0, min(hi, runs[i][1])-max(lo, runs[i][0]))
            i += 1
        return total

    def cluster(self, cluster, materials):
        c, lo, hi = (cluster['interval'][k] for k in ('contig', 'start', 'end'))
        if not materials:
            if cluster.get('resource_excluded'):
                reason='local_reconstruction_resource_limit'
                self.base_masks[c][reason].append((lo,hi))
                self.window_masks[c][reason].append((lo-22,hi))
            return
        exact_all = None
        for material in materials:
            exact, reasons, total, evaluated = material_windows(material)
            exact_all = exact if exact_all is None else intersect(exact_all, exact)
            counters = self.sample_counts[c]
            counters['model_window_delta'] += total-max(0, hi-lo-22)
            counters['evaluated_window_delta'] += evaluated-self._reference_windows(c, lo, max(lo, hi-22))
            counters['indel_path_windows_evaluated'] += evaluated-size(exact)
            for gap in material['excluded_intervals']:
                self.base_masks[c][gap['reason']].append((gap['start'], gap['end']))
                self.window_masks[c][gap['reason']].append((gap['start']-22, gap['end']))
            for reason, intervals in reasons.items():
                self.window_masks[c][reason].extend(intervals)
        # These are all reference-grid windows lacking a supported exact path on every copy.
        self.window_masks[c]['sample_path_or_local_status_unavailable'].extend(
            subtract([(lo, hi-22)], exact_all or []))

    def exclude_contig(self, contig, reason):
        length = self.reference[contig]['bases']
        self.base_masks[contig][reason].append((0, length))
        self.window_masks[contig][reason].append((0, max(0, length-22)))

    def finish(self, sample_by_contig, failed=False):
        rows, regions, windows = [], [], []
        for c, ref in self.reference.items():
            length, stop = ref['bases'], max(0, ref['bases']-22)
            ploidy = []
            for r in sorted((r for r in (self.sample or {}).get('ploidy_intervals', []) if r['contig'] == c), key=lambda r:r['start']):
                if ploidy and r['start'] == ploidy[-1]['end'] and r['copies'] == ploidy[-1]['copies']:
                    ploidy[-1]['end'] = r['end']
                else:
                    ploidy.append(dict(r))
            positive = [(r['start'], r['end']) for r in ploidy if r['copies'] in (1,2)]
            permitted = [(a, b-22) for a,b in positive if b-a >= 23]
            unknown = subtract([(0,length)], [(r['start'],r['end']) for r in ploidy])
            zero = [(r['start'],r['end']) for r in ploidy if r['copies']==0]
            for reason, masks in [('sample_ploidy_undeclared',unknown), ('zero_declared_copies',zero)]:
                self.base_masks[c][reason].extend(masks)
                self.window_masks[c][reason].extend((max(0,a-22),min(stop,b)) for a,b in masks)
            self.window_masks[c]['ploidy_boundary'].extend(subtract([(a,b-22) for a,b in merge(positive) if b-a>=23], permitted))
            for row in (self.sample or {}).get('excluded_contigs',[]):
                if row['contig']==c:
                    self.exclude_contig(c, 'sample_contig_excluded: '+row['reason'])
            if failed:
                self.exclude_contig(c,'sample_run_failed')
            base_masks = {k:intersect(v,[(0,length)]) for k,v in self.base_masks[c].items() if v}
            window_masks = {k:intersect(v,[(0,stop)]) for k,v in self.window_masks[c].items() if v}
            modeled = subtract([(0,length)], [r for v in base_masks.values() for r in v])
            comparable = subtract([(0,stop)], [r for v in window_masks.values() for r in v])
            base_excluded, window_excluded = length-size(modeled), stop-size(comparable)
            nominal = sum(max(0,r['end']-r['start']-22)*r['copies'] for r in ploidy if r['copies'])
            nominal_valid = sum(self._reference_windows(c,r['start'],max(r['start'],r['end']-22))*r['copies'] for r in ploidy if r['copies'])
            full_mask = any(v == [(0,length)] and k != 'reference_sequence_unknown' for k,v in base_masks.items())
            counters = self.sample_counts[c]
            resource_windows = window_masks.get('local_reconstruction_resource_limit', [])
            resource_excluded = sum(r['copies']*sum(self._reference_windows(c,a,b)
                for a,b in intersect(resource_windows, [(r['start'],max(r['start'],r['end']-22))]))
                for r in ploidy if r['copies'])
            confidence = None if self.confidence is None else size(intersect(self.confidence.get(c,[]),[(0,length)]))
            variant_counts = {k:0 for k in ('source_vcf_rows','variant_rows','variant_records',
                'applied_variant_records','reference_variant_records',
                'resolved_spanning_deletion_variant_records','unsupported_variant_records')}
            variant_counts.update(sample_by_contig.get(c,{}))
            row = dict(contig=c,reference_bases=length,reference_acgt_bases=length-ref['unknown_bases'],
                modeled_reference_bases=size(modeled),excluded_reference_bases=base_excluded,
                reference_oriented_windows_per_guide=2*stop,
                reference_evaluable_oriented_windows_per_guide=ref['evaluable_oriented_windows_per_guide'],
                comparable_oriented_windows_per_guide=2*size(comparable),
                excluded_comparison_oriented_windows_per_guide=2*window_excluded,
                nominal_reference_copy_oriented_windows_per_guide=2*nominal,
                conditional_model_oriented_windows_per_guide=2*(nominal+counters['model_window_delta']),
                sample_evaluated_copy_oriented_windows_per_guide=0 if full_mask else 2*(nominal_valid+counters['evaluated_window_delta']-resource_excluded),
                evaluated_indel_copy_oriented_windows_per_guide=0 if full_mask else 2*counters['indel_path_windows_evaluated'],
                benchmark_region_bases=confidence,confidently_genotyped_bases=None,
                exclusion_reasons={k:2*size(v) for k,v in window_masks.items() if v},
                **variant_counts)
            if full_mask:
                row['scope_excluded_variant_records'] = row.get('variant_rows',0)
                for key in ('applied_variant_records','reference_variant_records','resolved_spanning_deletion_variant_records'):
                    if key in row:
                        row[key] = 0
            rows.append(row)
            for kind, domain, masks, output in [('base',length,base_masks,regions),('window_start',stop,window_masks,windows)]:
                # Full partition; overlaps keep every reason, never duplicate denominators.
                events=defaultdict(list)
                events[0];events[domain]
                for reason, ranges in masks.items():
                    for a,b in merge(ranges):
                        events[a].append((reason,1));events[b].append((reason,-1))
                active=Counter();previous=0
                for pos in sorted(events):
                    if previous < pos:
                        reasons=sorted(k for k,v in active.items() if v)
                        output.append({'contig':c,'start':previous,'end':pos,
                            'status':'excluded' if reasons else 'evaluated','reasons':reasons})
                    for reason,delta in events[pos]:active[reason]+=delta
                    previous=pos
            self.window_masks[c] = window_masks
        return rows, regions, windows


def load_confidence(path, reference):
    """Read benchmark intervals for a separate coverage comparison."""
    from .common import PopulationError
    result=defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#','track','browser')):
                continue
            fields=line.split()
            if len(fields)<3:
                raise PopulationError('Confidence BED requires contig, start and end')
            c=fields[0]
            if c not in reference.aliases:
                raise PopulationError('Confidence BED contig is outside the reference')
            c=reference.aliases[c]
            a,b=map(int,fields[1:3])
            if not any(lo<=a<b<=hi for lo,hi in reference.bounds(c)):
                raise PopulationError('Confidence BED exceeds reference bounds')
            result[c].append((a,b))
    return {c:merge(r) for c,r in result.items()}
