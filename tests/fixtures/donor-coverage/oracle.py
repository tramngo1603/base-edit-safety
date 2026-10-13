#!/usr/bin/env python3
"""Enumerate 23-base NGG windows over the hand-declared fixture haplotypes.

Both strands, ungapped spacer <=6 mismatches vs the declared guide.
No VCF parsing, no genotype logic: haplotype strings, exclusion spans
(0-based half-open, forward coordinates) and excluded forward window
starts are fixture truth; start exclusions apply before orientation.
"""
import json
import sys

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def main():
    fx = json.load(open("fixtures.json"))
    guide = fx["guide"]
    assert len(guide) == 20
    tot_c = tot_h = 0
    ok = True
    for case in fx["cases"]:
        cand, hits = 0, []
        for hap in case["haplotypes"]:
            seq = hap["seq"]
            n = len(seq)
            exc = [tuple(x) for x in hap.get("exclude", [])]
            xst = [tuple(x) for x in hap.get("excluded_window_starts", [])]
            for strand, s in (("+", seq), ("-", revcomp(seq))):
                ex = exc if strand == "+" else [(n - e, n - b) for b, e in exc]
                for i in range(n - 22):
                    if any(i < e and i + 23 > b for b, e in ex):
                        continue
                    fstart = i if strand == "+" else n - 23 - i
                    if any(a <= fstart < b for a, b in xst):
                        continue
                    w = s[i:i + 23]
                    if w[21] != "G" or w[22] != "G":
                        continue
                    cand += 1
                    mm = sum(a != b for a, b in zip(w[:20], guide))
                    if mm <= 6:
                        hits.append((hap["name"], strand, fstart, mm, w))
        ec, eh = case["expected_candidates"], case["expected_hits"]
        good = cand == ec and len(hits) == eh
        ok = ok and good
        status = "OK" if good else "MISMATCH"
        print(f"{case['id']}: candidates={cand}/{ec} hits={len(hits)}/{eh} {status}")
        for name, strand, fstart, mm, w in hits:
            print(f"  hit {name} {strand} fwd_start={fstart} mm={mm} {w}")
        tot_c += cand
        tot_h += len(hits)
    print(f"TOTAL candidates={tot_c} hits={tot_h}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
