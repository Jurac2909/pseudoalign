#!/usr/bin/env python3
"""
Compare per-read transcript assignments between our SAM and kallisto's SAM.
Usage: python3 compare_sam.py <ours.sam> <kallisto.sam>
"""

import sys
from collections import defaultdict

def parse_sam(path):
    """Return dict: read_name -> frozenset of transcript names"""
    assignments = defaultdict(set)
    with open(path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            cols = line.rstrip('\n').split('\t')
            qname = cols[0]
            flag  = int(cols[1])
            rname = cols[2]
            if rname == '*':
                # explicitly unmapped
                assignments[qname]  # ensure key exists
                continue
            # skip supplementary (2048) but keep primary (0) and secondary (256)
            if flag & 2048:
                continue
            assignments[qname].add(rname)
    return {k: frozenset(v) for k, v in assignments.items()}

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} ours.sam kallisto.sam")
        sys.exit(1)

    ours = parse_sam(sys.argv[1])
    kall = parse_sam(sys.argv[2])

    all_reads = sorted(set(ours) | set(kall))

    exact = 0
    ours_subset = 0   # our EC ⊆ kallisto EC (more specific)
    kall_subset = 0   # kallisto EC ⊆ our EC (less specific)
    different   = 0
    only_ours   = 0   # we mapped, kallisto didn't
    only_kall   = 0   # kallisto mapped, we didn't
    both_unmapped = 0

    print(f"{'read':<25} {'ours_ec':>6} {'kall_ec':>6}  {'relation':<12}  notes")
    print('-' * 75)

    for r in all_reads:
        o = ours.get(r, frozenset())
        k = kall.get(r, frozenset())

        if not o and not k:
            both_unmapped += 1
            rel = 'both_unmapped'
        elif not o:
            only_kall += 1
            rel = 'only_kallisto'
        elif not k:
            only_ours += 1
            rel = 'only_ours'
        elif o == k:
            exact += 1
            rel = 'exact'
        elif o < k:
            ours_subset += 1
            rel = 'ours⊂kall'
        elif k < o:
            kall_subset += 1
            rel = 'kall⊂ours'
        else:
            different += 1
            rel = 'different'

        if rel != 'exact':
            notes = ''
            if o - k:
                notes += f'only_ours: {sorted(o-k)[:3]}{"…" if len(o-k)>3 else ""}  '
            if k - o:
                notes += f'only_kall: {sorted(k-o)[:3]}{"…" if len(k-o)>3 else ""}'
            print(f"{r:<25} {len(o):>6} {len(k):>6}  {rel:<12}  {notes}")

    total = len(all_reads)
    print()
    print(f"Total reads       : {total}")
    print(f"Exact match       : {exact}  ({100*exact/total:.1f}%)")
    print(f"ours ⊂ kallisto   : {ours_subset}")
    print(f"kallisto ⊂ ours   : {kall_subset}")
    print(f"Disjoint/diff     : {different}")
    print(f"Only ours mapped  : {only_ours}")
    print(f"Only kallisto     : {only_kall}")
    print(f"Both unmapped     : {both_unmapped}")

if __name__ == '__main__':
    main()
