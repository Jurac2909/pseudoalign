#!/usr/bin/env python3
"""
Per-read disagreement report.

One row per read where our mapping differs from kallisto's.
Columns: relation, read_name, read_seq, ours_transcripts, ours_seqs, kall_transcripts, kall_seqs

Usage:
  python3 error_report.py <ours.sam> <kall.sam> <reads.fastq> <transcripts.fa> [out.tsv]
"""
import sys
from collections import defaultdict


def parse_sam(path):
    reads = defaultdict(set)
    with open(path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.split('\t')
            if int(parts[1]) & 2048:
                continue
            ref = parts[2]
            if ref != '*':
                reads[parts[0]].add(ref)
    return {k: frozenset(v) for k, v in reads.items()}


def parse_fastq_selected(path, wanted):
    seqs = {}
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            name = h.strip().lstrip('@').split()[0]
            seq  = f.readline().strip()
            f.readline()
            f.readline()
            if name in wanted:
                seqs[name] = seq
    return seqs


def parse_fasta_selected(path, wanted):
    seqs = {}
    cur_name = None
    cur_seq  = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_name in wanted:
                    seqs[cur_name] = ''.join(cur_seq)
                cur_name = line[1:].split()[0]
                cur_seq  = []
            elif cur_name in wanted:
                cur_seq.append(line)
        if cur_name in wanted:
            seqs[cur_name] = ''.join(cur_seq)
    return seqs


def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} ours.sam kall.sam reads.fastq transcripts.fa [out.tsv]")
        sys.exit(1)

    ours_sam = sys.argv[1]
    kall_sam = sys.argv[2]
    fastq    = sys.argv[3]
    trans_fa = sys.argv[4]
    out_path = sys.argv[5] if len(sys.argv) > 5 else "error_report.tsv"

    print("Parsing SAMs...")
    ours = parse_sam(ours_sam)
    kall = parse_sam(kall_sam)

    print("Finding disagreements...")
    errors = []
    wanted_reads = set()
    wanted_trans = set()

    for rname in sorted(set(ours) | set(kall)):
        oe = ours.get(rname, frozenset())
        ke = kall.get(rname, frozenset())
        if oe == ke:
            continue

        if oe and ke and not oe <= ke and not ke <= oe:
            rel = 'different'
        elif oe < ke:
            rel = 'ours⊂kall'
        elif ke < oe:
            rel = 'kall⊂ours'
        elif not oe:
            rel = 'only_kall'
        else:
            rel = 'only_ours'

        errors.append((rname, rel, tuple(sorted(oe)), tuple(sorted(ke))))
        wanted_reads.add(rname)
        wanted_trans.update(oe)
        wanted_trans.update(ke)

    print(f"  {len(errors):,} disagreeing reads, {len(wanted_trans)} unique transcripts")

    print(f"Loading {len(wanted_reads):,} read sequences...")
    read_seqs = parse_fastq_selected(fastq, wanted_reads)

    print(f"Loading {len(wanted_trans)} transcript sequences...")
    trans_seqs = parse_fasta_selected(trans_fa, wanted_trans)

    print(f"Writing {len(errors):,} rows to {out_path}...")
    with open(out_path, 'w') as f:
        f.write('\t'.join(['relation', 'read_name', 'read_seq',
                           'ours_transcripts', 'ours_seqs',
                           'kall_transcripts', 'kall_seqs']) + '\n')
        for rname, rel, oe, ke in errors:
            def fmt(names):
                t_names = ','.join(names)
                t_seqs  = ';'.join(f"{t}={trans_seqs.get(t, 'N/A')}" for t in names)
                return t_names, t_seqs

            on, os_ = fmt(oe)
            kn, ks  = fmt(ke)
            f.write('\t'.join([rel, rname, read_seqs.get(rname, 'N/A'),
                                on, os_, kn, ks]) + '\n')

    print(f"Done — {out_path}")


if __name__ == '__main__':
    main()
