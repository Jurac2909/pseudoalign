#!/usr/bin/env python3
"""
Add k-mer match counts to error_report.tsv rows.

For each error pattern, counts how many 31-mers from the representative read
appear in each transcript in our EC vs kallisto's EC.  The result shows which
transcripts genuinely share sequence with the read and which are spurious.

Usage:
  python3 compare_errors.py <errors.tsv> <transcripts.fa> [out.tsv] [k]
"""
import sys

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

def kmer_set(seq, k):
    s = set()
    for i in range(len(seq) - k + 1):
        s.add(seq[i:i+k])
    return s

def rc(seq):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]

def canonical_kmers(seq, k):
    s = set()
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        s.add(min(km, rc(km)))
    return s

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} errors.tsv transcripts.fa [out.tsv] [k]")
        sys.exit(1)

    errors_path = sys.argv[1]
    trans_fa    = sys.argv[2]
    out_path    = sys.argv[3] if len(sys.argv) > 3 else "errors_compared.tsv"
    k           = int(sys.argv[4]) if len(sys.argv) > 4 else 31

    print("Reading error rows...")
    rows = []
    wanted_trans = set()
    with open(errors_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            row = dict(zip(header, parts))
            ours_names = [t for t in row['ours_ec'].split(',') if t]
            kall_names = [t for t in row['kall_ec'].split(',') if t]
            rows.append((row, ours_names, kall_names))
            wanted_trans.update(ours_names)
            wanted_trans.update(kall_names)

    print(f"Loading {len(wanted_trans)} full transcript sequences...")
    trans_seqs = parse_fasta_selected(trans_fa, wanted_trans)

    print(f"Computing k-mer overlaps (k={k}) for {len(rows)} patterns...")
    out_header = ['count', 'relation', 'rep_read', 'read_seq',
                  'ours_ec', 'ours_kmer_hits',
                  'kall_ec', 'kall_kmer_hits']

    with open(out_path, 'w') as f:
        f.write('\t'.join(out_header) + '\n')
        for row, ours_names, kall_names in rows:
            read_seq = row['read_seq']
            if len(read_seq) < k:
                read_kmers = set()
            else:
                read_kmers = canonical_kmers(read_seq, k)

            def hits(names):
                parts = []
                for t in names:
                    tseq = trans_seqs.get(t, '')
                    if tseq and read_kmers:
                        t_kmers = canonical_kmers(tseq, k)
                        n = len(read_kmers & t_kmers)
                        pct = 100.0 * n / len(read_kmers)
                        parts.append(f"{t}:{n}/{len(read_kmers)}({pct:.0f}%)")
                    else:
                        parts.append(f"{t}:N/A")
                return ','.join(parts)

            f.write('\t'.join([
                row['count'],
                row['relation'],
                row['rep_read'],
                read_seq,
                row['ours_ec'],   hits(ours_names),
                row['kall_ec'],   hits(kall_names),
            ]) + '\n')

    print(f"Done — {out_path}")

if __name__ == '__main__':
    main()
