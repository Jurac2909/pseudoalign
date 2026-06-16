#!/usr/bin/env python3
"""
Compare our full-transcriptome results against kallisto.
Usage: python3 compare_full.py <ours_results.tsv> <kallisto_abundance.tsv>

Outputs:
  - summary to stdout
  - comparison_full.tsv with per-transcript breakdown
"""
import sys
import math
from collections import defaultdict

def load_ours(path):
    data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            tid, counts, tpm = parts[0], float(parts[1]), float(parts[2])
            data[tid] = (counts, tpm)
    return data

def load_kallisto(path):
    data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            tid, length, eff_length, counts, tpm = parts[0], parts[1], parts[2], float(parts[3]), float(parts[4])
            data[tid] = (counts, tpm)
    return data

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} ours_results.tsv kallisto_abundance.tsv")
        sys.exit(1)

    ours = load_ours(sys.argv[1])
    kall = load_kallisto(sys.argv[2])

    all_tids = sorted(set(ours) | set(kall))

    both = ours_only = kall_only = neither = 0
    total_ours_counts = total_kall_counts = 0.0

    rows = []
    for tid in all_tids:
        oc, ot = ours.get(tid, (0.0, 0.0))
        kc, kt = kall.get(tid, (0.0, 0.0))
        total_ours_counts += oc
        total_kall_counts += kc

        if oc > 0 and kc > 0:
            status = 'both'
            both += 1
        elif oc > 0:
            status = 'ours_only'
            ours_only += 1
        elif kc > 0:
            status = 'kall_only'
            kall_only += 1
        else:
            status = 'neither'
            neither += 1

        rows.append((tid, oc, kc, ot, kt, status))

    # Write TSV
    out_path = 'comparison_full.tsv'
    with open(out_path, 'w') as f:
        f.write('transcript\tours_counts\tkall_counts\tours_tpm\tkall_tpm\tstatus\n')
        for tid, oc, kc, ot, kt, status in rows:
            f.write(f'{tid}\t{oc:.1f}\t{kc:.1f}\t{ot:.1f}\t{kt:.1f}\t{status}\n')

    total = len(all_tids)
    detected = both + ours_only + kall_only

    print(f"Transcripts total     : {total:,}")
    print(f"Both detected         : {both:,}  ({100*both/total:.1f}%)")
    print(f"Ours only             : {ours_only:,}  ({100*ours_only/total:.1f}%)")
    print(f"Kallisto only         : {kall_only:,}  ({100*kall_only/total:.1f}%)")
    print(f"Neither               : {neither:,}  ({100*neither/total:.1f}%)")
    print()
    print(f"Total est. counts (ours)    : {total_ours_counts:,.0f}")
    print(f"Total est. counts (kallisto): {total_kall_counts:,.0f}")
    print()

    # Correlation on transcripts detected by both
    both_rows = [(oc, kc) for _, oc, kc, _, _, s in rows if s == 'both']
    if both_rows:
        import math
        log_ours = [math.log2(oc + 1) for oc, _ in both_rows]
        log_kall = [math.log2(kc + 1) for _, kc in both_rows]
        n = len(log_ours)
        mean_o = sum(log_ours) / n
        mean_k = sum(log_kall) / n
        cov = sum((a - mean_o) * (b - mean_k) for a, b in zip(log_ours, log_kall)) / n
        std_o = math.sqrt(sum((a - mean_o)**2 for a in log_ours) / n)
        std_k = math.sqrt(sum((b - mean_k)**2 for b in log_kall) / n)
        r = cov / (std_o * std_k) if std_o * std_k > 0 else 0
        print(f"Pearson r (log2 counts, both detected): {r:.4f}")

    # Top discrepancies: ours >> kallisto and kallisto >> ours
    ours_over  = sorted([(oc-kc, tid, oc, kc) for tid, oc, kc, _, _, s in rows if s == 'both'], reverse=True)
    kall_over  = sorted([(kc-oc, tid, oc, kc) for tid, oc, kc, _, _, s in rows if s == 'both'], reverse=True)

    print()
    print("Top 10 transcripts where ours >> kallisto (both detected):")
    print(f"  {'transcript':<30} {'ours':>10} {'kall':>10} {'diff':>10}")
    for diff, tid, oc, kc in ours_over[:10]:
        print(f"  {tid:<30} {oc:>10.1f} {kc:>10.1f} {diff:>10.1f}")

    print()
    print("Top 10 transcripts where kallisto >> ours (both detected):")
    print(f"  {'transcript':<30} {'ours':>10} {'kall':>10} {'diff':>10}")
    for diff, tid, oc, kc in kall_over[:10]:
        print(f"  {tid:<30} {oc:>10.1f} {kc:>10.1f} {diff:>10.1f}")

    print()
    print(f"Comparison written to: {out_path}")

if __name__ == '__main__':
    main()
