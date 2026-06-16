import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

GTF     = "tests/data/hg38.gtf"
RESULTS = "results.tsv"
OUTDIR  = "plots"
BINS    = 800   # number of bins per chromosome

BG     = 'white'
FILL_C = '#1565c0'
EDGE_C = '#0288d1'
AXIS_C = '#444444'
TITLE_C = '#111111'

print("Parsing GTF...")
records = []
with open(GTF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 9 or parts[2] != 'transcript':
            continue
        m = re.search(r'transcript_id "([^"]+)"', parts[8])
        if m:
            records.append({
                'transcript': m.group(1),
                'chrom':      parts[0],
                'start':      int(parts[3]),
                'end':        int(parts[4]),
            })

gtf = pd.DataFrame(records)
print(f"  {len(gtf)} transcripts in GTF")

df = pd.read_csv(RESULTS, sep='\t')
df = df[df['est_counts'] > 0]
df = df.merge(gtf, on='transcript', how='inner')

standard = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
df = df[df['chrom'].isin(standard)].copy()
print(f"  {len(df)} transcripts with mapped reads on standard chromosomes")

df['mid'] = (df['start'] + df['end']) / 2

chrom_len   = gtf[gtf['chrom'].isin(standard)].groupby('chrom')['end'].max()
chrom_order = [c for c in standard if c in df['chrom'].unique()]

os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({
    'figure.facecolor':  BG,
    'axes.facecolor':    BG,
    'savefig.facecolor': BG,
    'text.color':        TITLE_C,
    'axes.labelcolor':   AXIS_C,
    'xtick.color':       AXIS_C,
    'ytick.color':       AXIS_C,
})

for chrom in chrom_order:
    sub  = df[df['chrom'] == chrom]
    clen = float(chrom_len.get(chrom, sub['mid'].max()))

    edges = np.linspace(0, clen, BINS + 1)
    mids  = (edges[:-1] + edges[1:]) / 2

    # sum TPM into bins
    tpm_sum, _ = np.histogram(sub['mid'].values,
                               bins=edges,
                               weights=sub['tpm'].values)

    # log-compress for display (0 stays 0)
    height = np.where(tpm_sum > 0, np.log10(tpm_sum + 1), 0)

    fig, ax = plt.subplots(figsize=(18, 3.5))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    bar_w = (edges[1] - edges[0]) / 1e6

    ax.bar(mids / 1e6, height, width=bar_w * 0.9,
           color=FILL_C, edgecolor='none', linewidth=0, alpha=0.75)

    # subtle glow on top edge
    ax.step(np.append(mids / 1e6, mids[-1] / 1e6 + bar_w),
            np.append(height, 0),
            where='mid', color=EDGE_C, lw=0.6, alpha=0.7)

    ax.set_xlim(0, clen / 1e6)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Genomic position (Mbp)', fontsize=9, color=AXIS_C)
    ax.set_ylabel('log₁₀(ΣTPM + 1)', fontsize=9, color=AXIS_C)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color('#cccccc')
    ax.spines['left'].set_color('#cccccc')
    ax.tick_params(colors=AXIS_C, labelsize=8)

    n_above = (sub['tpm'] >= 1).sum()
    ax.set_title(
        f'{chrom}  —  {len(sub):,} transcripts  '
        f'({n_above:,} with TPM ≥ 1)',
        fontsize=11, color=TITLE_C, pad=6)

    out = os.path.join(OUTDIR, f'{chrom}.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  {out}")

print(f"Done — {len(chrom_order)} images in {OUTDIR}/")
