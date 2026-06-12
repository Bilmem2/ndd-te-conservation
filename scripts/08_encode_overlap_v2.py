import subprocess
import pandas as pd
import numpy as np
from scipy import stats
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS = REPO_ROOT / "results"
DATA = REPO_ROOT / "data"

ctcf = f"{DATA}/hg38/ctcf_peaks.bed"
dnase = f"{DATA}/hg38/brain_dnase.bed.gz"

# NDD ve HK için Alu + CTCF + DNase birlikte analiz
all_results = []

for cat in ["HighConfNDD", "Housekeeping"]:
    prom = f"{DATA}/hg38/promoters/promoters_{cat}.bed"

    # Alu count
    alu_out = subprocess.run(
        ["bedtools", "intersect", "-a", prom,
         "-b", f"{DATA}/hg38/rmsk/Alu.bed", "-c"],
        capture_output=True, text=True).stdout
    df = pd.read_csv(pd.io.common.StringIO(alu_out),
                     sep='\t', header=None,
                     names=['chr','start','end','gene','score','strand','alu_count'])

    # CTCF count
    ctcf_out = subprocess.run(
        ["bedtools", "intersect", "-a", prom, "-b", ctcf, "-c"],
        capture_output=True, text=True).stdout
    df_ctcf = pd.read_csv(pd.io.common.StringIO(ctcf_out),
                           sep='\t', header=None,
                           names=['chr','start','end','gene','score','strand','ctcf_count'])

    # DNase count
    dnase_out = subprocess.run(
        ["bedtools", "intersect", "-a", prom, "-b", dnase, "-c"],
        capture_output=True, text=True).stdout
    df_dnase = pd.read_csv(pd.io.common.StringIO(dnase_out),
                            sep='\t', header=None,
                            names=['chr','start','end','gene','score','strand','dnase_count'])

    # Birleştir
    df = df.merge(df_ctcf[['gene','ctcf_count']], on='gene')
    df = df.merge(df_dnase[['gene','dnase_count']], on='gene')
    df['alu_free'] = df['alu_count'] == 0
    df['category'] = cat
    all_results.append(df)

merged = pd.concat(all_results)

print("=== Alu-free vs Alu-positive: CTCF ve DNase ===\n")
for cat in ["HighConfNDD", "Housekeeping"]:
    sub = merged[merged['category']==cat]
    alu_free = sub[sub['alu_free']]
    alu_pos  = sub[~sub['alu_free']]

    print(f"--- {cat} ---")
    print(f"  Alu-free: n={len(alu_free)}, CTCF median={alu_free['ctcf_count'].median():.1f}, DNase median={alu_free['dnase_count'].median():.1f}")
    print(f"  Alu-pos:  n={len(alu_pos)},  CTCF median={alu_pos['ctcf_count'].median():.1f}, DNase median={alu_pos['dnase_count'].median():.1f}")

    for metric in ['ctcf_count', 'dnase_count']:
        u, p = stats.mannwhitneyu(alu_free[metric], alu_pos[metric],
                                   alternative='greater')
        r = 1 - (2*u)/(len(alu_free)*len(alu_pos))
        sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
        print(f"  {metric}: Alu-free > Alu-pos | p={p:.2e} | r={r:.3f} | {sig}")
    print()

merged.to_csv(f"{RESULTS}/encode_overlap_v2.csv", index=False)
print("Kaydedildi.")
