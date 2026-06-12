import pandas as pd
import numpy as np
from scipy import stats
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS = REPO_ROOT / "results"
DATA = REPO_ROOT / "data"

# === 1. gnomAD pLI yükle ===
print("gnomAD constraint yükleniyor...")
gnomad = pd.read_csv(f"{DATA}/gnomad_constraint.tsv", sep='\t',
                     usecols=['gene', 'lof.pLI', 'canonical'])
gnomad = gnomad[gnomad['canonical'] == True].copy()
gnomad = gnomad.dropna(subset=['lof.pLI'])
gnomad = gnomad.rename(columns={'gene': 'gene_name', 'lof.pLI': 'pLI'})
gnomad = gnomad[['gene_name', 'pLI']].drop_duplicates('gene_name')
print(f"gnomAD: {len(gnomad)} gen")

# === 2. hg38 Alu density yükle ===
print("Alu density yükleniyor...")
dfs = []
for cat in ['HighConfNDD', 'Housekeeping']:
    f = f"{RESULTS}/hg38/{cat}_Alu.bed"
    df = pd.read_csv(f, sep='\t', header=None)
    df.columns = ['chr','start','end','gene_name','score','strand','alu_count']
    df['density'] = df['alu_count'] / ((df['end'] - df['start']) / 1000)
    df['category'] = cat
    dfs.append(df[['gene_name','density','category']])

alu_df = pd.concat(dfs)
print(f"Alu density: {len(alu_df)} gen")

# === 3. Birleştir ===
merged = alu_df.merge(gnomad, on='gene_name', how='inner')
print(f"Birleştirilen: {len(merged)} gen")

# === 4. Korelasyon ===
r, p = stats.spearmanr(merged['pLI'], merged['density'])
print(f"\nTüm genler Spearman r={r:.3f}, p={p:.2e}, n={len(merged)}")

ndd = merged[merged['category']=='HighConfNDD']
r_ndd, p_ndd = stats.spearmanr(ndd['pLI'], ndd['density'])
print(f"HighConfNDD Spearman r={r_ndd:.3f}, p={p_ndd:.2e}, n={len(ndd)}")

hk = merged[merged['category']=='Housekeeping']
r_hk, p_hk = stats.spearmanr(hk['pLI'], hk['density'])
print(f"Housekeeping Spearman r={r_hk:.3f}, p={p_hk:.2e}, n={len(hk)}")

# === 5. pLI gruplarına göre ===
merged['pLI_group'] = pd.cut(merged['pLI'],
    bins=[0, 0.1, 0.5, 0.9, 1.001],
    labels=['Tolerant\n(pLI<0.1)', 'Moderate\n(0.1-0.5)',
            'Constrained\n(0.5-0.9)', 'Highly Constrained\n(pLI>0.9)'])

print("\n=== pLI Gruplarına Göre Alu Density (median) ===")
print(merged.groupby('pLI_group', observed=True)['density'].agg(['median','count']))

merged.to_csv(f"{RESULTS}/pli_alu_merged.csv", index=False)
print("\nKaydedildi.")
