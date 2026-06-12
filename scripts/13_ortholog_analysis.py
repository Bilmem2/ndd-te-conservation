import pandas as pd
import numpy as np
from scipy import stats
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA    = REPO_ROOT / "data"
ORTHO   = f"{DATA}/orthologs"
RESULTS = REPO_ROOT / "results"

# Sembol tablosunu yükle
sym_df = pd.read_csv(f"{ORTHO}/ensembl_to_symbol.tsv", sep='\t',
                     header=0, names=['ensembl_id','hgnc_symbol'])
sym_df = sym_df.dropna(subset=['hgnc_symbol'])
sym_df = sym_df[sym_df['hgnc_symbol'] != '']
print(f"Sembol tablosu: {len(sym_df)} gen")

# Gen listeleri
ndd_genes = set(open(f"{DATA}/gene_lists/HighConfNDD_genes.txt").read().splitlines())
hk_genes  = set(open(f"{DATA}/gene_lists/Housekeeping_genes.txt").read().splitlines())

species = {
    "ponAbe3":  "pabelii",
    "nomLeu3":  "nleucogenys",
    "rheMac10": "mmulatta",
    "calJac4":  "cjacchus",
    "mm10":     "mmusculus",
    "canFam4":  "clfamiliaris",
}

sp_labels = {
    "ponAbe3": "Orangutan", "nomLeu3": "Gibbon",
    "rheMac10": "Macaque",  "calJac4": "Marmoset",
    "mm10": "Mouse",        "canFam4": "Dog",
}

te_map = {
    "ponAbe3": ["Alu", "LINE1"], "nomLeu3": ["Alu", "LINE1"],
    "rheMac10": ["Alu", "LINE1"], "calJac4": ["Alu", "LINE1"],
    "mm10": ["LINE1"], "canFam4": ["LINE1"],
}

def get_density(filepath):
    if not os.path.exists(filepath): return None
    df = pd.read_csv(filepath, sep='\t', header=None)
    return (df.iloc[:,6] / ((df.iloc[:,2]-df.iloc[:,1])/1000)).values

results = []

for sp_code, prefix in species.items():
    raw_file = f"{ORTHO}/{sp_code}_raw.tsv"
    col_names = ['ensembl_id', 'ortholog_id', 'orthology_type']
    ortho_df = pd.read_csv(raw_file, sep='\t', header=0, names=col_names)

    # 1:1 ortologları filtrele
    one2one = ortho_df[ortho_df['orthology_type']=='ortholog_one2one'].copy()
    one2one = one2one.merge(sym_df, on='ensembl_id', how='inner')

    # NDD ve HK kümelerini filtrele
    ndd_ortho = set(one2one[one2one['hgnc_symbol'].isin(ndd_genes)]['hgnc_symbol'])
    hk_ortho  = set(one2one[one2one['hgnc_symbol'].isin(hk_genes)]['hgnc_symbol'])

    print(f"\n{sp_labels[sp_code]}: 1:1 NDD={len(ndd_ortho)}, HK={len(hk_ortho)}")

    for te in te_map[sp_code]:
        # Orijinal intersect dosyalarını yükle ama sadece ortolog genleri tut
        hk_f  = f"{RESULTS}/{sp_code}/Housekeeping_{te}.bed"
        ndd_f = f"{RESULTS}/{sp_code}/HighConfNDD_{te}.bed"
        if not os.path.exists(hk_f) or not os.path.exists(ndd_f):
            continue

        hk_bed  = pd.read_csv(hk_f,  sep='\t', header=None)
        ndd_bed = pd.read_csv(ndd_f, sep='\t', header=None)

        # Sadece 1:1 ortolog olanları al
        hk_bed  = hk_bed[hk_bed.iloc[:,3].isin(hk_ortho)]
        ndd_bed = ndd_bed[ndd_bed.iloc[:,3].isin(ndd_ortho)]

        if sp_code == 'mm10':
            # Mouse büyük harf
            hk_bed  = hk_bed[hk_bed.iloc[:,3].str.upper().isin(
                {g.upper() for g in hk_ortho})]
            ndd_bed = ndd_bed[ndd_bed.iloc[:,3].str.upper().isin(
                {g.upper() for g in ndd_ortho})]

        hk_d  = (hk_bed.iloc[:,6] / ((hk_bed.iloc[:,2]-hk_bed.iloc[:,1])/1000)).values
        ndd_d = (ndd_bed.iloc[:,6] / ((ndd_bed.iloc[:,2]-ndd_bed.iloc[:,1])/1000)).values

        if len(hk_d) < 10 or len(ndd_d) < 10:
            print(f"  {te}: yetersiz n (HK={len(hk_d)}, NDD={len(ndd_d)})")
            continue

        u, p = stats.mannwhitneyu(hk_d, ndd_d, alternative='greater')
        r = 1 - (2*u)/(len(hk_d)*len(ndd_d))
        sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
        print(f"  {te}: n_HK={len(hk_d)}, n_NDD={len(ndd_d)} | "
              f"p={p:.2e} | r={r:.3f} | {sig}")
        results.append({
            "species": sp_labels[sp_code], "TE": te,
            "n_HK": len(hk_d), "n_NDD": len(ndd_d),
            "p_value": p, "r": round(r,3), "sig": sig,
            "method": "ortholog_1to1"
        })

df = pd.DataFrame(results)
df.to_csv(f"{RESULTS}/statistics_ortholog.csv", index=False)
print("\n=== ALU (ortolog-validated) ===")
print(df[df['TE']=='Alu'].to_string(index=False))
print("\n=== LINE-1 (ortolog-validated) ===")
print(df[df['TE']=='LINE1'].to_string(index=False))
