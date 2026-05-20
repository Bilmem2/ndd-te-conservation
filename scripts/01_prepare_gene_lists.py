import pandas as pd
import os
import re

sync = os.path.expanduser("~/synchronized_constraint/data")
out_dir = "data/gene_lists"
os.makedirs(out_dir, exist_ok=True)

# === 1. HIGH-CONFIDENCE NDD ===
sfari = pd.read_csv(f"{sync}/SFARI-Gene_genes_05-01-2026release_05-12-2026export.csv")
sfari_genes = set(sfari[sfari['gene-score'].isin([1,2])]['gene-symbol'].str.strip().dropna())
print(f"SFARI Tier 1+2: {len(sfari_genes)} gen")

clingen = pd.read_csv(f"{sync}/clingen_epilepsy.csv")
clingen_curated = clingen[clingen['Classification'].isin(['Definitive','Strong'])].copy()
clingen_curated['gene_clean'] = clingen_curated['Gene'].str.replace(r'HGNC:\d+', '', regex=True).str.strip()
clingen_genes = set(clingen_curated['gene_clean'].dropna())
print(f"ClinGen Definitive/Strong: {len(clingen_genes)} gen")

high_conf = sfari_genes | clingen_genes
print(f"High-confidence NDD toplam: {len(high_conf)} gen")

# === 2. BROAD NDD ===
hpo_files = {
    "ADHD":    "genes_for_HP_0007018.txt",
    "Seizure": "genes_for_HP_0001250.txt",
    "Autism":  "genes_for_HP_0000729.txt",
    "ID":      "genes_for_HP_0001249.txt",
}
broad = set()
for name, fname in hpo_files.items():
    df = pd.read_csv(f"{sync}/{fname}", sep='\t')
    genes = set(df.iloc[:,1].str.strip().dropna())
    broad |= genes
    print(f"  HPO {name}: {len(genes)} gen")
broad |= high_conf
print(f"Broad NDD toplam: {len(broad)} gen")

# === 3. HOUSEKEEPING ===
hk_df = pd.read_csv(f"{sync}/Housekeeping_GenesHuman.csv", sep=';')
hk_genes = set(hk_df['Gene.name'].str.strip().dropna())

brain_df = pd.read_csv(f"{sync}/rna_brain_gtex.tsv", sep='\t')
brain_max = brain_df.groupby("Gene name")["TPM"].max()
brain_expressed = set(brain_max[brain_max >= 1].index)

hk_final = (hk_genes & brain_expressed) - broad
print(f"\nHousekeeping (brain-expr, NDD-free): {len(hk_final)} gen")

# === OVERLAP KONTROLÜ ===
print(f"\nHigh-conf ∩ HK: {len(high_conf & hk_final)} (0 olmalı)")
print(f"Broad NDD ∩ HK: {len(broad & hk_final)} (0 olmalı)")

# === KAYDET ===
sets = {
    "HighConfNDD":  high_conf,
    "BroadNDD":     broad,
    "Housekeeping": hk_final,
}
for name, genes in sets.items():
    path = f"{out_dir}/{name}_genes.txt"
    with open(path, 'w') as f:
        f.write('\n'.join(sorted(genes)))
    print(f"Kaydedildi: {name}_genes.txt — {len(genes)} gen")

print("\nGen listeleri hazır.")
