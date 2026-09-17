import pandas as pd
import os
import re
import sys
from pathlib import Path

# the console is cp1254 on a Turkish Windows install, which cannot encode the
# set-intersection sign used below
sys.stdout.reconfigure(encoding="utf-8")

REPO_ROOT = Path(__file__).resolve().parent.parent
# Raw gene-set source files, downloaded by hand rather than by 00_download_data.sh
# because each comes from a portal with its own export button: SFARI Gene,
# ClinGen Epilepsy, the four HPO term lists, HRT Atlas and the Human Protein
# Atlas GTEx brain table. All are free; the README "Other datasets" table gives
# the URLs and the expected filenames. Place them under data/sources/.
sync = REPO_ROOT / "data" / "sources"
out_dir = REPO_ROOT / "data" / "gene_lists"
os.makedirs(out_dir, exist_ok=True)

# === 1. HIGH-CONFIDENCE NDD ===
# SFARI names its export after the release and export dates, so the filename
# changes every time it is downloaded. The published lists were built from the
# 01 May 2026 release; that file is used when present, and otherwise the newest
# export in data/sources/ is taken and its name printed, because the release
# determines the gene set and must be recorded with any result.
SFARI_PUBLISHED = "SFARI-Gene_genes_05-01-2026release_05-12-2026export.csv"
sfari_path = Path(sync) / SFARI_PUBLISHED
if not sfari_path.exists():
    found = sorted(Path(sync).glob("SFARI-Gene_genes_*.csv"))
    if not found:
        raise SystemExit(f"data/sources/ icinde SFARI export yok (beklenen: {SFARI_PUBLISHED})")
    sfari_path = found[-1]
    print(f"UYARI: yayinlanan SFARI surumu yok, kullanilan: {sfari_path.name}")
    print("       bu, data/gene_lists/ icindeki listelerden farkli bir kume uretir")
else:
    print(f"SFARI surumu: {sfari_path.name} (yayinlanan listelerin kaynagi)")

sfari = pd.read_csv(sfari_path)
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
# data/gene_lists/ holds the sets every reported analysis was run on. Writing
# them from a different source release would silently put the pipeline out of
# step with the manuscript, so that only happens on the published SFARI release
# or when the caller asks for it explicitly with --overwrite.
sets = {
    "HighConfNDD":  high_conf,
    "BroadNDD":     broad,
    "Housekeeping": hk_final,
}
on_published_release = sfari_path.name == SFARI_PUBLISHED
forced = "--overwrite" in sys.argv

if on_published_release or forced:
    for name, genes in sets.items():
        path = f"{out_dir}/{name}_genes.txt"
        with open(path, 'w') as f:
            f.write('\n'.join(sorted(genes)))
        print(f"Kaydedildi: {name}_genes.txt — {len(genes)} gen")
    print("\nGen listeleri hazır.")
else:
    # The lists are left alone, but the comparison is worth keeping: it says how
    # far the current source releases have moved from the ones behind the
    # reported results, gene by gene.
    print("\nYAZILMADI. Kaynak surumu yayinlanan listelerinkinden farkli:")
    diff_rows = []
    for name, genes in sets.items():
        cur_path = Path(out_dir) / f"{name}_genes.txt"
        if not cur_path.exists():
            # BroadNDD is an intermediate for the housekeeping exclusion and is
            # not committed, so there is nothing to compare it against
            print(f"  {name:<14} yayinlanan liste yok, atlandi ({len(genes)} gen)")
            continue
        cur = {l.strip() for l in open(cur_path, encoding="utf-8") if l.strip()}
        print(f"  {name:<14} mevcut {len(cur):>5}   bu kaynaklardan {len(genes):>5}"
              f"   ortak {len(cur & genes):>5}")
        for g in sorted(cur - genes):
            diff_rows.append((name, g, "only_in_published"))
        for g in sorted(genes - cur):
            diff_rows.append((name, g, "only_in_current_sources"))

    audit = Path(out_dir).parent.parent / "results" / "gene_set_source_diff.tsv"
    audit.parent.mkdir(parents=True, exist_ok=True)
    with open(audit, "w", encoding="utf-8", newline="") as f:
        f.write(f"# published SFARI: {SFARI_PUBLISHED}\n")
        f.write(f"# current  SFARI: {sfari_path.name}\n")
        f.write("gene_set\tgene\tstatus\n")
        for row in diff_rows:
            f.write("\t".join(row) + "\n")
    print(f"\nFark dokumu: {audit.relative_to(REPO_ROOT)} ({len(diff_rows)} satir)")
    print("data/gene_lists/ makalede raporlanan analizlerin kumelerini tutuyor.")
    print("Guncel kaynaklarla yeniden uretmek icin: --overwrite")
