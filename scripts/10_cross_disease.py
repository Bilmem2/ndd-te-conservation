import subprocess
import pandas as pd
import numpy as np
from scipy import stats
import os
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA = REPO_ROOT / "data"
RESULTS = REPO_ROOT / "results"
os.makedirs(f"{RESULTS}/cross_disease", exist_ok=True)

alu_bed = f"{DATA}/hg38/rmsk/Alu.bed"
gtf_gz = f"{DATA}/hg38/gtf/gencode.v47.gtf.gz"
clinvar_path = f"{DATA}/clinvar_variants.txt.gz"

# ClinVar oku — header'dan kolon isimlerini al
print("ClinVar işleniyor...")
clinvar = pd.read_csv(clinvar_path, sep='\t', low_memory=False,
                      usecols=['#AlleleID','GeneSymbol','ClinicalSignificance',
                               'PhenotypeList','Assembly'])
clinvar = clinvar.rename(columns={'#AlleleID':'AlleleID'})
clinvar = clinvar[clinvar['Assembly']=='GRCh38']
patho = clinvar[clinvar['ClinicalSignificance'].str.contains(
    'Pathogenic', na=False, case=False)]

# ClinVar'ın GeneSymbol alanı bir varyant için birden çok geni ";" ile ayırarak
# verir ve gen atanmamışsa düz "-" yazar. Alanı olduğu gibi kullanmak hem çoklu
# girdilerin içindeki genleri kaybettirir hem de "-" değerini bir "gen adı"na
# çevirir; bu değer aşağıdaki promotör seçiminde her eksi-strand geni eşleştirir.
SYMBOL_RE = re.compile(r'^[A-Za-z][A-Za-z0-9_.@-]*$')


def symbols(series):
    out = set()
    for entry in series.dropna():
        for tok in str(entry).split(';'):
            tok = tok.strip()
            if tok and tok != '-' and SYMBOL_RE.match(tok):
                out.add(tok)
    return out


# Kardiyovasküler
cardio_terms = ['cardiomyopathy','arrhythmia','channelopathy',
                'long QT','Brugada','heart failure']
cardio_mask = patho['PhenotypeList'].str.contains(
    '|'.join(cardio_terms), na=False, case=False)
cardio_genes = symbols(patho[cardio_mask]['GeneSymbol'])

# Mendelian
mendelian_genes = symbols(patho['GeneSymbol'])

ndd_genes = set(open(f"{DATA}/gene_lists/HighConfNDD_genes.txt").read().splitlines())
hk_genes  = set(open(f"{DATA}/gene_lists/Housekeeping_genes.txt").read().splitlines())

cardio_clean    = cardio_genes - ndd_genes - hk_genes
mendelian_clean = mendelian_genes - ndd_genes - hk_genes

print(f"Kardiyovasküler: {len(cardio_clean)} gen")
print(f"Mendelian: {len(mendelian_clean)} gen")

# Gen listelerini kaydet
for name, genes in [("Cardiovascular", cardio_clean),
                     ("Mendelian", mendelian_clean)]:
    with open(f"{DATA}/gene_lists/{name}_genes.txt", 'w') as f:
        f.write('\n'.join(sorted(genes)))

print("\nAlu analizi yapılıyor...")
results = []

gene_sets = ["HighConfNDD", "Cardiovascular", "Mendelian", "Housekeeping"]

for cat in gene_sets:
    if cat in ["HighConfNDD", "Housekeeping"]:
        prom_out = f"{DATA}/hg38/promoters/promoters_{cat}.bed"
    else:
        prom_out = f"{RESULTS}/cross_disease/{cat}_promoters.bed"
        gene_file = f"{DATA}/gene_lists/{cat}_genes.txt"
        # Gen sembolü ALANINA göre birebir eşleşme; satır bazlı `grep -Fw` değil
        # (o, sembolü koordinat/strand sütunlarıyla da eşleştirebiliyordu).
        awk_cmd = f"""zcat {gtf_gz} | awk -v W=2000 'BEGIN{{OFS="\\t"}}
        NR==FNR {{ want[$1]=1; next }}
        $3=="gene" && /protein_coding/ {{
            match($0, /gene_name "([^"]+)"/, a); gname=a[1];
            if (gname=="" || !(gname in want)) next
            if ($7=="+") {{ tss=$4; start=tss-W; end=tss+W }}
            else          {{ tss=$5; start=tss-W; end=tss+W }}
            if (start<1) start=1
            print $1, start, end, gname, ".", $7
        }}' {gene_file} - | sort -k1,1 -k2,2n > {prom_out}"""
        subprocess.run(awk_cmd, shell=True)

    int_out = subprocess.run(
        ["bedtools", "intersect", "-a", prom_out, "-b", alu_bed, "-c"],
        capture_output=True, text=True).stdout
    df = pd.read_csv(pd.io.common.StringIO(int_out), sep='\t', header=None)
    df['density'] = df.iloc[:,6] / ((df.iloc[:,2]-df.iloc[:,1])/1000)
    df['category'] = cat
    results.append(df[['density','category']])
    print(f"  {cat}: n={len(df)}, median={df['density'].median():.3f}")

df_all = pd.concat(results)
hk = df_all[df_all['category']=='Housekeeping']['density']

print("\n=== CROSS-DISEASE SPESİFİSİTE ===")
print(f"{'Grup':>20} {'n':>6} {'median':>8} {'p-value':>12} {'r':>8} {'sig':>5}")
for cat in ["HighConfNDD", "Cardiovascular", "Mendelian"]:
    grp = df_all[df_all['category']==cat]['density']
    u, p = stats.mannwhitneyu(hk, grp, alternative='greater')
    r = 1 - (2*u)/(len(hk)*len(grp))
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    print(f"{cat:>20} {len(grp):>6} {grp.median():>8.3f} {p:>12.2e} {r:>8.3f} {sig:>5}")

df_all.to_csv(f"{RESULTS}/cross_disease/cross_disease_results.csv", index=False)
print("\nKaydedildi.")
