import subprocess
import pandas as pd
import numpy as np
from scipy import stats
import os

DATA = os.path.expanduser("~/comparative_TE_study/data")
RESULTS = os.path.expanduser("~/comparative_TE_study/results")
os.makedirs(f"{RESULTS}/sensitivity", exist_ok=True)

alu_bed = f"{DATA}/hg38/rmsk/Alu.bed"
gtf_gz = f"{DATA}/hg38/gtf/gencode.v47.gtf.gz"
gene_lists = f"{DATA}/gene_lists"

windows = [500, 1000, 2000, 3000]
results = []

print("Pencere hassasiyeti testi (hg38, Alu)...")

for window in windows:
    print(f"\n  Window: ±{window}bp")

    # Her pencere için promotor BED oluştur
    for cat in ["HighConfNDD", "Housekeeping"]:
        prom_out = f"{RESULTS}/sensitivity/{cat}_w{window}.bed"

        # GTF'den promotor çıkar
        awk_cmd = f"""zcat {gtf_gz} | awk -v W={window} 'BEGIN{{OFS="\\t"}}
        $3=="gene" && /protein_coding/ {{
            match($0, /gene_name "([^"]+)"/, a); gname=a[1];
            if (gname=="") next
            if ($7=="+") {{ tss=$4; start=tss-W; end=tss+W }}
            else          {{ tss=$5; start=tss-W; end=tss+W }}
            if (start<1) start=1
            print $1, start, end, gname, ".", $7
        }}' | grep -Fw -f {gene_lists}/{cat}_genes.txt | sort -k1,1 -k2,2n > {prom_out}"""

        subprocess.run(awk_cmd, shell=True)

        # Alu intersect
        int_out = subprocess.run(
            ["bedtools", "intersect", "-a", prom_out, "-b", alu_bed, "-c"],
            capture_output=True, text=True).stdout
        df = pd.read_csv(pd.io.common.StringIO(int_out),
                         sep='\t', header=None)
        df['density'] = df.iloc[:,6] / ((df.iloc[:,2] - df.iloc[:,1]) / 1000)
        df['category'] = cat
        df['window'] = window
        results.append(df[['density','category','window']])

# İstatistikler
df_all = pd.concat(results)

print("\n=== PENCERE HASSASIYETI SONUÇLARI ===")
print(f"{'Window':>8} {'HK median':>12} {'NDD median':>12} {'p-value':>12} {'r':>8} {'sig':>5}")

for w in windows:
    hk = df_all[(df_all['window']==w) & (df_all['category']=='Housekeeping')]['density']
    ndd = df_all[(df_all['window']==w) & (df_all['category']=='HighConfNDD')]['density']
    u, p = stats.mannwhitneyu(hk, ndd, alternative='greater')
    r = 1 - (2*u)/(len(hk)*len(ndd))
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    print(f"{w:>8} {hk.median():>12.3f} {ndd.median():>12.3f} {p:>12.2e} {r:>8.3f} {sig:>5}")

df_all.to_csv(f"{RESULTS}/sensitivity/window_sensitivity.csv", index=False)
print("\nKaydedildi.")
