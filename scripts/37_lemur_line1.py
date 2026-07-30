"""
37_lemur_line1.py — LINE-1 at mouse lemur NDD promoters.

Method matches 20_lemur.py with LINE/L1 substituted for SINE/Alu: LINE-1 count
per kb at TSS +/-2 kb, HighConfNDD versus Housekeeping, one-sided
Mann-Whitney U (HK > NDD), rank-biserial r.

Inputs : data/mmur3/{mmur3.gtf.gz, chromAlias.txt, repeatMasker.out.gz}
Output : results/mmur3/{HighConfNDD_LINE1.bed, Housekeeping_LINE1.bed,
                       lemur_line1_stats.csv}
"""
import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
MM, OUT = ROOT / "data" / "mmur3", ROOT / "results" / "mmur3"
OUT.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("lem", ROOT / "scripts" / "20_lemur.py")
lem = importlib.util.module_from_spec(spec)
sys.modules["lem"] = lem
spec.loader.exec_module(lem)

L1_BED = MM / "line1_refseq.bed"

if not L1_BED.exists():
    print("extracting LINE/L1 from RepeatMasker output ...")
    n = 0
    with gzip.open(MM / "repeatMasker.out.gz", "rt", errors="replace") as fh, \
            open(L1_BED, "w", newline="") as out:
        for line in fh:
            p = line.split()
            if len(p) > 10 and p[0].isdigit() and p[10] == "LINE/L1":
                out.write(f"{p[4]}\t{int(p[5]) - 1}\t{p[6]}\n")
                n += 1
    print(f"  wrote {n:,} LINE-1 intervals")

l1 = pd.read_csv(L1_BED, sep="\t", header=None, names=["rs", "start", "end"])
l1["chrom"] = l1["rs"].map(lem.refseq_to_ensembl())
l1 = l1.dropna(subset=["chrom"])[["chrom", "start", "end"]]

ndd, hk = lem.gene_set("HighConfNDD_genes.txt"), lem.gene_set("Housekeeping_genes.txt")
prom = lem.build_promoters(ndd | hk)
prom["category"] = np.where(prom["gene"].isin(ndd), "HighConfNDD", "Housekeeping")
prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
prom["l1"] = lem.count_alu(prom, l1)          # generic interval counter
prom["density"] = prom["l1"] / prom["width_kb"]

for cat in ("HighConfNDD", "Housekeeping"):
    prom[prom["category"] == cat][["chrom", "start", "end", "gene", "category", "l1"]] \
        .to_csv(OUT / f"{cat}_LINE1.bed", sep="\t", header=False, index=False)

nd = prom[prom["category"] == "HighConfNDD"]["density"]
hkd = prom[prom["category"] == "Housekeeping"]["density"]
U, p = stats.mannwhitneyu(hkd, nd, alternative="greater")
u2, _ = stats.mannwhitneyu(hkd, nd, alternative="two-sided")
r = 1 - (2 * u2) / (len(hkd) * len(nd))
sig = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

row = dict(species="MouseLemur", assembly="Mmur_3.0", mya=70, TE="LINE1",
           n_HK=len(hkd), n_NDD=len(nd),
           median_HK=round(hkd.median(), 3), median_NDD=round(nd.median(), 3),
           mean_HK=round(hkd.mean(), 3), mean_NDD=round(nd.mean(), 3),
           p_value=p, r=round(r, 3), sig=sig)
pd.DataFrame([row]).to_csv(OUT / "lemur_line1_stats.csv", index=False)

print(f"\nLINE-1 elements mapped: {len(l1):,} on {l1['chrom'].nunique()} sequences")
print("\n=== Mouse lemur LINE-1 at NDD promoters ===")
for k, v in row.items():
    print(f"  {k}: {v}")
print("\nContext - LINE-1 effect sizes already in the manuscript:")
print("  Human -0.049** | Orangutan -0.044** | Macaque -0.048**")
print("  Gibbon, marmoset, mouse, dog: not significant")
