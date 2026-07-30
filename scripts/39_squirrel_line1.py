"""
39_squirrel_line1.py — LINE-1 at squirrel monkey NDD promoters.

Method matches 30_squirrel_monkey.py with LINE/L1 substituted for SINE/Alu.
Completes the LINE-1 panel across all nine genomes.

Inputs : data/saiBol1/{rmsk/rmsk.txt.gz, gtf/saiBol1.gtf.gz}
Output : results/saiBol1/{HighConfNDD_LINE1.bed, Housekeeping_LINE1.bed,
                          squirrel_line1_stats.csv}
"""
import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
SB, OUT = ROOT / "data" / "saiBol1", ROOT / "results" / "saiBol1"
OUT.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("sq", ROOT / "scripts" / "30_squirrel_monkey.py")
sq = importlib.util.module_from_spec(spec)
sys.modules["sq"] = sq
spec.loader.exec_module(sq)


def load_line1():
    rows = []
    with gzip.open(SB / "rmsk" / "rmsk.txt.gz", "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 12 and f[11] == "LINE" and f[12] == "L1":
                rows.append((f[5], int(f[6]), int(f[7])))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


ndd, hk = sq.gene_set("HighConfNDD_genes.txt"), sq.gene_set("Housekeeping_genes.txt")
l1 = load_line1()
prom = sq.build_promoters()
print(f"LINE-1 elements: {len(l1):,};  protein-coding promoters: {len(prom):,}")

prom["l1"] = sq.count_alu(prom, l1)
prom["density"] = prom["l1"] / ((prom["end"] - prom["start"]) / 1000.0)
prom["category"] = np.where(prom["gene"].isin(ndd), "HighConfNDD",
                            np.where(prom["gene"].isin(hk), "Housekeeping", "other"))

for cat in ("HighConfNDD", "Housekeeping"):
    prom[prom["category"] == cat][["chrom", "start", "end", "gene", "l1"]] \
        .assign(score=".", strand=".")[["chrom", "start", "end", "gene", "score", "strand", "l1"]] \
        .to_csv(OUT / f"{cat}_LINE1.bed", sep="\t", header=False, index=False)

nd = prom[prom["category"] == "HighConfNDD"]["density"]
hkd = prom[prom["category"] == "Housekeeping"]["density"]
U, p = stats.mannwhitneyu(hkd, nd, alternative="greater")
r = sq.rb(hkd, nd)
sig = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

row = dict(species="SquirrelMonkey", assembly="SaiBol1.0", mya=40, TE="LINE1",
           n_HK=len(hkd), n_NDD=len(nd),
           median_HK=round(hkd.median(), 3), median_NDD=round(nd.median(), 3),
           mean_HK=round(hkd.mean(), 3), mean_NDD=round(nd.mean(), 3),
           p_value=p, r=round(r, 3), sig=sig)
pd.DataFrame([row]).to_csv(OUT / "squirrel_line1_stats.csv", index=False)

print("\n=== Squirrel monkey LINE-1 at NDD promoters ===")
for k, v in row.items():
    print(f"  {k}: {v}")
print("\nContext - LINE-1 in the other eight genomes:")
print("  significant: mouse lemur -0.064, human -0.049, macaque -0.048, orangutan -0.044")
print("  not significant: gibbon, marmoset, mouse, dog")
