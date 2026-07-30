"""
38_squirrel_ortholog.py — 1:1 ortholog control for the squirrel monkey.

Repeats the Alu comparison restricted to genes with a one-to-one human ortholog,
which tests whether the effect depends on symbol-based gene matching or on how
completely the assembly is annotated.

Inputs : results/saiBol1/{HighConfNDD_Alu.bed, Housekeeping_Alu.bed}
         data/orthologs/saiBol1_raw.tsv  (Ensembl BioMart, human -> S. boliviensis)
Output : results/saiBol1/squirrel_ortholog_validation.csv
"""
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "saiBol1"
COLS = ["chrom", "start", "end", "gene", "score", "strand", "alu"]


def load(cat):
    df = pd.read_csv(RES / f"{cat}_Alu.bed", sep="\t", header=None, names=COLS,
                     lineterminator="\n")
    for c in ("start", "end", "alu"):
        df[c] = pd.to_numeric(df[c].astype(str).str.strip(), errors="coerce")
    df = df.dropna(subset=["start", "end", "alu"])
    df["density"] = df["alu"] / ((df["end"] - df["start"]) / 1000.0)
    return df


one2one = set()
for line in open(ROOT / "data" / "orthologs" / "saiBol1_raw.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 3 and f[2] == "ortholog_one2one" and f[0]:
        one2one.add(f[0])
print(f"1:1 orthologs from BioMart: {len(one2one):,}")

ndd, hk = load("HighConfNDD"), load("Housekeeping")


def test(h, n):
    U, p = stats.mannwhitneyu(h, n, alternative="greater")
    u2, _ = stats.mannwhitneyu(h, n, alternative="two-sided")
    return round(1 - (2 * u2) / (len(h) * len(n)), 3), p


rows = []
r_all, p_all = test(hk.density, ndd.density)
rows.append(dict(method="symbol-based (all mapped)", n_HK=len(hk), n_NDD=len(ndd),
                 p_value=p_all, r=r_all))

ndd_o = ndd[ndd.gene.isin(one2one)]
hk_o = hk[hk.gene.isin(one2one)]
r_o, p_o = test(hk_o.density, ndd_o.density)
rows.append(dict(method="ortholog one2one", n_HK=len(hk_o), n_NDD=len(ndd_o),
                 p_value=p_o, r=r_o))

out = pd.DataFrame(rows)
out.to_csv(RES / "squirrel_ortholog_validation.csv", index=False)

print()
print(out.to_string(index=False))
print(f"\n|delta r| = {abs(r_o - r_all):.3f}")
print(f"dropped by the ortholog restriction: {len(hk)-len(hk_o)} HK, {len(ndd)-len(ndd_o)} NDD")
print(f"\nwrote {RES / 'squirrel_ortholog_validation.csv'}")
