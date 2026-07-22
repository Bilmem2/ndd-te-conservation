#!/usr/bin/env python3
"""
20_lemur.py — Extend the primate Alu-depletion analysis to a STREPSIRRHINE
(gray mouse lemur, Microcebus murinus, Mmur_3.0), pushing the conservation test
beyond the Simian ancestor (~40 My) toward the base of primates (~65-70 My).

Mouse lemur carries genuine Alu (SINE/Alu, ~496k copies) that amplified
independently in the strepsirrhine lineage from the ancestral Alu — so this is
a direct, same-family test (unlike the mouse B1/B2 analog).

Method matches the manuscript: Alu count per kb at TSS +/-2 kb, HighConfNDD vs
Housekeeping, one-sided Mann-Whitney U (HK > NDD), rank-biserial r. Gene sets
mapped by HGNC symbol (as for the other non-human species).

Data:
  data/mmur3/mmur3.gtf.gz        Ensembl 112 GTF (seqnames '1','2',...)
  data/mmur3/alu_refseq.bed      SINE/Alu from GenArk RepeatMasker (RefSeq names)
  data/mmur3/chromAlias.txt      RefSeq <-> Ensembl seqname map
Outputs:
  results/mmur3/{HighConfNDD_Alu.bed, Housekeeping_Alu.bed, lemur_stats.csv}
"""
import gzip, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
MM = REPO / "data" / "mmur3"
GL = REPO / "data" / "gene_lists"
OUT = REPO / "results" / "mmur3"
OUT.mkdir(parents=True, exist_ok=True)
NAME_RE = re.compile(r'gene_name "([^"]+)"')


def gene_set(fn):
    return set(l.strip() for l in open(GL / fn) if l.strip())


def refseq_to_ensembl():
    """chromAlias: cols = refseq, assembly, genbank, ncbi, ucsc. GTF uses 'ncbi'."""
    m = {}
    for line in open(MM / "chromAlias.txt"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[3]:
            m[f[0]] = f[3]
    return m


def load_alu():
    a = pd.read_csv(MM / "alu_refseq.bed", sep="\t", header=None, names=["rs", "start", "end"])
    m = refseq_to_ensembl()
    a["chrom"] = a["rs"].map(m)
    return a.dropna(subset=["chrom"])[["chrom", "start", "end"]]


def build_promoters(wanted):
    rows, seen = [], set()
    for line in gzip.open(MM / "mmur3.gtf.gz", "rt"):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "gene":
            continue
        if 'protein_coding' not in f[8]:
            continue
        mm = NAME_RE.search(f[8])
        if not mm or mm.group(1) not in wanted or mm.group(1) in seen:
            continue
        seen.add(mm.group(1))
        start, end, strand = int(f[3]), int(f[4]), f[6]
        tss = start if strand == "+" else end
        rows.append((f[0], max(0, tss - 2000), tss + 2000, mm.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def count_alu(prom, alu):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in alu.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        s = np.sort(m["start"].values); e = m["end"].values[np.argsort(m["start"].values)]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def main():
    ndd, hk = gene_set("HighConfNDD_genes.txt"), gene_set("Housekeeping_genes.txt")
    alu = load_alu()
    prom = build_promoters(ndd | hk)
    prom["category"] = np.where(prom["gene"].isin(ndd), "HighConfNDD", "Housekeeping")
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    prom["alu"] = count_alu(prom, alu)
    prom["density"] = prom["alu"] / prom["width_kb"]

    for cat in ("HighConfNDD", "Housekeeping"):
        prom[prom["category"] == cat][["chrom", "start", "end", "gene", "category", "alu"]] \
            .to_csv(OUT / f"{cat}_Alu.bed", sep="\t", header=False, index=False)

    nd = prom[prom["category"] == "HighConfNDD"]["density"]
    hkd = prom[prom["category"] == "Housekeeping"]["density"]
    U, p = stats.mannwhitneyu(hkd, nd, alternative="greater")
    u2, _ = stats.mannwhitneyu(hkd, nd, alternative="two-sided")
    r = 1 - (2 * u2) / (len(hkd) * len(nd))
    sig = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
    row = dict(species="MouseLemur", assembly="Mmur_3.0", mya=70, TE="Alu",
               n_HK=len(hkd), n_NDD=len(nd),
               median_HK=round(hkd.median(), 3), median_NDD=round(nd.median(), 3),
               mean_HK=round(hkd.mean(), 3), mean_NDD=round(nd.mean(), 3),
               p_value=p, r=round(r, 3), sig=sig)
    pd.DataFrame([row]).to_csv(OUT / "lemur_stats.csv", index=False)

    print(f"Mouse lemur Alu: {len(alu)} elements on {alu['chrom'].nunique()} mapped seqs")
    print(f"Promoters: NDD n={len(nd)}, HK n={len(hkd)}")
    print("\n=== Mouse lemur (strepsirrhine, ~70 My) Alu depletion at NDD promoters ===")
    for k, v in row.items():
        print(f"  {k}: {v}")
    print("\nContext — primate Alu effect sizes from the manuscript:")
    print("  Macaque -0.384 | Orangutan -0.367 | Human -0.345 | Gibbon -0.291 | Marmoset -0.147")


if __name__ == "__main__":
    main()
