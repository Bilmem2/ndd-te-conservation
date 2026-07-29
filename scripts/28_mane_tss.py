#!/usr/bin/env python3
"""
28_mane_tss.py — does the human result depend on how the TSS is defined?

The main analysis places the promoter window on the outermost 5' coordinate of
the gene feature, which for a multi-promoter gene is the most distal annotated
start rather than the canonical one. Developmental regulators frequently have
alternative promoters, so this is worth checking. Here the window is instead
centred on the TSS of the MANE Select transcript (the single transcript agreed
on by NCBI and EMBL-EBI as the representative isoform) and the human Alu
comparison is repeated.

Output: results/sensitivity/mane_tss.csv
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "hg38"
GL = REPO / "data" / "gene_lists"
OUT = REPO / "results" / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}


def gene_set(fn):
    return {l.strip() for l in open(GL / fn) if l.strip()}


def build(mane_only):
    """Promoter windows; MANE Select transcript TSS, or gene-feature TSS."""
    rows, seen = [], set()
    with gzip.open(DATA / "gtf" / "gencode.v47.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[0] not in MAIN or 'gene_type "protein_coding"' not in f[8]:
                continue
            if mane_only:
                if f[2] != "transcript" or 'tag "MANE_Select"' not in f[8]:
                    continue
            elif f[2] != "gene":
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) in seen:
                continue
            seen.add(m.group(1))
            start, end, strand = int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            rows.append((f[0], max(0, tss - 2000), tss + 2000, m.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def count_alu(prom, alu):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in alu.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        o = np.argsort(m["start"].values)
        s, e = m["start"].values[o], m["end"].values[o]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    ndd, hk = gene_set("HighConfNDD_genes.txt"), gene_set("Housekeeping_genes.txt")
    alu = pd.read_csv(DATA / "alu_rmsk.bed", sep="\t", header=None,
                      usecols=[0, 1, 2], names=["chrom", "start", "end"])
    alu = alu[alu["chrom"].isin(MAIN)]

    rows = []
    for label, mane in [("gene-feature TSS (main analysis)", False),
                        ("MANE Select transcript TSS", True)]:
        prom = build(mane)
        prom["alu_d"] = count_alu(prom, alu) / ((prom["end"] - prom["start"]) / 1000.0)
        nd = prom["alu_d"].values[prom["gene"].isin(ndd).values]
        hkd = prom["alu_d"].values[prom["gene"].isin(hk).values]
        _, p = stats.mannwhitneyu(hkd, nd, alternative="greater")
        genome = prom["alu_d"].mean()
        rows.append(dict(tss_definition=label, n_all_PC=len(prom),
                         n_NDD=len(nd), n_HK=len(hkd),
                         median_alu_NDD=round(float(np.median(nd)), 3),
                         median_alu_HK=round(float(np.median(hkd)), 3),
                         NDD_over_genome=round(float(nd.mean()) / genome, 3),
                         r=round(rb(hkd, nd), 3), p_value=p))
        print(f"{label:34s} n_NDD={len(nd):5d} n_HK={len(hkd):5d} "
              f"r={rows[-1]['r']:+.3f} p={p:.2e} NDD/genome={rows[-1]['NDD_over_genome']:.3f}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "mane_tss.csv", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n=== TSS definition sensitivity (hg38, Alu) ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'mane_tss.csv'}")


if __name__ == "__main__":
    main()
