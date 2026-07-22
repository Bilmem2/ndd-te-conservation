#!/usr/bin/env python3
"""
15_context_controls.py — Confounder controls the MGG reviewer asked for:
local GENE DENSITY and RECOMBINATION RATE (hg38).

For each NDD / housekeeping promoter we compute:
  * gene density  = # protein-coding genes whose TSS lies within +/-500 kb of
                    this promoter's TSS (i.e. genes per Mb), from GENCODE v47.
  * recombination = local sex-averaged rate (cM/Mb) over a +/-50 kb window
                    around the TSS, interpolated from the Beagle GRCh38 genetic
                    map (deCODE/HapMap-derived).

Then two questions, mirroring the GC-content control already in the paper:
  1. Do NDD and HK promoters differ systematically in gene density / recombination?
     (two-sided Mann-Whitney U)
  2. Does Alu depletion (HK > NDD) persist WITHIN density- and recombination-
     matched strata?  If yes, neither is a sufficient explanation.

Inputs:
  results/hg38/{HighConfNDD,Housekeeping}_Alu.bed   (promoter windows + fixed Alu count)
  data/hg38/gtf/gencode.v47.gtf.gz                  (gene TSS universe)
  data/hg38/recomb/chr_in_chrom_field/plink.chrchr*.GRCh38.map  (Beagle map)
Outputs:
  results/context/{context_per_promoter.tsv, context_group_diff.csv,
                   alu_by_density_stratum.csv, alu_by_recomb_stratum.csv}
No bedtools dependency.
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
DATA = REPO / "data"
OUT = RESULTS / "context"
OUT.mkdir(parents=True, exist_ok=True)

DENSITY_FLANK = 500_000    # +/-500 kb -> count = genes per Mb
RECOMB_FLANK = 50_000      # +/-50 kb window for local cM/Mb


def load_promoters(cat):
    f = RESULTS / "hg38" / f"{cat}_Alu.bed"
    df = pd.read_csv(f, sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand", "fixed_alu"])
    df["category"] = cat
    df["tss"] = np.where(df["strand"] == "+", df["start"] + 2000, df["end"] - 2000)
    return df[["chrom", "start", "end", "gene", "strand", "fixed_alu", "category", "tss"]]


def gene_tss_from_gencode():
    """All protein-coding gene TSS (strand-aware) -> dict chrom -> sorted np.array of TSS bp."""
    path = DATA / "hg38" / "gtf" / "gencode.v47.gtf.gz"
    by_chrom = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            if 'gene_type "protein_coding"' not in f[8]:
                continue
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            by_chrom.setdefault(chrom, []).append(tss)
    return {c: np.sort(np.array(v)) for c, v in by_chrom.items()}


def recomb_map():
    """chrom -> (bp array, cM array) from Beagle GRCh38 plink maps."""
    d = DATA / "hg38" / "recomb" / "chr_in_chrom_field"
    maps = {}
    for f in sorted(d.glob("plink.chrchr*.GRCh38.map")):
        arr = np.loadtxt(f, usecols=(0, 2, 3), dtype=object)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        chrom = arr[0, 0]
        cm = arr[:, 1].astype(float)
        bp = arr[:, 2].astype(float)
        order = np.argsort(bp)
        maps[chrom] = (bp[order], cm[order])
    return maps


def local_gene_density(prom, tss_by_chrom):
    dens = np.zeros(len(prom))
    for i, row in enumerate(prom.itertuples()):
        arr = tss_by_chrom.get(row.chrom)
        if arr is None:
            dens[i] = np.nan
            continue
        lo = np.searchsorted(arr, row.tss - DENSITY_FLANK, "left")
        hi = np.searchsorted(arr, row.tss + DENSITY_FLANK, "right")
        dens[i] = (hi - lo) - 1  # exclude the gene itself; = genes per Mb
    return dens


def local_recomb(prom, maps):
    rate = np.full(len(prom), np.nan)
    for i, row in enumerate(prom.itertuples()):
        m = maps.get(row.chrom)
        if m is None:
            continue
        bp, cm = m
        a = np.interp(row.tss - RECOMB_FLANK, bp, cm)
        b = np.interp(row.tss + RECOMB_FLANK, bp, cm)
        rate[i] = (b - a) / (2 * RECOMB_FLANK / 1e6)  # cM per Mb
    return rate


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def strat_test(prom, value_col, nq=4):
    """Alu depletion (HK>NDD) within equal-count strata of value_col.

    Uses rank-based quantile binning so integer-valued metrics with ties still
    split into nq balanced strata (avoids pd.cut duplicate-edge failures)."""
    prom = prom.dropna(subset=[value_col]).copy()
    prom["alu_density"] = prom["fixed_alu"] / ((prom["end"] - prom["start"]) / 1000)
    ranks = prom[value_col].rank(method="first")
    prom["q"] = pd.qcut(ranks, nq, labels=[f"Q{i+1}" for i in range(nq)])
    rows = []
    for lab in [f"Q{i+1}" for i in range(nq)]:
        s = prom[prom["q"] == lab]
        hk = s[s["category"] == "Housekeeping"]["alu_density"]
        nd = s[s["category"] == "HighConfNDD"]["alu_density"]
        lo, hi = s[value_col].min(), s[value_col].max()
        if len(hk) < 15 or len(nd) < 15:
            rows.append(dict(stratum=lab, range=f"{lo:.2f}-{hi:.2f}", n_HK=len(hk),
                             n_NDD=len(nd), p_value=np.nan, r=np.nan))
            continue
        _, p = stats.mannwhitneyu(hk, nd, alternative="greater")
        rows.append(dict(stratum=lab, range=f"{lo:.2f}-{hi:.2f}", n_HK=len(hk), n_NDD=len(nd),
                         median_HK=round(hk.median(), 3), median_NDD=round(nd.median(), 3),
                         p_value=p, r=round(rb(hk, nd), 3)))
    return pd.DataFrame(rows)


def main():
    prom = pd.concat([load_promoters("HighConfNDD"), load_promoters("Housekeeping")],
                     ignore_index=True)
    print("Loading GENCODE gene TSS ..."); tss = gene_tss_from_gencode()
    print("Loading recombination map ..."); maps = recomb_map()
    prom["gene_density"] = local_gene_density(prom, tss)
    prom["recomb_cM_Mb"] = local_recomb(prom, maps)
    prom.to_csv(OUT / "context_per_promoter.tsv", sep="\t", index=False)

    # Q1: group differences
    diff_rows = []
    for col in ["gene_density", "recomb_cM_Mb"]:
        nd = prom[prom["category"] == "HighConfNDD"][col].dropna()
        hk = prom[prom["category"] == "Housekeeping"][col].dropna()
        _, p = stats.mannwhitneyu(nd, hk, alternative="two-sided")
        diff_rows.append(dict(metric=col, n_NDD=len(nd), n_HK=len(hk),
                              median_NDD=round(nd.median(), 3), median_HK=round(hk.median(), 3),
                              mean_NDD=round(nd.mean(), 3), mean_HK=round(hk.mean(), 3),
                              p_value=p, r=round(rb(nd, hk), 3)))
    diff = pd.DataFrame(diff_rows)

    # Q2: Alu depletion within matched strata (rank-based quartiles)
    dstrat = strat_test(prom, "gene_density")
    rstrat = strat_test(prom, "recomb_cM_Mb")

    diff.to_csv(OUT / "context_group_diff.csv", index=False)
    dstrat.to_csv(OUT / "alu_by_density_stratum.csv", index=False)
    rstrat.to_csv(OUT / "alu_by_recomb_stratum.csv", index=False)

    pd.set_option("display.width", 170); pd.set_option("display.max_columns", 20)
    print("\n=== Q1: do NDD vs HK promoters differ in context? (two-sided) ===")
    print(diff.to_string(index=False))
    print("\n=== Q2a: Alu depletion (HK>NDD) within GENE-DENSITY quartiles ===")
    print(dstrat.to_string(index=False))
    print("\n=== Q2b: Alu depletion (HK>NDD) within RECOMBINATION quartiles ===")
    print(rstrat.to_string(index=False))


if __name__ == "__main__":
    main()
