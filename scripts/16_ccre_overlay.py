#!/usr/bin/env python3
"""
16_ccre_overlay.py — Functional/epigenomic overlay (hg38) answering the
reviewer's "does not integrate functional genomic or epigenomic datasets to
establish the regulatory or biological consequences of the observed depletion."

Uses the ENCODE SCREEN Registry of candidate cis-Regulatory Elements (cCREs,
V4, GRCh38; 2,348,854 elements) with classes:
  PLS (promoter-like), pELS/dELS (enhancer-like), CA-CTCF, CA-H3K4me3, CA, TF, CA-TF.

Three questions:
  1. Are NDD promoters (which are Alu-depleted) themselves ENRICHED for active
     cis-regulatory elements vs housekeeping promoters?  i.e. the depleted
     windows are not inert — they are regulatory-dense.  -> per-promoter cCRE
     density, two-sided MWU, rank-biserial r.
  2. Is Alu ABSENCE linked to regulatory activity?  Compare cCRE density in
     Alu-free vs Alu-containing promoters (one-sided MWU, Alu-free > Alu-pos).
  3. Does Alu depletion (HK>NDD) hold — or strengthen — in the most
     regulatory-active promoters?  -> Alu depletion within cCRE-density strata.

Inputs:
  results/hg38/{HighConfNDD,Housekeeping}_Alu.bed
  data/hg38/GRCh38-cCREs.bed   (chrom,start,end,rDHS,cCRE,class)
Outputs:
  results/ccre/{ccre_per_promoter.tsv, ccre_group_diff.csv,
                alu_free_vs_pos.csv, alu_by_ccre_stratum.csv}
No bedtools dependency (numpy interval overlap).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
DATA = REPO / "data"
OUT = RESULTS / "ccre"
OUT.mkdir(parents=True, exist_ok=True)

# functional groupings of cCRE classes
GROUPS = {
    "PLS": {"PLS"},
    "ELS": {"pELS", "dELS"},
    "CTCF": {"CA-CTCF"},
    "active_all": {"PLS", "pELS", "dELS", "CA-H3K4me3", "CA-TF", "TF"},  # DNase+mark/TF-supported
}


def load_promoters(cat):
    f = RESULTS / "hg38" / f"{cat}_Alu.bed"
    df = pd.read_csv(f, sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand", "fixed_alu"])
    df["category"] = cat
    return df[["chrom", "start", "end", "gene", "fixed_alu", "category"]]


def load_ccre():
    df = pd.read_csv(DATA / "hg38" / "GRCh38-cCREs.bed", sep="\t", header=None,
                     usecols=[0, 1, 2, 5], names=["chrom", "start", "end", "cls"])
    return df


def count_overlaps(prom, feats):
    """count feats overlapping each promoter window (numpy, per chrom)."""
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in feats.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if len(idx) == 0:
            continue
        fs = np.sort(m["start"].values)
        fe = m["end"].values[np.argsort(m["start"].values)]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(fs, we, "left")
            if lo == 0:
                continue
            counts[i] = int(np.count_nonzero(fe[:lo] > ws))
    return counts


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    prom = pd.concat([load_promoters("HighConfNDD"), load_promoters("Housekeeping")],
                     ignore_index=True)
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    ccre = load_ccre()
    is_ndd = prom["category"].values == "HighConfNDD"

    for g, classes in GROUPS.items():
        prom[f"{g}_n"] = count_overlaps(prom, ccre[ccre["cls"].isin(classes)])
    prom["alu_free"] = prom["fixed_alu"] == 0
    prom.to_csv(OUT / "ccre_per_promoter.tsv", sep="\t", index=False)

    # Q1: NDD vs HK cCRE density
    q1 = []
    for g in GROUPS:
        dens = prom[f"{g}_n"] / prom["width_kb"]
        nd, hk = dens[is_ndd], dens[~is_ndd]
        _, p = stats.mannwhitneyu(nd, hk, alternative="two-sided")
        q1.append(dict(cCRE_group=g, median_NDD=round(nd.median(), 3), median_HK=round(hk.median(), 3),
                       mean_NDD=round(nd.mean(), 3), mean_HK=round(hk.mean(), 3),
                       p_value=p, r=round(rb(nd, hk), 3),
                       direction="NDD>HK" if nd.mean() > hk.mean() else "NDD<HK"))
    q1 = pd.DataFrame(q1)

    # Q2: Alu-free vs Alu-positive promoters — regulatory activity (one-sided free>pos)
    q2 = []
    for g in GROUPS:
        free = prom[prom["alu_free"]][f"{g}_n"]
        pos = prom[~prom["alu_free"]][f"{g}_n"]
        _, p = stats.mannwhitneyu(free, pos, alternative="greater")
        q2.append(dict(cCRE_group=g, n_alu_free=len(free), n_alu_pos=len(pos),
                       median_free=round(free.median(), 2), median_pos=round(pos.median(), 2),
                       p_free_gt_pos=p, r=round(rb(free, pos), 3)))
    q2 = pd.DataFrame(q2)

    # Q3: Alu depletion (HK>NDD) within active-cCRE-density quartiles (rank-based)
    prom["alu_density"] = prom["fixed_alu"] / prom["width_kb"]
    active_dens = prom["active_all_n"] / prom["width_kb"]
    prom["ccre_stratum"] = pd.qcut(active_dens.rank(method="first"), 4,
                                   labels=["Q1(low)", "Q2", "Q3", "Q4(high)"])
    q3 = []
    for lab in ["Q1(low)", "Q2", "Q3", "Q4(high)"]:
        s = prom[prom["ccre_stratum"] == lab]
        hk = s[s["category"] == "Housekeeping"]["alu_density"]
        nd = s[s["category"] == "HighConfNDD"]["alu_density"]
        adl, adh = active_dens[prom["ccre_stratum"] == lab].min(), active_dens[prom["ccre_stratum"] == lab].max()
        _, p = stats.mannwhitneyu(hk, nd, alternative="greater")
        q3.append(dict(active_cCRE_stratum=lab, active_range=f"{adl:.2f}-{adh:.2f}",
                       n_HK=len(hk), n_NDD=len(nd),
                       median_HK=round(hk.median(), 3), median_NDD=round(nd.median(), 3),
                       p_value=p, r=round(rb(hk, nd), 3)))
    q3 = pd.DataFrame(q3)

    q1.to_csv(OUT / "ccre_group_diff.csv", index=False)
    q2.to_csv(OUT / "alu_free_vs_pos.csv", index=False)
    q3.to_csv(OUT / "alu_by_ccre_stratum.csv", index=False)

    pd.set_option("display.width", 170); pd.set_option("display.max_columns", 20)
    print("\n=== Q1: cCRE density at NDD vs HK promoters (two-sided) ===")
    print(q1.to_string(index=False))
    print("\n=== Q2: regulatory activity in Alu-free vs Alu-positive promoters (free>pos) ===")
    print(q2.to_string(index=False))
    print("\n=== Q3: Alu depletion (HK>NDD) within active-cCRE-density quartiles ===")
    print(q3.to_string(index=False))


if __name__ == "__main__":
    main()
