#!/usr/bin/env python3
"""
probe_subfamily.py — FEASIBILITY PROBE: is there an age-dependent (purifying-
selection-consistent) signature in the Alu depletion at NDD promoters?

Within the comparative/genomic identity of the paper (no population data), we
split reference-genome Alu by subfamily AGE:
  AluJ (oldest, ~65 My), AluS (intermediate), AluY (youngest, still active).
milliDiv (divergence from consensus) is the per-element age proxy.

Two tests (NDD vs housekeeping promoters, TSS +/-2 kb, hg38):
  1. Per-subfamily depletion: is depletion (HK>NDD) stronger for YOUNGER
     subfamilies?  If the constraint is ongoing, recently-active AluY should be
     more depleted at NDD promoters than ancient AluJ.
  2. Age of survivors: are Alu that DID survive inside NDD promoters older
     (higher milliDiv) than those in HK promoters?  Consistent with recent
     insertions being purged at NDD loci.

Caveat (reported honestly): reference genome shows survivors only, so this is
selection-CONSISTENT, not a clean targeting-vs-selection separation.

Inputs: results/hg38/{HighConfNDD,Housekeeping}_Alu.bed  (promoter windows),
        data/hg38/alu_detailed.bed  (chrom,start,end,repName,milliDiv)
Output: results/functional/probe_subfamily.csv + console
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RES = REPO / "results"
DATA = REPO / "data" / "hg38"


def promoters(cat):
    df = pd.read_csv(RES / "hg38" / f"{cat}_Alu.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "sc", "strand", "c"])
    df["category"] = cat
    return df[["chrom", "start", "end", "gene", "category"]]


def subfam(name):
    for k in ("AluY", "AluS", "AluJ"):
        if name.startswith(k):
            return k
    return "other"


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    prom = pd.concat([promoters("HighConfNDD"), promoters("Housekeeping")], ignore_index=True)
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    is_ndd = prom["category"].values == "HighConfNDD"

    alu = pd.read_csv(DATA / "alu_detailed.bed", sep="\t", header=None,
                      names=["chrom", "start", "end", "repName", "milliDiv"])
    alu["fam"] = alu["repName"].map(subfam)
    alu = alu[alu["fam"] != "other"]
    print("Subfamily counts / mean age (milliDiv, higher=older):")
    for k in ("AluY", "AluS", "AluJ"):
        s = alu[alu["fam"] == k]
        print(f"  {k}: n={len(s)}  mean milliDiv={s['milliDiv'].mean():.1f}")

    # Test 1: per-subfamily depletion NDD vs HK, + collect divergences of overlaps
    rows = []
    div_by_cat = {"HighConfNDD": [], "Housekeeping": []}
    for fam in ("AluY", "AluS", "AluJ"):
        sub = alu[alu["fam"] == fam]
        counts = np.zeros(len(prom), dtype=int)
        for chrom, m in sub.groupby("chrom"):
            idx = np.where(prom["chrom"].values == chrom)[0]
            if not len(idx):
                continue
            o = np.argsort(m["start"].values)
            s = m["start"].values[o]; e = m["end"].values[o]; dv = m["milliDiv"].values[o]
            for i in idx:
                ws, we = prom.at[i, "start"], prom.at[i, "end"]
                lo = np.searchsorted(s, we, "left")
                if not lo:
                    continue
                hit = np.where(e[:lo] > ws)[0]
                counts[i] = len(hit)
                for h in hit:
                    div_by_cat[prom.at[i, "category"]].append(dv[h])
        dens = counts / prom["width_kb"].values
        nd, hk = dens[is_ndd], dens[~is_ndd]
        _, p = stats.mannwhitneyu(hk, nd, alternative="greater")
        rows.append(dict(subfamily=fam, mean_age_milliDiv=round(sub["milliDiv"].mean(), 1),
                         dens_HK=round(hk.mean(), 3), dens_NDD=round(nd.mean(), 3),
                         ratio_NDD_HK=round(nd.mean() / hk.mean(), 3) if hk.mean() else np.nan,
                         p_depletion=p, r=round(rb(hk, nd), 3)))
    res = pd.DataFrame(rows)

    # Test 2: age of surviving Alu inside NDD vs HK promoters (all subfamilies pooled)
    dn = np.array([d for c in ("HighConfNDD",) for d in div_by_cat[c]])
    dh = np.array(div_by_cat["Housekeeping"])
    _, p_age = stats.mannwhitneyu(dn, dh, alternative="greater")  # NDD older?

    res.to_csv(RES / "functional" / "probe_subfamily.csv", index=False)
    pd.set_option("display.width", 160)
    print("\n=== Test 1: per-subfamily Alu depletion at NDD promoters (HK>NDD) ===")
    print(res.to_string(index=False))
    print("\n=== Test 2: age of surviving promoter Alu (milliDiv; higher=older) ===")
    print(f"  NDD promoters: median milliDiv={np.median(dn):.0f} (n={len(dn)})")
    print(f"  HK  promoters: median milliDiv={np.median(dh):.0f} (n={len(dh)})")
    print(f"  MWU NDD>HK (survivors older at NDD?): p={p_age:.3g}")
    print("\nREAD: selection-consistent IF younger AluY depleted >= AluJ (|r| larger for AluY)")
    print("      and/or surviving NDD-promoter Alu are older (higher milliDiv).")


if __name__ == "__main__":
    main()
