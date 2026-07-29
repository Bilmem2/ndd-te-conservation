#!/usr/bin/env python3
"""
23_brain_overlay.py — Neurodevelopmental tissue-SPECIFICITY of the regulatory
overlay at NDD promoters (hg38), with multi-sample robustness.

ENCODE fetal DNase-seq, three brain donors vs three non-neural control tissues
(all mid-gestation, ~100-122 days):
  BRAIN   : fetal brain 122d (ENCFF955AQD), 117d (ENCFF631TDE), 101d (ENCFF670PXX)
  CONTROL : fetal liver 113d (ENCFF667IEN), lung 120d (ENCFF362PZG),
            stomach 110d (ENCFF016LYI)

Total peak counts differ widely across samples (~48k-454k), so ALL inference
uses the within-sample NDD-vs-housekeeping contrast (internally normalized).
Robustness question: does the NDD/HK enrichment go one way for every brain
donor and the other way for every control tissue?

Reports per-sample NDD/HK ratios, a brain-vs-control specificity test, the
Alu-free vs Alu-positive association (correlational, not causal), and Alu
depletion within brain-DNase strata.

Output: results/brain/{brain_overlay_persample.csv, brain_overlay_summary.csv}
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RES = REPO / "results"
BR = REPO / "data" / "hg38" / "brain"
OUT = RES / "brain"
OUT.mkdir(parents=True, exist_ok=True)

BRAIN = {"brain_122d": "fetal_brain_dnase.bed.gz", "brain_117d": "brain2_117d.bed.gz",
         "brain_101d": "brain3_101d.bed.gz"}
CTRL = {"liver_113d": "fetal_liver_dnase.bed.gz", "lung_120d": "lung_120d.bed.gz",
        "stomach_110d": "stomach_110d.bed.gz"}


def load_promoters(cat):
    df = pd.read_csv(RES / "hg38" / f"{cat}_Alu.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand", "fixed_alu"])
    df["category"] = cat
    return df[["chrom", "start", "end", "gene", "fixed_alu", "category"]]


def load_peaks(fn):
    rows = []
    with gzip.open(BR / fn, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            f = line.split("\t")
            rows.append((f[0], int(f[1]), int(f[2])))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def count_overlaps(prom, feats):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in feats.groupby("chrom"):
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


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    prom = pd.concat([load_promoters("HighConfNDD"), load_promoters("Housekeeping")],
                     ignore_index=True)
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    is_ndd = prom["category"].values == "HighConfNDD"

    rows = []
    dens = {}
    for grp, samples in [("brain", BRAIN), ("control", CTRL)]:
        for name, fn in samples.items():
            pk = load_peaks(fn)
            d = count_overlaps(prom, pk) / prom["width_kb"].values
            dens[name] = d
            nd, hk = d[is_ndd], d[~is_ndd]
            _, p = stats.mannwhitneyu(nd, hk, alternative="two-sided")
            rows.append(dict(sample=name, group=grp, n_peaks=len(pk),
                             mean_NDD=round(nd.mean(), 3), mean_HK=round(hk.mean(), 3),
                             ratio_NDD_HK=round(nd.mean() / hk.mean(), 3) if hk.mean() else np.nan,
                             r=round(rb(nd, hk), 3), p_value=p))
    persample = pd.DataFrame(rows)
    persample.to_csv(OUT / "brain_overlay_persample.csv", index=False)

    # composite brain / control signal (mean density across the samples in each group)
    brain_d = np.mean([dens[n] for n in BRAIN], axis=0)
    ctrl_d = np.mean([dens[n] for n in CTRL], axis=0)
    prom["alu_d"] = prom["fixed_alu"] / prom["width_kb"]

    # specificity: per-promoter brain minus control, NDD vs HK
    bias = brain_d - ctrl_d
    _, p_spec = stats.mannwhitneyu(bias[is_ndd], bias[~is_ndd], alternative="greater")

    # Q3: Alu-free vs Alu-positive, composite brain DNase (association, NOT causal)
    free = brain_d[prom["fixed_alu"].values == 0]
    pos = brain_d[prom["fixed_alu"].values > 0]
    _, p_free = stats.mannwhitneyu(free, pos, alternative="greater")

    # Q4: Alu depletion (HK>NDD) within composite-brain-DNase quartiles
    bq = pd.qcut(pd.Series(brain_d).rank(method="first"), 4,
                 labels=["Q1(low)", "Q2", "Q3", "Q4(high)"])
    q4 = []
    for lab in ["Q1(low)", "Q2", "Q3", "Q4(high)"]:
        s = prom[bq.values == lab]
        hk = s[s["category"] == "Housekeeping"]["alu_d"]; nd = s[s["category"] == "HighConfNDD"]["alu_d"]
        _, p = stats.mannwhitneyu(hk, nd, alternative="greater")
        q4.append(dict(brain_DNase_stratum=lab, n_HK=len(hk), n_NDD=len(nd),
                       median_alu_HK=round(hk.median(), 3), median_alu_NDD=round(nd.median(), 3),
                       p_value=p, r=round(rb(hk, nd), 3)))
    q4 = pd.DataFrame(q4)

    q4.to_csv(OUT / "alu_by_brain_dnase_stratum.csv", index=False)

    rmap = dict(zip(persample["sample"], persample["ratio_NDD_HK"]))
    summ = dict(
        brain_ratios=[rmap[n] for n in BRAIN],
        control_ratios=[rmap[n] for n in CTRL],
        p_NDD_more_brain_biased=round(p_spec, 5),
        p_alu_free_more_brain_DNase=round(p_free, 5),
        r_alu_free_vs_pos=round(rb(free, pos), 3))
    pd.DataFrame([summ]).to_csv(OUT / "brain_overlay_summary.csv", index=False)

    pd.set_option("display.width", 170); pd.set_option("display.max_columns", 20)
    print("=== Per-sample NDD/HK DNase enrichment (within-sample, normalized) ===")
    print(persample.to_string(index=False))
    print(f"\nBrain NDD/HK ratios : {summ['brain_ratios']}  (all > 1 = consistent brain enrichment)")
    print(f"Control NDD/HK ratios: {summ['control_ratios']}  (all <= ~1 = no control enrichment)")
    print(f"\nSpecificity: NDD more brain-biased than HK  p = {summ['p_NDD_more_brain_biased']}")
    print(f"Alu-free promoters have more composite-brain DNase (association)  p = "
          f"{summ['p_alu_free_more_brain_DNase']}, r = {summ['r_alu_free_vs_pos']}")
    print("\n=== Alu depletion (HK>NDD) within composite fetal-brain-DNase quartiles ===")
    print(q4.to_string(index=False))


if __name__ == "__main__":
    main()
