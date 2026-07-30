"""
41_dog_cansine.py — Can-SINE at dog NDD promoters.

Can-SINEs are the dominant SINE family of the dog lineage and are tRNA-derived
with AT-biased target-site preferences, unlike the GC-biased, L1-mobilised Alu
and B1/B2 families that carry the SINE-class comparisons elsewhere in this study.
They were therefore excluded from those comparisons. This script measures them
anyway, so that the boundary of the pattern is reported rather than assumed:
if Can-SINEs are not depleted, the constraint applies to GC-biased SINEs
specifically rather than to SINEs as a class.

Method matches the main analysis: Can-SINE count per kb at TSS +/-2 kb,
HighConfNDD versus Housekeeping, one-sided Mann-Whitney U (HK > NDD),
rank-biserial r. Promoter windows are taken from the committed dog BEDs.

Inputs : data/canFam4/rmsk/rmsk.txt.gz
         results/canFam4/{HighConfNDD_LINE1.bed, Housekeeping_LINE1.bed}
Output : results/canFam4/{HighConfNDD_CanSINE.bed, Housekeeping_CanSINE.bed,
                          dog_cansine_stats.csv}
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "canFam4"
COLS = ["chrom", "start", "end", "gene", "score", "strand", "n"]


def load_cansine():
    rows = []
    with gzip.open(ROOT / "data" / "canFam4" / "rmsk" / "rmsk.txt.gz", "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 12 and f[11] == "SINE" and f[12] == "tRNA":
                rows.append((f[5], int(f[6]), int(f[7])))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def load_promoters(cat):
    df = pd.read_csv(RES / f"{cat}_LINE1.bed", sep="\t", header=None, names=COLS,
                     lineterminator="\n")
    for c in ("start", "end"):
        df[c] = pd.to_numeric(df[c].astype(str).str.strip(), errors="coerce")
    return df.dropna(subset=["start", "end"]).reset_index(drop=True)[
        ["chrom", "start", "end", "gene", "score", "strand"]]


def count(prom, feats):
    out = np.zeros(len(prom), dtype=int)
    for chrom, m in feats.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        o = np.argsort(m["start"].values)
        s, e = m["start"].values[o], m["end"].values[o]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                out[i] = int(np.count_nonzero(e[:lo] > ws))
    return out


cs = load_cansine()
print(f"Can-SINE (SINE/tRNA) elements: {len(cs):,} on {cs['chrom'].nunique()} sequences")

frames = {}
for cat in ("HighConfNDD", "Housekeeping"):
    p = load_promoters(cat)
    p["n"] = count(p, cs)
    p["density"] = p["n"] / ((p["end"] - p["start"]) / 1000.0)
    p[["chrom", "start", "end", "gene", "score", "strand", "n"]].to_csv(
        RES / f"{cat}_CanSINE.bed", sep="\t", header=False, index=False)
    frames[cat] = p
    print(f"  {cat}: {len(p)} promoters, mean {p['density'].mean():.3f} per kb")

nd, hk = frames["HighConfNDD"]["density"], frames["Housekeeping"]["density"]
U, p_one = stats.mannwhitneyu(hk, nd, alternative="greater")
u2, _ = stats.mannwhitneyu(hk, nd, alternative="two-sided")
r = 1 - (2 * u2) / (len(hk) * len(nd))
sig = "***" if p_one < 1e-3 else "**" if p_one < 1e-2 else "*" if p_one < 0.05 else "ns"

row = dict(species="Dog", assembly="ROS_Cfam_1.0", mya=95, TE="Can-SINE",
           n_HK=len(hk), n_NDD=len(nd),
           median_HK=round(hk.median(), 3), median_NDD=round(nd.median(), 3),
           mean_HK=round(hk.mean(), 3), mean_NDD=round(nd.mean(), 3),
           p_value=p_one, r=round(r, 3), sig=sig)
pd.DataFrame([row]).to_csv(RES / "dog_cansine_stats.csv", index=False)

print("\n=== Dog Can-SINE at NDD promoters ===")
for k, v in row.items():
    print(f"  {k}: {v}")
print("\nContext - GC-biased SINE effects elsewhere in the panel:")
print("  mouse B1/B2 -0.413 | macaque Alu -0.384 | human Alu -0.345 | marmoset Alu -0.147")
