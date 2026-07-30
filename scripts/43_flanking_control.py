"""
43_flanking_control.py — is the Alu deficit local to the promoter or regional?

The effect grows with window size (r = -0.109 at TSS +/-500 bp, -0.381 at +/-3 kb),
which is the opposite of what a constraint acting on promoter architecture would
predict. That gradient is compatible with the deficit belonging to the wider
genomic neighbourhood of NDD loci rather than to their promoters.

This script separates the two. For each gene it scores the promoter window and a
series of equally wide windows placed at increasing distance from the TSS on both
sides, and repeats the NDD versus Housekeeping comparison at each distance. If the
deficit is promoter-borne it should decay as the windows move away; if it is
regional it should persist.

Gene sets are held constant across distances: a gene is used only when every one
of its windows, at every offset, falls inside the chromosome, so the curve is not
confounded by which genes drop out at the far offsets.

Method otherwise matches the main analysis: 4 kb windows, Alu records per kb,
one-sided Mann-Whitney U (HK > NDD), rank-biserial r. Offset 0 reproduces the
published promoter figure.

Inputs : data/hg38/rmsk/Alu.bed
         results/hg38/{HighConfNDD,Housekeeping}_Alu.bed
Output : results/hg38/flanking_control.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "hg38" / "flanking_control.csv"

WIDTH = 4000                                     # same width as the promoter window
OFFSETS = [0, 10_000, 25_000, 50_000, 100_000, 250_000]

alu = pd.read_csv(ROOT / "data/hg38/rmsk/Alu.bed", sep="\t", header=None,
                  usecols=[0, 1, 2], names=["chr", "start", "end"])

idx = {}
for chrom, sub in alu.groupby("chr"):
    idx[chrom] = (np.sort(sub.start.to_numpy()), np.sort(sub.end.to_numpy()))
chrom_max = {c: int(e[-1]) for c, (s, e) in idx.items()}


def count(chrom, qs, qe):
    """Records overlapping [qs, qe). Matches `bedtools intersect -c`."""
    if chrom not in idx:
        return 0
    starts, ends = idx[chrom]
    return (np.searchsorted(starts, qe, side="left")
            - np.searchsorted(ends, qs, side="right"))


def load(name):
    df = pd.read_csv(ROOT / f"results/hg38/{name}_Alu.bed", sep="\t", header=None,
                     names=["chr", "start", "end", "gene", "score", "strand", "n"])
    df["tss"] = (df.start + df.end) // 2
    return df


ndd, hk = load("HighConfNDD"), load("Housekeeping")
half, far = WIDTH // 2, max(OFFSETS)


def in_bounds(df):
    """Keep genes whose windows all fit, so the gene set is fixed across offsets."""
    lim = df.chr.map(chrom_max)
    return (df.tss - far - half >= 0) & (df.tss + far + half <= lim) & lim.notna()


ndd, hk = ndd[in_bounds(ndd)].copy(), hk[in_bounds(hk)].copy()
print(f"windows {WIDTH // 1000} kb wide, offsets up to {far // 1000} kb")
print(f"genes retained: NDD {len(ndd)}, HK {len(hk)}\n")


def density(df, off, side):
    """Alu per kb in one window at `off` from the TSS. side -1 is 5', +1 is 3'.

    Flanks are kept separate and oriented by transcription: NDD genes are longer,
    so a 3' window falls inside the locus more often for them than for the shorter
    housekeeping genes, and averaging the two sides would fold that difference in.
    """
    out = np.zeros(len(df))
    for i, (c, t, st) in enumerate(zip(df.chr.to_numpy(), df.tss.to_numpy(),
                                       df.strand.to_numpy())):
        if off == 0:
            out[i] = count(c, t - half, t + half)
        else:
            d = off * side * (1 if st == "+" else -1)
            out[i] = count(c, t + d - half, t + d + half)
    return out / (WIDTH / 1000)


rows = []
for off in OFFSETS:
    for side, label in ((-1, "upstream"), (1, "downstream")):
        if off == 0 and side == 1:
            continue
        d_ndd, d_hk = density(ndd, off, side), density(hk, off, side)
        u = stats.mannwhitneyu(d_hk, d_ndd, alternative="greater")
        r = 1 - (2 * u.statistic) / (len(d_hk) * len(d_ndd))
        rows.append(dict(offset_kb=off // 1000,
                         side="promoter" if off == 0 else label,
                         mean_HK=round(d_hk.mean(), 3),
                         mean_NDD=round(d_ndd.mean(), 3),
                         ratio=round(d_ndd.mean() / d_hk.mean(), 3),
                         r=round(r, 4), p=u.pvalue))

res = pd.DataFrame(rows)
print(f"{'offset':>8}  {'side':<11}{'HK/kb':>9}{'NDD/kb':>9}{'NDD/HK':>9}"
      f"{'r':>9}{'p':>12}")
print("-" * 68)
for _, x in res.iterrows():
    print(f"{int(x.offset_kb):>6} kb  {x.side:<11}{x.mean_HK:>9.3f}{x.mean_NDD:>9.3f}"
          f"{x.ratio:>9.3f}{x.r:>9.3f}{x.p:>12.2e}")
print("-" * 68)

prom = res[res.side == "promoter"].iloc[0]
up = res[res.side == "upstream"]
print(f"\npromoter           : r = {prom.r:.3f}, NDD/HK = {prom.ratio:.3f}"
      f"  ({100 * (1 - prom.ratio):.0f} % depleted)")
print(f"upstream 10-250 kb : r = {up.r.min():.3f} to {up.r.max():.3f}, "
      f"NDD/HK = {up.ratio.min():.3f} to {up.ratio.max():.3f}"
      f"  ({100 * (1 - up.ratio.max()):.0f}-{100 * (1 - up.ratio.min()):.0f} % depleted)")

OUT.parent.mkdir(parents=True, exist_ok=True)
res.to_csv(OUT, index=False)
print(f"\nwrote {OUT.relative_to(ROOT)}")
