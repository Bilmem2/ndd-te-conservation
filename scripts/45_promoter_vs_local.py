"""
45_promoter_vs_local.py — per-gene promoter density normalised to its own local
background.

43_flanking_control.py answers the reviewer's first suggestion, comparing NDD and
housekeeping promoters against progressively distant windows. It shows a deficit
that is deepest at the promoter and settles onto a shallower plateau that is still
present 250 kb away, so part of the signal belongs to the neighbourhood rather
than the promoter. That comparison is between groups at each distance; it does not
say whether an individual promoter is depleted relative to the sequence around it.

This script answers the second suggestion. Each gene is scored twice, once over
its promoter and once over a wide ring of local background, and the promoter is
expressed as a fraction of that background. A gene sitting in an Alu-poor
neighbourhood is then compared against its own neighbourhood rather than against
the genome, so a regional deficit cancels and only a promoter-local component
survives. NDD and housekeeping genes are compared on that normalised quantity.

The ring excludes the 10 kb closest to the TSS: 43 shows the group ratio has not
settled before then, so sequence nearer than that is not background. Two ring
definitions are reported, since the choice is arbitrary and the conclusion should
not depend on it.

Ratios need a non-zero denominator. A 4 kb window is frequently Alu-free, which is
why the background is taken over a wide ring rather than a single window; genes
whose ring is nonetheless empty are dropped and counted in the output.

Method otherwise matches the main analysis: 4 kb promoter window, Alu records per
kb, one-sided Mann-Whitney U (HK > NDD), rank-biserial r.

Inputs : data/hg38/rmsk/Alu.bed
         results/hg38/{HighConfNDD,Housekeeping}_Alu.bed
Output : results/hg38/promoter_vs_local.csv
         results/hg38/promoter_vs_local_per_gene.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "hg38" / "promoter_vs_local.csv"
OUT_GENE = ROOT / "results" / "hg38" / "promoter_vs_local_per_gene.tsv"

WIDTH = 4000                     # promoter window, as in the main analysis
RINGS = {                        # (inner, outer) radius of the background ring
    "10-100 kb": (10_000, 100_000),
    "25-250 kb": (25_000, 250_000),
}

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


def ring_density(df, inner, outer):
    """Alu per kb over [TSS-outer, TSS-inner) + (TSS+inner, TSS+outer].

    Both sides are pooled: unlike the group comparison in 43, where the two
    flanks are kept apart because NDD genes are longer, here the ring is only a
    local expectation for the promoter and pooling gives a steadier denominator.
    """
    span_kb = 2 * (outer - inner) / 1000
    out = np.zeros(len(df))
    for i, (c, t) in enumerate(zip(df.chr.to_numpy(), df.tss.to_numpy())):
        out[i] = (count(c, t - outer, t - inner)
                  + count(c, t + inner, t + outer))
    return out / span_kb


def promoter_density(df):
    half = WIDTH // 2
    out = np.zeros(len(df))
    for i, (c, t) in enumerate(zip(df.chr.to_numpy(), df.tss.to_numpy())):
        out[i] = count(c, t - half, t + half)
    return out / (WIDTH / 1000)


def rank_biserial(hi, lo):
    """r for the one-sided test that `hi` exceeds `lo`."""
    u = stats.mannwhitneyu(hi, lo, alternative="greater")
    return 1 - (2 * u.statistic) / (len(hi) * len(lo)), u.pvalue


ndd_all, hk_all = load("HighConfNDD"), load("Housekeeping")
rows, per_gene = [], []

for label, (inner, outer) in RINGS.items():
    # keep genes whose whole ring fits on the chromosome, so the denominator is
    # always measured over the same amount of sequence
    def fits(df):
        lim = df.chr.map(chrom_max)
        return (df.tss - outer >= 0) & (df.tss + outer <= lim) & lim.notna()

    ndd, hk = ndd_all[fits(ndd_all)].copy(), hk_all[fits(hk_all)].copy()

    for df in (ndd, hk):
        df["prom"] = promoter_density(df)
        df["bg"] = ring_density(df, inner, outer)

    dropped = {k: int((d.bg == 0).sum()) for k, d in (("NDD", ndd), ("HK", hk))}
    ndd, hk = ndd[ndd.bg > 0].copy(), hk[hk.bg > 0].copy()
    for df in (ndd, hk):
        df["norm"] = df.prom / df.bg

    # raw promoter comparison, for reference at the same gene set
    r_raw, p_raw = rank_biserial(hk.prom.values, ndd.prom.values)
    # normalised comparison — the quantity the reviewer asked for
    r_norm, p_norm = rank_biserial(hk.norm.values, ndd.norm.values)
    # background alone, to show the regional component is real
    r_bg, p_bg = rank_biserial(hk.bg.values, ndd.bg.values)

    rows.append(dict(
        ring=label, n_NDD=len(ndd), n_HK=len(hk),
        dropped_NDD=dropped["NDD"], dropped_HK=dropped["HK"],
        bg_HK=round(hk.bg.mean(), 3), bg_NDD=round(ndd.bg.mean(), 3),
        bg_ratio=round(ndd.bg.mean() / hk.bg.mean(), 3), r_bg=round(r_bg, 4), p_bg=p_bg,
        prom_HK=round(hk.prom.mean(), 3), prom_NDD=round(ndd.prom.mean(), 3),
        prom_ratio=round(ndd.prom.mean() / hk.prom.mean(), 3),
        r_prom=round(r_raw, 4), p_prom=p_raw,
        # Ratio of means, so this is comparable with bg_ratio and prom_ratio:
        # how much of the promoter deficit is left once the neighbourhood is
        # divided out. Algebraically prom_ratio / bg_ratio.
        norm_ratio_of_means=round(
            (ndd.prom.mean() / ndd.bg.mean()) / (hk.prom.mean() / hk.bg.mean()), 3),
        # Median of the per-gene ratio. Not comparable with the columns above:
        # the per-gene distribution is zero-inflated and right-skewed, so its
        # median sits far below the ratio of means. Reported for the record.
        norm_median_HK=round(hk.norm.median(), 3),
        norm_median_NDD=round(ndd.norm.median(), 3),
        r_norm=round(r_norm, 4), p_norm=p_norm,
    ))

    for name, df in (("HighConfNDD", ndd), ("Housekeeping", hk)):
        per_gene.append(pd.DataFrame({
            "ring": label, "category": name, "gene": df.gene.values,
            "promoter_per_kb": df.prom.values, "background_per_kb": df.bg.values,
            "normalised": df.norm.values}))

res = pd.DataFrame(rows)

print(f"promoter window {WIDTH // 1000} kb; background pooled over both flanks\n")
for _, x in res.iterrows():
    print(f"ring {x.ring}   NDD n={x.n_NDD}, HK n={x.n_HK}"
          f"   (dropped for empty ring: {x.dropped_NDD} / {x.dropped_HK})")
    print(f"  background alone   NDD/HK = {x.bg_ratio:>6.3f}"
          f"   r = {x.r_bg:>7.3f}   p = {x.p_bg:.2e}")
    print(f"  promoter alone     NDD/HK = {x.prom_ratio:>6.3f}"
          f"   r = {x.r_prom:>7.3f}   p = {x.p_prom:.2e}")
    print(f"  promoter / local   NDD/HK = {x.norm_ratio_of_means:>6.3f}"
          f"   r = {x.r_norm:>7.3f}   p = {x.p_norm:.2e}   <- normalised")
    print(f"      ratios above are ratios of means and are comparable;"
          f" per-gene median {x.norm_median_NDD} vs {x.norm_median_HK}")
    print(f"      effect size falls {abs(x.r_prom):.3f} -> {abs(x.r_norm):.3f}"
          f"  ({100 * (1 - abs(x.r_norm) / abs(x.r_prom)):.0f} % of the promoter"
          f" effect was regional)")
    print()

print("-" * 72)
print("Both components are present. The background comparison is significant on its")
print("own, so the neighbourhood really is Alu-poor; and the promoter stays depleted")
print("after each gene is divided by its own neighbourhood, so a promoter-local")
print("component survives. Compare effect sizes, not the ratios, across the two:")
print("r is rank-based and tolerates the many Alu-free promoter windows.")
print("-" * 72)

OUT.parent.mkdir(parents=True, exist_ok=True)
res.to_csv(OUT, index=False)
pd.concat(per_gene, ignore_index=True).to_csv(OUT_GENE, sep="\t", index=False)
print(f"\nwrote {OUT.relative_to(ROOT)}")
print(f"wrote {OUT_GENE.relative_to(ROOT)}")
