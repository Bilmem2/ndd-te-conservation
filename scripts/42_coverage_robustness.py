"""
42_coverage_robustness.py — Alu depletion measured as base-pair coverage.

RepeatMasker splits an insertion interrupted by a later one into several records,
so a density expressed as records per kb counts a fragmented element more than
once. This script repeats the human Alu comparison with that dependence removed:
overlapping and adjacent Alu intervals are merged, and each promoter is scored by
the fraction of its window covered by Alu sequence rather than by how many records
fall inside it. A coverage measure is invariant to how an element is split.

Method otherwise matches the main analysis: TSS +/-2 kb, HighConfNDD versus
Housekeeping, one-sided Mann-Whitney U (HK > NDD), rank-biserial r. Promoter
windows and their record counts are taken from the committed human BEDs, so the
count-based row reproduces the published figure alongside the coverage-based one.

Inputs : data/hg38/rmsk/Alu.bed
         results/hg38/{HighConfNDD,Housekeeping}_Alu.bed
Output : results/hg38/coverage_robustness.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "hg38" / "coverage_robustness.csv"


def merge_by_chrom(bed):
    """Collapse overlapping and book-ended intervals into disjoint runs."""
    out = {}
    for chrom, sub in bed.groupby("chr"):
        iv = sub[["start", "end"]].to_numpy()
        iv = iv[np.argsort(iv[:, 0])]
        s, e = [iv[0, 0]], [iv[0, 1]]
        for a, b in iv[1:]:
            if a <= e[-1]:
                if b > e[-1]:
                    e[-1] = b
            else:
                s.append(a)
                e.append(b)
        out[chrom] = (np.asarray(s), np.asarray(e))
    return out


def covered_bp(df, merged):
    """Base pairs of each window overlapped by the merged intervals."""
    cov = np.zeros(len(df))
    pos = {k: i for i, k in enumerate(df.index)}
    for chrom, idx in df.groupby("chr").groups.items():
        if chrom not in merged:
            continue
        ms, me = merged[chrom]
        sub = df.loc[idx]
        for k, qs, qe in zip(idx, sub.start.to_numpy(), sub.end.to_numpy()):
            i = np.searchsorted(me, qs, side="right")
            j = np.searchsorted(ms, qe, side="left")
            if j > i:
                cov[pos[k]] = np.sum(np.minimum(me[i:j], qe)
                                     - np.maximum(ms[i:j], qs))
    return cov


def effect(hk_vals, ndd_vals):
    """One-sided Mann-Whitney U (HK > NDD) with rank-biserial r."""
    u = stats.mannwhitneyu(hk_vals, ndd_vals, alternative="greater")
    return 1 - (2 * u.statistic) / (len(hk_vals) * len(ndd_vals)), u.pvalue


alu = pd.read_csv(ROOT / "data/hg38/rmsk/Alu.bed", sep="\t", header=None,
                  usecols=[0, 1, 2], names=["chr", "start", "end"])
widths = (alu.end - alu.start).to_numpy()
merged = merge_by_chrom(alu)
n_merged = sum(len(v[0]) for v in merged.values())

print(f"Alu records          : {len(alu):,}")
print(f"  median length      : {np.median(widths):.0f} bp")
print(f"  fragments < 100 bp : {(widths < 100).mean() * 100:.1f} %")
print(f"  after merging      : {n_merged:,} "
      f"({100 * (1 - n_merged / len(alu)):.1f} % collapsed)")

sets = {}
for name in ("HighConfNDD", "Housekeeping"):
    df = pd.read_csv(ROOT / f"results/hg38/{name}_Alu.bed", sep="\t", header=None,
                     names=["chr", "start", "end", "gene", "score", "strand", "count"])
    df["width"] = df.end - df.start
    df["density_per_kb"] = df["count"] / (df.width / 1000)
    df["coverage_fraction"] = covered_bp(df, merged) / df.width
    sets[name] = df

ndd, hk = sets["HighConfNDD"], sets["Housekeeping"]
rows = []
for label, col in (("records per kb", "density_per_kb"),
                   ("merged bp coverage", "coverage_fraction")):
    r, p = effect(hk[col], ndd[col])
    rows.append(dict(measure=label, n_HK=len(hk), n_NDD=len(ndd),
                     median_HK=round(hk[col].median(), 4),
                     median_NDD=round(ndd[col].median(), 4),
                     r=round(r, 4), p=p))

res = pd.DataFrame(rows)
print()
print(res.to_string(index=False))
print(f"\n|r| difference: {abs(res.r.iloc[1]) - abs(res.r.iloc[0]):+.4f}")

OUT.parent.mkdir(parents=True, exist_ok=True)
res.to_csv(OUT, index=False)
print(f"\nwrote {OUT.relative_to(ROOT)}")
