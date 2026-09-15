"""
46_paired_matched_tests.py — re-analyse both matched-control comparisons with
procedures that respect the pairing.

18_matched_control.py and 27_constraint_matched.py each build a 1:1 matched
control by nearest-neighbour matching and then compare the two groups with an
independent-samples Mann-Whitney U test. That discards the matching: it treats
the control promoters as an unrelated sample when each one was chosen to
resemble a specific NDD promoter. The reviewer is right that the pairing should
be preserved, and this script does that for both analyses.

Three paired procedures are reported.

  Wilcoxon signed-rank on the within-pair difference (NDD minus control),
  one-sided against the alternative that NDD is lower. Pairs whose difference
  is exactly zero carry no information about direction and are excluded, which
  is the standard treatment; because Alu density in a 4 kb window is coarse and
  frequently zero, the number of such pairs is large and is reported rather
  than buried.

  A within-pair permutation test, which makes no distributional assumption at
  all: the sign of each observed difference is flipped at random, the mean
  difference is recomputed, and the observed mean is placed in that null. This
  is the exact analogue of the matching design, since under the null the label
  within a matched pair is exchangeable.

  The matched-pairs rank-biserial correlation (Kerby 2014) as the paired effect
  size: the difference between the sums of positive and negative signed ranks,
  divided by their total. It is on the same [-1, 1] scale as the rank-biserial r
  used elsewhere in the paper, and the sign convention is kept so that negative
  means NDD promoters are the depleted member of the pair.

The unpaired statistic is recomputed alongside on the identical pairs, so the
two can be compared directly rather than across different gene sets.

Inputs : results/matched/matched_pairs.tsv             (18_matched_control.py)
         results/matched/constraint_matched_pairs.tsv  (27_constraint_matched.py)
         data/hg38/alu_rmsk.bed
Output : results/matched/paired_tests.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
REPO = Path(__file__).resolve().parents[1]
MATCH = REPO / "results" / "matched"
OUT = MATCH / "paired_tests.csv"

N_PERM = 20_000
RNG = np.random.default_rng(20260915)


def alu_index():
    alu = pd.read_csv(REPO / "data/hg38/alu_rmsk.bed", sep="\t", header=None,
                      usecols=[0, 1, 2], names=["chrom", "start", "end"])
    idx = {}
    for chrom, sub in alu.groupby("chrom"):
        idx[chrom] = (np.sort(sub.start.to_numpy()), np.sort(sub.end.to_numpy()))
    return idx


def density(idx, chrom, start, end):
    """Alu records per kb in [start, end), matching the main analysis."""
    out = np.zeros(len(chrom))
    for i, (c, s, e) in enumerate(zip(chrom, start, end)):
        if c in idx:
            st, en = idx[c]
            out[i] = (np.searchsorted(st, e, side="left")
                      - np.searchsorted(en, s, side="right"))
    return out / ((np.asarray(end) - np.asarray(start)) / 1000.0)


def matched_pairs_rank_biserial(d):
    """Kerby's matched-pairs r from the signed ranks of the differences d.

    Zero differences are dropped first, as in the Wilcoxon test. Returns r with
    the sign of the mean rank direction: negative when d tends to be negative,
    i.e. when the NDD member of the pair is the depleted one.
    """
    nz = d[d != 0]
    if nz.size == 0:
        return np.nan, 0
    ranks = stats.rankdata(np.abs(nz))
    w_pos = ranks[nz > 0].sum()
    w_neg = ranks[nz < 0].sum()
    return (w_pos - w_neg) / (w_pos + w_neg), nz.size


def permutation_p(d, n_perm=N_PERM):
    """One-sided p for mean(d) < 0 under within-pair sign exchangeability.

    Returns (p, at_floor). With no permutation at least as extreme as the
    observed mean the estimate is 1/(n+1), which is the resolution of the test
    rather than a measured value; at_floor marks that case so it is reported as
    a bound.
    """
    obs = d.mean()
    flips = RNG.choice([-1.0, 1.0], size=(n_perm, d.size))
    null = (flips * d).mean(axis=1)
    n_extreme = int(np.sum(null <= obs))
    return (n_extreme + 1) / (n_perm + 1), n_extreme == 0


def analyse(label, ndd, ctrl, n_note=""):
    d = ndd - ctrl
    ties = int((d == 0).sum())

    w = stats.wilcoxon(ndd, ctrl, alternative="less", zero_method="wilcox")
    r_paired, n_eff = matched_pairs_rank_biserial(d)
    p_perm, perm_floor = permutation_p(d)

    # sign test on the discordant pairs only
    n_lower = int((d < 0).sum())
    n_higher = int((d > 0).sum())
    p_sign = stats.binomtest(n_lower, n_lower + n_higher, 0.5,
                             alternative="greater").pvalue

    # the unpaired statistic, recomputed on exactly these pairs
    u = stats.mannwhitneyu(ctrl, ndd, alternative="greater")
    r_unpaired = 1 - (2 * u.statistic) / (len(ctrl) * len(ndd))

    row = dict(
        analysis=label, n_pairs=len(d), tied_pairs=ties, informative_pairs=n_eff,
        mean_NDD=round(ndd.mean(), 4), mean_control=round(ctrl.mean(), 4),
        median_diff=round(float(np.median(d)), 4), mean_diff=round(d.mean(), 4),
        pairs_NDD_lower=n_lower, pairs_NDD_higher=n_higher,
        pct_lower_of_discordant=round(100 * n_lower / max(n_lower + n_higher, 1), 1),
        r_paired=round(r_paired, 4), p_wilcoxon=w.pvalue,
        p_permutation=p_perm, p_permutation_at_floor=perm_floor,
        p_sign=p_sign,
        r_unpaired=round(r_unpaired, 4), p_unpaired=u.pvalue,
    )

    print(f"=== {label} ===")
    if n_note:
        print(f"  {n_note}")
    print(f"  pairs {len(d)}   tied {ties} ({100 * ties / len(d):.1f} %)"
          f"   informative {n_eff}")
    print(f"  mean Alu/kb   NDD {ndd.mean():.4f}   control {ctrl.mean():.4f}"
          f"   mean diff {d.mean():+.4f}")
    print(f"  discordant pairs: NDD lower in {n_lower}, higher in {n_higher}"
          f"  ({row['pct_lower_of_discordant']:.1f} % lower)")
    print(f"  PAIRED    r = {r_paired:+.4f}"
          f"   Wilcoxon p = {w.pvalue:.3e}"
          f"   permutation p {'<' if perm_floor else '='} {p_perm:.1e}"
          f"   sign p = {p_sign:.3e}")
    print(f"  unpaired  r = {r_unpaired:+.4f}   p = {u.pvalue:.3e}"
          f"   (same pairs, pairing ignored)")
    print()
    return row


rows = []

# ── 1. constraint + brain expression + GC matching (script 27) ──────────────
cp = pd.read_csv(MATCH / "constraint_matched_pairs.tsv", sep="\t")
idx = alu_index()
cp["ndd_alu"] = density(idx, cp.ndd_chrom, cp.ndd_start, cp.ndd_end)
cp["ctrl_alu"] = density(idx, cp.control_chrom, cp.control_start, cp.control_end)
rows.append(analyse("constraint + brain expression + GC (script 27)",
                    cp.ndd_alu.values, cp.ctrl_alu.values,
                    "densities recomputed from the exported pair coordinates"))

# ── 2. GC + gene density + recombination matching (script 18) ───────────────
mp = pd.read_csv(MATCH / "matched_pairs.tsv", sep="\t")
rows.append(analyse("GC + gene density + recombination (script 18)",
                    mp.ndd_alu_density.values, mp.hk_alu_density.values,
                    "densities taken from the exported pair table"))

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print("-" * 74)
print("A paired effect size is not on the same footing as the unpaired one: it")
print("summarises the direction of within-pair differences, not the separation of")
print("two distributions, so the two are reported side by side rather than one")
print("being presented as a correction of the other.")
print("-" * 74)
print(f"\nwrote {OUT.relative_to(REPO)}")
