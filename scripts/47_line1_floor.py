"""
47_line1_floor.py — is the LINE-1 null a real null, or too little signal to see?

The argument that the deficit is specific to a class of elements rests on LINE-1
showing no comparable depletion. That argument is only worth making if a deficit
the size of the Alu one could have been detected in the LINE-1 data at all. The
manuscript already makes this check for dog Can-SINEs, where it is settled by
noting that Can-SINE density at promoters is comparable to Alu density, so there
was room for a deficit to show. LINE-1 is sparser than Alu at promoters, so the
same argument cannot simply be repeated and a stronger test is needed.

Two things are measured.

First the observed quantities: LINE-1 frequency at NDD and housekeeping
promoters and across all protein-coding promoters, with the depletion test, so
the null is stated with its effect size rather than only its p value.

Then detectability, by thinning. Each Alu record is retained with the
probability that brings mean Alu frequency down to the observed mean LINE-1
frequency, which preserves the between-promoter structure of the Alu data while
reducing its density to LINE-1 levels. The NDD versus housekeeping comparison is
repeated on the thinned data. If the Alu deficit is still detected after
thinning, sparsity is not why LINE-1 shows no deficit, and the class-selectivity
argument stands. If it is not, the LINE-1 null is uninformative and the argument
has to be made another way.

Method matches the main analysis: TSS +/-2 kb, records per kb, one-sided
Mann-Whitney U (HK > NDD), rank-biserial r.

Inputs : results/hg38/{HighConfNDD,Housekeeping}_{Alu,LINE1}.bed
         data/hg38/promoters/promoters_all.bed, data/hg38/rmsk/LINE1.bed
Output : results/hg38/line1_floor.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "hg38" / "line1_floor.csv"

N_THIN = 400
RNG = np.random.default_rng(20260915)
COLS = ["chr", "start", "end", "gene", "score", "strand", "n"]


def load(cat, te):
    df = pd.read_csv(REPO / f"results/hg38/{cat}_{te}.bed", sep="\t",
                     header=None, names=COLS)
    df["kb"] = (df.end - df.start) / 1000.0
    df["density"] = df.n / df.kb
    return df


def test(hk, ndd):
    """One-sided HK > NDD, with rank-biserial r. Matches the main analysis."""
    u = stats.mannwhitneyu(hk, ndd, alternative="greater")
    r = 1 - (2 * u.statistic) / (len(hk) * len(ndd))
    return r, u.pvalue


# ── observed ────────────────────────────────────────────────────────────────
obs = {}
for te in ("Alu", "LINE1"):
    nd, hk = load("HighConfNDD", te), load("Housekeeping", te)
    r, p = test(hk.density.values, nd.density.values)
    obs[te] = dict(te=te, n_NDD=len(nd), n_HK=len(hk),
                   mean_NDD=nd.density.mean(), mean_HK=hk.density.mean(),
                   ratio=nd.density.mean() / hk.density.mean(), r=r, p=p,
                   nd=nd, hk=hk)

# genome-wide promoter baseline for each element, for context
prom = pd.read_csv(REPO / "data/hg38/promoters/promoters_all.bed", sep="\t",
                   header=None, names=["chr", "start", "end", "gene", "score", "strand"])
baseline = {}
for te, bed in (("Alu", "Alu.bed"), ("LINE1", "LINE1.bed")):
    rec = pd.read_csv(REPO / f"data/hg38/rmsk/{bed}", sep="\t", header=None,
                      usecols=[0, 1, 2], names=["chr", "start", "end"])
    idx = {c: (np.sort(s.start.to_numpy()), np.sort(s.end.to_numpy()))
           for c, s in rec.groupby("chr")}
    cnt = np.zeros(len(prom))
    for i, (c, s, e) in enumerate(zip(prom.chr.to_numpy(), prom.start.to_numpy(),
                                      prom.end.to_numpy())):
        if c in idx:
            st, en = idx[c]
            cnt[i] = (np.searchsorted(st, e, side="left")
                      - np.searchsorted(en, s, side="right"))
    baseline[te] = (cnt / ((prom.end - prom.start) / 1000.0)).mean()

print("=== observed frequency at TSS +/- 2 kb (records per kb) ===")
print(f"{'element':<8}{'NDD':>9}{'HK':>9}{'NDD/HK':>9}{'genome':>9}"
      f"{'r':>9}{'p':>12}")
print("-" * 65)
for te in ("Alu", "LINE1"):
    o = obs[te]
    print(f"{te:<8}{o['mean_NDD']:>9.3f}{o['mean_HK']:>9.3f}{o['ratio']:>9.3f}"
          f"{baseline[te]:>9.3f}{o['r']:>9.3f}{o['p']:>12.2e}")
print("-" * 65)

ratio_density = obs["LINE1"]["mean_HK"] / obs["Alu"]["mean_HK"]
print(f"\nLINE-1 is {1 / ratio_density:.1f}x sparser than Alu at housekeeping"
      f" promoters, so the density-comparability argument used for Can-SINE"
      f"\ndoes not transfer and detectability has to be tested directly.\n")

# ── thinning ────────────────────────────────────────────────────────────────
# retain each Alu record with probability p so that mean Alu frequency falls to
# the observed mean LINE-1 frequency
pooled_alu = np.concatenate([obs["Alu"]["nd"].n.values, obs["Alu"]["hk"].n.values])
pooled_kb = np.concatenate([obs["Alu"]["nd"].kb.values, obs["Alu"]["hk"].kb.values])
pooled_l1 = np.concatenate([obs["LINE1"]["nd"].n.values, obs["LINE1"]["hk"].n.values])
p_keep = pooled_l1.sum() / pooled_alu.sum()

nd_n, nd_kb = obs["Alu"]["nd"].n.values, obs["Alu"]["nd"].kb.values
hk_n, hk_kb = obs["Alu"]["hk"].n.values, obs["Alu"]["hk"].kb.values

rs, ps = np.empty(N_THIN), np.empty(N_THIN)
for i in range(N_THIN):
    nd_t = RNG.binomial(nd_n, p_keep) / nd_kb
    hk_t = RNG.binomial(hk_n, p_keep) / hk_kb
    rs[i], ps[i] = test(hk_t, nd_t)

sig = float((ps < 0.05).mean())
print(f"=== Alu thinned to LINE-1 density ({N_THIN} draws, keep p = {p_keep:.3f}) ===")
print(f"  mean thinned Alu frequency: "
      f"{(RNG.binomial(hk_n, p_keep) / hk_kb).mean():.3f} per kb"
      f"   (observed LINE-1 at HK: {obs['LINE1']['mean_HK']:.3f})")
print(f"  effect size r : median {np.median(rs):.3f}"
      f"   IQR [{np.percentile(rs, 25):.3f}, {np.percentile(rs, 75):.3f}]")
print(f"  p < 0.05 in   : {100 * sig:.1f} % of draws")
print(f"  observed LINE-1 r = {obs['LINE1']['r']:.3f}, and the thinned Alu r")
print(f"  falls below it in {100 * float((rs < obs['LINE1']['r']).mean()):.1f} % of draws")
print()

verdict = ("detectable: the Alu deficit survives thinning to LINE-1 density, so"
           " sparsity\n  does not explain the LINE-1 null and the"
           " class-selectivity argument stands"
           if sig >= 0.8 else
           "NOT detectable at this density: the LINE-1 null is uninformative and"
           " the\n  class-selectivity argument must be made another way"
           " (see the 4 kb argument)")
print(f"VERDICT: {verdict}")

rows = []
for te in ("Alu", "LINE1"):
    o = obs[te]
    rows.append(dict(analysis=f"observed {te}", n_NDD=o["n_NDD"], n_HK=o["n_HK"],
                     mean_NDD=round(o["mean_NDD"], 4), mean_HK=round(o["mean_HK"], 4),
                     ratio=round(o["ratio"], 4), genome_mean=round(baseline[te], 4),
                     r=round(o["r"], 4), p=o["p"]))
rows.append(dict(analysis="Alu thinned to LINE-1 density", n_NDD=len(nd_n),
                 n_HK=len(hk_n), mean_NDD=np.nan, mean_HK=np.nan,
                 ratio=np.nan, genome_mean=np.nan,
                 r=round(float(np.median(rs)), 4), p=float(np.median(ps)),
                 thin_keep_p=round(float(p_keep), 4), n_draws=N_THIN,
                 r_iqr_low=round(float(np.percentile(rs, 25)), 4),
                 r_iqr_high=round(float(np.percentile(rs, 75)), 4),
                 frac_p_below_05=round(sig, 4)))

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"\nwrote {OUT.relative_to(REPO)}")
