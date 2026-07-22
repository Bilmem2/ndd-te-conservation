#!/usr/bin/env python3
"""
18_matched_control.py — Joint (multivariate) confounder control for Alu
depletion at NDD promoters (hg38).

The paper controls GC, then (script 15) gene density and recombination, each
marginally.  A skeptic can argue NDD and HK promoters differ on several
correlated properties at once, so one-variable-at-a-time stratification is not
enough.  Here we build a covariate-MATCHED housekeeping control by 1:1
nearest-neighbour Mahalanobis matching on three covariates simultaneously:

    GC content, local gene density (genes/Mb), local recombination (cM/Mb)

then test whether Alu depletion (HK > NDD) survives against the matched set,
and report covariate balance (standardized mean differences) before vs after
matching.

GC is computed directly from hg38.2bit for each promoter window (twobitreader).
Gene density and recombination are reused from results/context/.

Inputs:
  results/context/context_per_promoter.tsv   (from 15_context_controls.py)
  data/hg38/hg38.2bit
Outputs:
  results/matched/{matched_covariate_balance.csv, matched_alu_test.csv,
                   matched_pairs.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from twobitreader import TwoBitFile

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
DATA = REPO / "data"
OUT = RESULTS / "matched"
OUT.mkdir(parents=True, exist_ok=True)
COVS = ["gc", "gene_density", "recomb_cM_Mb"]


def gc_content(tb, chrom, start, end):
    try:
        seq = tb[chrom][int(start):int(end)].upper()
    except Exception:
        return np.nan
    g = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    n = g + at
    return g / n if n else np.nan


def smd(a, b):
    """standardized mean difference (pooled sd)."""
    sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    return (np.mean(a) - np.mean(b)) / sd if sd else np.nan


def greedy_match(nd, hk, cov_cols):
    """1:1 nearest-neighbour Mahalanobis matching, HK without replacement.

    NDD promoters matched hardest-first (largest min-distance) is overkill;
    we use natural order with greedy nearest available HK. Returns HK indices."""
    X = pd.concat([nd[cov_cols], hk[cov_cols]]).values
    cov = np.cov(X, rowvar=False)
    VI = np.linalg.pinv(cov)
    ndX = nd[cov_cols].values
    hkX = hk[cov_cols].values
    used = np.zeros(len(hk), dtype=bool)
    chosen = np.full(len(nd), -1, dtype=int)
    for i in range(len(nd)):
        diff = hkX - ndX[i]
        d2 = np.einsum("ij,jk,ik->i", diff, VI, diff)
        d2[used] = np.inf
        j = int(np.argmin(d2))
        if np.isfinite(d2[j]):
            chosen[i] = j
            used[j] = True
    return chosen


def main():
    prom = pd.read_csv(RESULTS / "context" / "context_per_promoter.tsv", sep="\t")
    tb = TwoBitFile(str(DATA / "hg38" / "hg38.2bit"))
    prom["gc"] = [gc_content(tb, c, s, e) for c, s, e in
                  zip(prom["chrom"], prom["start"], prom["end"])]
    prom["alu_density"] = prom["fixed_alu"] / ((prom["end"] - prom["start"]) / 1000)
    prom = prom.dropna(subset=COVS + ["alu_density"]).reset_index(drop=True)

    nd = prom[prom["category"] == "HighConfNDD"].reset_index(drop=True)
    hk = prom[prom["category"] == "Housekeeping"].reset_index(drop=True)

    chosen = greedy_match(nd, hk, COVS)
    ok = chosen >= 0
    matched_hk = hk.iloc[chosen[ok]].reset_index(drop=True)
    nd_m = nd[ok].reset_index(drop=True)

    # covariate balance before (all HK) vs after (matched HK)
    bal = []
    for c in COVS:
        bal.append(dict(covariate=c,
                        NDD_mean=round(nd_m[c].mean(), 4),
                        HK_all_mean=round(hk[c].mean(), 4),
                        HK_matched_mean=round(matched_hk[c].mean(), 4),
                        SMD_before=round(smd(nd_m[c].values, hk[c].values), 3),
                        SMD_after=round(smd(nd_m[c].values, matched_hk[c].values), 3)))
    balance = pd.DataFrame(bal)

    # Alu depletion: before (all HK) vs after (matched HK)
    def alu_test(hk_set, label):
        _, p = stats.mannwhitneyu(hk_set["alu_density"], nd_m["alu_density"], alternative="greater")
        u, _ = stats.mannwhitneyu(hk_set["alu_density"], nd_m["alu_density"], alternative="two-sided")
        r = 1 - (2 * u) / (len(hk_set) * len(nd_m))
        return dict(comparison=label, n_HK=len(hk_set), n_NDD=len(nd_m),
                    median_alu_HK=round(hk_set["alu_density"].median(), 3),
                    median_alu_NDD=round(nd_m["alu_density"].median(), 3),
                    mean_alu_HK=round(hk_set["alu_density"].mean(), 3),
                    mean_alu_NDD=round(nd_m["alu_density"].mean(), 3),
                    p_value=p, r=round(r, 3))
    tests = pd.DataFrame([alu_test(hk, "NDD vs ALL housekeeping"),
                          alu_test(matched_hk, "NDD vs MATCHED housekeeping")])

    pairs = pd.DataFrame({
        "ndd_gene": nd_m["gene"].values, "matched_hk_gene": matched_hk["gene"].values,
        "ndd_gc": nd_m["gc"].round(3).values, "hk_gc": matched_hk["gc"].round(3).values,
        "ndd_genedens": nd_m["gene_density"].values, "hk_genedens": matched_hk["gene_density"].values,
        "ndd_recomb": nd_m["recomb_cM_Mb"].round(3).values, "hk_recomb": matched_hk["recomb_cM_Mb"].round(3).values,
        "ndd_alu_density": nd_m["alu_density"].round(3).values,
        "hk_alu_density": matched_hk["alu_density"].round(3).values})

    balance.to_csv(OUT / "matched_covariate_balance.csv", index=False)
    tests.to_csv(OUT / "matched_alu_test.csv", index=False)
    pairs.to_csv(OUT / "matched_pairs.tsv", sep="\t", index=False)

    pd.set_option("display.width", 180); pd.set_option("display.max_columns", 25)
    print(f"\nMatched {ok.sum()} of {len(nd)} NDD promoters to unique HK controls "
          f"(HK pool = {len(hk)}).")
    print("\n=== Covariate balance (|SMD|<0.1 = well matched) ===")
    print(balance.to_string(index=False))
    print("\n=== Alu depletion (HK > NDD): before vs after joint matching ===")
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
