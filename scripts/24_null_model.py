#!/usr/bin/env python3
"""
24_null_model.py — Consolidated permutation null model for every species-TE
combination reported in the manuscript (Results, "Depletion Is Robust to Null
Model Permutation"; Fig. 4).

Earlier runs left the permutation evidence split across two places: the panels
of Fig. 4 (computed inline by 12_figures_final.py) and a standalone
null_model_new.csv. Neither covered every combination the manuscript reports,
and the standalone file had gone stale - its observed p-values predated later
pipeline fixes and no longer matched statistics_final.csv, so it was removed.
This script recomputes all five combinations from the committed promoter BEDs
and writes a single table whose observed p-values reproduce the final
statistics exactly.

For each combination the observed one-sided Mann-Whitney U p-value (HK > NDD)
is compared with n = 10,000 random partitions of the pooled promoter set into
groups of the observed sizes. The empirical p-value is the fraction of permuted
p-values at least as extreme as the observed one; the null false-positive rate
is the fraction of permuted p-values below 0.05 (expected ~0.05 if calibrated).

No raw genome data required - reads results/<assembly>/<category>_<TE>.bed.
Output: results/null_model_full.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
N_PERM = 10000
SEED = 42

COMBOS = [
    ("hg38", "Human", "Alu"),
    ("ponAbe3", "Orangutan", "Alu"),
    ("rheMac10", "Macaque", "Alu"),
    ("calJac4", "Marmoset", "Alu"),
    ("mm10", "Mouse", "B1B2"),
]


def density(path):
    """TE count per kb for each promoter window (col 6 = count)."""
    df = pd.read_csv(path, sep="\t", header=None)
    return (df.iloc[:, 6] / ((df.iloc[:, 2] - df.iloc[:, 1]) / 1000)).values


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for assembly, label, te in COMBOS:
        hk = density(RESULTS / assembly / f"Housekeeping_{te}.bed")
        ndd = density(RESULTS / assembly / f"HighConfNDD_{te}.bed")
        _, obs_p = stats.mannwhitneyu(hk, ndd, alternative="greater")

        pooled = np.concatenate([hk, ndd])
        n_hk = len(hk)
        perm_p = np.empty(N_PERM)
        for i in range(N_PERM):
            s = rng.permutation(pooled)
            _, perm_p[i] = stats.mannwhitneyu(
                s[:n_hk], s[n_hk:], alternative="greater")

        rows.append(dict(
            assembly=assembly, species=label, TE=te,
            n_HK=n_hk, n_NDD=len(ndd),
            observed_p=obs_p,
            empirical_p=float(np.mean(perm_p <= obs_p)),
            null_fpr=float(np.mean(perm_p < 0.05)),
            n_perm=N_PERM))
        print(f"{label:10s} {te:5s}  observed p={obs_p:.3e}  "
              f"empirical p={rows[-1]['empirical_p']:.4f}  "
              f"null FPR={rows[-1]['null_fpr']:.4f}")

    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "null_model_full.csv", index=False)
    print(f"\nWrote {RESULTS / 'null_model_full.csv'}")


if __name__ == "__main__":
    main()
