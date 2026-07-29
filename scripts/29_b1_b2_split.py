#!/usr/bin/env python3
"""
29_b1_b2_split.py — is the mouse SINE depletion carried by B1, by B2, or both?

The main analysis pools mouse B1 and B2 as the functional SINE analogues of
primate Alu. They have different ancestries: B1 derives from 7SL/SRP RNA, as Alu
does, whereas B2 derives from tRNA. The convergent-exclusion argument is stronger
if the depletion holds for each family separately rather than being carried by
the 7SL-derived family alone. This script repeats the mouse comparison for B1 and
B2 individually, reusing the committed promoter windows.

Also reports Benjamini-Hochberg q-values for the primary species-TE tests.

Output: results/mm10/b1_b2_split.csv, results/consolidated/fdr_table.csv
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RMSK = REPO / "data" / "mm10" / "rmsk" / "rmsk.txt.gz"
RES = REPO / "results"


def load_promoters(cat):
    df = pd.read_csv(RES / "mm10" / f"{cat}_B1B2.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand", "n"])
    return df[["chrom", "start", "end", "gene"]]


def load_family(family):
    """Stream rmsk.txt.gz and keep one SINE family (repFamily column 12)."""
    rows = []
    with gzip.open(RMSK, "rt") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 12 and f[11] == "SINE" and f[12] == family:
                rows.append((f[5], int(f[6]), int(f[7])))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def count(prom, feats):
    counts = np.zeros(len(prom), dtype=int)
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
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        r = n - rank
        prev = min(prev, p[idx] * n / r)
        q[idx] = prev
    return q


def main():
    nd = load_promoters("HighConfNDD")
    hk = load_promoters("Housekeeping")
    prom = pd.concat([nd.assign(cat="NDD"), hk.assign(cat="HK")], ignore_index=True)
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    is_nd = prom["cat"].values == "NDD"

    rows = []
    for family, label in [("Alu", "B1 (SINE/Alu, 7SL-derived)"),
                          ("B2", "B2 (SINE/B2, tRNA-derived)")]:
        feats = load_family(family)
        d = count(prom, feats) / prom["width_kb"].values
        a, b = d[~is_nd], d[is_nd]
        _, p = stats.mannwhitneyu(a, b, alternative="greater")
        rows.append(dict(family=label, n_elements=len(feats),
                         n_HK=int((~is_nd).sum()), n_NDD=int(is_nd.sum()),
                         median_HK=round(float(np.median(a)), 3),
                         median_NDD=round(float(np.median(b)), 3),
                         r=round(rb(a, b), 3), p_value=p))
        print(f"{label:32s} {len(feats):8d} elements  r={rows[-1]['r']:+.3f}  p={p:.3e}")

    out = pd.DataFrame(rows)
    out.to_csv(RES / "mm10" / "b1_b2_split.csv", index=False)

    # BH correction over the primary species-TE tests reported in Table 1
    tab = pd.read_csv(RES / "consolidated" / "cross_species.csv")
    tab["q_value_BH"] = bh(tab["p_value"].values)
    tab[["species", "TE", "p_value", "q_value_BH", "r"]].to_csv(
        RES / "consolidated" / "fdr_table.csv", index=False)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n=== Mouse SINE depletion by family ===")
    print(out.to_string(index=False))
    print("\n=== Benjamini-Hochberg q-values over the 14 primary tests ===")
    print(tab[["species", "TE", "p_value", "q_value_BH", "r"]].to_string(index=False))


if __name__ == "__main__":
    main()
