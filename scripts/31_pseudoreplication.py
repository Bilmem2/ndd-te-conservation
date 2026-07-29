#!/usr/bin/env python3
"""
31_pseudoreplication.py — how much of the human signal survives when
non-independent promoters are removed?

Promoters are treated as independent observations throughout, but they are not:
disease and housekeeping genes both cluster in the genome (NDD genes notably at
recurrent CNV loci such as 16p11.2 and 22q11.2), and both sets contain paralogue
families. The nominal sample size therefore overstates the effective one. This
script re-runs the human Alu comparison under two deliberately aggressive
thinning schemes:

  1. genomic thinning - within each gene set, retain a promoter only if it lies
     at least D bp from the previously retained promoter on the same chromosome
     (D = 1 Mb and 5 Mb), which removes local clusters;
  2. one-per-family - retain one gene per symbol root (trailing digits stripped,
     e.g. SOX1/SOX2/SOX3 -> SOX), a crude but conservative proxy for paralogue
     families in the absence of a curated family table.

Both schemes discard real data and so lose power; the question is whether the
effect size, not the p-value, is stable.

Output: results/sensitivity/pseudoreplication.csv
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RES = REPO / "results"
OUT = RES / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

ROOT_RE = re.compile(r"^([A-Za-z]+?)\d*[A-Z]?$")


def load(cat):
    df = pd.read_csv(RES / "hg38" / f"{cat}_Alu.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand", "n"])
    df["alu_d"] = df["n"] / ((df["end"] - df["start"]) / 1000.0)
    return df


def thin_by_distance(df, min_dist):
    """Keep promoters at least min_dist apart on each chromosome (greedy, left to right)."""
    keep = []
    for chrom, sub in df.sort_values(["chrom", "start"]).groupby("chrom"):
        last = -np.inf
        for row in sub.itertuples():
            if row.start - last >= min_dist:
                keep.append(row.Index)
                last = row.start
    return df.loc[keep]


def thin_by_family(df):
    """Keep one gene per symbol root, as a rough paralogue-family proxy."""
    roots = []
    for g in df["gene"]:
        m = ROOT_RE.match(str(g))
        roots.append(m.group(1) if m else str(g))
    return df.assign(root=roots).drop_duplicates("root")


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    nd_all, hk_all = load("HighConfNDD"), load("Housekeeping")

    schemes = [("none (full set)", lambda d: d)]
    for mb in (1, 5):
        schemes.append((f"genomic thinning, {mb} Mb", lambda d, m=mb: thin_by_distance(d, m * 1_000_000)))
    schemes.append(("one per symbol root", thin_by_family))

    rows = []
    for label, fn in schemes:
        nd, hk = fn(nd_all), fn(hk_all)
        a, b = hk["alu_d"].values, nd["alu_d"].values
        _, p = stats.mannwhitneyu(a, b, alternative="greater")
        rows.append(dict(scheme=label, n_HK=len(hk), n_NDD=len(nd),
                         pct_retained_NDD=round(100 * len(nd) / len(nd_all), 1),
                         median_HK=round(float(np.median(a)), 3),
                         median_NDD=round(float(np.median(b)), 3),
                         r=round(rb(a, b), 3), p_value=p))
        print(f"{label:26s} n={len(hk):5d}/{len(nd):5d}  r={rows[-1]['r']:+.3f}  p={p:.2e}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "pseudoreplication.csv", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n=== Alu depletion under thinning of non-independent promoters (hg38) ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'pseudoreplication.csv'}")


if __name__ == "__main__":
    main()
