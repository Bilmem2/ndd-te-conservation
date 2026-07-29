#!/usr/bin/env python3
"""
32_matching_sensitivity.py — is the constraint-matched result an artefact of one
matching run?

The NDD-specific residual (r = -0.106 against controls matched on LOEUF, brain
expression and GC) is the number the "extreme of a gradient" reading rests on,
and it comes from a single greedy nearest-neighbour pass in natural gene order.
Greedy matching is order-dependent: genes matched early get their closest
partners, and later ones take what is left. This script varies the things that
could plausibly drive the result:

  * matching order      - 25 random permutations of the NDD set
  * distance metric     - Mahalanobis vs standardised Euclidean
  * caliper             - discarding pairs beyond the 90th percentile distance
  * covariate set       - LOEUF alone; LOEUF + brain expression; all three

If the residual is real, the effect size should sit in the same range across all
of these.

Output: results/matched/matching_sensitivity.csv
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from twobitreader import TwoBitFile

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "hg38"
GL = REPO / "data" / "gene_lists"
OUT = REPO / "results" / "matched"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}
N_ORDERS = 25


def gene_set(fn):
    return {l.strip() for l in open(GL / fn) if l.strip()}


def promoters():
    rows, seen = [], set()
    with gzip.open(DATA / "gtf" / "gencode.v47.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or 'gene_type "protein_coding"' not in f[8] or f[0] not in MAIN:
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) in seen:
                continue
            seen.add(m.group(1))
            start, end, strand = int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            rows.append((f[0], max(0, tss - 2000), tss + 2000, m.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def count_alu(prom, alu):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in alu.groupby("chrom"):
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


def gc_content(prom):
    tb = TwoBitFile(str(DATA / "hg38.2bit"))
    out = np.full(len(prom), np.nan)
    for i, r in enumerate(prom.itertuples()):
        try:
            s = tb[r.chrom][int(r.start):int(r.end)].upper()
        except Exception:
            continue
        g, at = s.count("G") + s.count("C"), s.count("A") + s.count("T")
        if g + at:
            out[i] = g / (g + at)
    return out


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def match(ndX, poX, metric, order, caliper=None):
    """Greedy 1:1 nearest neighbour; returns (ndd_idx, pool_idx) arrays."""
    if metric == "mahalanobis":
        # np.cov collapses to a scalar for a single covariate
        cov = np.atleast_2d(np.cov(np.vstack([ndX, poX]), rowvar=False))
        VI = np.linalg.pinv(cov)
    else:  # standardised euclidean
        sd = np.vstack([ndX, poX]).std(axis=0)
        VI = np.diag(1.0 / np.where(sd > 0, sd, 1.0) ** 2)
    used = np.zeros(len(poX), bool)
    nd_keep, po_keep, dists = [], [], []
    for i in order:
        d = poX - ndX[i]
        d2 = np.einsum("ij,jk,ik->i", d, VI, d)
        d2[used] = np.inf
        j = int(np.argmin(d2))
        if np.isfinite(d2[j]):
            used[j] = True
            nd_keep.append(i)
            po_keep.append(j)
            dists.append(d2[j])
    nd_keep, po_keep, dists = np.array(nd_keep), np.array(po_keep), np.array(dists)
    if caliper is not None and len(dists):
        thr = np.quantile(dists, caliper)
        ok = dists <= thr
        nd_keep, po_keep = nd_keep[ok], po_keep[ok]
    return nd_keep, po_keep


def main():
    ndd = gene_set("HighConfNDD_genes.txt")
    con = pd.read_csv(DATA / "gnomad_constraint.tsv", sep="\t",
                      usecols=["gene", "mane_select", "lof.oe_ci.upper"])
    con = con[con["mane_select"] == True]
    con = (con.groupby("gene")["lof.oe_ci.upper"].mean()
           .rename("loeuf").reset_index())
    gx = pd.read_csv(DATA / "gtex_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
    bcols = [c for c in gx.columns if c.startswith("Brain")]
    gx = pd.DataFrame({"gene": gx["Description"],
                       "brain_tpm": gx[bcols].max(axis=1)}).drop_duplicates("gene")

    prom = promoters()
    alu = pd.read_csv(DATA / "alu_rmsk.bed", sep="\t", header=None,
                      usecols=[0, 1, 2], names=["chrom", "start", "end"])
    alu = alu[alu["chrom"].isin(MAIN)]
    prom["alu_d"] = count_alu(prom, alu) / ((prom["end"] - prom["start"]) / 1000.0)
    prom["gc"] = gc_content(prom)
    prom = prom.merge(con, on="gene", how="left").merge(gx, on="gene", how="left")
    prom["brain_log"] = np.log10(prom["brain_tpm"] + 1)

    usable = prom.dropna(subset=["loeuf", "brain_log", "gc"]).reset_index(drop=True)
    nd = usable[usable["gene"].isin(ndd)].reset_index(drop=True)
    pool = usable[~usable["gene"].isin(ndd)].reset_index(drop=True)
    nd_alu, pool_alu = nd["alu_d"].values, pool["alu_d"].values

    rows = []

    def run(label, covs, metric, orders, caliper=None):
        ndX, poX = nd[covs].values, pool[covs].values
        rs, ps, ns = [], [], []
        for seed in orders:
            order = (np.arange(len(ndX)) if seed is None
                     else np.random.default_rng(seed).permutation(len(ndX)))
            i, j = match(ndX, poX, metric, order, caliper)
            a, b = pool_alu[j], nd_alu[i]
            _, p = stats.mannwhitneyu(a, b, alternative="greater")
            rs.append(rb(a, b)); ps.append(p); ns.append(len(i))
        rows.append(dict(scheme=label, covariates="+".join(covs), metric=metric,
                         n_runs=len(orders), n_pairs=int(np.median(ns)),
                         r_median=round(float(np.median(rs)), 3),
                         r_min=round(float(np.min(rs)), 3),
                         r_max=round(float(np.max(rs)), 3),
                         p_max=float(np.max(ps))))
        print(f"{label:38s} r={rows[-1]['r_median']:+.3f} "
              f"[{rows[-1]['r_min']:+.3f},{rows[-1]['r_max']:+.3f}]  "
              f"n={rows[-1]['n_pairs']}  p_max={rows[-1]['p_max']:.1e}")

    ALL = ["loeuf", "brain_log", "gc"]
    run("original (natural order)", ALL, "mahalanobis", [None])
    run("25 random matching orders", ALL, "mahalanobis", list(range(25)))
    run("standardised Euclidean", ALL, "euclidean", list(range(10)))
    run("caliper, closest 90%", ALL, "mahalanobis", list(range(10)), caliper=0.90)
    run("LOEUF only", ["loeuf"], "mahalanobis", list(range(10)))
    run("LOEUF + brain expression", ["loeuf", "brain_log"], "mahalanobis", list(range(10)))

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "matching_sensitivity.csv", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n=== Constraint-matched residual under alternative matching schemes ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'matching_sensitivity.csv'}")


if __name__ == "__main__":
    main()
