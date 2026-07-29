#!/usr/bin/env python3
"""
33_insertion_opportunity.py — is the deficit explained by where Alu *can* insert?

Extant Alu density is the product of insertion opportunity and post-insertion
retention, and the two cannot be separated by counting fixed elements. Alu is
non-autonomous and retrotransposes using the L1 ORF2p endonuclease, which nicks
DNA at a degenerate A/T-rich recognition site (canonically 5'-TTAAAA-3'). Alu can
therefore only integrate where such a site exists, and the density of these sites
is a sequence-intrinsic proxy for insertion opportunity that does not depend on
what happened after insertion.

The comparison is unusually well controlled here: the recognition site is A/T
rich, and NDD and housekeeping promoters have already been shown to have nearly
identical GC content, so the base-composition confound that would normally
dominate such a test is absent.

  fewer sites at NDD promoters  -> lower insertion opportunity; targeting bias
                                   contributes to the deficit
  equal site density            -> equal opportunity; the deficit arises after
                                   insertion, which favours selection

The script first validates the proxy by asking whether site density predicts Alu
density across all protein-coding promoters. If it does not, the test is
uninformative and is reported as such.

Output: results/mechanism/insertion_opportunity.csv
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
OUT = REPO / "results" / "mechanism"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}

# L1 ORF2p endonuclease recognition site. STRICT is the canonical consensus;
# RELAXED admits the commonly reported degenerate variants. Both strands are
# counted, since the nick can occur on either.
STRICT = [re.compile("TTAAAA"), re.compile("TTTTAA")]
RELAXED = [re.compile("TT[AG]AAA"), re.compile("TTT[CT]AA")]


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


def n_overlapping(pattern, seq):
    """Count matches including overlaps."""
    n, pos = 0, 0
    while True:
        m = pattern.search(seq, pos)
        if not m:
            return n
        n += 1
        pos = m.start() + 1


def scan(prom):
    tb = TwoBitFile(str(DATA / "hg38.2bit"))
    strict = np.zeros(len(prom))
    relaxed = np.zeros(len(prom))
    gc = np.full(len(prom), np.nan)
    for i, r in enumerate(prom.itertuples()):
        try:
            s = tb[r.chrom][int(r.start):int(r.end)].upper()
        except Exception:
            continue
        if not s or s.count("N") > 0.5 * len(s):
            continue
        strict[i] = sum(n_overlapping(p, s) for p in STRICT)
        relaxed[i] = sum(n_overlapping(p, s) for p in RELAXED)
        g, at = s.count("G") + s.count("C"), s.count("A") + s.count("T")
        if g + at:
            gc[i] = g / (g + at)
    return strict, relaxed, gc


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    ndd, hk = gene_set("HighConfNDD_genes.txt"), gene_set("Housekeeping_genes.txt")
    prom = promoters()
    alu = pd.read_csv(DATA / "alu_rmsk.bed", sep="\t", header=None,
                      usecols=[0, 1, 2], names=["chrom", "start", "end"])
    alu = alu[alu["chrom"].isin(MAIN)]

    kb = (prom["end"] - prom["start"]) / 1000.0
    prom["alu_d"] = count_alu(prom, alu) / kb
    print(f"scanning {len(prom)} promoter windows for L1 endonuclease sites ...")
    strict, relaxed, gc = scan(prom)
    prom["en_strict"] = strict / kb
    prom["en_relaxed"] = relaxed / kb
    prom["gc"] = gc
    prom = prom.dropna(subset=["gc"]).reset_index(drop=True)

    # --- validity check: does site density predict Alu density at all? ---
    rho_s, p_s = stats.spearmanr(prom["en_strict"], prom["alu_d"])
    rho_r, p_r = stats.spearmanr(prom["en_relaxed"], prom["alu_d"])
    print(f"\nProxy validity (all {len(prom)} promoters):")
    print(f"  Spearman rho, strict site density vs Alu density : {rho_s:+.3f} (p={p_s:.2e})")
    print(f"  Spearman rho, relaxed site density vs Alu density: {rho_r:+.3f} (p={p_r:.2e})")

    is_nd = prom["gene"].isin(ndd).values
    is_hk = prom["gene"].isin(hk).values

    rows = []
    for col, label in [("en_strict", "L1 EN site (TTAAAA)"),
                       ("en_relaxed", "L1 EN site (degenerate)"),
                       ("gc", "GC content"),
                       ("alu_d", "Alu density")]:
        nd, hkv = prom[col].values[is_nd], prom[col].values[is_hk]
        _, p_two = stats.mannwhitneyu(nd, hkv, alternative="two-sided")
        rows.append(dict(measure=label, n_NDD=len(nd), n_HK=len(hkv),
                         mean_NDD=round(float(nd.mean()), 4),
                         mean_HK=round(float(hkv.mean()), 4),
                         NDD_over_HK=round(float(nd.mean() / hkv.mean()), 3),
                         genome_mean=round(float(prom[col].mean()), 4),
                         NDD_over_genome=round(float(nd.mean() / prom[col].mean()), 3),
                         r=round(rb(nd, hkv), 3), p_two_sided=p_two))
        print(f"  {label:26s} NDD={nd.mean():.4f}  HK={hkv.mean():.4f}  "
              f"NDD/HK={rows[-1]['NDD_over_HK']:.3f}  p={p_two:.2e}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "insertion_opportunity.csv", index=False)
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 20)
    print("\n=== Insertion opportunity vs realised Alu content ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'insertion_opportunity.csv'}")


if __name__ == "__main__":
    main()
