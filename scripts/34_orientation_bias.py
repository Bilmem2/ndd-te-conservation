#!/usr/bin/env python3
"""
34_orientation_bias.py — a post-insertion signature that does not depend on
insertion opportunity.

L1-mediated retrotransposition is agnostic to the orientation of the host gene:
ORF2p nicks whichever strand carries the recognition site, so a newly inserted
Alu is equally likely to end up sense or antisense with respect to the
transcript. Anything that skews the orientation of the elements that survive
must therefore have acted after insertion. Antisense Alu in promoters and 5'
regions is the more disruptive configuration, since the antisense strand carries
cryptic polyadenylation and splice-acceptor signals, so selection is expected to
remove it preferentially.

Two questions follow:

  1. Do surviving Alu at NDD promoters show a stronger sense bias than those at
     housekeeping promoters? A difference implies stronger post-insertion
     selection at NDD loci.
  2. Does the bias deepen with element age? AluJ and AluS have had far longer to
     be purged than AluY, so progressive removal predicts an age gradient,
     whereas an insertion-level explanation predicts none.

Output: results/mechanism/orientation_bias.csv
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "hg38"
GL = REPO / "data" / "gene_lists"
OUT = REPO / "results" / "mechanism"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}


def gene_set(fn):
    return {l.strip() for l in open(GL / fn) if l.strip()}


def promoters():
    """TSS +/-2 kb windows, keeping the gene's own strand."""
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
            rows.append((f[0], max(0, tss - 2000), tss + 2000, m.group(1), strand))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "gene_strand"])


def load_alu():
    """Alu elements with strand and subfamily, from the raw RepeatMasker table."""
    rows = []
    with gzip.open(DATA / "rmsk" / "rmsk.txt.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 12 and f[11] == "SINE" and f[12] == "Alu" and f[5] in MAIN:
                name = f[10]
                sub = ("AluY" if name.startswith("AluY") else
                       "AluS" if name.startswith("AluS") else
                       "AluJ" if name.startswith("AluJ") else "other")
                rows.append((f[5], int(f[6]), int(f[7]), f[9], sub))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "alu_strand", "subfamily"])


def assign(prom, alu):
    """For each Alu, find the promoter it falls in; return per-element records."""
    recs = []
    pg = {c: g for c, g in prom.groupby("chrom")}
    for chrom, sub in alu.groupby("chrom"):
        if chrom not in pg:
            continue
        p = pg[chrom].sort_values("start")
        ps, pe = p["start"].values, p["end"].values
        gs, gene = p["gene_strand"].values, p["gene"].values
        order = np.argsort(ps)
        ps, pe, gs, gene = ps[order], pe[order], gs[order], gene[order]
        for r in sub.itertuples():
            lo = np.searchsorted(ps, r.end, "left")
            if not lo:
                continue
            hit = np.where(pe[:lo] > r.start)[0]
            for h in hit:
                recs.append((gene[h], gs[h], r.alu_strand, r.subfamily))
    return pd.DataFrame(recs, columns=["gene", "gene_strand", "alu_strand", "subfamily"])


def sense_stats(df):
    same = (df["gene_strand"].values == df["alu_strand"].values)
    return int(same.sum()), int((~same).sum())


def main():
    ndd, hk = gene_set("HighConfNDD_genes.txt"), gene_set("Housekeeping_genes.txt")
    prom = promoters()
    prom = prom[prom["gene"].isin(ndd | hk)].reset_index(drop=True)
    prom["category"] = np.where(prom["gene"].isin(ndd), "NDD", "HK")
    print(f"{len(prom)} promoters; loading Alu with strand ...")
    alu = load_alu()
    print(f"{len(alu)} Alu elements on main chromosomes")

    hits = assign(prom, alu)
    hits = hits.merge(prom[["gene", "category"]], on="gene", how="left")
    print(f"{len(hits)} Alu elements inside promoter windows\n")

    rows = []

    def compare(label, sub):
        a = sub[sub["category"] == "NDD"]
        b = sub[sub["category"] == "HK"]
        if len(a) < 30 or len(b) < 30:
            return
        a_s, a_a = sense_stats(a)
        b_s, b_a = sense_stats(b)
        odds, p = stats.fisher_exact([[a_s, a_a], [b_s, b_a]])
        rows.append(dict(stratum=label,
                         NDD_sense=a_s, NDD_antisense=a_a,
                         NDD_pct_sense=round(100 * a_s / (a_s + a_a), 1),
                         HK_sense=b_s, HK_antisense=b_a,
                         HK_pct_sense=round(100 * b_s / (b_s + b_a), 1),
                         odds_ratio=round(odds, 3), p_fisher=p))
        print(f"{label:16s} NDD {a_s:5d}/{a_a:5d} ({rows[-1]['NDD_pct_sense']:.1f}% sense)   "
              f"HK {b_s:5d}/{b_a:5d} ({rows[-1]['HK_pct_sense']:.1f}%)   "
              f"OR={odds:.3f}  p={p:.3f}")

    compare("all Alu", hits)
    for fam in ["AluY", "AluS", "AluJ"]:
        compare(fam, hits[hits["subfamily"] == fam])

    # is either set different from the 50:50 expectation of orientation-neutral insertion?
    print()
    for cat in ["NDD", "HK"]:
        s, a = sense_stats(hits[hits["category"] == cat])
        p = stats.binomtest(s, s + a, 0.5).pvalue
        print(f"{cat} vs 50:50 expectation: {s}/{s+a} sense ({100*s/(s+a):.1f}%), p={p:.2e}")
        rows.append(dict(stratum=f"{cat} vs 50:50", NDD_sense=s, NDD_antisense=a,
                         NDD_pct_sense=round(100 * s / (s + a), 1),
                         HK_sense="", HK_antisense="", HK_pct_sense="",
                         odds_ratio="", p_fisher=p))

    # --- per-subfamily depletion with bootstrap CIs on the NDD/HK ratio ---
    print("\n=== NDD/HK density ratio by subfamily (bootstrap 95% CI) ===")
    kb = (prom["end"] - prom["start"]) / 1000.0
    is_nd = (prom["category"] == "NDD").values
    sub_rows = []
    rng = np.random.default_rng(0)
    for fam in ["AluY", "AluS", "AluJ"]:
        fam_alu = alu[alu["subfamily"] == fam]
        cnt = np.zeros(len(prom), dtype=int)
        pg = {c: g for c, g in fam_alu.groupby("chrom")}
        for i, r in enumerate(prom.itertuples()):
            g = pg.get(r.chrom)
            if g is None:
                continue
            s = g["start"].values
            e = g["end"].values
            o = np.argsort(s)
            s, e = s[o], e[o]
            lo = np.searchsorted(s, r.end, "left")
            if lo:
                cnt[i] = int(np.count_nonzero(e[:lo] > r.start))
        d = cnt / kb.values
        nd_d, hk_d = d[is_nd], d[~is_nd]
        ratio = nd_d.mean() / hk_d.mean()
        b = np.empty(4000)
        for k in range(4000):
            b[k] = (rng.choice(nd_d, len(nd_d), True).mean() /
                    rng.choice(hk_d, len(hk_d), True).mean())
        lo_ci, hi_ci = np.nanpercentile(b, [2.5, 97.5])
        sub_rows.append(dict(subfamily=fam, dens_NDD=round(float(nd_d.mean()), 4),
                             dens_HK=round(float(hk_d.mean()), 4),
                             ratio_NDD_HK=round(float(ratio), 3),
                             ci_low=round(float(lo_ci), 3), ci_high=round(float(hi_ci), 3)))
        print(f"  {fam}: ratio={ratio:.3f}  95% CI [{lo_ci:.3f}, {hi_ci:.3f}]")
    pd.DataFrame(sub_rows).to_csv(OUT / "subfamily_ratio_ci.csv", index=False)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "orientation_bias.csv", index=False)
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 20)
    print("\n=== Alu orientation relative to the host gene ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'orientation_bias.csv'}")


if __name__ == "__main__":
    main()
