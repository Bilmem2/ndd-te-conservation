#!/usr/bin/env python3
"""
26_cross_disease_recompute.py — repair the cross-disease gene sets and redo the
comparison.

The original 10_cross_disease.py took ClinVar's GeneSymbol column verbatim.
That column holds semicolon-separated multi-gene strings and a literal "-" for
variants with no assigned gene, and the promoter set was then selected with
`grep -Fw -f <gene list>` against whole BED lines. Two things went wrong:

  * composite entries such as "BGLT3;HBE1;HBG1" never matched any promoter, so
    the genes inside them were silently lost; and
  * the bare "-" pattern matched the strand column of every minus-strand
    promoter in the genome, injecting ~9,900 arbitrary genes into the
    "Mendelian" set (9,912 minus- vs 2,704 plus-strand promoters, where a real
    gene set is near 50/50).

The Mendelian comparison was therefore computed largely on a random half of the
genome. This script rebuilds both sets by splitting on ";", dropping non-symbol
tokens, and re-removing the NDD and housekeeping genes, then matches promoters
on exact gene symbol rather than by line grep. Results are reported both against
housekeeping promoters (as before) and against the within-species genome
baseline, since housekeeping promoters are themselves Alu-rich.

Inputs : data/gene_lists/*.txt, data/hg38/gtf/gencode.v47.gtf.gz,
         data/hg38/alu_rmsk.bed
Outputs: corrected data/gene_lists/{Cardiovascular,Mendelian}_genes.txt
         results/cross_disease/cross_disease_recomputed.csv
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
GL = DATA / "gene_lists"
OUT = REPO / "results" / "cross_disease"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.@-]*$")
N_BOOT = 4000


def read_list(fn):
    return [l.strip() for l in (GL / fn).read_text(encoding="utf-8").splitlines()]


def clean_symbols(raw_lines):
    """Split ClinVar composite entries and keep plausible gene symbols."""
    out = set()
    for line in raw_lines:
        for tok in line.replace("\r", "").split(";"):
            tok = tok.strip()
            if tok and tok != "-" and SYMBOL_RE.match(tok):
                out.add(tok)
    return out


def build_promoters():
    """TSS +/-2 kb for every protein-coding gene in GENCODE v47."""
    rows, seen = [], set()
    with gzip.open(DATA / "hg38" / "gtf" / "gencode.v47.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or 'gene_type "protein_coding"' not in f[8]:
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) in seen:
                continue
            seen.add(m.group(1))
            start, end, strand = int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            rows.append((f[0], max(0, tss - 2000), tss + 2000, m.group(1), strand))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene", "strand"])


def count_alu(prom, alu):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in alu.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        order = np.argsort(m["start"].values)
        s = m["start"].values[order]
        e = m["end"].values[order]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def ratio_ci(group, baseline, seed):
    rng = np.random.default_rng(seed)
    b = np.empty(N_BOOT)
    for k in range(N_BOOT):
        b[k] = rng.choice(group, len(group), True).mean() / \
            rng.choice(baseline, len(baseline), True).mean()
    return np.nanpercentile(b, [2.5, 97.5])


def main():
    ndd = set(read_list("HighConfNDD_genes.txt"))
    hk = set(read_list("Housekeeping_genes.txt"))

    sets = {}
    for name in ("Cardiovascular", "Mendelian"):
        raw = read_list(f"{name}_genes.txt")
        cleaned = clean_symbols(raw) - ndd - hk
        sets[name] = cleaned
        (GL / f"{name}_genes.txt").write_text("\n".join(sorted(cleaned)) + "\n",
                                              encoding="utf-8")
        print(f"{name}: {len(raw)} raw lines -> {len(cleaned)} clean symbols")

    prom = build_promoters()
    alu = pd.read_csv(DATA / "hg38" / "alu_rmsk.bed", sep="\t", header=None,
                      usecols=[0, 1, 2], names=["chrom", "start", "end"])
    prom["alu_d"] = count_alu(prom, alu) / ((prom["end"] - prom["start"]) / 1000.0)
    genome_mean = prom["alu_d"].mean()
    print(f"\nGenome baseline: {genome_mean:.3f} Alu/kb over {len(prom)} promoters")

    hk_d = prom["alu_d"].values[prom["gene"].isin(hk).values]
    rows = []
    for label, genes in [("HighConfNDD", ndd), ("Cardiovascular", sets["Cardiovascular"]),
                         ("Mendelian", sets["Mendelian"])]:
        mask = prom["gene"].isin(genes).values
        d = prom["alu_d"].values[mask]
        strand = prom["strand"].values[mask]
        _, p_hk = stats.mannwhitneyu(hk_d, d, alternative="greater")
        rest = prom["alu_d"].values[~mask]
        _, p_gen = stats.mannwhitneyu(rest, d, alternative="greater")
        ci = ratio_ci(d, prom["alu_d"].values, 7)
        rows.append(dict(
            gene_set=label, n_promoters=int(mask.sum()),
            pct_plus_strand=round(100 * float((strand == "+").mean()), 1),
            median_alu=round(float(np.median(d)), 3),
            r_vs_housekeeping=round(rb(hk_d, d), 3), p_vs_housekeeping=p_hk,
            ratio_vs_genome=round(float(d.mean()) / genome_mean, 3),
            genome_CI=f"[{ci[0]:.3f}, {ci[1]:.3f}]", p_vs_genome=p_gen))
        print(f"  {label:16s} n={mask.sum():5d}  +strand={rows[-1]['pct_plus_strand']:.1f}%  "
              f"r_vs_HK={rows[-1]['r_vs_housekeeping']:+.3f}  "
              f"genome={rows[-1]['ratio_vs_genome']:.3f} {rows[-1]['genome_CI']}  "
              f"p_gen={p_gen:.2e}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "cross_disease_recomputed.csv", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n=== Cross-disease comparison, corrected gene sets ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'cross_disease_recomputed.csv'}")


if __name__ == "__main__":
    main()
