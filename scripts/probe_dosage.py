#!/usr/bin/env python3
"""
probe_dosage.py — FEASIBILITY PROBE (not final analysis).

Question: does the targeting-vs-selection signal SCALE with gene dosage
sensitivity across ALL protein-coding genes (not just NDD vs HK)?

Design (gene-level, genome-wide, for power):
  For every protein-coding gene with a gnomAD constraint score, take its
  promoter (TSS +/- 2 kb) and measure:
    * FIXED Alu density      (RepeatMasker, ~1.28M Alu)          -> post-selection
    * POLYMORPHIC Alu density(gnomAD v4.1 INS:ME:ALU)            -> pre-selection snapshot
    * AFS of promoter MEIs    (singleton fraction, median AF)     -> selection signature
  Then stratify genes by constraint (pLI bins; LOEUF deciles) and ask:
    - does FIXED density DECLINE with constraint?
    - does POLYMORPHIC density stay ~FLAT? (=> deficit is post-insertion selection)
    - does singleton fraction RISE with constraint? (=> stronger ongoing selection)

  POSITIVE result = fixed declines + poly flat + singletons rise with dosage
  sensitivity => the selection-decomposition generalizes into a continuum,
  and NDD (extreme pLI) is its endpoint. NEGATIVE = signal is NDD-only / weak.

Inputs (all local): gencode.v47.gtf.gz, gnomad_constraint.tsv, alu_rmsk.bed,
gnomad_mei.tsv.  No bedtools.
Output: results/functional/probe_dosage_strata.csv (+ console summary)
"""
import gzip, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "hg38"
OUT = REPO / "results" / "functional"
OUT.mkdir(parents=True, exist_ok=True)
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}
NAME_RE = re.compile(r'gene_name "([^"]+)"')


def load_constraint():
    c = pd.read_csv(DATA / "gnomad_constraint.tsv", sep="\t",
                    usecols=["gene", "mane_select", "lof.pLI", "lof.oe"], dtype={"mane_select": str})
    c = c[c["mane_select"].str.lower() == "true"].copy()
    c["lof.pLI"] = pd.to_numeric(c["lof.pLI"], errors="coerce")
    c["lof.oe"] = pd.to_numeric(c["lof.oe"], errors="coerce")
    return c.dropna(subset=["lof.pLI"]).drop_duplicates("gene")[["gene", "lof.pLI", "lof.oe"]]


def build_promoters(genes_wanted):
    rows = []
    seen = set()
    with gzip.open(DATA / "gtf" / "gencode.v47.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or 'gene_type "protein_coding"' not in f[8]:
                continue
            chrom = f[0]
            if chrom not in MAIN:
                continue
            m = NAME_RE.search(f[8])
            if not m:
                continue
            g = m.group(1)
            if g not in genes_wanted or g in seen:
                continue
            seen.add(g)
            start, end, strand = int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            rows.append((chrom, max(0, tss - 2000), tss + 2000, g))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def load_bed_alu():
    df = pd.read_csv(DATA / "alu_rmsk.bed", sep="\t", header=None, names=["chrom", "start", "end"])
    return df[df["chrom"].isin(MAIN)]


def load_mei_alu():
    df = pd.read_csv(DATA / "gnomad_mei.tsv", sep="\t").rename(columns={"#chrom": "chrom"})
    df = df[(df["svtype"] == "INS:ME:ALU") & (df["chrom"].isin(MAIN))].copy()
    df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
    df["AC"] = pd.to_numeric(df["AC"], errors="coerce")
    return df.dropna(subset=["AF"])


def count_fixed(prom, feats):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in feats.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        s = np.sort(m["start"].values)
        e = m["end"].values[np.argsort(m["start"].values)]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def count_poly(prom, mei):
    counts = np.zeros(len(prom), dtype=int)
    af_by_row = [[] for _ in range(len(prom))]
    ac_by_row = [[] for _ in range(len(prom))]
    for chrom, m in mei.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        order = np.argsort(m["start"].values)
        s = m["start"].values[order]; e = m["end"].values[order]
        af = m["AF"].values[order]; ac = m["AC"].values[order]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if not lo:
                continue
            hit = np.where(e[:lo] > ws)[0]
            counts[i] = len(hit)
            for h in hit:
                af_by_row[i].append(af[h]); ac_by_row[i].append(ac[h])
    return counts, af_by_row, ac_by_row


def main():
    con = load_constraint()
    prom = build_promoters(set(con["gene"]))
    prom = prom.merge(con, on="gene").reset_index(drop=True)
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    print(f"Genes with promoter + constraint: {len(prom)}")

    prom["fixed_alu"] = count_fixed(prom, load_bed_alu())
    pc, af_rows, ac_rows = count_poly(prom, load_mei_alu())
    prom["poly_alu"] = pc
    prom["fixed_density"] = prom["fixed_alu"] / prom["width_kb"]
    prom["poly_density"] = prom["poly_alu"] / prom["width_kb"]

    # continuous correlations
    print("\n=== Continuous (Spearman vs constraint) ===")
    for lab, x in [("pLI", prom["lof.pLI"]), ("LOEUF(lof.oe)", prom["lof.oe"])]:
        rf, pf = stats.spearmanr(x, prom["fixed_density"], nan_policy="omit")
        rp, pp = stats.spearmanr(x, prom["poly_density"], nan_policy="omit")
        sign = "+" if lab == "pLI" else "-"  # pLI up=constrained; LOEUF down=constrained
        print(f"  {lab:14s} vs FIXED: r={rf:+.3f} (p={pf:.1e})   vs POLY: r={rp:+.3f} (p={pp:.1e})"
              f"   [constraint {'increases' if lab=='pLI' else 'decreases'} along axis]")

    # pLI strata
    prom["pLI_bin"] = pd.cut(prom["lof.pLI"], [0, .1, .5, .9, 1.0001],
                             labels=["<0.1", "0.1-0.5", "0.5-0.9", ">0.9"], include_lowest=True)
    rows = []
    base_fixed = base_poly = None
    for lab in ["<0.1", "0.1-0.5", "0.5-0.9", ">0.9"]:
        sub = prom[prom["pLI_bin"] == lab]
        ii = sub.index
        afs = [a for i in ii for a in af_rows[i]]
        acs = [a for i in ii for a in ac_rows[i]]
        mf, mp = sub["fixed_density"].mean(), sub["poly_density"].mean()
        if base_fixed is None:
            base_fixed, base_poly = mf, mp
        rows.append(dict(pLI_bin=lab, n_genes=len(sub),
                         fixed_density=round(mf, 3), poly_density=round(mp, 4),
                         fixed_rel=round(mf / base_fixed, 3), poly_rel=round(mp / base_poly, 3),
                         n_promoter_MEI=len(afs),
                         singleton_frac=round(np.mean(np.array(acs) == 1), 3) if acs else np.nan,
                         median_AF=f"{np.median(afs):.2e}" if afs else "NA"))
    strata = pd.DataFrame(rows)
    strata.to_csv(OUT / "probe_dosage_strata.csv", index=False)
    pd.set_option("display.width", 170); pd.set_option("display.max_columns", 20)
    print("\n=== Decomposition across pLI strata (baseline = least-constrained <0.1) ===")
    print(strata.to_string(index=False))
    print("\nREAD: selection scales with dosage IF  fixed_rel falls  &  poly_rel ~1  &  singleton_frac rises.")


if __name__ == "__main__":
    main()
