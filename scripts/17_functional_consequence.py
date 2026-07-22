#!/usr/bin/env python3
"""
17_functional_consequence.py — Link promoter Alu occupancy to gene-level
functional properties (hg38), addressing the reviewer's request for the
"regulatory or biological consequences" of the depletion.

Two public-data readouts, reported honestly (including confounds):

  A. DOSAGE SENSITIVITY (gnomAD v4.1 constraint).  Do genes under stronger
     dosage constraint carry fewer promoter Alu?  If Alu exclusion tracks
     functional importance, promoter Alu density should be inversely related
     to pLI / directly related to LoF o/e.  This connects the evolutionary
     signal to a gene-level fitness consequence (haploinsufficiency intolerance)
     and to the paper's "extreme of a broader constraint spectrum" framing.

  B. BRAIN EXPRESSION (GTEx v8).  Do genes with a promoter Alu differ in brain
     expression from those without?  NOTE the built-in confound: active
     promoters recruit Alu (open-chromatin insertion bias), so a naive genome-
     wide test is expected to be positive/neutral and does NOT test repression.
     We therefore test WITHIN each gene set and report the result as-is,
     without forcing a causal reading.

Inputs:
  results/hg38/{HighConfNDD,Housekeeping}_Alu.bed
  data/hg38/gnomad_constraint.tsv          (gene, lof.pLI, lof.oe, mis.z_score, mane_select)
  data/hg38/gtex_gene_median_tpm.gct.gz    (GTEx v8 median TPM by tissue)
Outputs:
  results/functional/{constraint_vs_alu.csv, constraint_by_pLI_bin.csv,
                      expression_vs_alu.csv}
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
DATA = REPO / "data"
OUT = RESULTS / "functional"
OUT.mkdir(parents=True, exist_ok=True)


def load_alu():
    frames = []
    for cat in ["HighConfNDD", "Housekeeping"]:
        f = RESULTS / "hg38" / f"{cat}_Alu.bed"
        df = pd.read_csv(f, sep="\t", header=None,
                         names=["chrom", "start", "end", "gene", "score", "strand", "fixed_alu"])
        df["category"] = cat
        df["alu_density"] = df["fixed_alu"] / ((df["end"] - df["start"]) / 1000)
        frames.append(df[["gene", "category", "fixed_alu", "alu_density"]])
    return pd.concat(frames, ignore_index=True)


def load_constraint():
    c = pd.read_csv(DATA / "hg38" / "gnomad_constraint.tsv", sep="\t",
                    usecols=["gene", "mane_select", "lof.pLI", "lof.oe", "mis.z_score"],
                    dtype={"mane_select": str})
    c = c[c["mane_select"].str.lower() == "true"].copy()
    for col in ["lof.pLI", "lof.oe", "mis.z_score"]:
        c[col] = pd.to_numeric(c[col], errors="coerce")
    return c.dropna(subset=["lof.pLI"]).drop_duplicates("gene")[
        ["gene", "lof.pLI", "lof.oe", "mis.z_score"]]


def load_brain_expression():
    g = pd.read_csv(DATA / "hg38" / "gtex_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
    brain_cols = [c for c in g.columns if c.startswith("Brain -")]
    g["brain_tpm"] = g[brain_cols].median(axis=1)
    return g[["Description", "brain_tpm"]].rename(columns={"Description": "gene"}).drop_duplicates("gene")


def spearman_row(label, x, y):
    m = pd.notna(x) & pd.notna(y)
    r, p = stats.spearmanr(x[m], y[m])
    return dict(subset=label, n=int(m.sum()), spearman_r=round(r, 3), p_value=p)


def main():
    alu = load_alu()
    con = load_constraint()
    expr = load_brain_expression()

    # ---- A. constraint vs Alu ----
    d = alu.merge(con, on="gene", how="inner")
    rows = [spearman_row("ALL", d["lof.pLI"], d["alu_density"]),
            spearman_row("HighConfNDD", d[d.category == "HighConfNDD"]["lof.pLI"],
                         d[d.category == "HighConfNDD"]["alu_density"]),
            spearman_row("Housekeeping", d[d.category == "Housekeeping"]["lof.pLI"],
                         d[d.category == "Housekeeping"]["alu_density"])]
    # robustness with LoF o/e (lower=more constrained) and missense Z (higher=more constrained)
    rows.append(spearman_row("ALL_lof.oe_vs_alu", d["lof.oe"], d["alu_density"]))
    rows.append(spearman_row("ALL_mis.z_vs_alu", d["mis.z_score"], d["alu_density"]))
    con_vs = pd.DataFrame(rows)

    # pLI bins -> median Alu density (monotone trend check), pooled across all genes
    d["pLI_bin"] = pd.cut(d["lof.pLI"], [0, .1, .5, .9, 1.0001],
                          labels=["<0.1 (tolerant)", "0.1-0.5", "0.5-0.9", ">0.9 (constrained)"],
                          include_lowest=True)
    binned = d.groupby("pLI_bin", observed=True).agg(
        n=("alu_density", "size"),
        median_alu_density=("alu_density", "median"),
        mean_alu_density=("alu_density", "mean")).reset_index()

    # ---- B. brain expression vs Alu ----
    e = alu.merge(expr, on="gene", how="inner")
    e["log_tpm"] = np.log10(e["brain_tpm"] + 1)
    e["has_alu"] = e["fixed_alu"] > 0
    erows = []
    for cat in ["HighConfNDD", "Housekeeping", "ALL"]:
        s = e if cat == "ALL" else e[e.category == cat]
        with_a = s[s.has_alu]["brain_tpm"]
        no_a = s[~s.has_alu]["brain_tpm"]
        u, p = stats.mannwhitneyu(with_a, no_a, alternative="two-sided") if len(with_a) and len(no_a) else (np.nan, np.nan)
        rr, pp = stats.spearmanr(s["alu_density"], s["log_tpm"])
        erows.append(dict(subset=cat, n=len(s),
                          n_with_alu=len(with_a), n_no_alu=len(no_a),
                          median_tpm_withAlu=round(with_a.median(), 2) if len(with_a) else np.nan,
                          median_tpm_noAlu=round(no_a.median(), 2) if len(no_a) else np.nan,
                          p_withVSno=p, spearman_aluDensity_logTPM=round(rr, 3), p_spearman=pp))
    expr_vs = pd.DataFrame(erows)

    con_vs.to_csv(OUT / "constraint_vs_alu.csv", index=False)
    binned.to_csv(OUT / "constraint_by_pLI_bin.csv", index=False)
    expr_vs.to_csv(OUT / "expression_vs_alu.csv", index=False)

    pd.set_option("display.width", 180); pd.set_option("display.max_columns", 25)
    print("\n=== A. Dosage constraint (gnomAD pLI) vs promoter Alu density (Spearman) ===")
    print(con_vs.to_string(index=False))
    print("\n--- promoter Alu density by pLI bin (all merged genes) ---")
    print(binned.to_string(index=False))
    print("\n=== B. Brain expression (GTEx) vs promoter Alu  [confounded — see header] ===")
    print(expr_vs.to_string(index=False))


if __name__ == "__main__":
    main()
