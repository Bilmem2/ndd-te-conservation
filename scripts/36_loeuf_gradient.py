"""
36_loeuf_gradient.py — promoter Alu density across the loss-of-function
constraint spectrum.

Tests whether promoter Alu content varies with gnomAD LOEUF across
protein-coding genes, and whether NDD promoters remain Alu-poor at fixed
constraint:

  * Alu density at TSS +/-2 kb for every protein-coding promoter in hg38
  * genes binned into LOEUF deciles (decile 1 = most constrained)
  * Spearman correlation across genes and across the ten decile means
  * NDD versus all other genes within each decile

Inputs : data/hg38/gtf/gencode.v47.gtf.gz, data/hg38/alu_rmsk.bed,
         data/hg38/gnomad_constraint.tsv, data/gene_lists/HighConfNDD_genes.txt
Output : results/matched/loeuf_gradient.csv
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parent.parent
DATA, GL, OUT = ROOT / "data", ROOT / "data" / "gene_lists", ROOT / "results" / "matched"
OUT.mkdir(parents=True, exist_ok=True)

# promoter construction and Alu counting are shared with the genome-baseline script
spec = importlib.util.spec_from_file_location("gb", ROOT / "scripts" / "25_genome_baseline.py")
gb = importlib.util.module_from_spec(spec)
sys.modules["gb"] = gb
spec.loader.exec_module(gb)

WINDOW_KB = 4.0

print("building all protein-coding promoters (hg38) ...")
prom = gb.build_all_pc_promoters(DATA / "hg38" / "gtf" / "gencode.v47.gtf.gz",
                                 chrom_filter=None, add_chr=False, upper=False)
print(f"  {len(prom):,} promoters")

print("counting Alu ...")
alu = gb.load_alu(DATA / "hg38" / "alu_rmsk.bed", None)
prom["alu"] = gb.count_alu(prom, alu)
prom["density"] = prom["alu"] / WINDOW_KB

print("merging gnomAD LOEUF ...")
con = (pd.read_csv(DATA / "hg38" / "gnomad_constraint.tsv", sep="\t",
                   usecols=["gene", "lof.oe_ci.upper"])
       .rename(columns={"lof.oe_ci.upper": "loeuf"})
       .dropna(subset=["loeuf"])
       .groupby("gene", as_index=False)["loeuf"].min())

df = prom.merge(con, on="gene", how="inner")
ndd = {l.strip() for l in open(GL / "HighConfNDD_genes.txt") if l.strip()}
df["ndd"] = df["gene"].isin(ndd)
print(f"  {len(df):,} promoters with LOEUF; {int(df['ndd'].sum())} are NDD genes")

# ---- the gradient -------------------------------------------------------
df["decile"] = pd.qcut(df["loeuf"], 10, labels=False) + 1   # 1 = most constrained

rows = []
for d, g in df.groupby("decile"):
    nd, rest = g[g.ndd], g[~g.ndd]
    row = dict(decile=int(d), n=len(g),
               loeuf_median=round(g.loeuf.median(), 3),
               alu_mean=round(g.density.mean(), 4),
               alu_median=round(g.density.median(), 4),
               n_ndd=len(nd),
               alu_mean_ndd=round(nd.density.mean(), 4) if len(nd) else np.nan,
               alu_mean_other=round(rest.density.mean(), 4))
    if len(nd) >= 20:
        u = mannwhitneyu(rest.density, nd.density, alternative="greater")
        row["p_ndd_vs_rest"] = u.pvalue
        row["r_ndd_vs_rest"] = round(-(2 * u.statistic / (len(rest) * len(nd)) - 1), 3)
    rows.append(row)

grad = pd.DataFrame(rows)
grad.to_csv(OUT / "loeuf_gradient.csv", index=False)

rho, p_rho = spearmanr(df.loeuf, df.density)
rho_dec, p_dec = spearmanr(grad.decile, grad.alu_mean)

print("\n=== Alu density by LOEUF decile (1 = most constrained) ===")
print(grad.to_string(index=False))

print(f"\nacross all {len(df):,} genes : Spearman rho(LOEUF, Alu density) = {rho:+.3f}  p = {p_rho:.3e}")
print(f"across the 10 decile means : Spearman rho = {rho_dec:+.3f}  p = {p_dec:.3e}")

# where do NDD genes sit?
pct = 100 * (df.loeuf.rank(pct=True)[df.ndd]).mean()
print(f"\nmean LOEUF percentile of NDD genes: {pct:.1f} (lower = more constrained)")

# the key test: are NDD genes still depleted inside the most-constrained tenth?
d1 = df[df.decile == 1]
u = mannwhitneyu(d1[~d1.ndd].density, d1[d1.ndd].density, alternative="greater")
r1 = -(2 * u.statistic / ((~d1.ndd).sum() * d1.ndd.sum()) - 1)
print(f"\nwithin decile 1 (most constrained): NDD n={int(d1.ndd.sum())}, other n={int((~d1.ndd).sum())}")
print(f"   Alu density  NDD {d1[d1.ndd].density.mean():.4f}  vs other {d1[~d1.ndd].density.mean():.4f}")
print(f"   r = {r1:+.3f}   p = {u.pvalue:.3e}")
print(f"\nwrote {OUT / 'loeuf_gradient.csv'}")
