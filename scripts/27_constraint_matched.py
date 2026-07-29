#!/usr/bin/env python3
"""
27_constraint_matched.py — is the promoter SINE deficit specific to NDD genes,
or a property of constrained, brain-expressed genes in general?

The housekeeping contrast establishes that NDD promoters are Alu-poor relative
to a non-disease control, and the genome baseline shows the deficit is not an
artefact of that control. Neither addresses the harder question: NDD genes are
overwhelmingly haploinsufficiency-intolerant regulators expressed in the
developing brain, so a control that does not match those properties cannot
separate "NDD-associated" from "constrained and brain-expressed".

Here each NDD promoter is matched 1:1, by nearest-neighbour Mahalanobis
distance, to a non-NDD protein-coding promoter with similar
  * LOEUF          (gnomAD loss-of-function observed/expected upper bound)
  * brain expression (max GTEx brain median TPM, log10)
  * promoter GC content
and Alu density is then compared between the two. If the deficit disappears
against this control, the pattern belongs to constrained developmental genes
generally; if it persists, it is specific to NDD-associated loci beyond what
constraint and brain expression explain.

Inputs : data/hg38/{gnomad_constraint.tsv, gtex_gene_median_tpm.gct.gz,
         hg38.2bit, alu_rmsk.bed, gtf/gencode.v47.gtf.gz}, data/gene_lists/
Output : results/matched/constraint_matched.csv
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
COVS = ["loeuf", "brain_log", "gc"]
N_BOOT = 4000


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


def smd(a, b):
    sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    return (np.mean(a) - np.mean(b)) / sd if sd else np.nan


def match(nd, pool, cols):
    X = pd.concat([nd[cols], pool[cols]]).values
    VI = np.linalg.pinv(np.cov(X, rowvar=False))
    ndX, poX = nd[cols].values, pool[cols].values
    used = np.zeros(len(pool), bool)
    chosen = np.full(len(nd), -1, int)
    for i in range(len(nd)):
        d = poX - ndX[i]
        d2 = np.einsum("ij,jk,ik->i", d, VI, d)
        d2[used] = np.inf
        j = int(np.argmin(d2))
        if np.isfinite(d2[j]):
            chosen[i] = j
            used[j] = True
    return chosen


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    ndd, hk = gene_set("HighConfNDD_genes.txt"), gene_set("Housekeeping_genes.txt")

    con = pd.read_csv(DATA / "gnomad_constraint.tsv", sep="\t",
                      usecols=["gene", "mane_select", "lof.oe_ci.upper", "lof.pLI"])
    con = con[con["mane_select"] == True]
    con = (con.groupby("gene")[["lof.oe_ci.upper", "lof.pLI"]].mean()
           .rename(columns={"lof.oe_ci.upper": "loeuf", "lof.pLI": "pli"}).reset_index())

    gx = pd.read_csv(DATA / "gtex_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
    brain_cols = [c for c in gx.columns if c.startswith("Brain")]
    gx = pd.DataFrame({"gene": gx["Description"],
                       "brain_tpm": gx[brain_cols].max(axis=1)}).drop_duplicates("gene")

    prom = promoters()
    alu = pd.read_csv(DATA / "alu_rmsk.bed", sep="\t", header=None,
                      usecols=[0, 1, 2], names=["chrom", "start", "end"])
    alu = alu[alu["chrom"].isin(MAIN)]
    prom["alu_d"] = count_alu(prom, alu) / ((prom["end"] - prom["start"]) / 1000.0)
    prom["gc"] = gc_content(prom)
    prom = prom.merge(con, on="gene", how="left").merge(gx, on="gene", how="left")
    prom["brain_log"] = np.log10(prom["brain_tpm"] + 1)
    genome_mean = prom["alu_d"].mean()

    usable = prom.dropna(subset=COVS).reset_index(drop=True)
    nd = usable[usable["gene"].isin(ndd)].reset_index(drop=True)
    pool = usable[~usable["gene"].isin(ndd)].reset_index(drop=True)
    print(f"promoters with LOEUF + GTEx + GC: {len(usable)}  "
          f"(NDD {len(nd)}, pool {len(pool)})")

    chosen = match(nd, pool, COVS)
    ok = chosen >= 0
    ctrl = pool.iloc[chosen[ok]].reset_index(drop=True)
    ndm = nd[ok].reset_index(drop=True)

    print("\n=== Covariate balance (|SMD| < 0.1 = well matched) ===")
    bal = []
    for c in COVS + ["pli"]:
        if c not in ndm:
            continue
        bal.append(dict(covariate=c, NDD_mean=round(ndm[c].mean(), 4),
                        pool_mean=round(pool[c].mean(), 4),
                        matched_mean=round(ctrl[c].mean(), 4),
                        SMD_before=round(smd(ndm[c].values, pool[c].dropna().values), 3),
                        SMD_after=round(smd(ndm[c].values, ctrl[c].values), 3)))
    bal = pd.DataFrame(bal)
    print(bal.to_string(index=False))

    # primary test: NDD vs its own matched controls (not mediated by the genome)
    _, p = stats.mannwhitneyu(ctrl["alu_d"], ndm["alu_d"], alternative="greater")
    r = rb(ctrl["alu_d"].values, ndm["alu_d"].values)
    n_hk_in_ctrl = int(ctrl["gene"].isin(hk).sum())

    # is the matched control itself below the genome baseline?
    ctrl_genes = set(ctrl["gene"])
    rest = prom["alu_d"].values[~prom["gene"].isin(ctrl_genes).values]
    rest = rest[~np.isnan(rest)]
    _, p_ctrl_vs_genome = stats.mannwhitneyu(rest, ctrl["alu_d"].values,
                                             alternative="greater")

    def boot_ci(vals, seed):
        rng = np.random.default_rng(seed)
        allv = prom["alu_d"].dropna().values
        b = np.empty(N_BOOT)
        for k in range(N_BOOT):
            b[k] = (rng.choice(vals, len(vals), True).mean() /
                    rng.choice(allv, len(allv), True).mean())
        return np.nanpercentile(b, [2.5, 97.5])

    ci_nd = boot_ci(ndm["alu_d"].values, 11)
    ci_ct = boot_ci(ctrl["alu_d"].values, 12)

    rows = [dict(comparison="NDD vs constraint+brain+GC matched control",
                 n_NDD=len(ndm), n_control=len(ctrl),
                 median_alu_NDD=round(float(ndm["alu_d"].median()), 3),
                 median_alu_control=round(float(ctrl["alu_d"].median()), 3),
                 mean_alu_NDD=round(float(ndm["alu_d"].mean()), 3),
                 mean_alu_control=round(float(ctrl["alu_d"].mean()), 3),
                 NDD_over_genome=round(float(ndm["alu_d"].mean()) / genome_mean, 3),
                 NDD_genome_CI=f"[{ci_nd[0]:.3f}, {ci_nd[1]:.3f}]",
                 control_over_genome=round(float(ctrl["alu_d"].mean()) / genome_mean, 3),
                 control_genome_CI=f"[{ci_ct[0]:.3f}, {ci_ct[1]:.3f}]",
                 p_control_below_genome=p_ctrl_vs_genome,
                 r=round(r, 3), p_value=p,
                 housekeeping_genes_in_control=n_hk_in_ctrl)]
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "constraint_matched.csv", index=False)
    bal.to_csv(OUT / "constraint_matched_balance.csv", index=False)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n=== Alu depletion against a constraint- and brain-expression-matched control ===")
    print(res.to_string(index=False))
    print(f"\nGenome baseline {genome_mean:.3f} Alu/kb; "
          f"{n_hk_in_ctrl} of {len(ctrl)} matched controls are housekeeping genes.")
    print(f"\nWrote {OUT / 'constraint_matched.csv'}")


if __name__ == "__main__":
    main()
