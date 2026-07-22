#!/usr/bin/env python3
"""
19_rebaseline.py — DECISIVE honest re-baselining of the NDD Alu decomposition.

The submitted paper contrasted NDD vs brain-expressed HOUSEKEEPING promoters.
But housekeeping/broadly-expressed genes are Alu-ENRICHED (Eller et al. 2007),
so NDD/HK = 0.57 conflates "housekeeping enrichment" with "NDD depletion".

Here we re-baseline against honest references:
  (1) GENOME baseline  — all protein-coding gene promoters.
  (2) MATCHED control  — non-disease genes matched to NDD on EXPRESSION BREADTH
                         (# GTEx tissues with TPM>=1) and promoter GC, the two
                         properties that most strongly predict Alu content.

For each reference we compute the targeting-vs-selection decomposition:
  * FIXED Alu density   (RepeatMasker)          -> post-selection
  * POLYMORPHIC Alu     (gnomAD v4.1 INS:ME:ALU) -> pre-selection targeting proxy
  * AFS of promoter MEIs (singleton frac, median AF)
Question: does NDD-specific depletion + the poly~flat / fixed~down signature
SURVIVE once housekeeping-enrichment / expression-breadth is removed?

Inputs (local): gencode.v47.gtf.gz, alu_rmsk.bed, gnomad_mei.tsv,
gtex_gene_median_tpm.gct.gz, hg38.2bit, data/gene_lists/*.txt
Output: results/matched/rebaseline_decomposition.csv + console
"""
import gzip, re
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
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}
NAME_RE = re.compile(r'gene_name "([^"]+)"')


def gene_list(fn):
    return set(l.strip() for l in open(GL / fn) if l.strip())


def build_promoters(genes_wanted):
    rows, seen = [], set()
    with gzip.open(DATA / "gtf" / "gencode.v47.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or 'gene_type "protein_coding"' not in f[8] or f[0] not in MAIN:
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) not in genes_wanted or m.group(1) in seen:
                continue
            seen.add(m.group(1))
            start, end, strand = int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            rows.append((f[0], max(0, tss - 2000), tss + 2000, m.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def expression_breadth():
    g = pd.read_csv(DATA / "gtex_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
    tissues = [c for c in g.columns if c not in ("Name", "Description")]
    br = (g[tissues] >= 1).sum(axis=1)
    return pd.DataFrame({"gene": g["Description"], "breadth": br}).drop_duplicates("gene")


def gc_of(prom):
    tb = TwoBitFile(str(DATA / "hg38.2bit"))
    out = np.full(len(prom), np.nan)
    for i, r in enumerate(prom.itertuples()):
        try:
            s = tb[r.chrom][int(r.start):int(r.end)].upper()
        except Exception:
            continue
        gc = s.count("G") + s.count("C"); at = s.count("A") + s.count("T")
        if gc + at:
            out[i] = gc / (gc + at)
    return out


def count_fixed(prom, feats):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in feats.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        s = np.sort(m["start"].values); e = m["end"].values[np.argsort(m["start"].values)]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def count_poly(prom, mei):
    counts = np.zeros(len(prom), dtype=int)
    afr = [[] for _ in range(len(prom))]; acr = [[] for _ in range(len(prom))]
    for chrom, m in mei.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        o = np.argsort(m["start"].values)
        s = m["start"].values[o]; e = m["end"].values[o]; af = m["AF"].values[o]; ac = m["AC"].values[o]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if not lo:
                continue
            hit = np.where(e[:lo] > ws)[0]
            counts[i] = len(hit)
            for h in hit:
                afr[i].append(af[h]); acr[i].append(ac[h])
    return counts, afr, acr


def smd(a, b):
    sd = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    return (np.mean(a) - np.mean(b)) / sd if sd else np.nan


def match(nd, pool, cols):
    X = pd.concat([nd[cols], pool[cols]]).values
    VI = np.linalg.pinv(np.cov(X, rowvar=False))
    ndX, poX = nd[cols].values, pool[cols].values
    used = np.zeros(len(pool), bool); chosen = np.full(len(nd), -1, int)
    for i in range(len(nd)):
        d = poX - ndX[i]
        d2 = np.einsum("ij,jk,ik->i", d, VI, d); d2[used] = np.inf
        j = int(np.argmin(d2))
        if np.isfinite(d2[j]):
            chosen[i] = j; used[j] = True
    return chosen


def ratio_ci(dens, ndmask, comask, seed):
    rng = np.random.default_rng(seed)
    ni, ci = np.where(ndmask)[0], np.where(comask)[0]
    b = np.empty(4000)
    for k in range(4000):
        n = dens[rng.choice(ni, len(ni), True)].mean(); c = dens[rng.choice(ci, len(ci), True)].mean()
        b[k] = n / c if c else np.nan
    return np.nanpercentile(b, [2.5, 97.5])


def main():
    ndd = gene_list("HighConfNDD_genes.txt")
    disease = ndd | gene_list("Cardiovascular_genes.txt") | gene_list("Mendelian_genes.txt")

    prom = build_promoters(_all_pc_genes())
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    prom = prom.merge(expression_breadth(), on="gene", how="left")
    prom["gc"] = gc_of(prom)
    # finalize row set BEFORE overlap counting so af/ac lists stay aligned
    prom = prom.dropna(subset=["breadth", "gc"]).reset_index(drop=True)
    prom["fixed_alu"] = count_fixed(prom, _load_alu())
    pc, afr, acr = count_poly(prom, _load_mei())
    prom["poly_alu"] = pc
    prom["fixed_d"] = prom["fixed_alu"] / prom["width_kb"]
    prom["poly_d"] = prom["poly_alu"] / prom["width_kb"]
    prom["is_ndd"] = prom["gene"].isin(ndd)
    prom["is_disease"] = prom["gene"].isin(disease)
    gpos = {g: i for i, g in enumerate(prom["gene"].values)}  # gene -> row index for af/ac lookup
    print(f"Promoters: {len(prom)} ; NDD in set: {int(prom['is_ndd'].sum())}")

    nd = prom[prom["is_ndd"]].reset_index(drop=True)
    pool = prom[~prom["is_disease"]].reset_index(drop=True)  # non-disease control pool
    chosen = match(nd, pool, ["breadth", "gc"])
    ok = chosen >= 0
    ctrl = pool.iloc[chosen[ok]].reset_index(drop=True)
    ndm = nd[ok].reset_index(drop=True)

    # balance
    print("\n=== Covariate balance (NDD vs matched non-disease control) ===")
    for c in ["breadth", "gc"]:
        print(f"  {c:8s}  NDD={ndm[c].mean():.3f}  ctrl={ctrl[c].mean():.3f}  "
              f"SMD_after={smd(ndm[c].values, ctrl[c].values):+.3f}")

    # decomposition vs genome and vs matched control
    genome_fixed = prom["fixed_d"].mean(); genome_poly = prom["poly_d"].mean()
    rows = []
    for label, ref in [("GENOME (all PC genes)", None), ("MATCHED non-disease control", ctrl)]:
        if ref is None:
            rf, rp = ndm["fixed_d"].mean() / genome_fixed, ndm["poly_d"].mean() / genome_poly
            ci_f = ci_p = ("—", "—")
        else:
            rf = ndm["fixed_d"].mean() / ref["fixed_d"].mean()
            rp = ndm["poly_d"].mean() / ref["poly_d"].mean()
            # bootstrap needs joint array + masks
            allf = np.concatenate([ndm["fixed_d"].values, ref["fixed_d"].values])
            allp = np.concatenate([ndm["poly_d"].values, ref["poly_d"].values])
            m_nd = np.r_[np.ones(len(ndm), bool), np.zeros(len(ref), bool)]
            m_co = ~m_nd
            ci_f = tuple(round(x, 3) for x in ratio_ci(allf, m_nd, m_co, 1))
            ci_p = tuple(round(x, 3) for x in ratio_ci(allp, m_nd, m_co, 2))
        rows.append(dict(reference=label,
                         fixed_NDD=round(ndm["fixed_d"].mean(), 3),
                         fixed_ratio=round(rf, 3), fixed_CI=str(ci_f),
                         poly_NDD=round(ndm["poly_d"].mean(), 4),
                         poly_ratio=round(rp, 3), poly_CI=str(ci_p)))
    dec = pd.DataFrame(rows)

    # AFS: NDD promoter MEIs vs matched-control promoter MEIs
    def gather(sub):
        af, ac = [], []
        for g in sub["gene"]:
            i = gpos.get(g)
            if i is not None:
                af += afr[i]; ac += acr[i]
        return np.array(af), np.array(ac)
    af_nd, ac_nd = gather(ndm); af_co, ac_co = gather(ctrl)
    afs = dict(n_MEI_NDD=len(af_nd), n_MEI_ctrl=len(af_co),
               singleton_NDD=round(np.mean(ac_nd == 1), 3) if len(ac_nd) else np.nan,
               singleton_ctrl=round(np.mean(ac_co == 1), 3) if len(ac_co) else np.nan,
               medAF_NDD=f"{np.median(af_nd):.2e}" if len(af_nd) else "NA",
               medAF_ctrl=f"{np.median(af_co):.2e}" if len(af_co) else "NA")
    if len(af_nd) >= 10 and len(af_co) >= 10:
        afs["p_AF_lower"] = round(stats.mannwhitneyu(af_nd, af_co, alternative="less")[1], 4)
        _, afs["p_singleton_higher"] = stats.fisher_exact(
            [[int((ac_nd == 1).sum()), int((ac_nd != 1).sum())],
             [int((ac_co == 1).sum()), int((ac_co != 1).sum())]], alternative="greater")
        afs["p_singleton_higher"] = round(afs["p_singleton_higher"], 4)

    dec.to_csv(OUT / "rebaseline_decomposition.csv", index=False)
    pd.set_option("display.width", 180); pd.set_option("display.max_columns", 20)
    print(f"\nGenome baseline: fixed={genome_fixed:.3f}/kb  poly={genome_poly:.4f}/kb")
    print("\n=== NDD Alu decomposition vs honest references ===")
    print(dec.to_string(index=False))
    print("\n=== AFS (selection signature): NDD vs matched-control promoter MEIs ===")
    for k, v in afs.items():
        print(f"  {k}: {v}")
    print("\nREAD: NDD-specific selection survives IF fixed_ratio<1 (CI excludes 1) "
          "AND poly_ratio~1 AND NDD MEIs rarer/more singleton.")


def _all_pc_genes():
    genes = set()
    with gzip.open(DATA / "gtf" / "gencode.v47.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] == "gene" and 'gene_type "protein_coding"' in f[8] and f[0] in MAIN:
                m = NAME_RE.search(f[8])
                if m:
                    genes.add(m.group(1))
    return genes


def _load_alu():
    df = pd.read_csv(DATA / "alu_rmsk.bed", sep="\t", header=None, names=["chrom", "start", "end"])
    return df[df["chrom"].isin(MAIN)]


def _load_mei():
    df = pd.read_csv(DATA / "gnomad_mei.tsv", sep="\t").rename(columns={"#chrom": "chrom"})
    df = df[(df["svtype"] == "INS:ME:ALU") & (df["chrom"].isin(MAIN))].copy()
    df["AF"] = pd.to_numeric(df["AF"], errors="coerce"); df["AC"] = pd.to_numeric(df["AC"], errors="coerce")
    return df.dropna(subset=["AF"])


if __name__ == "__main__":
    main()
