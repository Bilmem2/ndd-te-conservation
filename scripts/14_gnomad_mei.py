#!/usr/bin/env python3
"""
14_gnomad_mei.py — Population-level test of purifying selection against
polymorphic mobile-element insertions (MEIs) at NDD gene promoters (hg38).

Complements the cross-species (fixed/ancient TE density) analysis with a
within-species population-genetic layer using gnomAD v4.1 SV MEI calls.

Core logic — separating insertion TARGETING bias from post-insertion SELECTION
(the two alternatives the reviewer asked us to disentangle):

  * FIXED Alu density (RepeatMasker, ancient+fixed): strongly depleted at NDD
    promoters (main analysis, r = -0.345).
  * POLYMORPHIC Alu density (gnomAD SV, currently segregating): if the deficit
    were caused by biased *targeting* away from NDD promoters, recent insertions
    would be depleted too.  If instead it is caused by purifying *selection*,
    recent insertions should land at ~equal density but be held at lower allele
    frequency (excess of rare variants / singletons).

  Test 1 (density):    per-promoter polymorphic MEI density, one-sided
                       Mann-Whitney U (HK > NDD), rank-biserial r.
  Test 2 (SFS/select): among insertions that fall in promoters, compare the
                       allele-frequency spectrum of NDD- vs HK-promoter
                       insertions (median AF, singleton fraction, fraction rare)
                       with a one-sided MWU (NDD AF < HK) and a one-sided
                       Fisher test on singleton counts (NDD > HK).
  Decomposition:       NDD/HK density ratio for FIXED vs POLYMORPHIC insertions.

Inputs:
  results/hg38/{HighConfNDD,Housekeeping}_Alu.bed   (promoter windows + fixed Alu count)
  data/hg38/gnomad_mei.tsv   (chrom,start,end,name,svtype,AN,AC,AF; INS:ME:* rows of
                              gnomad.v4.1.sv.sites.bed.gz)
Outputs:
  results/gnomad_mei/{mei_depletion.csv, mei_af_selection.csv,
                      decomposition.csv, promoter_mei_counts.tsv}
No bedtools dependency (numpy interval overlap).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
DATA = REPO / "data"
OUT = RESULTS / "gnomad_mei"
OUT.mkdir(parents=True, exist_ok=True)

# INS:ME classes, plus a pooled SINE-related group (Alu + SVA are both SINE-derived)
MEI_CLASSES = {"ALU": ["INS:ME:ALU"], "LINE1": ["INS:ME:LINE1"],
               "SVA": ["INS:ME:SVA"], "SINE(ALU+SVA)": ["INS:ME:ALU", "INS:ME:SVA"]}
RARE_AF = 1e-3  # "rare" threshold for fraction-rare metric


def load_promoters(cat):
    f = RESULTS / "hg38" / f"{cat}_Alu.bed"
    df = pd.read_csv(f, sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand", "fixed_alu"])
    df = df[["chrom", "start", "end", "gene", "fixed_alu"]].copy()
    df["category"] = cat
    return df


def load_mei():
    df = pd.read_csv(DATA / "hg38" / "gnomad_mei.tsv", sep="\t")
    df = df.rename(columns={"#chrom": "chrom"})
    df = df[df["svtype"].str.startswith("INS:ME:")].copy()
    df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
    df["AC"] = pd.to_numeric(df["AC"], errors="coerce")
    return df.dropna(subset=["AF"])


def overlap_counts_and_af(prom, mei):
    """counts aligned to prom rows; af_df of every (gene,category,AF,AC) overlap."""
    counts = np.zeros(len(prom), dtype=int)
    af_records = []
    for chrom, m in mei.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if len(idx) == 0:
            continue
        order = np.argsort(m["start"].values)
        ms_s = m["start"].values[order]
        me_s = m["end"].values[order]
        maf_s = m["AF"].values[order]
        mac_s = m["AC"].values[order]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(ms_s, we, side="left")
            if lo == 0:
                continue
            hit = np.where(me_s[:lo] > ws)[0]
            counts[i] = len(hit)
            for h in hit:
                af_records.append((prom.at[i, "gene"], prom.at[i, "category"],
                                   maf_s[h], mac_s[h]))
    af_df = pd.DataFrame(af_records, columns=["gene", "category", "AF", "AC"])
    return counts, af_df


def rb(a, b):
    """rank-biserial for MWU(a=HK,b=NDD); negative => lower in NDD."""
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def main():
    prom = pd.concat([load_promoters("HighConfNDD"), load_promoters("Housekeeping")],
                     ignore_index=True)
    prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
    is_ndd = prom["category"].values == "HighConfNDD"
    is_hk = ~is_ndd
    mei_all = load_mei()

    dep_rows, af_rows = [], []
    counts_out = prom[["chrom", "start", "end", "gene", "category", "width_kb", "fixed_alu"]].copy()

    for cls, svts in MEI_CLASSES.items():
        mei = mei_all[mei_all["svtype"].isin(svts)]
        counts, af_df = overlap_counts_and_af(prom, mei)
        if cls in ("ALU", "LINE1", "SVA"):
            counts_out[f"{cls}_count"] = counts
        dens = counts / prom["width_kb"].values
        ndd_d, hk_d = dens[is_ndd], dens[is_hk]
        _, p = stats.mannwhitneyu(hk_d, ndd_d, alternative="greater")
        dep_rows.append(dict(
            TE=cls, n_HK=int(is_hk.sum()), n_NDD=int(is_ndd.sum()),
            ins_in_HK=int(counts[is_hk].sum()), ins_in_NDD=int(counts[is_ndd].sum()),
            dens_HK=round(float(hk_d.mean()), 4), dens_NDD=round(float(ndd_d.mean()), 4),
            ratio_NDD_HK=round(float(ndd_d.mean() / hk_d.mean()), 3) if hk_d.mean() else np.nan,
            p_depletion=p, r=round(rb(hk_d, ndd_d), 3)))

        # Test 2 — allele-frequency spectrum among promoter insertions
        a_ndd = af_df[af_df["category"] == "HighConfNDD"]
        a_hk = af_df[af_df["category"] == "Housekeeping"]
        if len(a_ndd) >= 10 and len(a_hk) >= 10:
            _, p_af = stats.mannwhitneyu(a_ndd["AF"], a_hk["AF"], alternative="less")
            s_ndd, s_hk = int((a_ndd["AC"] == 1).sum()), int((a_hk["AC"] == 1).sum())
            # one-sided Fisher: NDD singleton fraction > HK
            _, p_single = stats.fisher_exact(
                [[s_ndd, len(a_ndd) - s_ndd], [s_hk, len(a_hk) - s_hk]], alternative="greater")
            af_rows.append(dict(
                TE=cls, n_ins_NDD=len(a_ndd), n_ins_HK=len(a_hk),
                med_AF_NDD=f"{a_ndd['AF'].median():.2e}", med_AF_HK=f"{a_hk['AF'].median():.2e}",
                singleton_NDD=round(s_ndd / len(a_ndd), 3), singleton_HK=round(s_hk / len(a_hk), 3),
                rare_NDD=round(float((a_ndd["AF"] < RARE_AF).mean()), 3),
                rare_HK=round(float((a_hk["AF"] < RARE_AF).mean()), 3),
                p_AF_lower=round(p_af, 4), p_singleton_higher=round(p_single, 4)))

    # Decomposition: fixed vs polymorphic Alu density ratio (targeting vs selection)
    fixed_dens = prom["fixed_alu"].values / prom["width_kb"].values
    poly_alu = counts_out["ALU_count"].values / prom["width_kb"].values

    def ratio_ci(dens, n_boot=10000, seed=0):
        """Bootstrap 95% CI of NDD/HK mean-density ratio (resample promoters within group)."""
        rng = np.random.default_rng(seed)
        h_idx = np.where(is_hk)[0]
        n_idx = np.where(is_ndd)[0]
        boots = np.empty(n_boot)
        for b in range(n_boot):
            h = dens[rng.choice(h_idx, len(h_idx), replace=True)].mean()
            n = dens[rng.choice(n_idx, len(n_idx), replace=True)].mean()
            boots[b] = n / h if h else np.nan
        return np.nanpercentile(boots, [2.5, 97.5]), boots

    (f_lo, f_hi), f_boot = ratio_ci(fixed_dens, seed=1)
    (p_lo, p_hi), p_boot = ratio_ci(poly_alu, seed=2)
    p_diff = float(np.mean(p_boot > f_boot))  # P(polymorphic ratio > fixed ratio)
    decomp = pd.DataFrame([
        dict(insertion_class="FIXED Alu (RepeatMasker)",
             dens_HK=round(fixed_dens[is_hk].mean(), 3), dens_NDD=round(fixed_dens[is_ndd].mean(), 3),
             ratio_NDD_HK=round(fixed_dens[is_ndd].mean() / fixed_dens[is_hk].mean(), 3),
             ci95=f"[{f_lo:.3f}, {f_hi:.3f}]"),
        dict(insertion_class="POLYMORPHIC Alu (gnomAD SV)",
             dens_HK=round(poly_alu[is_hk].mean(), 4), dens_NDD=round(poly_alu[is_ndd].mean(), 4),
             ratio_NDD_HK=round(poly_alu[is_ndd].mean() / poly_alu[is_hk].mean(), 3),
             ci95=f"[{p_lo:.3f}, {p_hi:.3f}]"),
    ])

    dep = pd.DataFrame(dep_rows)
    af = pd.DataFrame(af_rows)
    counts_out.to_csv(OUT / "promoter_mei_counts.tsv", sep="\t", index=False)
    dep.to_csv(OUT / "mei_depletion.csv", index=False)
    af.to_csv(OUT / "mei_af_selection.csv", index=False)
    decomp.to_csv(OUT / "decomposition.csv", index=False)

    pd.set_option("display.width", 170); pd.set_option("display.max_columns", 25)
    print("\n=== TEST 1: polymorphic MEI density, NDD vs HK (one-sided HK>NDD) ===")
    print(dep.to_string(index=False))
    print("\n=== TEST 2: allele-frequency spectrum of promoter insertions ===")
    print(af.to_string(index=False) if len(af) else "  (insufficient)")
    print("\n=== DECOMPOSITION: targeting (polymorphic) vs selection (fixed) ===")
    print(decomp.to_string(index=False))
    print(f"\nBootstrap P(polymorphic ratio > fixed ratio) = {p_diff:.4f}  (10,000 reps)")
    print("Interpretation: FIXED ratio << 1 but POLYMORPHIC ratio ~ 1 (CIs separated) => deficit")
    print("arises AFTER insertion (purifying selection), not from biased targeting.")


if __name__ == "__main__":
    main()
