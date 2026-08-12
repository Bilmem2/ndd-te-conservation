"""
44_mei_matched.py — polymorphic mobile-element insertions at NDD promoters,
tested against matched controls and with detection bias addressed directly.

The existing population-genetic analysis compares NDD promoters with the whole
housekeeping set. That contrast carries the same confounds as the cross-species
one: housekeeping promoters differ from NDD promoters in GC content, gene density
and recombination rate, and all three affect both the true insertion rate and the
probability that a short-read structural-variant caller detects an insertion.

This script repeats the comparison against the context-matched control set, adds
a GC-stratified test, compares the full allele-frequency spectrum rather than two
summary statistics, and separates elements by length so that the size dependence
the selection model predicts can be read directly.

Inputs : data/hg38/gnomad_mei.tsv
         results/hg38/{HighConfNDD,Housekeeping}_Alu.bed
         results/matched/matched_pairs.tsv
Output : results/gnomad_mei/matched_mei.csv
         results/gnomad_mei/afs_bins.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "gnomad_mei"
OUT.mkdir(parents=True, exist_ok=True)

CLASSES = {"INS:ME:ALU": "Alu", "INS:ME:SVA": "SVA", "INS:ME:LINE1": "LINE1"}
APPROX_BP = {"Alu": 300, "SVA": 2000, "LINE1": 6000}


def rank_biserial(hk, ndd):
    """One-sided Mann-Whitney U (HK > NDD) with rank-biserial r."""
    u = stats.mannwhitneyu(hk, ndd, alternative="greater")
    return 1 - (2 * u.statistic) / (len(hk) * len(ndd)), u.pvalue


# ── inputs ──────────────────────────────────────────────────────────────────
mei = pd.read_csv(ROOT / "data/hg38/gnomad_mei.tsv", sep="\t")
mei = mei.rename(columns={"#chrom": "chr", "chrom": "chr"})
mei = mei[mei.svtype.isin(CLASSES)].copy()
mei["cls"] = mei.svtype.map(CLASSES)
mei["pos"] = (mei.start + mei.end) // 2

prom = {}
for name in ("HighConfNDD", "Housekeeping"):
    df = pd.read_csv(ROOT / f"results/hg38/{name}_Alu.bed", sep="\t", header=None,
                     names=["chr", "start", "end", "gene", "score", "strand", "n"])
    df["kb"] = (df.end - df.start) / 1000
    prom[name] = df

pairs = pd.read_csv(ROOT / "results/matched/matched_pairs.tsv", sep="\t")
cpairs = pd.read_csv(ROOT / "results/matched/constraint_matched_pairs.tsv", sep="\t")
print(f"bağlam-eşleştirilmiş çift    : {len(pairs)}   (GC, gen yoğunluğu, rekombinasyon)")
print(f"kısıtlılık-eşleştirilmiş çift : {len(cpairs)}   (LOEUF, beyin ifadesi, GC)")

# ── count insertions per promoter, per element class ───────────────────────
idx = {c: {k: np.sort(g.pos.to_numpy()) for k, g in sub.groupby("chr")}
       for c, sub in mei.groupby("cls")}


def counts(df, cls):
    out = np.zeros(len(df), dtype=int)
    for i, (c, s, e) in enumerate(zip(df.chr, df.start, df.end)):
        arr = idx[cls].get(c)
        if arr is not None:
            out[i] = np.searchsorted(arr, e) - np.searchsorted(arr, s)
    return out


for name, df in prom.items():
    for cls in CLASSES.values():
        df[cls] = counts(df, cls)
    df["SINE"] = df["Alu"] + df["SVA"]

ndd_all, hk_all = prom["HighConfNDD"], prom["Housekeeping"]

# matched subsets, keyed by the pair table
ndd_m = ndd_all[ndd_all.gene.isin(pairs.ndd_gene)].copy()
hk_m = hk_all[hk_all.gene.isin(pairs.matched_hk_gene)].copy()
gc = dict(zip(pairs.ndd_gene, pairs.ndd_gc)) | dict(zip(pairs.matched_hk_gene, pairs.hk_gc))
for d in (ndd_m, hk_m):
    d["gc"] = d.gene.map(gc)

print(f"bağlam-eşleştirilmiş NDD {len(ndd_m)}, kontrol {len(hk_m)}")

# constraint-matched set: controls are drawn from all protein-coding genes,
# so their promoter windows come from the pair table rather than the two
# curated promoter BEDs, which between them cover only NDD and housekeeping.
def frame(prefix):
    d = cpairs[[f"{prefix}_gene", f"{prefix}_chrom", f"{prefix}_start",
                f"{prefix}_end", f"{prefix}_gc"]].copy()
    d.columns = ["gene", "chr", "start", "end", "gc"]
    d["kb"] = (d.end - d.start) / 1000
    for cls in CLASSES.values():
        d[cls] = counts(d, cls)
    d["SINE"] = d["Alu"] + d["SVA"]
    return d


ndd_c, ctl_c = frame("ndd"), frame("control")
print(f"kısıtlılık-eşleştirilmiş NDD {len(ndd_c)}, kontrol {len(ctl_c)}"
      f"   (kontrollerin {ctl_c.gene.isin(hk_all.gene).sum()} tanesi housekeeping)")

# ── test 1: density under three control definitions ────────────────────────
rows = []
for cls in ("Alu", "SVA", "LINE1", "SINE"):
    for label, a, b in (("tüm housekeeping", hk_all, ndd_all),
                        ("bağlam-eşleştirilmiş", hk_m, ndd_m),
                        ("kısıtlılık-eşleştirilmiş", ctl_c, ndd_c)):
        dh, dn = (a[cls] / a.kb).to_numpy(), (b[cls] / b.kb).to_numpy()
        r, p = rank_biserial(dh, dn)
        rows.append(dict(karsilastirma=label, TE=cls, n_HK=len(a), n_NDD=len(b),
                         ins_HK=int(a[cls].sum()), ins_NDD=int(b[cls].sum()),
                         dens_HK=round(dh.mean(), 4), dens_NDD=round(dn.mean(), 4),
                         oran=round(dn.mean() / dh.mean(), 3) if dh.mean() else np.nan,
                         r=round(r, 4), p=p))

res = pd.DataFrame(rows)
print("\n=== YOĞUNLUK: tüm housekeeping vs eşleştirilmiş kontrol ===")
print(f"{'karşılaştırma':<27}{'TE':<7}{'HK/kb':>9}{'NDD/kb':>9}{'NDD/HK':>9}{'r':>9}{'p':>11}")
print("-" * 81)
for _, x in res.iterrows():
    print(f"{x.karsilastirma:<27}{x.TE:<7}{x.dens_HK:>9.4f}{x.dens_NDD:>9.4f}"
          f"{x.oran:>9.3f}{x.r:>9.3f}{x.p:>11.2g}")

# ── test 2: GC-stratified, the detection-bias control ──────────────────────
print("\n=== GC KATMANLI (saptama yanlılığı kontrolü) ===")
both = pd.concat([ndd_m.assign(set="NDD"), hk_m.assign(set="HK")])
both = both[both.gc.notna()]
edges = both.gc.quantile([0, .25, .5, .75, 1]).to_numpy()
edges[-1] += 1e-9
gc_rows = []
print(f"{'GC çeyreği':<16}{'n NDD':>7}{'n HK':>7}{'NDD/kb':>9}{'HK/kb':>9}{'oran':>8}{'r':>8}{'p':>10}")
print("-" * 74)
for q in range(4):
    lo, hi = edges[q], edges[q + 1]
    sub = both[(both.gc >= lo) & (both.gc < hi)]
    a = sub[sub.set == "HK"]; b = sub[sub.set == "NDD"]
    if len(a) < 20 or len(b) < 20:
        continue
    dh, dn = (a.SINE / a.kb).to_numpy(), (b.SINE / b.kb).to_numpy()
    r, p = rank_biserial(dh, dn)
    gc_rows.append(dict(stratum=f"{lo:.3f}-{hi:.3f}", n_NDD=len(b), n_HK=len(a),
                        dens_NDD=round(dn.mean(), 4), dens_HK=round(dh.mean(), 4),
                        oran=round(dn.mean() / dh.mean(), 3) if dh.mean() else np.nan,
                        r=round(r, 4), p=p))
    print(f"{lo:.3f}–{hi:.3f}   {len(b):>7}{len(a):>7}{dn.mean():>9.4f}{dh.mean():>9.4f}"
          f"{dn.mean()/dh.mean() if dh.mean() else np.nan:>8.3f}{r:>8.3f}{p:>10.2g}")

# ── test 3: full allele-frequency spectrum ─────────────────────────────────
print("\n=== ALLEL FREKANS SPEKTRUMU (kısıtlılık-eşleştirilmiş kontrole karşı) ===")
def insertions_in(df):
    keep = []
    for c, s, e in zip(df.chr, df.start, df.end):
        sub = mei[(mei.chr == c) & (mei.pos >= s) & (mei.pos < e)]
        if len(sub):
            keep.append(sub)
    return pd.concat(keep) if keep else mei.iloc[:0]


ins_ndd, ins_hk = insertions_in(ndd_c), insertions_in(ctl_c)
BINS = [0, 1e-5, 1e-4, 1e-3, 1e-2, 1.0]
LBL = ["≤1e-5", "1e-5–1e-4", "1e-4–1e-3", "1e-3–1e-2", ">1e-2"]
afs = []
print(f"{'AF aralığı':<14}{'NDD n':>8}{'NDD %':>8}{'HK n':>8}{'HK %':>8}")
print("-" * 46)
for k in range(5):
    a = ((ins_ndd.AF > BINS[k]) & (ins_ndd.AF <= BINS[k + 1])).sum()
    b = ((ins_hk.AF > BINS[k]) & (ins_hk.AF <= BINS[k + 1])).sum()
    afs.append(dict(bin=LBL[k], n_NDD=int(a), n_HK=int(b),
                    pct_NDD=round(100 * a / max(len(ins_ndd), 1), 1),
                    pct_HK=round(100 * b / max(len(ins_hk), 1), 1)))
    print(f"{LBL[k]:<14}{a:>8}{100*a/max(len(ins_ndd),1):>8.1f}{b:>8}"
          f"{100*b/max(len(ins_hk),1):>8.1f}")

u = stats.mannwhitneyu(ins_ndd.AF, ins_hk.AF, alternative="less")
print(f"\ntüm AF dağılımı, tek yönlü MWU (NDD < HK) : p = {u.pvalue:.4g}")
print(f"  n insersiyon: NDD {len(ins_ndd)}, HK {len(ins_hk)}")
print(f"  medyan AF   : NDD {ins_ndd.AF.median():.3g}, HK {ins_hk.AF.median():.3g}")

# ── test 4: size axis ──────────────────────────────────────────────────────
print("\n=== BOYUT EKSENİ (seçilim modeli büyükte daha güçlü tükenme öngörür) ===")
print(f"{'TE':<7}{'yaklaşık bp':>12}{'NDD/HK oranı':>15}{'r':>9}{'p':>11}")
print("-" * 56)
for cls in ("Alu", "SVA", "LINE1"):
    x = res[(res.TE == cls) & (res.karsilastirma == "kısıtlılık-eşleştirilmiş")].iloc[0]
    print(f"{cls:<7}{APPROX_BP[cls]:>12}{x.oran:>15.3f}{x.r:>9.3f}{x.p:>11.2g}")

res.to_csv(OUT / "matched_mei.csv", index=False)
pd.DataFrame(afs).to_csv(OUT / "afs_bins.csv", index=False)
pd.DataFrame(gc_rows).to_csv(OUT / "gc_stratified_mei.csv", index=False)
print(f"\nyazıldı: {OUT.relative_to(ROOT)}/{{matched_mei,afs_bins,gc_stratified_mei}}.csv")
