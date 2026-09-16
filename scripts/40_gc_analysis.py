"""
40_gc_analysis.py — promoter GC content and GC-stratified Alu depletion.

Compares promoter GC composition between the NDD and housekeeping sets, then
repeats the Alu comparison within GC strata to test whether the depletion is
explained by GC composition. Produces the supplementary GC figure.

Inputs : data/hg38/hg38.2bit
         results/hg38/{HighConfNDD_Alu.bed, Housekeeping_Alu.bed}
Outputs: results/hg38/gc_analysis.csv
         figures/ESM4_GC_Analysis.{pdf,png}
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import twobitreader
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT / "results" / "hg38", ROOT / "figures"
COLORS = {"Housekeeping": "#4C72B0", "HighConfNDD": "#DD8452"}
STRATA = [(0.0, 0.40, "<40%"), (0.40, 0.50, "40-50%"),
          (0.50, 0.60, "50-60%"), (0.60, 0.70, "60-70%")]


def load(cat):
    df = pd.read_csv(RES / f"{cat}_Alu.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "gene", "score", "strand", "alu"],
                     lineterminator="\n")
    for c in ("start", "end", "alu"):
        df[c] = pd.to_numeric(df[c].astype(str).str.strip(), errors="coerce")
    df = df.dropna(subset=["start", "end", "alu"]).reset_index(drop=True)
    df["density"] = df["alu"] / ((df["end"] - df["start"]) / 1000.0)
    df["category"] = cat
    return df


def gc_of(genome, df):
    out = np.full(len(df), np.nan)
    for i, r in df.iterrows():
        try:
            seq = genome[r["chrom"]][int(r["start"]):int(r["end"])].upper()
        except Exception:
            continue
        acgt = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")
        if acgt:
            out[i] = (seq.count("G") + seq.count("C")) / acgt
    return out


def rb(hk, nd):
    """rank-biserial, one-sided HK > NDD"""
    u = stats.mannwhitneyu(hk, nd, alternative="greater")
    return -(2 * u.statistic / (len(hk) * len(nd)) - 1), u.pvalue


print("reading hg38.2bit ...")
genome = twobitreader.TwoBitFile(str(ROOT / "data" / "hg38" / "hg38.2bit"))

frames = []
for cat in ("Housekeeping", "HighConfNDD"):
    d = load(cat)
    print(f"  {cat}: {len(d)} promoters -> GC")
    d["gc"] = gc_of(genome, d)
    frames.append(d.dropna(subset=["gc"]))
hk, nd = frames

u2 = stats.mannwhitneyu(nd.gc, hk.gc, alternative="two-sided")
r_gc = 2 * u2.statistic / (len(nd) * len(hk)) - 1

rows = [dict(stratum="all", n_HK=len(hk), n_NDD=len(nd),
             median_gc_HK=round(hk.gc.median(), 3), median_gc_NDD=round(nd.gc.median(), 3),
             gc_p=u2.pvalue, gc_r=round(r_gc, 3), alu_r=np.nan, alu_p=np.nan)]
for lo, hi, lab in STRATA:
    h = hk[(hk.gc >= lo) & (hk.gc < hi)].density
    n = nd[(nd.gc >= lo) & (nd.gc < hi)].density
    if len(h) < 10 or len(n) < 10:
        continue
    r, p = rb(h, n)
    rows.append(dict(stratum=lab, n_HK=len(h), n_NDD=len(n),
                     median_gc_HK=np.nan, median_gc_NDD=np.nan, gc_p=np.nan, gc_r=np.nan,
                     alu_r=round(r, 3), alu_p=p))

out = pd.DataFrame(rows)
out.to_csv(RES / "gc_analysis.csv", index=False)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for d, lab in ((hk, "Housekeeping"), (nd, "HighConfNDD")):
    ax[0].hist(d.gc, bins=40, density=True, alpha=0.6, color=COLORS[lab],
               label=f"{lab.replace('HighConfNDD','NDD')} (n={len(d)})")
ax[0].set_xlabel("Promoter GC content")
ax[0].set_ylabel("Density")
ax[0].legend(frameon=False)
ax[0].set_title("a   Promoter GC content", loc="left", fontweight="bold", fontsize=11)

# The right panel used to be a scatter of Alu frequency against GC content. Alu
# frequency in a 4 kb window takes about a dozen distinct values, so the points
# fell into horizontal stripes and the two groups overplotted into one mass: the
# claim the panel exists to support, that the deficit holds at every GC level,
# could not be read off it. The stratified comparison was already being computed
# for the CSV, so it is drawn here instead.
strata = [row for row in rows if row["stratum"] != "all"]
xs = np.arange(len(strata))
W = 0.36
for off, (frame, lab) in ((-W / 2, (hk, "Housekeeping")), (+W / 2, (nd, "HighConfNDD"))):
    means, errs = [], []
    for lo, hi, lab_s in STRATA:
        if not any(r["stratum"] == lab_s for r in strata):
            continue
        v = frame[(frame.gc >= lo) & (frame.gc < hi)].density
        means.append(v.mean())
        errs.append(1.96 * v.std(ddof=1) / np.sqrt(len(v)))
    ax[1].bar(xs + off, means, W, yerr=errs, capsize=3, color=COLORS[lab],
              edgecolor="black", linewidth=0.6,
              label=lab.replace("HighConfNDD", "NDD"), zorder=2)

top = max(b.get_height() for b in ax[1].patches)
for x, row in zip(xs, strata):
    ax[1].text(x, top * 1.23, f"$r$ = {row['alu_r']:+.3f}", ha="center", fontsize=9)
    ax[1].text(x, top * 1.14, f"n = {row['n_HK']} / {row['n_NDD']}", ha="center",
               fontsize=8, color="#666")
ax[1].set_xticks(xs)
ax[1].set_xticklabels([r["stratum"] for r in strata])
ax[1].set_xlabel("Promoter GC content stratum")
ax[1].set_ylabel("Mean Alu frequency (count per kb)")
ax[1].set_title("b   Alu frequency within GC stratum", loc="left",
                fontweight="bold", fontsize=11)
ax[1].set_ylim(0, top * 1.34)
ax[1].spines[["top", "right"]].set_visible(False)
# The left panel already carries the colour key for both groups.
ax[0].spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(FIGS / "ESM4_GC_Analysis.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGS / "ESM4_GC_Analysis.png", dpi=150, bbox_inches="tight")
plt.close()

print("\n=== recomputed ===")
print(out.to_string(index=False))
