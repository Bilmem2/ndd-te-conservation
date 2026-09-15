#!/usr/bin/env python3
"""Figure 1 — Alu depletion at NDD promoters across the seven primates.

Replaces a seven-panel boxplot. Alu count in a 4 kb window is an integer, so
promoter Alu frequency takes only the values 0, 0.25, 0.5 ... and quartiles of
it are coarse: the fourteen boxes of that figure resolved to six distinct
(Q1, median, Q3) triples, four species drew identical boxes, and the two
marmoset medians coincided at 0.5 despite a highly significant rank difference.

Panel A shows the human distribution at full resolution, where the deficit reads
as what it is, a shift of mass toward zero. Panel B carries all seven species as
effect sizes on one axis, ordered by divergence time, which is the comparison
the section is about. Every value is recomputed from the committed promoter BEDs
and reproduces Table 1.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parent.parent
R, FIG = ROOT / "results", ROOT / "figures"

HK_C, NDD_C = "#5B8FD4", "#E08B4F"
WIN_KB = 4.0
RNG = np.random.default_rng(20260915)

SPECIES = [("Human", "hg38", 0), ("Orangutan", "ponAbe3", 16),
           ("Gibbon", "nomLeu3", 20), ("Macaque", "rheMac10", 25),
           ("Marmoset", "calJac4", 40), ("Squirrel monkey", "saiBol1", 40),
           ("Mouse lemur", "mmur3", 70)]


def counts(asm, grp):
    df = pd.read_csv(R / asm / f"{grp}_Alu.bed", sep="\t", header=None)
    return df.iloc[:, -1].to_numpy(float)


def rank_biserial(a, b):
    """r for a (NDD) lying below b (HK); the Mann-Whitney U effect size."""
    n1, n2 = len(a), len(b)
    ranks = pd.Series(np.concatenate([a, b])).rank().to_numpy()
    u1 = ranks[:n1].sum() - n1 * (n1 + 1) / 2
    return 2 * u1 / (n1 * n2) - 1


def boot_ci(a, b, n=2000):
    vals = np.empty(n)
    for i in range(n):
        vals[i] = rank_biserial(RNG.choice(a, len(a), replace=True),
                                RNG.choice(b, len(b), replace=True))
    return np.percentile(vals, [2.5, 97.5])


fig = plt.figure(figsize=(13.2, 5.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.22], wspace=0.28)

# ── A. what the depletion looks like, in the anchor genome ─────────────────
axA = fig.add_subplot(gs[0, 0])
hk, ndd = counts("hg38", "Housekeeping"), counts("hg38", "HighConfNDD")
TOP = 5
bins = np.arange(TOP + 2)


def frac(v):
    c = np.bincount(np.minimum(v, TOP).astype(int), minlength=TOP + 1)
    return c / c.sum()


fh, fn = frac(hk), frac(ndd)
w = 0.38
axA.bar(bins[:-1] - w / 2, fh, w, color=HK_C, edgecolor="black", linewidth=0.6,
        label="Housekeeping")
axA.bar(bins[:-1] + w / 2, fn, w, color=NDD_C, edgecolor="black", linewidth=0.6,
        label="NDD")
for x, (p, q) in enumerate(zip(fh, fn)):
    axA.text(x - w / 2, p + .012, f"{100*p:.0f}", ha="center", fontsize=8.5,
             color="#2A4D7A")
    axA.text(x + w / 2, q + .012, f"{100*q:.0f}", ha="center", fontsize=8.5,
             color="#8A4A12")

axA.set_xticks(bins[:-1])
axA.set_xticklabels([str(i) for i in range(TOP)] + [f"{TOP}+"])
axA.set_xlabel("Alu elements per promoter window (TSS $\\pm$2 kb)")
axA.set_ylabel("Fraction of promoters")
axA.set_ylim(0, max(fn.max(), fh.max()) * 1.20)
axA.legend(frameon=False, loc="upper right", fontsize=10, bbox_to_anchor=(1.0, 0.94))
axA.set_title("A   Human promoters", loc="left", fontweight="bold", fontsize=12)
axA.spines[["top", "right"]].set_visible(False)
axA.text(0.30, 0.97, f"mean {hk.mean()/WIN_KB:.2f} vs {ndd.mean()/WIN_KB:.2f} per kb",
         transform=axA.transAxes, fontsize=9.5, color="#444", ha="center")

# ── B. the cross-species claim on one axis ────────────────────────────────
axB = fig.add_subplot(gs[0, 1])
rows = []
for name, asm, mya in SPECIES:
    h, d = counts(asm, "Housekeeping"), counts(asm, "HighConfNDD")
    r = rank_biserial(d, h)
    lo, hi = boot_ci(d, h)
    rows.append((name, mya, r, lo, hi, h.mean() / WIN_KB, d.mean() / WIN_KB))
    print(f"  {name:<16} r={r:+.3f}  CI[{lo:+.3f},{hi:+.3f}]  "
          f"HK={h.mean()/WIN_KB:.3f} NDD={d.mean()/WIN_KB:.3f}")

ys = np.arange(len(rows))[::-1]
for y, (name, mya, r, lo, hi, mh, md) in zip(ys, rows):
    axB.plot([lo, hi], [y, y], color="#333", linewidth=1.4, zorder=2)
    for edge in (lo, hi):
        axB.plot([edge, edge], [y - .12, y + .12], color="#333", linewidth=1.4)
    axB.scatter([r], [y], s=86, color=NDD_C, edgecolor="black", linewidth=0.8,
                zorder=3)
    axB.text(0.965, y, f"{mh:.2f} / {md:.2f}", transform=axB.get_yaxis_transform(),
             va="center", ha="right", fontsize=9, color="#555")

axB.axvline(0, color="#999", linestyle="--", linewidth=1)
axB.set_yticks(ys)
axB.set_yticklabels([f"{n}\n$\\sim${m} Mya" for n, m, *_ in rows], fontsize=9.5)
axB.set_xlabel("Rank-biserial $r$ (NDD vs housekeeping)")
axB.set_xlim(-0.46, 0.22)
axB.set_xticks([-0.4, -0.3, -0.2, -0.1, 0.0])
axB.set_ylim(-0.7, len(rows) - 0.25)
axB.set_title("B   Effect size across the primate panel", loc="left",
              fontweight="bold", fontsize=12)
axB.spines[["top", "right", "left"]].set_visible(False)
axB.tick_params(axis="y", length=0)
axB.text(0.965, len(rows) - 0.48, "mean per kb\nHK / NDD",
         transform=axB.get_yaxis_transform(), va="center", ha="right",
         fontsize=8.5, color="#555", style="italic")
axB.grid(axis="x", color="#EEE", linewidth=0.8)
axB.set_axisbelow(True)

fig.savefig(FIG / "Fig1_Alu_Primates.pdf", bbox_inches="tight")
fig.savefig(FIG / "Fig1_Alu_Primates.png", dpi=200, bbox_inches="tight")
print("\nyazildi: Fig1_Alu_Primates.pdf / .png")
