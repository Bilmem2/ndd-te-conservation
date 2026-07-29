#!/usr/bin/env python3
"""Regenerate Fig1 (Alu across primates, now incl. mouse lemur) and Fig3
(significance heatmap, now incl. mouse lemur) for the revised manuscript.
Focused update; does not touch bedtools-dependent figures."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RES = REPO / "results"
FIGS = REPO / "figures"
COLORS = {"Housekeeping": "#4878CF", "HighConfNDD": "#E07B39"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})


def dens(path):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", header=None)
    return (df.iloc[:, -1] / ((df.iloc[:, 2] - df.iloc[:, 1]) / 1000)).values  # last col = count


def sigstars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def add_sig(ax, x1, x2, y, p):
    ax.plot([x1, x1, x2, x2], [y * 0.97, y, y, y * 0.97], "k-", lw=0.8)
    ax.text((x1 + x2) / 2, y * 1.01, sigstars(p), ha="center", va="bottom", fontsize=10)


# ------- FIG 1: Alu across 7 primates (two Platyrrhini + a strepsirrhine) -------
primates = [
    ("hg38", "Human\n(hg38)", 0),
    ("ponAbe3", "Orangutan\n(ponAbe3)", 16),
    ("nomLeu3", "Gibbon\n(nomLeu3)", 20),
    ("rheMac10", "Macaque\n(rheMac10)", 25),
    ("calJac4", "Marmoset\n(calJac4)", 40),
    ("saiBol1", "Squirrel monkey\n(SaiBol1.0)", 40),
    ("mmur3", "Mouse lemur\n(Mmur_3.0)", 70),
]
fig, axes = plt.subplots(1, len(primates), figsize=(3.5 * len(primates), 6))
for ax, (sp, label, mya) in zip(axes, primates):
    hk = dens(f"{RES}/{sp}/Housekeeping_Alu.bed")
    nd = dens(f"{RES}/{sp}/HighConfNDD_Alu.bed")
    bp = ax.boxplot([hk, nd], patch_artist=True,
                    medianprops=dict(color="black", linewidth=2),
                    flierprops=dict(marker="o", markersize=2, alpha=0.3), widths=0.6)
    for patch, cat in zip(bp["boxes"], ["Housekeeping", "HighConfNDD"]):
        patch.set_facecolor(COLORS[cat]); patch.set_alpha(0.85)
    _, p = stats.mannwhitneyu(hk, nd, alternative="greater")
    ymax = max(np.quantile(hk, 0.95), np.quantile(nd, 0.95))
    add_sig(ax, 1, 2, (ymax if ymax > 0 else 1) * 1.18, p)
    ax.set_title(f"{label}\n~{mya} Mya", fontweight="bold", fontsize=10)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["HK", "NDD"], fontsize=10)
    if ax is axes[0]:
        ax.set_ylabel("Alu Frequency (count per kb)", fontsize=11)
    ax.text(0.5, -0.18, f"n={len(hk)}, {len(nd)}", transform=ax.transAxes,
            ha="center", fontsize=8, color="gray")

patches = [mpatches.Patch(color=COLORS[c], label=c.replace("HighConfNDD", "NDD"), alpha=0.85)
           for c in ["Housekeeping", "HighConfNDD"]]
fig.legend(handles=patches, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02),
           frameon=False, fontsize=10)
plt.suptitle("Alu Depletion at NDD Gene Promoters Across Primates (incl. strepsirrhine)\n"
             "TSS ± 2 kb | HighConfNDD vs Housekeeping | Mann-Whitney U",
             fontweight="bold", fontsize=12)
plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.savefig(f"{FIGS}/Fig1_Alu_Primates.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIGS}/Fig1_Alu_Primates.png", dpi=150, bbox_inches="tight")
plt.close()
print("Fig1 (7 primates incl. two Platyrrhini and a strepsirrhine) saved.")

# ---------------- FIG 3: heatmap with lemur row ----------------
sp_info = [
    ("hg38", "Human"), ("ponAbe3", "Orangutan"), ("nomLeu3", "Gibbon"),
    ("rheMac10", "Macaque"), ("calJac4", "Marmoset"),
    ("saiBol1", "Squirrel monkey"), ("mmur3", "Mouse lemur"),
    ("mm10", "Mouse"), ("canFam4", "Dog"),
]
n = len(sp_info)
pmat = np.full((n, 2), np.nan)
sig_mat = []
for i, (sp, label) in enumerate(sp_info):
    row = []
    for j, te in enumerate(["Alu", "LINE1"]):
        te_file = "B1B2" if (sp == "mm10" and te == "Alu") else te
        fn = f"{RES}/{sp}/HighConfNDD_{te_file}.bed"
        hn = f"{RES}/{sp}/Housekeeping_{te_file}.bed"
        if not (os.path.exists(fn) and os.path.exists(hn)):
            row.append("—"); continue
        hk, nd = dens(hn), dens(fn)
        _, p = stats.mannwhitneyu(hk, nd, alternative="greater")
        pmat[i, j] = min(-np.log10(p), 60) if p > 0 else 60
        row.append(f"{sigstars(p)}\np={p:.1e}")
    sig_mat.append(row)

fig, ax = plt.subplots(figsize=(7, 7.5))
im = ax.imshow(np.ma.masked_invalid(pmat), cmap="RdYlGn", aspect="auto", vmin=0, vmax=20)
plt.colorbar(im, ax=ax, label="-log₁₀(p-value)  HK > NDD")
ax.set_xticks([0, 1]); ax.set_xticklabels(["Alu / B1-B2\n(SINE)", "LINE-1"], fontsize=11)
ax.set_yticks(range(n)); ax.set_yticklabels([s[1] for s in sp_info], fontsize=12, fontweight="bold")
for i in range(n):
    for j in range(2):
        v = pmat[i, j]
        if np.isnan(v):
            ax.text(j, i, "—", ha="center", va="center", fontsize=11, color="white")
        else:
            ax.text(j, i, sig_mat[i][j], ha="center", va="center", fontsize=8,
                    color="white" if v > 10 else "black", fontweight="bold")
ax.set_title("TE Depletion Significance: HK vs NDD Gene Promoters\n"
             "Across Species and TE Classes", fontweight="bold", fontsize=11)
plt.tight_layout()
plt.savefig(f"{FIGS}/Fig3_Heatmap.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIGS}/Fig3_Heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("Fig3 (9 species) saved.")
