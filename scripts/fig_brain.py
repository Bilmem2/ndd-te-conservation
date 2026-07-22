#!/usr/bin/env python3
"""Figure 9 — neurodevelopmental preferential association.
Panel A: per-sample NDD/HK fetal-DNase enrichment ratio (3 brain donors vs 3
non-neural control tissues). Panel B: Alu depletion effect size across
composite fetal-brain-DNase quartiles. Reads results/brain/."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
R = REPO / "results" / "brain"
FIGS = REPO / "figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
BRAINC, CTRLC, NDDC = "#5b53b0", "#9aa1ab", "#E07B39"

ps = pd.read_csv(R / "brain_overlay_persample.csv")
label = {"brain_122d": "Brain 122d", "brain_117d": "Brain 117d", "brain_101d": "Brain 101d",
         "liver_113d": "Liver 113d", "lung_120d": "Lung 120d", "stomach_110d": "Stomach 110d"}
order = ["brain_122d", "brain_117d", "brain_101d", "liver_113d", "lung_120d", "stomach_110d"]
ps = ps.set_index("sample").loc[order].reset_index()

# Panel B values (from 23_brain_overlay.py, composite brain-DNase quartiles)
QLAB = ["Q1(low)", "Q2", "Q3", "Q4(high)"]
QR = [-0.294, -0.265, -0.275, -0.391]

fig, ax = plt.subplots(1, 2, figsize=(13, 5))

# Panel A
cols = [BRAINC if g == "brain" else CTRLC for g in ps["group"]]
bars = ax[0].bar(range(len(ps)), ps["ratio_NDD_HK"], color=cols, edgecolor="black", lw=0.6, alpha=0.9)
ax[0].axhline(1.0, color="k", lw=1, ls="--")
for i, (r, p) in enumerate(zip(ps["ratio_NDD_HK"], ps["p_value"])):
    sig = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
    ax[0].text(i, r + 0.015, f"{r:.2f}\n{sig}", ha="center", va="bottom", fontsize=8)
ax[0].set_xticks(range(len(ps)))
ax[0].set_xticklabels([label[s] for s in ps["sample"]], rotation=30, ha="right", fontsize=9)
ax[0].set_ylabel("NDD / housekeeping DNase enrichment")
ax[0].set_ylim(0, 1.45)
ax[0].set_title("A  Fetal DNase enrichment at NDD promoters\n(3 brain donors vs 3 non-neural tissues)",
                fontweight="bold", fontsize=10, loc="left")
from matplotlib.patches import Patch
ax[0].legend(handles=[Patch(color=BRAINC, label="fetal brain"), Patch(color=CTRLC, label="non-neural control")],
             fontsize=9, loc="upper right")

# Panel B
b = ax[1].bar(QLAB, QR, color=NDDC, edgecolor="black", lw=0.6, alpha=0.9)
for i, r in enumerate(QR):
    ax[1].text(i, r - 0.012, f"{r:.3f}", ha="center", va="top", fontsize=8)
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set_ylim(min(QR) * 1.25, 0.03)
ax[1].set_ylabel("Alu depletion  (rank-biserial r, HK > NDD)")
ax[1].set_xlabel("composite fetal-brain DNase density")
ax[1].set_title("B  Alu depletion is strongest in the most\nfetal-brain-active promoters",
                fontweight="bold", fontsize=10, loc="left")
for i in range(4):
    ax[1].text(i, 0.008, "***", ha="center", fontsize=9)

plt.suptitle("Alu depletion at NDD promoters preferentially associates with fetal-brain regulatory elements",
             fontweight="bold", fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(f"{FIGS}/Fig9_BrainSpecificity.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIGS}/Fig9_BrainSpecificity.png", dpi=150, bbox_inches="tight")
plt.close()
print("Fig9 saved.")
