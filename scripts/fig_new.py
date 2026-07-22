#!/usr/bin/env python3
"""Generate the two NEW manuscript figures for the revised paper:
  Fig 7 — Alu depletion is robust to genomic-context confounders
          (gene density, recombination, joint matched control).
  Fig 8 — Depleted NDD promoters are regulatory-active (ENCODE cCRE overlay).
Reads the already-computed result CSVs (scripts 15/16/18)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
R = REPO / "results"
FIGS = REPO / "figures"
NDDC, HKC = "#E07B39", "#4878CF"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})


def barlabel(ax, bars, vals, fmt="{:.3f}"):
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() - 0.012 if v < 0 else b.get_height() + 0.005,
                fmt.format(v), ha="center", va="top" if v < 0 else "bottom", fontsize=8)


# ================= FIG 7: confounder robustness =================
ds = pd.read_csv(R / "context" / "alu_by_density_stratum.csv")
rs = pd.read_csv(R / "context" / "alu_by_recomb_stratum.csv")
mc = pd.read_csv(R / "matched" / "matched_alu_test.csv")

fig, ax = plt.subplots(1, 3, figsize=(16, 5))

for a, tab, title in [(ax[0], ds, "Gene-density quartile"),
                      (ax[1], rs, "Recombination quartile")]:
    labs = tab["stratum"].astype(str).tolist()
    rvals = tab["r"].tolist()
    bars = a.bar(labs, rvals, color=NDDC, alpha=0.85, edgecolor="black", lw=0.6)
    barlabel(a, bars, rvals)
    a.axhline(0, color="k", lw=0.8)
    a.set_title(f"Alu depletion within\n{title}", fontweight="bold", fontsize=11)
    a.set_ylabel("rank-biserial r (HK vs NDD)")
    a.set_ylim(min(rvals) * 1.35, 0.05)
    for i, p in enumerate(tab["p_value"]):
        a.text(i, 0.01, "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns",
               ha="center", fontsize=9)

# matched control
lab = ["NDD vs\nALL housekeeping", "NDD vs\nMATCHED control\n(GC+density+recomb)"]
rvals = mc["r"].tolist()
bars = ax[2].bar(lab, rvals, color=[HKC, NDDC], alpha=0.85, edgecolor="black", lw=0.6)
barlabel(ax[2], bars, rvals)
ax[2].axhline(0, color="k", lw=0.8)
ax[2].set_title("Joint multivariate matching", fontweight="bold", fontsize=11)
ax[2].set_ylabel("rank-biserial r")
ax[2].set_ylim(min(rvals) * 1.35, 0.05)

plt.suptitle("Alu Depletion at NDD Promoters Is Robust to Genomic-Context Confounders\n"
             "NDD and housekeeping promoters differ in gene density and recombination, "
             "yet depletion persists in every stratum and after joint matching",
             fontweight="bold", fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig(f"{FIGS}/Fig7_ContextControls.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIGS}/Fig7_ContextControls.png", dpi=150, bbox_inches="tight")
plt.close()
print("Fig7 saved.")

# ================= FIG 8: cCRE functional overlay =================
gd = pd.read_csv(R / "ccre" / "ccre_group_diff.csv")
fp = pd.read_csv(R / "ccre" / "alu_free_vs_pos.csv")
cs = pd.read_csv(R / "ccre" / "alu_by_ccre_stratum.csv")

fig, ax = plt.subplots(1, 3, figsize=(16, 5))

# A: cCRE density NDD vs HK by class
groups = ["PLS", "ELS", "active_all"]
gsub = gd[gd["cCRE_group"].isin(groups)].set_index("cCRE_group").loc[groups]
x = np.arange(len(groups)); w = 0.38
ax[0].bar(x - w / 2, gsub["mean_HK"], w, label="Housekeeping", color=HKC, alpha=0.85, edgecolor="black", lw=0.6)
ax[0].bar(x + w / 2, gsub["mean_NDD"], w, label="NDD", color=NDDC, alpha=0.85, edgecolor="black", lw=0.6)
ax[0].set_xticks(x); ax[0].set_xticklabels(["Promoter-like\n(PLS)", "Enhancer-like\n(ELS)", "Active\n(all)"], fontsize=9)
ax[0].set_ylabel("cCRE density (per kb)")
ax[0].set_title("NDD promoters are regulatory-active\n(comparable cCRE density to HK)", fontweight="bold", fontsize=10)
ax[0].legend(fontsize=9)

# B: Alu-free vs Alu-positive active cCRE density
fa = fp.set_index("cCRE_group")
gg = ["PLS", "ELS", "active_all"]
mf = [fa.loc[g, "median_free"] for g in gg]
mp = [fa.loc[g, "median_pos"] for g in gg]
x = np.arange(len(gg))
ax[1].bar(x - w / 2, mf, w, label="Alu-free promoters", color="#5aa469", alpha=0.9, edgecolor="black", lw=0.6)
ax[1].bar(x + w / 2, mp, w, label="Alu-positive promoters", color="#b0763a", alpha=0.9, edgecolor="black", lw=0.6)
ax[1].set_xticks(x); ax[1].set_xticklabels(["PLS", "ELS", "Active"], fontsize=10)
ax[1].set_ylabel("cCRE count (median)")
ax[1].set_title("Alu-free promoters are MORE\nregulatory-dense", fontweight="bold", fontsize=10)
ax[1].legend(fontsize=9)

# C: Alu depletion within active-cCRE quartiles
labs = cs["active_cCRE_stratum"].astype(str).tolist()
rvals = cs["r"].tolist()
bars = ax[2].bar(labs, rvals, color=NDDC, alpha=0.85, edgecolor="black", lw=0.6)
barlabel(ax[2], bars, rvals)
ax[2].axhline(0, color="k", lw=0.8)
ax[2].set_ylim(min(rvals) * 1.3, 0.05)
ax[2].set_ylabel("rank-biserial r (HK vs NDD)")
ax[2].set_title("Depletion persists across\nregulatory-density quartiles", fontweight="bold", fontsize=10)
for i, p in enumerate(cs["p_value"]):
    ax[2].text(i, 0.01, "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns", ha="center", fontsize=9)

plt.suptitle("Alu-Depleted NDD Promoters Are Functionally Active Regulatory Regions (ENCODE cCREs)",
             fontweight="bold", fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(f"{FIGS}/Fig8_cCRE.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIGS}/Fig8_cCRE.png", dpi=150, bbox_inches="tight")
plt.close()
print("Fig8 saved.")
