#!/usr/bin/env python3
"""Figures 3 and 4 of the manuscript.

  Figure 3 — the deficit against genomic-context confounders
             (gene density, recombination rate, joint covariate matching)
  Figure 4 — the depleted promoters against ENCODE cCRE annotation

Both read the effect sizes already computed by scripts 15, 16 and 18; nothing is
recalculated here, so the numbers match the tables.

Two things this file used to get wrong. Every panel scaled its own y-axis, so
three panels showing the same quantity could not be compared by bar length --
the recombination panel reached -0.5 while the one beside it stopped at -0.4.
And significance was marked on some panels but not others. Both are fixed: the
panels of a figure share one axis, and every comparison carries its marker and
its sample sizes.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
R, FIGS = REPO / "results", REPO / "figures"
NDDC, HKC = "#E07B39", "#4878CF"
FREE_C, POS_C = "#5AA469", "#B0763A"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def panel(ax, letter, title):
    ax.set_title(f"{letter}   {title}", loc="left", fontweight="bold", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)


def rbars(ax, labels, rvals, pvals, ns, lo):
    """One panel of rank-biserial effect sizes, marked and labelled."""
    bars = ax.bar(labels, rvals, color=NDDC, edgecolor="black", lw=0.6, zorder=2)
    for b, v, p, (nh, nd) in zip(bars, rvals, pvals, ns):
        x = b.get_x() + b.get_width() / 2
        ax.text(x, v - abs(lo) * 0.035, f"{v:.3f}", ha="center", va="top", fontsize=8.5)
        ax.text(x, abs(lo) * 0.045, stars(p), ha="center", fontsize=9)
        ax.text(x, abs(lo) * 0.115, f"{nh} / {nd}", ha="center", fontsize=7.5,
                color="#777")
    ax.axhline(0, color="k", lw=0.8)
    ax.grid(axis="y", color="#EEE", lw=0.8)
    ax.set_axisbelow(True)


# ── Figure 3: genomic-context confounders ────────────────────────────────────
ds = pd.read_csv(R / "context" / "alu_by_density_stratum.csv")
rs = pd.read_csv(R / "context" / "alu_by_recomb_stratum.csv")
mc = pd.read_csv(R / "matched" / "matched_alu_test.csv")

allr = ds["r"].tolist() + rs["r"].tolist() + mc["r"].tolist()
LO, HI = min(allr) * 1.30, abs(min(allr)) * 0.20     # one scale for all three

fig, ax = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
for a, tab, letter, title in [(ax[0], ds, "a", "Within gene-density quartile"),
                              (ax[1], rs, "b", "Within recombination quartile")]:
    rbars(a, tab["stratum"].astype(str).tolist(), tab["r"].tolist(),
          tab["p_value"].tolist(), list(zip(tab["n_HK"], tab["n_NDD"])), LO)
    panel(a, letter, title)

rbars(ax[2], ["vs all\nhousekeeping", "vs matched control\n(GC + density + recomb)"],
      mc["r"].tolist(), mc["p_value"].tolist(),
      list(zip(mc["n_HK"], mc["n_NDD"])), LO)
panel(ax[2], "c", "Joint covariate matching")

ax[0].set_ylabel("Rank-biserial $r$ (NDD vs control)")
ax[0].set_ylim(LO, HI)
fig.text(0.5, -0.02, "Bar labels give $r$; above each bar, $n$ as housekeeping / NDD. "
         "*** $p<0.001$", ha="center", fontsize=9, color="#555")
plt.tight_layout()
plt.savefig(FIGS / "Fig3_ContextControls.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGS / "Fig3_ContextControls.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Fig3 saved  (ortak eksen {LO:.3f}..{HI:.3f})")

# ── Figure 4: cCRE overlay ───────────────────────────────────────────────────
gd = pd.read_csv(R / "ccre" / "ccre_group_diff.csv")
fp = pd.read_csv(R / "ccre" / "alu_free_vs_pos.csv")
cs = pd.read_csv(R / "ccre" / "alu_by_ccre_stratum.csv")

GROUPS = ["PLS", "ELS", "active_all"]
NICE = ["Promoter-like\n(PLS)", "Enhancer-like\n(ELS)", "Active\n(all)"]
fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
w, x = 0.38, np.arange(len(GROUPS))

g = gd[gd["cCRE_group"].isin(GROUPS)].set_index("cCRE_group").loc[GROUPS]
ax[0].bar(x - w / 2, g["mean_HK"], w, label="Housekeeping", color=HKC,
          edgecolor="black", lw=0.6, zorder=2)
ax[0].bar(x + w / 2, g["mean_NDD"], w, label="NDD", color=NDDC,
          edgecolor="black", lw=0.6, zorder=2)
top = max(g["mean_HK"].max(), g["mean_NDD"].max())
for i, p in enumerate(g["p_value"]):
    ax[0].text(i, max(g["mean_HK"].iloc[i], g["mean_NDD"].iloc[i]) + top * 0.03,
               stars(p), ha="center", fontsize=9)
ax[0].set_xticks(x, NICE, fontsize=9)
ax[0].set_ylabel("cCRE density (per kb)")
ax[0].set_ylim(0, top * 1.18)
ax[0].legend(fontsize=9, frameon=False)
panel(ax[0], "a", "cCRE density by class")

f = fp.set_index("cCRE_group").loc[GROUPS]
ax[1].bar(x - w / 2, f["median_free"], w, label="Alu-free promoters", color=FREE_C,
          edgecolor="black", lw=0.6, zorder=2)
ax[1].bar(x + w / 2, f["median_pos"], w, label="Alu-positive promoters", color=POS_C,
          edgecolor="black", lw=0.6, zorder=2)
top = max(f["median_free"].max(), f["median_pos"].max())
for i, p in enumerate(f["p_free_gt_pos"]):
    ax[1].text(i, max(f["median_free"].iloc[i], f["median_pos"].iloc[i]) + top * 0.03,
               stars(p), ha="center", fontsize=9)
ax[1].set_xticks(x, ["PLS", "ELS", "Active"], fontsize=10)
ax[1].set_ylabel("cCRE count (median)")
ax[1].set_ylim(0, top * 1.18)
ax[1].legend(fontsize=9, frameon=False)
panel(ax[1], "b", "Alu-free versus Alu-positive promoters")

LO4 = min(cs["r"]) * 1.30
rbars(ax[2], cs["active_cCRE_stratum"].astype(str).tolist(), cs["r"].tolist(),
      cs["p_value"].tolist(), list(zip(cs["n_HK"], cs["n_NDD"])), LO4)
ax[2].set_ylabel("Rank-biserial $r$ (NDD vs housekeeping)")
ax[2].set_ylim(LO4, abs(LO4) * 0.20)
panel(ax[2], "c", "Within active-cCRE density quartile")

fig.text(0.5, -0.02, "*** $p<0.001$; ** $p<0.01$; * $p<0.05$; ns, not significant. "
         "In c, $n$ is given as housekeeping / NDD.", ha="center", fontsize=9,
         color="#555")
plt.tight_layout()
plt.savefig(FIGS / "Fig4_cCRE.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGS / "Fig4_cCRE.png", dpi=150, bbox_inches="tight")
plt.close()
print("Fig4 saved")
