#!/usr/bin/env python3
"""Regenerate Fig 4 (permutation null model) with all five species-TE
combinations reported in the Results, adding the marmoset panel that the
earlier four-panel version omitted. Mirrors 24_null_model.py (same seed and
combinations) so the figure and results/null_model_full.csv agree."""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGS = REPO / "figures"
N_PERM = 10000
SEED = 42

COMBOS = [
    ("hg38", "Human", "Alu"),
    ("ponAbe3", "Orangutan", "Alu"),
    ("rheMac10", "Macaque", "Alu"),
    ("calJac4", "Marmoset", "Alu"),
    ("mm10", "Mouse", "B1B2"),
]
TE_LABEL = {"Alu": "Alu", "B1B2": "B1/B2 (SINE)"}


def density(path):
    df = pd.read_csv(path, sep="\t", header=None)
    return (df.iloc[:, 6] / ((df.iloc[:, 2] - df.iloc[:, 1]) / 1000)).values


rng = np.random.default_rng(SEED)
NCOL3 = 3
NROW3 = -(-len(COMBOS) // NCOL3)
fig, axgrid3 = plt.subplots(NROW3, NCOL3, figsize=(5.2 * NCOL3, 4.2 * NROW3))
axes = axgrid3.ravel()
# Five panels in a 3x2 grid leave one cell over. It carries the legend, so the
# key is stated once instead of being redrawn inside every panel.
spare = axes[len(COMBOS):]
for extra in spare:
    extra.axis("off")

for ax, letter, (assembly, label, te) in zip(axes, "abcde", COMBOS):
    hk = density(RESULTS / assembly / f"Housekeeping_{te}.bed")
    ndd = density(RESULTS / assembly / f"HighConfNDD_{te}.bed")
    _, obs_p = stats.mannwhitneyu(hk, ndd, alternative="greater")

    pooled = np.concatenate([hk, ndd])
    n_hk = len(hk)
    perm_p = np.empty(N_PERM)
    for i in range(N_PERM):
        s = rng.permutation(pooled)
        _, perm_p[i] = stats.mannwhitneyu(s[:n_hk], s[n_hk:], alternative="greater")

    emp_p = float(np.mean(perm_p <= obs_p))
    null_fpr = float(np.mean(perm_p < 0.05))

    ax.hist(perm_p, bins=50, color="#4878CF", alpha=0.75, label="Random gene sets")
    ax.axvline(obs_p, color="red", linestyle="--", linewidth=2,
               label="Observed NDD set")
    ax.set_xlabel("$p$-value", fontsize=10)
    if ax in (axes[0], axes[NCOL3]):
        ax.set_ylabel("Frequency", fontsize=10)
    ax.set_title(f"{letter}   {label} | {TE_LABEL[te]}", loc="left",
                 fontweight="bold", fontsize=11)
    ax.text(0.97, 0.97, f"empirical $p$ = {emp_p:.4f}\nnull FPR = {null_fpr:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="#FFFDF0", edgecolor="#DDD"))
    ax.spines[["top", "right"]].set_visible(False)
    print(f"{label:10s} {te:5s} observed p={obs_p:.3e}  emp p={emp_p:.4f}  FPR={null_fpr:.3f}")

if len(spare):
    h, l = axes[0].get_legend_handles_labels()
    spare[0].legend(h, l, frameon=False, fontsize=11, loc="center",
                    title=f"{N_PERM:,} permutations\nper species",
                    title_fontproperties={"weight": "bold", "size": 11})

plt.tight_layout()
plt.savefig(FIGS / "ESM2_NullModel.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGS / "ESM2_NullModel.png", dpi=150, bbox_inches="tight")
plt.close()
print("Fig4 regenerated with 5 panels.")
