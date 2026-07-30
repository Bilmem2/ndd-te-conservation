#!/usr/bin/env python3
"""Supplementary figure: SINE-class effect size against divergence time.

Not a phylogenetic regression - with one measurement per species and seven taxa
there are too few degrees of freedom for that. The point of the figure is the
argument made in the Discussion: the depletion does not track divergence time
(marmoset at 40 Mya is the weakest, the mouse lemur at 70 Mya is not), and the
one lineage whose SINEs amplified independently of primate Alu, the mouse, shows
the deficit as strongly as any primate."""
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
FIGS = REPO / "figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

df = pd.read_csv(REPO / "results" / "consolidated" / "cross_species.csv")
sine = df[df["TE"].isin(["Alu", "B1B2"])].copy()
sine["label"] = sine["species"].replace({"MouseLemur": "Mouse lemur", "SquirrelMonkey": "Squirrel monkey"})

prim = sine[sine["TE"] == "Alu"]
rod = sine[sine["TE"] == "B1B2"]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(prim["mya"], prim["r"], "o-", color="#E07B39", lw=1.2, ms=9,
        markeredgecolor="black", markeredgewidth=0.6,
        label="Primates (Alu, shared ancestry)")
ax.plot(rod["mya"], rod["r"], "D", color="#5b53b0", ms=11,
        markeredgecolor="black", markeredgewidth=0.6,
        label="Mouse (B1/B2, independently amplified)")

for _, r in sine.iterrows():
    ax.annotate(r["label"], (r["mya"], r["r"]), textcoords="offset points",
                xytext=(0, 11), ha="center", fontsize=9)

ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("Divergence time from human (Mya)")
ax.set_ylabel("SINE depletion  (rank-biserial $r$, HK > NDD)")
ax.set_ylim(-0.48, 0.06)
ax.set_xlim(-5, 100)
ax.legend(fontsize=9, loc="lower left")
ax.grid(alpha=0.25, ls=":")

plt.tight_layout()
plt.savefig(FIGS / "FigS5_PhyloEffect.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGS / "FigS5_PhyloEffect.png", dpi=150, bbox_inches="tight")
plt.close()
print("FigS5 saved:", sine[["label", "mya", "TE", "r"]].to_dict("records"))
