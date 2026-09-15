"""
fig_heatmap_effect.py — the cross-species heatmap, coloured by effect size.

The published version of this figure colours each cell by $-\\log_{10}(p)$. That
conflates two different things: with these sample sizes a negligible LINE-1
difference can reach a very small p, so the LINE-1 column reads as strongly as
the SINE column even though its effect sizes are an order of magnitude smaller.
Colouring by rank-biserial r instead puts the magnitude of the difference on the
colour axis and leaves significance to the markers printed in each cell, which
are retained.

The scale is diverging and symmetric about zero, so equal magnitudes carry equal
visual weight whichever side of zero they fall. Stretching the narrow positive
range to fill the blue half would make an r of +0.02 look as strong as a
depletion of -0.2, which is the confusion this figure is being redrawn to avoid.
The blue half therefore goes largely unused, and that is itself the result:
nothing in the panel is enriched.

Sign convention matches the rest of the paper: r is computed from the one-sided
Mann-Whitney U test of HK > NDD, so negative means depleted at NDD promoters.

Inputs : results/<species>/{HighConfNDD,Housekeeping}_{Alu,B1B2,LINE1}.bed
Output : figures/Fig2_Heatmap.{pdf,png}
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results"
FIGS = REPO / "figures"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

SPECIES = [
    ("hg38", "Human"), ("ponAbe3", "Orangutan"), ("nomLeu3", "Gibbon"),
    ("rheMac10", "Macaque"), ("calJac4", "Marmoset"),
    ("saiBol1", "Squirrel monkey"), ("mmur3", "Mouse lemur"),
    ("mm10", "Mouse"), ("canFam6", "Dog"),
]
def dens(path):
    """Element frequency per kb.

    Two bed layouts coexist in results/: the primate and rodent files carry
    BED6 plus a count (chrom, start, end, gene, score, strand, n) while the
    mouse lemur files carry chrom, start, end, gene, category, n. Reading with
    a fixed column list silently misaligns the second form, so the count is
    taken as the last column and the interval from the first three.
    """
    d = pd.read_csv(path, sep="\t", header=None)
    start, end, count = d.iloc[:, 1], d.iloc[:, 2], d.iloc[:, d.shape[1] - 1]
    return (count / ((end - start) / 1000.0)).to_numpy()


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


n = len(SPECIES)
rmat = np.full((n, 2), np.nan)
cells = []
for i, (sp, label) in enumerate(SPECIES):
    row = []
    for j, te in enumerate(["Alu", "LINE1"]):
        te_file = "B1B2" if (sp == "mm10" and te == "Alu") else te
        fn = RES / sp / f"HighConfNDD_{te_file}.bed"
        hn = RES / sp / f"Housekeeping_{te_file}.bed"
        if not (fn.exists() and hn.exists()):
            row.append("—")
            continue
        hk, nd = dens(hn), dens(fn)
        u = stats.mannwhitneyu(hk, nd, alternative="greater")
        r = 1 - (2 * u.statistic) / (len(hk) * len(nd))
        rmat[i, j] = r
        row.append(f"{r:+.3f}\n{stars(u.pvalue)}  p={u.pvalue:.1e}")
    cells.append(row)

# Symmetric about zero. A two-slope norm would stretch the thin positive range
# over the whole blue half and make an r of +0.02 look as strong as a depletion
# of -0.2, which is the very confusion this figure is being redrawn to avoid.
# The unused blue half is itself informative: nothing here is enriched.
lim = np.ceil(np.nanmax(np.abs(rmat)) * 20) / 20
vmin, vmax = -lim, lim
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

fig, ax = plt.subplots(figsize=(7.2, 7.8))
im = ax.imshow(np.ma.masked_invalid(rmat), cmap="RdBu", norm=norm, aspect="auto")
cb = plt.colorbar(im, ax=ax, extend="neither")
cb.set_label("Rank-biserial $r$   (negative = depleted at NDD promoters)",
             fontsize=10)

ax.set_xticks([0, 1])
ax.set_xticklabels(["Alu / B1-B2\n(SINE)", "LINE-1"], fontsize=11)
ax.set_yticks(range(n))
ax.set_yticklabels([s[1] for s in SPECIES], fontsize=12, fontweight="bold")

for i in range(n):
    for j in range(2):
        v = rmat[i, j]
        if np.isnan(v):
            ax.text(j, i, "—", ha="center", va="center", fontsize=14, color="black")
        else:
            # white text only on the saturated end of the scale
            col = "white" if v <= vmin * 0.55 else "black"
            ax.text(j, i, cells[i][j], ha="center", va="center", fontsize=8,
                    color=col, fontweight="bold", linespacing=1.35)

ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.4)
ax.tick_params(which="minor", length=0)

plt.tight_layout()
fig.savefig(FIGS / "Fig2_Heatmap.pdf", dpi=300, bbox_inches="tight")
fig.savefig(FIGS / "Fig2_Heatmap.png", dpi=150, bbox_inches="tight")

out = pd.DataFrame({"species": [s[1] for s in SPECIES],
                    "r_SINE": rmat[:, 0].round(4),
                    "r_LINE1": rmat[:, 1].round(4)})
print(out.to_string(index=False))
print(f"\ncolour scale: {vmin:.2f} to {vmax:.2f}, centred on 0")
print(f"SINE column  |r| range: {np.nanmin(np.abs(rmat[:, 0])):.3f}"
      f" to {np.nanmax(np.abs(rmat[:, 0])):.3f}")
print(f"LINE-1 column |r| range: {np.nanmin(np.abs(rmat[:, 1])):.3f}"
      f" to {np.nanmax(np.abs(rmat[:, 1])):.3f}")
print(f"\nwrote {(FIGS / 'Fig2_Heatmap.pdf').relative_to(REPO)}")
