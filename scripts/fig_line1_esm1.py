r"""ESM 1: the two element classes side by side, same species, same construction.

The point of showing LINE-1 in the supplement is the contrast with the SINEs, so
the figure shows both and lets the reader see it: a visible gap between the two
bars in every SINE row, and no gap in any LINE-1 row.

Two earlier attempts are worth recording as things not to do. The nine-panel
boxplot gave each species its own y-axis, so a figure whose only purpose is a
comparison between species could not support one, and because most promoters
carry no LINE-1 the quartiles collapsed to flat lines. Replacing it with LINE-1
bars plus SINE tick marks on a shared axis was no better: the axis then had to
reach the SINE densities, which squeezed the LINE-1 bars into the left third.
Each class gets its own axis here, and the element scored is named on every row.
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
WIN_KB, RNG = 4.0, np.random.default_rng(20260915)
HK_C, NDD_C = "#5B8FD4", "#E08B4F"

PANEL = [("Human", "hg38", 0, "Alu"), ("Orangutan", "ponAbe3", 16, "Alu"),
         ("Gibbon", "nomLeu3", 20, "Alu"), ("Macaque", "rheMac10", 25, "Alu"),
         ("Marmoset", "calJac4", 40, "Alu"), ("Squirrel monkey", "saiBol1", 40, "Alu"),
         ("Mouse lemur", "mmur3", 70, "Alu"), ("Mouse", "mm10", 90, "B1B2"),
         ("Dog", "canFam6", 95, "CanSINE")]
PRETTY = {"Alu": "Alu", "B1B2": "B1/B2", "CanSINE": "Can-SINE", "LINE1": "LINE-1"}


def dens(asm, grp, te):
    df = pd.read_csv(R / asm / f"{grp}_{te}.bed", sep="\t", header=None)
    return df.iloc[:, -1].to_numpy(float) / WIN_KB


def ci(v, n=2000):
    return np.percentile(RNG.choice(v, (n, len(v)), replace=True).mean(axis=1),
                         [2.5, 97.5])


rows = []
for name, asm, mya, sine in PANEL:
    rec = {"name": name, "mya": mya, "sine_lbl": PRETTY[sine]}
    for key, te in (("sine", sine), ("l1", "LINE1")):
        h, d = dens(asm, "Housekeeping", te), dens(asm, "HighConfNDD", te)
        rec[key] = (h.mean(), ci(h), d.mean(), ci(d))
    rows.append(rec)
    print(f"  {name:<16} {PRETTY[sine]:<8} HK={rec['sine'][0]:.3f} NDD={rec['sine'][2]:.3f}"
          f"   LINE-1 HK={rec['l1'][0]:.3f} NDD={rec['l1'][2]:.3f}")

fig, (axS, axL) = plt.subplots(1, 2, figsize=(12.4, 5.6), sharey=True,
                               gridspec_kw=dict(wspace=0.08))
ys = np.arange(len(rows))[::-1]
H = 0.36


def draw(ax, key, title, xmax):
    for y, r in zip(ys, rows):
        hm, hci, dm, dci = r[key]
        ax.barh(y + H / 2, hm, H, color=HK_C, edgecolor="black", linewidth=0.6, zorder=2)
        ax.barh(y - H / 2, dm, H, color=NDD_C, edgecolor="black", linewidth=0.6, zorder=2)
        ax.plot(hci, [y + H / 2] * 2, color="#222", linewidth=1.1, zorder=3)
        ax.plot(dci, [y - H / 2] * 2, color="#222", linewidth=1.1, zorder=3)
        # Only the SINE panel needs a per-row label: the family scored changes
        # between lineages. Panel B is LINE-1 on every row and says so once.
        if key == "sine":
            ax.text(0.98, y, r["sine_lbl"], transform=ax.get_yaxis_transform(),
                    va="center", ha="right", fontsize=8.5, color="#888",
                    style="italic")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=12, pad=8)
    ax.set_xlabel("Mean density (elements per kb)")
    ax.set_xlim(0, xmax)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#EEE", linewidth=0.8)
    ax.set_axisbelow(True)


top = max(max(r["sine"][0], r["sine"][2]) for r in rows)
draw(axS, "sine", "a   SINE-class elements", top * 1.30)
draw(axL, "l1", "b   LINE-1", top * 1.30)

axS.set_yticks(ys)
axS.set_yticklabels([f"{r['name']}\n$\\sim${r['mya']} Mya" for r in rows], fontsize=9.5)
axS.set_ylim(-0.7, len(rows) - 0.3)
fig.legend(handles=[
    plt.Rectangle((0, 0), 1, 1, fc=HK_C, ec="black", lw=0.6, label="Housekeeping"),
    plt.Rectangle((0, 0), 1, 1, fc=NDD_C, ec="black", lw=0.6, label="NDD")],
    frameon=False, fontsize=10, ncol=2, loc="lower center",
    bbox_to_anchor=(0.5, -0.06))

fig.savefig(FIG / "ESM1_LINE1_Mammals.pdf", bbox_inches="tight")
fig.savefig(FIG / "ESM1_LINE1_Mammals.png", dpi=200, bbox_inches="tight")
print("\nyazildi: ESM1_LINE1_Mammals.pdf / .png")
