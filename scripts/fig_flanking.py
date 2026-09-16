"""
fig_flanking.py — the promoter-versus-region figure.

Panel A plots the NDD/housekeeping Alu ratio against distance from the TSS. The
promoter sits well below the flanks, and the flanks settle onto a plateau that
never reaches parity out to 250 kb: the deficit has a regional part and a
promoter-local part, and the figure shows both at once.

Panel B states the split. The promoter ratio is exactly the product of the
background ratio and the promoter-to-background ratio, so the 43 % promoter
deficit factors into a 27 % regional component and a 22 % promoter-local one.

Inputs : results/hg38/flanking_control.csv      (43_flanking_control.py)
         results/hg38/promoter_vs_local.csv     (45_promoter_vs_local.py)
Output : figures/Fig6_Flanking.{pdf,png}
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
REPO = Path(__file__).resolve().parents[1]
FIGS = REPO / "figures"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})

UP, DOWN, PROM = "#4878CF", "#E07B39", "#333333"
RING = "10-100 kb"

fl = pd.read_csv(REPO / "results/hg38/flanking_control.csv")
pv = pd.read_csv(REPO / "results/hg38/promoter_vs_local.csv")
row = pv[pv.ring == RING].iloc[0]

prom = fl[fl.side == "promoter"].iloc[0]
up = fl[fl.side == "upstream"].sort_values("offset_kb")
dn = fl[fl.side == "downstream"].sort_values("offset_kb")

fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(12.4, 5.0), gridspec_kw={"width_ratios": [1.65, 1]})

# ── Panel A: ratio against distance ─────────────────────────────────────────
# Offsets are plotted at even spacing rather than to scale: the sampled
# distances span two orders of magnitude and a linear axis collapses everything
# below 100 kb into the left margin, hiding the plateau that is the point.
OFFS = list(up.offset_kb)
pos = {o: i + 1 for i, o in enumerate(OFFS)}

axA.axhline(1.0, color="black", lw=0.9, ls="--", zorder=1)
axA.text(len(OFFS) + 0.42, 1.006, "no depletion", fontsize=8.5, color="black",
         ha="right", va="bottom")

# the band the flanks settle into, taken from the flank points themselves
lo = min(up.ratio.min(), dn.ratio.min())
hi = max(up.ratio.max(), dn.ratio.max())
axA.axhspan(lo, hi, color="#BBBBBB", alpha=0.22, zorder=0)
axA.text(0.30, (lo + hi) / 2, "regional\nplateau", fontsize=8.5, color="#555555",
         ha="center", va="center", style="italic", linespacing=1.3)

axA.plot([pos[o] for o in up.offset_kb], up.ratio, "o-", color=UP, lw=1.6, ms=8,
         markeredgecolor="black", markeredgewidth=0.6,
         label="upstream flank", zorder=3)
axA.plot([pos[o] for o in dn.offset_kb], dn.ratio, "s-", color=DOWN, lw=1.6, ms=8,
         markeredgecolor="black", markeredgewidth=0.6,
         label="downstream flank", zorder=3)
axA.plot([0], [prom.ratio], "D", color=PROM, ms=11,
         markeredgecolor="black", markeredgewidth=0.7,
         label="promoter (TSS ± 2 kb)", zorder=4)

axA.annotate(f"{prom.ratio:.3f}", (0, prom.ratio), textcoords="offset points",
             xytext=(11, -3), fontsize=9.5, fontweight="bold", color=PROM)
for s, d in ((up, 11), (dn, -13)):
    axA.annotate(f"{s.ratio.iloc[-1]:.3f}", (len(OFFS), s.ratio.iloc[-1]),
                 textcoords="offset points", xytext=(4, d), fontsize=9)

axA.set_xlabel("Distance of 4 kb window from TSS (kb, not to scale)")
axA.set_ylabel("Alu density ratio, NDD / housekeeping")
axA.set_xlim(-0.55, len(OFFS) + 0.55)
axA.set_ylim(0.50, 1.05)
axA.set_xticks([0] + [pos[o] for o in OFFS])
axA.set_xticklabels(["0"] + [str(int(o)) for o in OFFS])
axA.legend(fontsize=9.5, loc="lower right", framealpha=0.95)
axA.set_title("a   Depletion by distance from the TSS", loc="left",
              fontsize=12, fontweight="bold")
axA.grid(axis="y", color="#DDDDDD", lw=0.6, zorder=0)
axA.set_axisbelow(True)

# ── Panel B: the multiplicative split ───────────────────────────────────────
total = 1 - row.prom_ratio
regional = 1 - row.bg_ratio
local = 1 - row.norm_ratio_of_means

axB.barh([2], [regional], color="#9AA7B5", edgecolor="black", lw=0.7, height=0.52)
axB.barh([1], [local], color="#C98B5E", edgecolor="black", lw=0.7, height=0.52)
axB.barh([0], [total], color=PROM, edgecolor="black", lw=0.7, height=0.52)

for y, v, lab in ((2, regional, f"{regional * 100:.0f} %"),
                  (1, local, f"{local * 100:.0f} %"),
                  (0, total, f"{total * 100:.0f} %")):
    axB.text(v + 0.012, y, lab, va="center", fontsize=10.5, fontweight="bold")

axB.set_yticks([2, 1, 0])
axB.set_yticklabels(["Regional\nbackground", "Promoter-local\n(normalised)",
                     "Promoter\n(observed)"], fontsize=10)
axB.set_xlabel("Alu depletion relative to housekeeping genes")
axB.set_xlim(0, 0.56)
axB.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
axB.set_xticklabels(["0", "10 %", "20 %", "30 %", "40 %", "50 %"])
axB.set_title("b   The promoter deficit factors in two", loc="left",
              fontsize=12, fontweight="bold")
axB.grid(axis="x", color="#DDDDDD", lw=0.6)
axB.set_axisbelow(True)
axB.invert_yaxis()

axB.text(0.5, -0.19,
         f"{row.bg_ratio:.3f} × {row.norm_ratio_of_means:.3f}"
         f" = {row.prom_ratio:.3f}\nbackground ring {RING};"
         f" rank-biserial r falls {row.r_prom:.3f} → {row.r_norm:.3f}",
         transform=axB.transAxes, ha="center", va="top", fontsize=9.5,
         color="#444444")

fig.tight_layout()
fig.subplots_adjust(bottom=0.24)
fig.savefig(FIGS / "Fig6_Flanking.pdf", dpi=300, bbox_inches="tight")
fig.savefig(FIGS / "Fig6_Flanking.png", dpi=150, bbox_inches="tight")
print(f"wrote {(FIGS / 'Fig6_Flanking.pdf').relative_to(REPO)}")
print(f"wrote {(FIGS / 'Fig6_Flanking.png').relative_to(REPO)}")
print(f"\nsplit: {row.bg_ratio:.3f} x {row.norm_ratio_of_means:.3f}"
      f" = {row.prom_ratio:.3f}")
