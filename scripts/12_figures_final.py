import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
import os, subprocess
from io import StringIO
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS = REPO_ROOT / "results"
DATA = REPO_ROOT / "data"
FIGS = REPO_ROOT / "figures"
os.makedirs(FIGS, exist_ok=True)

# Renk paleti
COLORS = {'Housekeeping': '#4878CF', 'HighConfNDD': '#E07B39'}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})

def get_density(filepath):
    if not os.path.exists(filepath): return None
    df = pd.read_csv(filepath, sep='\t', header=None)
    return (df.iloc[:,6] / ((df.iloc[:,2]-df.iloc[:,1])/1000)).values

def add_significance(ax, x1, x2, y, p):
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    ax.plot([x1,x1,x2,x2],[y*0.97,y,y,y*0.97],'k-',lw=0.8)
    ax.text((x1+x2)/2, y*1.01, sig, ha='center', va='bottom', fontsize=10)

# ══════════════════════════════════════════════════════
# FIG 1: Alu depletion — 5 primat
# ══════════════════════════════════════════════════════
print("Fig1 çiziliyor...")
primates = [
    ("hg38",     "Human\n(hg38)",       0),
    ("ponAbe3",  "Orangutan\n(ponAbe3)", 16),
    ("nomLeu3",  "Gibbon\n(nomLeu3)",    20),
    ("rheMac10", "Macaque\n(rheMac10)",  25),
    ("calJac4",  "Marmoset\n(calJac4)",  40),
]

fig, axes = plt.subplots(1, 5, figsize=(18, 6), sharey=False)
for ax, (sp, label, mya) in zip(axes, primates):
    hk  = get_density(f"{RESULTS}/{sp}/Housekeeping_Alu.bed")
    ndd = get_density(f"{RESULTS}/{sp}/HighConfNDD_Alu.bed")
    bp = ax.boxplot([hk, ndd], patch_artist=True,
                    medianprops=dict(color='black', linewidth=2),
                    flierprops=dict(marker='o', markersize=2, alpha=0.3),
                    widths=0.6)
    for patch, cat in zip(bp['boxes'], ['Housekeeping','HighConfNDD']):
        patch.set_facecolor(COLORS[cat]); patch.set_alpha(0.85)
    _, p = stats.mannwhitneyu(hk, ndd, alternative='greater')
    ymax = max(np.quantile(hk, 0.95), np.quantile(ndd, 0.95))
    add_significance(ax, 1, 2, ymax*1.18, p)
    ax.set_title(f"{label}\n~{mya} Mya", fontweight='bold', fontsize=10)
    ax.set_xticks([1,2])
    ax.set_xticklabels(['HK','NDD'], fontsize=10)
    if ax == axes[0]: ax.set_ylabel('Alu Frequency (count per kb)', fontsize=11)
    n_hk, n_ndd = len(hk), len(ndd)
    ax.text(0.5, -0.18, f'n={n_hk}, {n_ndd}', transform=ax.transAxes,
            ha='center', fontsize=8, color='gray')

patches = [mpatches.Patch(color=COLORS[c], label=c.replace('HighConfNDD','NDD'), alpha=0.85)
           for c in ['Housekeeping','HighConfNDD']]
fig.legend(handles=patches, loc='lower center', ncol=2,
           bbox_to_anchor=(0.5,-0.02), frameon=False, fontsize=10)
plt.suptitle('Alu Depletion at NDD Gene Promoters Across Primates\n'
             'TSS ± 2 kb | HighConfNDD vs Housekeeping | Mann-Whitney U',
             fontweight='bold', fontsize=12)
plt.tight_layout(rect=[0,0.05,1,0.95])
plt.savefig(f"{FIGS}/Fig1_Alu_Primates.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIGS}/Fig1_Alu_Primates.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Fig1 kaydedildi.")

# ══════════════════════════════════════════════════════
# FIG 2: LINE-1 depletion — 7 tür
# ══════════════════════════════════════════════════════
print("Fig2 çiziliyor...")
all_species = [
    ("hg38",     "Human",     0),
    ("ponAbe3",  "Orangutan", 16),
    ("nomLeu3",  "Gibbon",    20),
    ("rheMac10", "Macaque",   25),
    ("calJac4",  "Marmoset",  40),
    ("mm10",     "Mouse",     90),
    ("canFam4",  "Dog",       95),
]

fig, axes = plt.subplots(1, 7, figsize=(22, 6), sharey=False)
for ax, (sp, label, mya) in zip(axes, all_species):
    hk  = get_density(f"{RESULTS}/{sp}/Housekeeping_LINE1.bed")
    ndd = get_density(f"{RESULTS}/{sp}/HighConfNDD_LINE1.bed")
    bp = ax.boxplot([hk, ndd], patch_artist=True,
                    medianprops=dict(color='black', linewidth=2),
                    flierprops=dict(marker='o', markersize=2, alpha=0.3),
                    widths=0.6)
    for patch, cat in zip(bp['boxes'], ['Housekeeping','HighConfNDD']):
        patch.set_facecolor(COLORS[cat]); patch.set_alpha(0.85)
    _, p = stats.mannwhitneyu(hk, ndd, alternative='greater')
    ymax = max(np.quantile(hk, 0.95), np.quantile(ndd, 0.95))
    if ymax > 0:
        add_significance(ax, 1, 2, ymax*1.18, p)
    ax.set_title(f"{label}\n~{mya} Mya", fontweight='bold', fontsize=9)
    ax.set_xticks([1,2])
    ax.set_xticklabels(['HK','NDD'], fontsize=9)
    if ax == axes[0]: ax.set_ylabel('LINE-1 Frequency (count per kb)', fontsize=11)
    ax.text(0.5,-0.18, f'n={len(hk)},{len(ndd)}', transform=ax.transAxes,
            ha='center', fontsize=7, color='gray')

fig.legend(handles=patches, loc='lower center', ncol=2,
           bbox_to_anchor=(0.5,-0.02), frameon=False, fontsize=10)
plt.suptitle('LINE-1 Depletion at NDD Gene Promoters Across Mammals\n'
             'TSS ± 2 kb | HighConfNDD vs Housekeeping | Mann-Whitney U',
             fontweight='bold', fontsize=12)
plt.tight_layout(rect=[0,0.05,1,0.95])
plt.savefig(f"{FIGS}/Fig2_LINE1_Mammals.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIGS}/Fig2_LINE1_Mammals.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Fig2 kaydedildi.")

# ══════════════════════════════════════════════════════
# FIG 3: Significance Heatmap
# ══════════════════════════════════════════════════════
print("Fig3 çiziliyor...")
sp_info = [
    ("hg38",     "Human",     "Alu",  "LINE-1"),
    ("ponAbe3",  "Orangutan", "Alu",  "LINE-1"),
    ("nomLeu3",  "Gibbon",    "Alu",  "LINE-1"),
    ("rheMac10", "Macaque",   "Alu",  "LINE-1"),
    ("calJac4",  "Marmoset",  "Alu",  "LINE-1"),
    ("mm10",     "Mouse",     None,   "LINE-1"),
    ("canFam4",  "Dog",       None,   "LINE-1"),
]
te_cols = ["Alu", "LINE-1"]
sp_labels = [s[1] for s in sp_info]
pmat = np.full((7,2), np.nan)
sig_mat = []

for i, (sp, label, alu_te, l1_te) in enumerate(sp_info):
    row = []
    for j, te in enumerate(["Alu", "LINE1"]):
        # Mouse uses B1/B2 as the SINE-class analog of primate Alu
        te_file = te
        if sp == "mm10" and te == "Alu":
            te_file = "B1B2"
        fname = f"{RESULTS}/{sp}/HighConfNDD_{te_file}.bed"
        hname = f"{RESULTS}/{sp}/Housekeeping_{te_file}.bed"
        if not os.path.exists(fname) or not os.path.exists(hname):
            row.append("—"); continue
        hk  = get_density(hname)
        ndd = get_density(fname)
        _, p = stats.mannwhitneyu(hk, ndd, alternative='greater')
        pmat[i,j] = min(-np.log10(p), 60) if p > 0 else 60
        sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
        row.append(f"{sig}\np={p:.2e}")
    sig_mat.append(row)

fig, ax = plt.subplots(figsize=(7, 7))
masked = np.ma.masked_invalid(pmat)
im = ax.imshow(masked, cmap='RdYlGn', aspect='auto', vmin=0, vmax=20)
plt.colorbar(im, ax=ax, label='-log₁₀(p-value)\nHK > NDD')
ax.set_xticks([0,1]); ax.set_xticklabels(['Alu / B1-B2\n(SINE)', 'LINE-1'], fontsize=11)
ax.set_yticks(range(7)); ax.set_yticklabels(sp_labels, fontsize=12, fontweight='bold')

for i in range(7):
    for j in range(2):
        val = pmat[i,j]
        if np.isnan(val):
            ax.text(j, i, '—', ha='center', va='center', fontsize=11, color='white')
        else:
            color = 'white' if val > 10 else 'black'
            ax.text(j, i, sig_mat[i][j], ha='center', va='center',
                    fontsize=8, color=color, fontweight='bold')

ax.set_title('TE Depletion Significance: HK vs NDD Gene Promoters\n'
             'Across Species and TE Classes', fontweight='bold', fontsize=11)
plt.tight_layout()
plt.savefig(f"{FIGS}/Fig3_Heatmap.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIGS}/Fig3_Heatmap.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Fig3 kaydedildi.")

# ══════════════════════════════════════════════════════
# FIG 4: Null Model — Alu (4 primat: Human, Orangutan, Macaque, Mouse LINE-1)
# ══════════════════════════════════════════════════════
print("Fig4 çiziliyor...")
np.random.seed(42)
N_PERM = 10000

null_species = [
    ("hg38",     "Human",    "Alu"),
    ("ponAbe3",  "Orangutan","Alu"),
    ("rheMac10", "Macaque",  "Alu"),
    ("mm10",     "Mouse",    "B1B2"),
]

fig, axes = plt.subplots(1, 4, figsize=(16,5))
for ax, (sp, label, te) in zip(axes, null_species):
    hk  = get_density(f"{RESULTS}/{sp}/Housekeeping_{te}.bed")
    ndd = get_density(f"{RESULTS}/{sp}/HighConfNDD_{te}.bed")
    _, obs_p = stats.mannwhitneyu(hk, ndd, alternative='greater')
    bg = np.concatenate([hk, ndd])
    perm_ps = []
    for _ in range(N_PERM):
        s = np.random.permutation(bg)
        _, p = stats.mannwhitneyu(s[:len(hk)], s[len(hk):len(hk)+len(ndd)],
                                   alternative='greater')
        perm_ps.append(p)
    perm_ps = np.array(perm_ps)
    emp_p = np.mean(perm_ps <= obs_p)
    null_fpr = np.mean(perm_ps < 0.05)
    ax.hist(perm_ps, bins=50, color='#4878CF', alpha=0.7, label='Random sets')
    ax.axvline(obs_p, color='red', linestyle='--', linewidth=2, label='Observed NDD')
    ax.set_xlabel('p-value', fontsize=10)
    if ax == axes[0]: ax.set_ylabel('Frequency', fontsize=10)
    te_label = {"Alu": "Alu", "B1B2": "B1/B2 (SINE)", "LINE1": "LINE-1"}[te]
    ax.set_title(f"{label} | {te_label}", fontweight='bold', fontsize=11)
    info = f"Emp. p={emp_p:.4f}\nNull FPR={null_fpr:.3f}"
    ax.text(0.97, 0.97, info, transform=ax.transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax.legend(fontsize=8, loc='upper left')

plt.suptitle(f'Permutation Test: Observed NDD Depletion vs. Random Gene Sets\n'
             f'n={N_PERM:,} permutations',
             fontweight='bold', fontsize=12)
plt.tight_layout()
plt.savefig(f"{FIGS}/Fig4_NullModel.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIGS}/Fig4_NullModel.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Fig4 kaydedildi.")

# ══════════════════════════════════════════════════════
# FIG 5: CpG Confounder
# ══════════════════════════════════════════════════════
print("Fig5 çiziliyor...")
cpg_bed = f"{DATA}/hg38/cpg_islands.bed"
alu_bed = f"{DATA}/hg38/rmsk/Alu.bed"

fig, axes = plt.subplots(1, 2, figsize=(12,6))
for ax, (cpg_status, cpg_label) in zip(axes,
    [('with','CpG Island Present'), ('without','No CpG Island')]):
    plot_data = []
    ns = []
    for cat in ['Housekeeping','HighConfNDD']:
        # Use existing intersect results (already contains TE counts)
        alu_file = f"{RESULTS}/hg38/{cat}_Alu.bed"
        df_alu = pd.read_csv(alu_file, sep='\t', header=None)
        df_alu.columns = ['chr','start','end','gene','score','strand','alu_count']
        df_alu['density'] = df_alu['alu_count'] / ((df_alu['end']-df_alu['start'])/1000)
        # Add CpG overlap
        alu_file_path = alu_file
        cpg_out = subprocess.run(
            ["bedtools","intersect","-a",alu_file_path,"-b",cpg_bed,"-c"],
            capture_output=True, text=True).stdout
        df_cpg = pd.read_csv(StringIO(cpg_out), sep='\t', header=None)
        df_cpg.columns = ['chr','start','end','gene','score','strand','alu_count','cpg_n']
        df_cpg['density'] = df_cpg['alu_count'] / ((df_cpg['end']-df_cpg['start'])/1000)
        df_cpg['gene'] = df_cpg['gene']
        df = df_cpg[['gene','density','cpg_n']]
        if cpg_status == 'with':
            subset = df[df['cpg_n']>0]['density'].values
        else:
            subset = df[df['cpg_n']==0]['density'].values
        plot_data.append(subset)
        ns.append(len(subset))

    bp = ax.boxplot(plot_data, patch_artist=True,
                    medianprops=dict(color='black',linewidth=2),
                    flierprops=dict(marker='o',markersize=2,alpha=0.3),
                    widths=0.6)
    for patch, cat in zip(bp['boxes'],['Housekeeping','HighConfNDD']):
        patch.set_facecolor(COLORS[cat]); patch.set_alpha(0.85)
    _, p = stats.mannwhitneyu(plot_data[0], plot_data[1], alternative='greater')
    ymax = max(np.quantile(plot_data[0],0.95), np.quantile(plot_data[1],0.95))
    add_significance(ax, 1, 2, ymax*1.18, p)
    ax.set_title(f'{cpg_label}\nn = {ns[0]}, {ns[1]}', fontweight='bold')
    ax.set_xticks([1,2])
    ax.set_xticklabels(['Housekeeping','NDD'], fontsize=10)
    ax.set_ylabel('Alu Frequency (count per kb)' if ax==axes[0] else '')

fig.legend(handles=patches, loc='lower center', ncol=2,
           bbox_to_anchor=(0.5,-0.02), frameon=False, fontsize=10)
plt.suptitle('Alu Depletion: CpG Island Confounder Test\n'
             'hg38 | TSS ± 2 kb | HighConfNDD vs Housekeeping',
             fontweight='bold', fontsize=12)
plt.tight_layout(rect=[0,0.05,1,1])
plt.savefig(f"{FIGS}/Fig5_CpG.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{FIGS}/Fig5_CpG.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Fig5 kaydedildi.")

print("\n=== TÜM FİGÜRLER TAMAMLANDI ===")
print(f"Konum: {FIGS}")
