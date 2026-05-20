import pandas as pd
import numpy as np
from scipy import stats
import os

RESULTS = os.path.expanduser("~/comparative_TE_study/results")

species = {
    "hg38":     {"label": "Human",     "mya": 0,  "alu": True},
    "ponAbe3":  {"label": "Orangutan", "mya": 16, "alu": True},
    "nomLeu3":  {"label": "Gibbon",    "mya": 20, "alu": True},
    "rheMac10": {"label": "Macaque",   "mya": 25, "alu": True},
    "calJac4":  {"label": "Marmoset",  "mya": 40, "alu": True},
    "mm10":     {"label": "Mouse",     "mya": 90, "alu": False},
    "canFam4":  {"label": "Dog",       "mya": 95, "alu": False},
}

def get_density(filepath):
    if not os.path.exists(filepath):
        return None
    df = pd.read_csv(filepath, sep='\t', header=None)
    return (df.iloc[:,6] / ((df.iloc[:,2] - df.iloc[:,1]) / 1000)).values

results = []

for sp, info in species.items():
    for te in ["Alu", "LINE1"]:
        if te == "Alu" and not info["alu"]:
            continue
        hk  = get_density(f"{RESULTS}/{sp}/Housekeeping_{te}.bed")
        ndd = get_density(f"{RESULTS}/{sp}/HighConfNDD_{te}.bed")
        if hk is None or ndd is None:
            continue
        u, p = stats.mannwhitneyu(hk, ndd, alternative='greater')
        r = 1 - (2*u)/(len(hk)*len(ndd))
        sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
        results.append({
            "species": info["label"],
            "mya": info["mya"],
            "TE": te,
            "n_HK": len(hk),
            "n_NDD": len(ndd),
            "median_HK": round(np.median(hk), 4),
            "median_NDD": round(np.median(ndd), 4),
            "p_value": p,
            "r": round(r, 3),
            "sig": sig
        })

df = pd.DataFrame(results)
df.to_csv(f"{RESULTS}/statistics_final.csv", index=False)

print("=== ALU (5 primat) ===")
alu = df[df['TE']=='Alu'].sort_values('mya')
print(alu[['species','mya','n_HK','n_NDD','median_HK','median_NDD','p_value','r','sig']].to_string(index=False))

print("\n=== LINE-1 (7 tür) ===")
l1 = df[df['TE']=='LINE1'].sort_values('mya')
print(l1[['species','mya','n_HK','n_NDD','median_HK','median_NDD','p_value','r','sig']].to_string(index=False))
