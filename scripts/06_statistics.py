import pandas as pd
import numpy as np
from scipy import stats
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS = REPO_ROOT / "results"

species = {
    "hg38":     {"label": "Human",     "line1": True, "alu": True},
    "rheMac10": {"label": "Macaque",   "line1": True, "alu": True},
    "calJac4":  {"label": "Marmoset",  "line1": True, "alu": True},
    "ponAbe3":  {"label": "Orangutan", "line1": True, "alu": True},
    "mm10":     {"label": "Mouse",     "line1": True, "alu": False},
    "canFam4":  {"label": "Dog",       "line1": True, "alu": False},
}

categories = ["HighConfNDD", "Housekeeping"]
results = []

def get_density(filepath):
    if not os.path.exists(filepath):
        return None
    df = pd.read_csv(filepath, sep='\t', header=None)
    density = df.iloc[:,6] / ((df.iloc[:,2] - df.iloc[:,1]) / 1000)
    return density.values

for sp, info in species.items():
    for te in ["LINE1", "Alu"]:
        if te == "Alu" and not info["alu"]:
            continue

        hk = get_density(f"{RESULTS}/{sp}/Housekeeping_{te}.bed")
        ndd = get_density(f"{RESULTS}/{sp}/HighConfNDD_{te}.bed")

        if hk is None or ndd is None:
            continue

        u, p = stats.mannwhitneyu(hk, ndd, alternative='greater')
        r = 1 - (2*u) / (len(hk) * len(ndd))
        sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"

        results.append({
            "species": info["label"],
            "sp_code": sp,
            "TE": te,
            "n_HK": len(hk),
            "n_NDD": len(ndd),
            "median_HK": round(np.median(hk), 5),
            "median_NDD": round(np.median(ndd), 5),
            "p_value": round(p, 6),
            "r": round(r, 3),
            "sig": sig
        })

df = pd.DataFrame(results)
df.to_csv(f"{RESULTS}/statistics_main.csv", index=False)

print("=== LINE-1 SONUÇLAR ===")
l1 = df[df['TE']=='LINE1'][['species','n_HK','n_NDD','median_HK','median_NDD','p_value','r','sig']]
print(l1.to_string(index=False))

print("\n=== Alu SONUÇLAR (primatlar) ===")
alu = df[df['TE']=='Alu'][['species','n_HK','n_NDD','median_HK','median_NDD','p_value','r','sig']]
print(alu.to_string(index=False))
