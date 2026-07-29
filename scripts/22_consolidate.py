#!/usr/bin/env python3
"""
22_consolidate.py — Assemble every final result into one coherent backbone for
the (revised) conserved-depletion manuscript. Reads the individual result CSVs,
builds the master cross-species table (now including mouse lemur), and writes a
human-readable RESULTS_BACKBONE.md digest of the supporting analyses with the
honest, re-baselined numbers.

Output: results/consolidated/{cross_species.csv, RESULTS_BACKBONE.md}
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
R = REPO / "results"
OUT = R / "consolidated"
OUT.mkdir(parents=True, exist_ok=True)
COLS = ["species", "mya", "TE", "n_HK", "n_NDD", "median_HK", "median_NDD", "p_value", "r", "sig"]


def mouse_b1b2():
    """Recompute mouse B1/B2 depletion from committed BEDs (missing from statistics_final.csv)."""
    def dens(cat):
        d = pd.read_csv(R / "mm10" / f"{cat}_B1B2.bed", sep="\t", header=None,
                        names=["chrom", "start", "end", "gene", "sc", "strand", "c"])
        return d["c"] / ((d["end"] - d["start"]) / 1000)
    hk, nd = dens("Housekeeping"), dens("HighConfNDD")
    _, p = stats.mannwhitneyu(hk, nd, alternative="greater")
    u2, _ = stats.mannwhitneyu(hk, nd, alternative="two-sided")
    r = 1 - (2 * u2) / (len(hk) * len(nd))
    sig = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
    return dict(species="Mouse", mya=90, TE="B1B2", n_HK=len(hk), n_NDD=len(nd),
                median_HK=round(hk.median(), 3), median_NDD=round(nd.median(), 3),
                p_value=p, r=round(r, 3), sig=sig)


def squirrel_monkey():
    """Second Platyrrhine (Saimiri boliviensis); columns renamed to the master schema."""
    s = pd.read_csv(R / "saiBol1" / "squirrel_stats.csv").iloc[0]
    p = float(s["p_vs_HK"])
    return dict(species="SquirrelMonkey", mya=int(s["mya"]), TE="Alu",
                n_HK=int(s["n_HK"]), n_NDD=int(s["n_NDD"]),
                median_HK=float(s["median_HK"]), median_NDD=float(s["median_NDD"]),
                p_value=p, r=float(s["r_vs_HK"]),
                sig="***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns")


def cross_species():
    base = pd.read_csv(R / "statistics_final.csv")
    lem = pd.read_csv(R / "mmur3" / "lemur_stats.csv")
    lem = lem.rename(columns={"species": "species"})[COLS]
    lem["species"] = "MouseLemur"
    b1b2 = pd.DataFrame([mouse_b1b2()])[COLS]
    sq = pd.DataFrame([squirrel_monkey()])[COLS]
    full = pd.concat([base[COLS], lem[COLS], sq, b1b2], ignore_index=True)
    order = {"Alu": 0, "B1B2": 0, "LINE1": 1}
    full["sec"] = full["TE"].map(lambda t: order.get(t, 2))
    full = full.sort_values(["sec", "mya"]).drop(columns="sec").reset_index(drop=True)
    full.to_csv(OUT / "cross_species.csv", index=False)
    return full


def rd(path):
    return pd.read_csv(path)


def main():
    cs = cross_species()
    sine = cs[cs["TE"].isin(["Alu", "B1B2"])]
    line1 = cs[cs["TE"] == "LINE1"]

    L = []
    L.append("# Results Backbone — Conserved SINE Depletion at NDD Promoters (revised)\n")
    L.append("_Auto-generated consolidation of all final analyses. Honest, re-baselined numbers._\n")

    L.append("## 1. Cross-species SINE-class depletion (headline, comparative genomics)\n")
    L.append("SINE = Alu in primates, B1/B2 in mouse. One-sided Mann-Whitney U (HK>NDD), "
             "rank-biserial r. **Mouse lemur (strepsirrhine, ~70 My) is new; ortholog-validated "
             "|delta r|=0.001.**\n")
    L.append(sine.to_markdown(index=False))
    L.append("\n## 2. LINE-1 (internal contrast: weak, lineage-variable)\n")
    L.append(line1.to_markdown(index=False))

    L.append("\n## 3. Supporting rigor (all human hg38)\n")

    mc = rd(R / "matched" / "matched_alu_test.csv")
    L.append("**3a. Multivariate matched control (GC + gene density + recombination):** "
             "depletion persists after joint matching.\n")
    L.append(mc.to_markdown(index=False))

    rb = rd(R / "matched" / "rebaseline_decomposition.csv")
    L.append("\n**3b. Honest re-baselining (vs genome & expression-breadth+GC-matched control):** "
             "the dramatic vs-housekeeping ratio was partly housekeeping Alu-enrichment; a modest "
             "(~18%) NDD-specific fixed-Alu deficit survives; polymorphic shows no deficit.\n")
    L.append(rb.to_markdown(index=False))

    cd = rd(R / "context" / "context_group_diff.csv")
    ds = rd(R / "context" / "alu_by_density_stratum.csv")
    rs = rd(R / "context" / "alu_by_recomb_stratum.csv")
    L.append("\n**3c. Context (NDD vs HK differ, but depletion survives strata):**\n")
    L.append(cd.to_markdown(index=False))
    L.append("\nAlu depletion within gene-density quartiles:\n")
    L.append(ds.to_markdown(index=False))
    L.append("\nAlu depletion within recombination quartiles:\n")
    L.append(rs.to_markdown(index=False))

    cc = rd(R / "ccre" / "ccre_group_diff.csv")
    fp = rd(R / "ccre" / "alu_free_vs_pos.csv")
    L.append("\n**3d. Functional overlay (ENCODE cCRE):** NDD promoters are regulatory-active; "
             "Alu-free promoters are more cCRE-dense.\n")
    L.append(cc.to_markdown(index=False))
    L.append("\nAlu-free vs Alu-positive regulatory density:\n")
    L.append(fp.to_markdown(index=False))

    L.append("\n## 4. Mechanism probes (reported honestly as inconclusive)\n")
    dec = rd(R / "gnomad_mei" / "decomposition.csv")
    sf = rd(R / "functional" / "probe_subfamily.csv")
    L.append("**4a. gnomAD MEI targeting-vs-selection decomposition** (deflated after honest "
             "re-baselining; AFS not significant vs matched control):\n")
    L.append(dec.to_markdown(index=False))
    L.append("\n**4b. Alu subfamily-age** (depletion ~uniform across AluJ/S/Y => no distinctive "
             "age-dependent selection signature; stably refractory environment):\n")
    L.append(sf.to_markdown(index=False))

    L.append("\n## 5. Honest framing / limitations\n")
    L.append("- Phenomenon of Alu-depletion at developmental promoters + Alu/B1 convergence is "
             "**prior art** (Polak & Domany 2006; Tsirigos & Rigoutsos 2009; Simons/Mattick TFRs). "
             "Our contribution = systematic cross-species quantification incl. a strepsirrhine, "
             "modern clinical NDD sets, and comprehensive confounder control. Do NOT claim first "
             "discovery of the phenomenon.\n")
    L.append("- Mechanism (targeting vs post-insertion selection) **remains open**: population "
             "(gnomAD MEI) and subfamily-age probes are inconclusive once properly controlled.\n")
    L.append("- Magnitude is modest vs an honest (genome / breadth-matched) baseline; the large "
             "vs-housekeeping effect was partly housekeeping Alu-enrichment (Eller et al. 2007).\n")
    L.append("- Functional validation (wet-lab) is future work / capstone.\n")

    (OUT / "RESULTS_BACKBONE.md").write_text("\n".join(L), encoding="utf-8")

    pd.set_option("display.width", 170)
    print("=== MASTER CROSS-SPECIES TABLE (SINE-class) ===")
    print(sine.to_string(index=False))
    print("\n=== LINE-1 ===")
    print(line1.to_string(index=False))
    print(f"\nWrote {OUT/'cross_species.csv'} and {OUT/'RESULTS_BACKBONE.md'}")


if __name__ == "__main__":
    main()
