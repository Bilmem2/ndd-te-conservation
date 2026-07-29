#!/usr/bin/env python3
"""
30_squirrel_monkey.py — a second Platyrrhine, to test whether the marmoset
exception is species-specific or shared across New World monkeys.

Marmoset is the one species in the panel whose NDD promoters are not depleted
relative to its own genome baseline. With a single Platyrrhine analysed that
observation is ambiguous: it could be a marmoset idiosyncrasy or a property of
the whole infraorder. Adding the Bolivian squirrel monkey (Saimiri boliviensis,
SaiBol1.0 / UCSC saiBol1), which diverged from marmoset roughly 20 Mya and from
human at about the same time as marmoset, distinguishes the two.

Sequence naming differs between the sources: the Ensembl GTF carries versioned
INSDC accessions (JH378105.1) while the UCSC RepeatMasker track uses the
unversioned form (JH378105), so the version suffix is stripped before matching
(362 of 366 annotated sequences map; the remainder are tiny contigs and MT).

Output: results/saiBol1/{HighConfNDD_Alu.bed, Housekeeping_Alu.bed,
        squirrel_stats.csv}
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
SB = REPO / "data" / "saiBol1"
GL = REPO / "data" / "gene_lists"
OUT = REPO / "results" / "saiBol1"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
N_BOOT = 4000


def gene_set(fn):
    return {l.strip() for l in open(GL / fn) if l.strip()}


def load_alu():
    rows = []
    with gzip.open(SB / "rmsk" / "rmsk.txt.gz", "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 12 and f[11] == "SINE" and f[12] == "Alu":
                rows.append((f[5], int(f[6]), int(f[7])))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def build_promoters():
    """All protein-coding promoters; seq names stripped of the version suffix."""
    rows, seen = [], set()
    with gzip.open(SB / "gtf" / "saiBol1.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or "protein_coding" not in f[8]:
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) in seen:
                continue
            seen.add(m.group(1))
            start, end, strand = int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            rows.append((f[0].split(".")[0], max(0, tss - 2000), tss + 2000, m.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def count_alu(prom, alu):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in alu.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        o = np.argsort(m["start"].values)
        s, e = m["start"].values[o], m["end"].values[o]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def rb(a, b):
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 1 - (2 * u) / (len(a) * len(b))


def ratio_ci(group, baseline, seed):
    rng = np.random.default_rng(seed)
    b = np.empty(N_BOOT)
    for k in range(N_BOOT):
        b[k] = (rng.choice(group, len(group), True).mean() /
                rng.choice(baseline, len(baseline), True).mean())
    return np.nanpercentile(b, [2.5, 97.5])


def main():
    ndd, hk = gene_set("HighConfNDD_genes.txt"), gene_set("Housekeeping_genes.txt")
    alu = load_alu()
    prom = build_promoters()
    print(f"Alu elements: {len(alu)};  protein-coding promoters: {len(prom)}")

    prom["alu_d"] = count_alu(prom, alu) / ((prom["end"] - prom["start"]) / 1000.0)
    genome_mean = prom["alu_d"].mean()

    is_nd = prom["gene"].isin(ndd).values
    is_hk = prom["gene"].isin(hk).values
    nd, hkd = prom["alu_d"].values[is_nd], prom["alu_d"].values[is_hk]

    # committed BEDs, in the same layout as the other species
    for mask, cat in [(is_nd, "HighConfNDD"), (is_hk, "Housekeeping")]:
        sub = prom[mask].copy()
        sub["score"], sub["strand"] = ".", "."
        sub["n"] = (sub["alu_d"] * ((sub["end"] - sub["start"]) / 1000.0)).round().astype(int)
        sub[["chrom", "start", "end", "gene", "score", "strand", "n"]].to_csv(
            OUT / f"{cat}_Alu.bed", sep="\t", header=False, index=False)

    _, p_hk = stats.mannwhitneyu(hkd, nd, alternative="greater")
    rest = prom["alu_d"].values[~is_nd]
    _, p_gen = stats.mannwhitneyu(rest, nd, alternative="greater")
    ci_nd = ratio_ci(nd, prom["alu_d"].values, 1)
    ci_hk = ratio_ci(hkd, prom["alu_d"].values, 2)

    row = dict(species="SquirrelMonkey", assembly="SaiBol1.0", mya=40, TE="Alu",
               n_all_PC=len(prom), n_HK=len(hkd), n_NDD=len(nd),
               median_HK=round(float(np.median(hkd)), 3),
               median_NDD=round(float(np.median(nd)), 3),
               genome_mean_alu=round(genome_mean, 3),
               r_vs_HK=round(rb(hkd, nd), 3), p_vs_HK=p_hk,
               NDD_over_genome=round(float(nd.mean()) / genome_mean, 3),
               NDD_CI=f"[{ci_nd[0]:.3f}, {ci_nd[1]:.3f}]",
               HK_over_genome=round(float(hkd.mean()) / genome_mean, 3),
               HK_CI=f"[{ci_hk[0]:.3f}, {ci_hk[1]:.3f}]",
               p_NDD_below_genome=p_gen)
    pd.DataFrame([row]).to_csv(OUT / "squirrel_stats.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n=== Bolivian squirrel monkey (Saimiri boliviensis), Alu ===")
    for k, v in row.items():
        print(f"  {k:20s} {v}")


if __name__ == "__main__":
    main()
