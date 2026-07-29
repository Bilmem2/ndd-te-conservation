#!/usr/bin/env python3
"""
25_genome_baseline.py — Express Alu depletion against a within-species GENOME
baseline (all protein-coding promoters) instead of the housekeeping contrast.

Housekeeping promoters are themselves Alu-rich, so NDD-vs-housekeeping mixes
"housekeeping enrichment" with "NDD depletion". 19_rebaseline.py established
this for human; that script also uses GTEx expression breadth and gnomAD
polymorphic insertions, neither of which exists outside human. The genome
baseline, by contrast, needs only a GTF and an Alu BED, so it ports to any
species whose raw annotation is available locally.

For each species this reports, at TSS +/-2 kb:
  * mean fixed Alu density over ALL protein-coding promoters (the baseline)
  * NDD / genome and housekeeping / genome density ratios, with bootstrap CIs
  * one-sided Mann-Whitney U for NDD vs all other protein-coding promoters

Species are included only when their local annotation is present; missing ones
are skipped with a notice rather than failing, so the table grows as raw data
is re-downloaded.

Output: results/matched/genome_baseline_by_species.csv
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
GL = DATA / "gene_lists"
OUT = REPO / "results" / "matched"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
ID_RE = re.compile(r'gene_id "([^"]+)"')
HUMAN_MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}
# unplaced scaffolds excluded for the Ensembl assemblies, matching 03_get_promoters.sh
SCAFFOLD_PREFIXES = ("chrMT", "chrKB", "chrKZ", "chrML", "chrJH", "chrGL")
N_BOOT = 4000

# label, GTF, TE BED, chromAlias, chrom filter, add "chr" prefix, uppercase symbols
SPECIES = {
    "hg38": dict(label="Human", gtf=DATA / "hg38" / "gtf" / "gencode.v47.gtf.gz",
                 bed=DATA / "hg38" / "alu_rmsk.bed", alias=None,
                 chroms=HUMAN_MAIN, add_chr=False, upper=False, te="Alu"),
    "ponAbe3": dict(label="Orangutan", gtf=None, bed=DATA / "ponAbe3" / "rmsk" / "Alu.bed",
                    alias=None, chroms=None, add_chr=True, upper=False, te="Alu"),
    "nomLeu3": dict(label="Gibbon", gtf=None, bed=DATA / "nomLeu3" / "rmsk" / "Alu.bed",
                    alias=None, chroms=None, add_chr=True, upper=False, te="Alu"),
    "rheMac10": dict(label="Macaque", gtf=None, bed=DATA / "rheMac10" / "rmsk" / "Alu.bed",
                     alias=None, chroms=None, add_chr=True, upper=False, te="Alu"),
    "calJac4": dict(label="Marmoset", gtf=None, bed=DATA / "calJac4" / "rmsk" / "Alu.bed",
                    alias=None, chroms=None, add_chr=True, upper=False, te="Alu"),
    "mmur3": dict(label="Mouse lemur", gtf=DATA / "mmur3" / "mmur3.gtf.gz",
                  bed=DATA / "mmur3" / "alu_refseq.bed",
                  alias=DATA / "mmur3" / "chromAlias.txt",
                  chroms=None, add_chr=False, upper=False, te="Alu"),
    "mm10": dict(label="Mouse", gtf=None, bed=DATA / "mm10" / "rmsk" / "B1B2.bed",
                 alias=None, chroms=None, add_chr=False, upper=True, te="B1/B2"),
}


def find_gtf(key, configured):
    """Use the configured path when given, else the single *.gtf.gz in data/<key>/gtf/."""
    if configured is not None:
        return configured if configured.exists() else None
    hits = sorted((DATA / key / "gtf").glob("*.gtf.gz")) if (DATA / key / "gtf").is_dir() else []
    return hits[0] if hits else None


def gene_set(fn):
    return {l.strip() for l in open(GL / fn) if l.strip()}


def load_alu(bed, alias_path):
    if alias_path is None:
        a = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2],
                        names=["chrom", "start", "end"])
        return a
    alias = {}
    for line in open(alias_path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[3]:
            alias[f[0]] = f[3]
    a = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["rs", "start", "end"])
    a["chrom"] = a["rs"].map(alias)
    return a.dropna(subset=["chrom"])[["chrom", "start", "end"]]


def build_all_pc_promoters(gtf, chrom_filter, add_chr, upper):
    """TSS +/-2 kb for every protein-coding gene (one window per gene symbol).

    Mirrors 03_get_promoters.sh: Ensembl contigs gain a 'chr' prefix and
    unplaced scaffolds are dropped; mouse symbols are upper-cased so they match
    the human-derived gene lists (as 04_split_promoters.sh does)."""
    rows, seen = [], set()
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or "protein_coding" not in f[8]:
                continue
            chrom = f"chr{f[0]}" if add_chr else f[0]
            if chrom_filter is not None and chrom not in chrom_filter:
                continue
            if add_chr and chrom.startswith(SCAFFOLD_PREFIXES):
                continue
            m = NAME_RE.search(f[8]) or ID_RE.search(f[8])
            if not m:
                continue
            gene = m.group(1).upper() if upper else m.group(1)
            if gene in seen:
                continue
            seen.add(gene)
            start, end, strand = int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end
            rows.append((chrom, max(0, tss - 2000), tss + 2000, gene))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def count_alu(prom, alu):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in alu.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        order = np.argsort(m["start"].values)
        s = m["start"].values[order]
        e = m["end"].values[order]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def ratio_ci(group, baseline, seed):
    """Bootstrap CI for mean(group)/mean(baseline); baseline resampled too."""
    rng = np.random.default_rng(seed)
    b = np.empty(N_BOOT)
    for k in range(N_BOOT):
        g = rng.choice(group, len(group), True).mean()
        r = rng.choice(baseline, len(baseline), True).mean()
        b[k] = g / r if r else np.nan
    return np.nanpercentile(b, [2.5, 97.5])


def main():
    ndd, hk = gene_set("HighConfNDD_genes.txt"), gene_set("Housekeeping_genes.txt")
    rows = []

    for key, cfg in SPECIES.items():
        label, bed = cfg["label"], cfg["bed"]
        gtf = find_gtf(key, cfg["gtf"])
        if gtf is None or not bed.exists():
            print(f"[skip] {label}: local annotation not present "
                  f"({'GTF' if gtf is None else 'TE BED'} missing)")
            continue

        prom = build_all_pc_promoters(gtf, cfg["chroms"], cfg["add_chr"], cfg["upper"])
        alu = load_alu(bed, cfg["alias"])
        prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
        prom["alu_d"] = count_alu(prom, alu) / prom["width_kb"].values

        is_ndd = prom["gene"].isin(ndd).values
        is_hk = prom["gene"].isin(hk).values
        genome_mean = prom["alu_d"].mean()

        nd_d = prom["alu_d"].values[is_ndd]
        hk_d = prom["alu_d"].values[is_hk]
        rest = prom["alu_d"].values[~is_ndd]
        _, p_nd = stats.mannwhitneyu(rest, nd_d, alternative="greater")

        ci_nd = ratio_ci(nd_d, prom["alu_d"].values, 1)
        ci_hk = ratio_ci(hk_d, prom["alu_d"].values, 2)

        rows.append(dict(
            species=label, n_all_PC=len(prom), n_NDD=int(is_ndd.sum()),
            n_HK=int(is_hk.sum()),
            genome_mean_alu=round(genome_mean, 3),
            NDD_mean_alu=round(nd_d.mean(), 3), HK_mean_alu=round(hk_d.mean(), 3),
            NDD_over_genome=round(nd_d.mean() / genome_mean, 3),
            NDD_CI=f"[{ci_nd[0]:.3f}, {ci_nd[1]:.3f}]",
            HK_over_genome=round(hk_d.mean() / genome_mean, 3),
            HK_CI=f"[{ci_hk[0]:.3f}, {ci_hk[1]:.3f}]",
            p_NDD_below_genome=p_nd))
        print(f"{label:12s} genome={genome_mean:.3f}/kb  "
              f"NDD/genome={rows[-1]['NDD_over_genome']:.3f} {rows[-1]['NDD_CI']}  "
              f"HK/genome={rows[-1]['HK_over_genome']:.3f} {rows[-1]['HK_CI']}  "
              f"p={p_nd:.2e}")

    if not rows:
        print("No species processed.")
        return
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "genome_baseline_by_species.csv", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print("\n=== Alu density relative to the within-species genome baseline ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'genome_baseline_by_species.csv'}")


if __name__ == "__main__":
    main()
