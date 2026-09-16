#!/usr/bin/env python3
"""
50_genome_baseline_by_class.py — the genome baseline of 25_genome_baseline.py,
run for LINE-1 and for dog Can-SINE as well as the GC-biased SINEs.

Why this is needed. The manuscript argues that the deficit distinguishes between
element families, and offers that as a reason to doubt indiscriminate removal of
sequence. That argument is made on the housekeeping contrast and on the rank
scale. But housekeeping promoters are Alu-rich and are not LINE-1-rich, so the
two element classes are not being measured against comparable references. The
question this script answers is what the classes look like against the same
neutral reference: all protein-coding promoters of the same species.

Promoter construction, chromosome filtering and the bootstrap follow
25_genome_baseline.py exactly, so the LINE-1 rows are comparable with the Alu
rows already reported. Two species need their own convention to stay comparable:
the squirrel monkey rows follow 30_squirrel_monkey.py (gene_name required, no
gene_id fallback) and the mouse lemur rows follow 20_lemur.py (RefSeq sequence
names mapped through chromAlias), because those are the scripts that produced
the published rows for those species.

Output: results/matched/genome_baseline_by_class.csv
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
GL = DATA / "gene_lists"
OUT = REPO / "results" / "matched"
OUT.mkdir(parents=True, exist_ok=True)

NAME_RE = re.compile(r'gene_name "([^"]+)"')
ID_RE = re.compile(r'gene_id "([^"]+)"')
HUMAN_MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}
SCAFFOLD_PREFIXES = ("chrMT", "chrKB", "chrKZ", "chrML", "chrJH", "chrGL")
N_BOOT = 4000

# label, GTF (None = autodiscover), element -> BED, chrom filter, chr prefix, upper
SPECIES = {
    "hg38":     dict(label="Human", gtf=DATA / "hg38" / "gtf" / "gencode.v47.gtf.gz",
                     chroms=HUMAN_MAIN, add_chr=False, upper=False,
                     beds={"Alu": DATA / "hg38" / "alu_rmsk.bed",
                           "LINE-1": DATA / "hg38" / "rmsk" / "LINE1.bed"}),
    "ponAbe3":  dict(label="Orangutan", gtf=None, chroms=None, add_chr=True, upper=False,
                     beds={"Alu": DATA / "ponAbe3" / "rmsk" / "Alu.bed",
                           "LINE-1": DATA / "ponAbe3" / "rmsk" / "LINE1.bed"}),
    "nomLeu3":  dict(label="Gibbon", gtf=None, chroms=None, add_chr=True, upper=False,
                     beds={"Alu": DATA / "nomLeu3" / "rmsk" / "Alu.bed",
                           "LINE-1": DATA / "nomLeu3" / "rmsk" / "LINE1.bed"}),
    "rheMac10": dict(label="Macaque", gtf=None, chroms=None, add_chr=True, upper=False,
                     beds={"Alu": DATA / "rheMac10" / "rmsk" / "Alu.bed",
                           "LINE-1": DATA / "rheMac10" / "rmsk" / "LINE1.bed"}),
    "calJac4":  dict(label="Marmoset", gtf=None, chroms=None, add_chr=True, upper=False,
                     beds={"Alu": DATA / "calJac4" / "rmsk" / "Alu.bed",
                           "LINE-1": DATA / "calJac4" / "rmsk" / "LINE1.bed"}),
    # saiBol1 scaffolds carry an INSDC version suffix in the GTF (JH378105.1) but
    # not in the RepeatMasker BED (JH378105); 30_squirrel_monkey.py strips it the
    # same way, requires gene_name (no gene_id fallback, which would add ~4000
    # unnamed predictions and shift the baseline), and adds no chr prefix.
    "saiBol1":  dict(label="Squirrel monkey", gtf=None, chroms=None, add_chr=False,
                     upper=False, strip_version=True, name_only=True,
                     beds={"Alu": DATA / "saiBol1" / "rmsk" / "Alu.bed",
                           "LINE-1": DATA / "saiBol1" / "rmsk" / "LINE1.bed"}),
    # mmur3 has a flat data directory and RepeatMasker coordinates on RefSeq
    # sequence names, so the BEDs are translated through chromAlias as in
    # 20_lemur.py. Both BEDs are written by 20_lemur.py / 37_lemur_line1.py.
    "mmur3":    dict(label="Mouse lemur", gtf=DATA / "mmur3" / "mmur3.gtf.gz",
                     chroms=None, add_chr=False, upper=False,
                     alias=DATA / "mmur3" / "chromAlias.txt",
                     beds={"Alu": DATA / "mmur3" / "alu_refseq.bed",
                           "LINE-1": DATA / "mmur3" / "line1_refseq.bed"}),
    "mm10":     dict(label="Mouse", gtf=None, chroms=None, add_chr=False, upper=True,
                     beds={"B1/B2": DATA / "mm10" / "rmsk" / "B1B2.bed",
                           "LINE-1": DATA / "mm10" / "rmsk" / "LINE1.bed"}),
    "canFam6":  dict(label="Dog", gtf=None, chroms=None, add_chr=True, upper=False,
                     beds={"Can-SINE": DATA / "canFam6" / "rmsk" / "CanSINE.bed",
                           "LINE-1": DATA / "canFam6" / "rmsk" / "LINE1.bed"}),
}


def find_gtf(key, configured):
    if configured is not None:
        return configured if configured.exists() else None
    d = DATA / key / "gtf"
    hits = sorted(d.glob("*.gtf.gz")) if d.is_dir() else []
    return hits[0] if hits else None


def gene_set(fn):
    return {l.strip() for l in open(GL / fn) if l.strip()}


def build_all_pc_promoters(gtf, chrom_filter, add_chr, upper, strip_version=False,
                           name_only=False):
    rows, seen = [], set()
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or "protein_coding" not in f[8]:
                continue
            chrom = f[0].split(".")[0] if strip_version else f[0]
            if add_chr:
                chrom = f"chr{chrom}"
            if chrom_filter is not None and chrom not in chrom_filter:
                continue
            if add_chr and chrom.startswith(SCAFFOLD_PREFIXES):
                continue
            m = NAME_RE.search(f[8]) if name_only \
                else (NAME_RE.search(f[8]) or ID_RE.search(f[8]))
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


def load_te(bed, alias_path):
    """TE intervals as chrom/start/end, translated through chromAlias if needed.

    chromAlias columns are refseq, assembly, genbank, ncbi, ucsc; the mouse lemur
    GTF uses the 'ncbi' names, which is what 20_lemur.py maps onto."""
    te = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"])
    if alias_path is None:
        return te
    alias = {}
    for line in open(alias_path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[3]:
            alias[f[0]] = f[3]
    te["chrom"] = te["chrom"].map(alias)
    return te.dropna(subset=["chrom"])


def count_te(prom, te):
    counts = np.zeros(len(prom), dtype=int)
    for chrom, m in te.groupby("chrom"):
        idx = np.where(prom["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        order = np.argsort(m["start"].values)
        s, e = m["start"].values[order], m["end"].values[order]
        for i in idx:
            ws, we = prom.at[i, "start"], prom.at[i, "end"]
            lo = np.searchsorted(s, we, "left")
            if lo:
                counts[i] = int(np.count_nonzero(e[:lo] > ws))
    return counts


def ratio_ci(group, baseline, seed):
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
        gtf = find_gtf(key, cfg["gtf"])
        if gtf is None:
            print(f"[atla] {cfg['label']}: GTF yok")
            continue
        prom = build_all_pc_promoters(gtf, cfg["chroms"], cfg["add_chr"], cfg["upper"],
                                      cfg.get("strip_version", False),
                                      cfg.get("name_only", False))
        prom["width_kb"] = (prom["end"] - prom["start"]) / 1000.0
        is_ndd = prom["gene"].isin(ndd).values
        is_hk = prom["gene"].isin(hk).values

        for element, bed in cfg["beds"].items():
            if not bed.exists():
                print(f"[atla] {cfg['label']} {element}: {bed.name} yok")
                continue
            te = load_te(bed, cfg.get("alias"))
            d = count_te(prom, te) / prom["width_kb"].values
            gm = d.mean()
            nd_d, hk_d = d[is_ndd], d[is_hk]
            _, p_nd = stats.mannwhitneyu(d[~is_ndd], nd_d, alternative="greater")
            ci_nd = ratio_ci(nd_d, d, 1)
            ci_hk = ratio_ci(hk_d, d, 2)
            rows.append(dict(
                species=cfg["label"], element=element, n_all_PC=len(prom),
                n_NDD=int(is_ndd.sum()), n_HK=int(is_hk.sum()),
                genome_mean=round(gm, 4),
                NDD_mean=round(nd_d.mean(), 4), HK_mean=round(hk_d.mean(), 4),
                NDD_over_genome=round(nd_d.mean() / gm, 3),
                NDD_CI=f"[{ci_nd[0]:.3f}, {ci_nd[1]:.3f}]",
                HK_over_genome=round(hk_d.mean() / gm, 3),
                HK_CI=f"[{ci_hk[0]:.3f}, {ci_hk[1]:.3f}]",
                p_NDD_below_genome=p_nd))
            print(f"  {cfg['label']:<16}{element:<9} "
                  f"NDD/genome {rows[-1]['NDD_over_genome']:.3f} {rows[-1]['NDD_CI']}   "
                  f"HK/genome {rows[-1]['HK_over_genome']:.3f}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "genome_baseline_by_class.csv", index=False)
    print(f"\nyazildi: results/matched/genome_baseline_by_class.csv ({len(out)} satir)")


if __name__ == "__main__":
    main()
