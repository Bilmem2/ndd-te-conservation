"""
48_cross_disease_strict.py — rebuild the cross-disease gene sets from an
unambiguous pathogenic filter, and measure what changes.

10_cross_disease.py selects ClinVar records with

    ClinicalSignificance.str.contains('Pathogenic', case=False)

which is a substring test. ClinVar's vocabulary includes "Conflicting
classifications of pathogenicity", and the word pathogenicity contains the word
pathogenic, so that value is selected too. The manuscript describes the sets as
pathogenic/likely pathogenic, and records on which submitters disagree are not
that. The rest of the paper takes the high-confidence branch of every resource it
uses, SFARI Tier 1 and 2 and ClinGen Definitive and Strong, so admitting
contested calls here is also out of keeping with its own design.

The rule applied instead: strip any secondary annotation after a semicolon, and
keep the record only if what remains is an unambiguous assertion of
pathogenicity, that is Pathogenic, Likely pathogenic or Pathogenic/Likely
pathogenic, optionally qualified as low penetrance. Anything containing
Conflicting or risk allele is excluded. The kept and dropped values are printed
so the rule can be checked rather than taken on trust, and so the list can be
quoted in the methods.

The ClinVar release is pinned to the dated monthly archive matching the vintage
of the other gene-set resources, so this run differs from the published one in
the filter alone and not in the underlying data.

Inputs : data/clinvar_variants.txt.gz, data/gene_lists/{HighConfNDD,Housekeeping}_genes.txt
         data/hg38/gtf/gencode.v47.gtf.gz, data/hg38/alu_rmsk.bed
Output : results/cross_disease/cross_disease_strict.csv
         data/gene_lists/{Cardiovascular,Mendelian}_genes_strict.txt
"""
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
REPO = Path(__file__).resolve().parents[1]
DATA, GL = REPO / "data", REPO / "data" / "gene_lists"
OUT = REPO / "results" / "cross_disease"
OUT.mkdir(parents=True, exist_ok=True)

# release is a command-line argument so the same code can be run against the
# dated archive and against the current download, and the two compared
CLINVAR = DATA / (sys.argv[1] if len(sys.argv) > 1 else "clinvar_variants.txt.gz")
CARDIO_TERMS = ["cardiomyopathy", "arrhythmia", "channelopathy",
                "long QT", "Brugada", "heart failure"]
MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}
NAME_RE = re.compile(r'gene_name "([^"]+)"')
WIN = 2000

KEEP = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
LOW_PEN = re.compile(r",\s*low penetrance\s*$", re.I)


def is_pathogenic(value):
    """True for an unambiguous P/LP assertion, ignoring secondary annotations."""
    v = str(value).split(";")[0].strip()
    if not v or "conflicting" in v.lower() or "risk allele" in v.lower():
        return False
    v = LOW_PEN.sub("", v).strip()
    # a combined call such as Pathogenic/Likely pathogenic/Pathogenic
    parts = {p.strip().lower() for p in v.split("/")}
    return bool(parts) and parts <= {"pathogenic", "likely pathogenic"}


# identical to SYMBOL_RE in 10_cross_disease.py and 26_cross_disease_recompute.py,
# so that the loose branch below reproduces the published sets exactly and the
# only difference between the two branches is the pathogenicity filter
SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.@-]*$")


def symbols(series):
    out = set()
    for cell in series.dropna():
        for tok in str(cell).replace("\r", "").split(";"):
            tok = tok.strip()
            if tok and tok != "-" and SYMBOL_RE.match(tok):
                out.add(tok)
    return out


def promoters():
    rows, seen = [], set()
    with gzip.open(DATA / "hg38" / "gtf" / "gencode.v47.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if (f[2] != "gene" or 'gene_type "protein_coding"' not in f[8]
                    or f[0] not in MAIN):
                continue
            m = NAME_RE.search(f[8])
            if not m or m.group(1) in seen:
                continue
            seen.add(m.group(1))
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            rows.append((f[0], max(0, tss - WIN), tss + WIN, m.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


print(f"reading {CLINVAR.name}")
cv = pd.read_csv(CLINVAR, sep="\t", low_memory=False,
                 usecols=["GeneSymbol", "ClinicalSignificance",
                          "PhenotypeList", "Assembly"])
cv = cv[cv["Assembly"] == "GRCh38"]
sig = cv["ClinicalSignificance"].fillna("")

loose = sig.str.contains("Pathogenic", case=False, na=False)
strict = sig.map(is_pathogenic)

print(f"\nGRCh38 records: {len(cv):,}")
print(f"  loose  (published filter) : {int(loose.sum()):,}")
print(f"  strict (this script)      : {int(strict.sum()):,}")
print(f"  removed                   : {int(loose.sum() - strict.sum()):,}"
      f"  ({100 * (loose.sum() - strict.sum()) / loose.sum():.1f} %)")

dropped = Counter(sig[loose & ~strict])
print("\nvalues the published filter accepted and this one drops:")
for v, n in dropped.most_common(8):
    print(f"  {n:>9,}  {v[:66]}")
kept_vals = Counter(sig[strict])
print("\nvalues kept (top 8):")
for v, n in kept_vals.most_common(8):
    print(f"  {n:>9,}  {v[:66]}")

ndd = set((GL / "HighConfNDD_genes.txt").read_text().split())
hk = set((GL / "Housekeeping_genes.txt").read_text().split())

sets = {}
for label, mask in (("loose", loose), ("strict", strict)):
    sub = cv[mask]
    cardio_mask = sub["PhenotypeList"].fillna("").str.contains(
        "|".join(CARDIO_TERMS), case=False, na=False)
    sets[label] = {
        "Cardiovascular": symbols(sub[cardio_mask]["GeneSymbol"]) - ndd - hk,
        "Mendelian": symbols(sub["GeneSymbol"]) - ndd - hk,
    }

print("\ngene sets")
for name in ("Cardiovascular", "Mendelian"):
    lo, st = len(sets["loose"][name]), len(sets["strict"][name])
    print(f"  {name:<16} loose {lo:>6,}   strict {st:>6,}"
          f"   removed {lo - st:>5,} ({100 * (lo - st) / lo:.1f} %)")

# ── promoter Alu frequency for each set, against the genome baseline ────────
prom = promoters()
alu = pd.read_csv(DATA / "hg38" / "alu_rmsk.bed", sep="\t", header=None,
                  usecols=[0, 1, 2], names=["chrom", "start", "end"])
idx = {c: (np.sort(s.start.to_numpy()), np.sort(s.end.to_numpy()))
       for c, s in alu.groupby("chrom")}
cnt = np.zeros(len(prom))
for i, (c, s, e) in enumerate(zip(prom.chrom, prom.start, prom.end)):
    if c in idx:
        st, en = idx[c]
        cnt[i] = (np.searchsorted(st, e, side="left")
                  - np.searchsorted(en, s, side="right"))
prom["alu_d"] = cnt / ((prom.end - prom.start) / 1000.0)
genome_mean = prom.alu_d.mean()
print(f"\ngenome baseline {genome_mean:.3f} Alu per kb"
      f" over {len(prom):,} protein-coding promoters")

RNG = np.random.default_rng(20260915)


def boot_ci(v, n=4000):
    m = RNG.choice(v, size=(n, v.size), replace=True).mean(axis=1) / genome_mean
    return np.percentile(m, [2.5, 97.5])


rows = []
hk_d = prom.alu_d.values[prom.gene.isin(hk).values]
for label in ("loose", "strict"):
    for name in ("HighConfNDD", "Cardiovascular", "Mendelian"):
        genes = ndd if name == "HighConfNDD" else sets[label][name]
        d = prom.alu_d.values[prom.gene.isin(genes).values]
        if d.size == 0:
            continue
        lo, hi = boot_ci(d)
        u = stats.mannwhitneyu(hk_d, d, alternative="greater")
        r = 1 - (2 * u.statistic) / (len(hk_d) * len(d))
        rows.append(dict(filter=label, set=name, n_genes=len(genes),
                         n_promoters=int(d.size),
                         mean_alu=round(float(d.mean()), 4),
                         over_genome=round(float(d.mean() / genome_mean), 3),
                         ci_low=round(float(lo), 3), ci_high=round(float(hi), 3),
                         r_vs_HK=round(float(r), 4), p_vs_HK=float(u.pvalue)))

res = pd.DataFrame(rows)
print(f"\n{'filter':<8}{'set':<16}{'genes':>7}{'prom':>7}{'Alu/kb':>9}"
      f"{'/genome':>9}{'95% CI':>18}")
print("-" * 76)
for _, x in res.iterrows():
    print(f"{x['filter']:<8}{x['set']:<16}{x.n_genes:>7,}{x.n_promoters:>7,}"
          f"{x.mean_alu:>9.3f}{x.over_genome:>9.3f}"
          f"   [{x.ci_low:.3f}, {x.ci_high:.3f}]")
print("-" * 76)

res.to_csv(OUT / "cross_disease_strict.csv", index=False)
for name in ("Cardiovascular", "Mendelian"):
    (GL / f"{name}_genes_strict.txt").write_text(
        "\n".join(sorted(sets["strict"][name])) + "\n", encoding="utf-8")
print(f"\nwrote {(OUT / 'cross_disease_strict.csv').relative_to(REPO)}")
print("wrote data/gene_lists/{Cardiovascular,Mendelian}_genes_strict.txt")
