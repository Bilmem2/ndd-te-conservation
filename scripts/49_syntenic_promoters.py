"""
49_syntenic_promoters.py — repeat the cross-species comparison on promoters
defined by alignment rather than by annotation.

The published analysis matches genes across species by 1:1 orthology and then
scores each species' own annotated TSS. That establishes depletion at the
annotated promoters of orthologous genes, but not that the same piece of
sequence is being scored in each species: annotation completeness differs
between assemblies, and a gene may use a different primary promoter in one
species than in another.

Here the human promoter windows themselves are carried into each other genome
with UCSC liftOver, and the SINE frequency of the resulting interval is scored
in that genome. Any promoter that survives the lift is the same sequence in both
species by alignment, so a deficit measured this way cannot come from
annotation differences.

Two things are watched rather than assumed. Lift success is reported separately
for the NDD and housekeeping sets, because if NDD promoters are the more
conserved of the two a fixed threshold would retain them preferentially and bias
the comparison; the rates are printed so that the reader can see whether they
diverge. And the whole analysis is repeated across a range of minMatch values,
since that threshold is the main free parameter and the conclusion should not
depend on it.

Lifted intervals vary in length, so frequency is expressed per kb of the lifted
interval. Intervals that shrink or expand beyond a factor of two relative to the
4 kb source window are dropped as unreliable and counted.

Mouse lemur is not included: its RepeatMasker annotation is distributed in
RepeatMasker .out format under RefSeq sequence names, which do not correspond to
the UCSC assembly the chain file targets.

Inputs : results/hg38/{HighConfNDD,Housekeeping}_Alu.bed
         data/chains/hg38To*.over.chain.gz, tools/liftOver
         data/<species>/rmsk/{Alu,B1B2,LINE1}.bed
Output : results/synteny/syntenic_promoters.csv
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "synteny"
OUT.mkdir(parents=True, exist_ok=True)
LIFTOVER = REPO / "tools" / "liftOver"

SRC_WIDTH = 4000
LEN_RANGE = (SRC_WIDTH / 2, SRC_WIDTH * 2)      # drop lifts outside this span
MIN_MATCH = [0.10, 0.30, 0.50]                  # the free parameter, swept

# species -> (chain suffix, rmsk folder, element bed, element label)
TARGETS = [
    ("Orangutan",       "PonAbe3",  "ponAbe3",  "Alu.bed",   "Alu"),
    ("Gibbon",          "NomLeu3",  "nomLeu3",  "Alu.bed",   "Alu"),
    ("Macaque",         "RheMac10", "rheMac10", "Alu.bed",   "Alu"),
    ("Marmoset",        "CalJac4",  "calJac4",  "Alu.bed",   "Alu"),
    ("Squirrel monkey", "SaiBol1",  "saiBol1",  "Alu.bed",   "Alu"),
    ("Mouse",           "Mm10",     "mm10",     "B1B2.bed",  "B1/B2"),
    ("Dog",             "CanFam6",  "canFam6",  "LINE1.bed", "LINE-1"),
]


def load_human(cat):
    df = pd.read_csv(REPO / f"results/hg38/{cat}_Alu.bed", sep="\t", header=None,
                     names=["chr", "start", "end", "gene", "score", "strand", "n"])
    return df[["chr", "start", "end", "gene"]]


def to_wsl(p):
    """C:\\a\\b -> /mnt/c/a/b. liftOver is a Linux binary, so it is invoked
    through WSL while the analysis itself runs under Windows Python, which is
    where pandas and scipy are installed."""
    p = Path(p).resolve()
    drive = p.drive.rstrip(":").lower()
    return "/mnt/" + drive + "/" + p.as_posix()[len(p.drive) + 1:]


def lift(bed_df, chain, min_match, workdir):
    """Run liftOver and return the mapped intervals as a frame."""
    src = workdir / "in.bed"
    dst = workdir / "out.bed"
    unm = workdir / "unmapped.bed"
    bed_df.to_csv(src, sep="\t", header=False, index=False, lineterminator="\n")
    inner = (f"{to_wsl(LIFTOVER)} -minMatch={min_match} "
             f"{to_wsl(src)} {to_wsl(chain)} {to_wsl(dst)} {to_wsl(unm)}")
    res = subprocess.run(["wsl", "-e", "bash", "-lc", inner],
                         capture_output=True, text=True)
    if not dst.exists():
        raise RuntimeError((res.stderr or res.stdout).strip()[:300])
    if dst.stat().st_size == 0:
        return pd.DataFrame(columns=["chr", "start", "end", "gene"])
    return pd.read_csv(dst, sep="\t", header=None,
                       names=["chr", "start", "end", "gene"])


def index_bed(path):
    rec = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chr", "start", "end"])
    return {c: (np.sort(s.start.to_numpy()), np.sort(s.end.to_numpy()))
            for c, s in rec.groupby("chr")}


def density(idx, df):
    out = np.zeros(len(df))
    for i, (c, s, e) in enumerate(zip(df.chr.to_numpy(), df.start.to_numpy(),
                                      df.end.to_numpy())):
        if c in idx:
            st, en = idx[c]
            out[i] = (np.searchsorted(st, e, side="left")
                      - np.searchsorted(en, s, side="right"))
    return out / ((df.end.to_numpy() - df.start.to_numpy()) / 1000.0)


def test(hk, ndd):
    u = stats.mannwhitneyu(hk, ndd, alternative="greater")
    r = 1 - (2 * u.statistic) / (len(hk) * len(ndd))
    return r, u.pvalue


ndd_h, hk_h = load_human("HighConfNDD"), load_human("Housekeeping")
print(f"human source windows: NDD {len(ndd_h)}, HK {len(hk_h)}"
      f"   ({SRC_WIDTH // 1000} kb each)\n")

rows = []
with tempfile.TemporaryDirectory() as td:
    work = Path(td)
    for label, chain_tag, folder, bed, elem in TARGETS:
        chain = REPO / "data" / "chains" / f"hg38To{chain_tag}.over.chain.gz"
        rmsk = REPO / "data" / folder / "rmsk" / bed
        if not chain.exists() or not rmsk.exists():
            print(f"{label:<17} SKIP (missing chain or rmsk)")
            continue
        idx = index_bed(rmsk)

        for mm in MIN_MATCH:
            lifted = {}
            dropped = {}
            for cat, src in (("NDD", ndd_h), ("HK", hk_h)):
                d = lift(src, chain, mm, work)
                span = d.end - d.start
                keep = span.between(*LEN_RANGE)
                dropped[cat] = int((~keep).sum())
                lifted[cat] = d[keep].reset_index(drop=True)

            n_nd, n_hk = len(lifted["NDD"]), len(lifted["HK"])
            if n_nd < 50 or n_hk < 50:
                print(f"{label:<17} minMatch {mm}: too few lifted "
                      f"({n_nd}/{n_hk}), skipped")
                continue

            d_nd = density(idx, lifted["NDD"])
            d_hk = density(idx, lifted["HK"])
            r, p = test(d_hk, d_nd)
            rate_nd = 100 * n_nd / len(ndd_h)
            rate_hk = 100 * n_hk / len(hk_h)

            rows.append(dict(
                species=label, element=elem, min_match=mm,
                n_NDD=n_nd, n_HK=n_hk,
                lift_rate_NDD=round(rate_nd, 1), lift_rate_HK=round(rate_hk, 1),
                rate_gap=round(rate_nd - rate_hk, 1),
                dropped_NDD=dropped["NDD"], dropped_HK=dropped["HK"],
                mean_NDD=round(float(d_nd.mean()), 4),
                mean_HK=round(float(d_hk.mean()), 4),
                ratio=round(float(d_nd.mean() / d_hk.mean()), 3),
                r=round(float(r), 4), p=float(p)))

res = pd.DataFrame(rows)

print(f"{'species':<17}{'elem':<8}{'mM':>5}{'nNDD':>6}{'nHK':>6}"
      f"{'lift% N':>9}{'lift% H':>9}{'gap':>7}{'NDD/HK':>9}{'r':>8}{'p':>11}")
print("-" * 95)
for _, x in res.iterrows():
    print(f"{x.species:<17}{x.element:<8}{x.min_match:>5.2f}{x.n_NDD:>6}{x.n_HK:>6}"
          f"{x.lift_rate_NDD:>9.1f}{x.lift_rate_HK:>9.1f}{x.rate_gap:>7.1f}"
          f"{x.ratio:>9.3f}{x.r:>8.3f}{x.p:>11.1e}")
print("-" * 95)

mid = res[res.min_match == 0.30]
print(f"\nAt minMatch 0.30, depletion is significant in "
      f"{int((mid.p < 0.05).sum())} of {len(mid)} species tested "
      f"(r {mid.r.min():.3f} to {mid.r.max():.3f}).")
gap = res.rate_gap.abs().max()
print(f"Largest NDD-minus-HK lift-rate gap across all runs: {gap:.1f} points"
      f"  ({'no strong retention bias' if gap < 10 else 'CHECK: retention bias'}).")

res.to_csv(OUT / "syntenic_promoters.csv", index=False)
print(f"\nwrote {(OUT / 'syntenic_promoters.csv').relative_to(REPO)}")
