#!/usr/bin/env python3
"""
35_alu_age_by_species.py — is the marmoset Alu complement unusually young?

The marmoset is the one species whose NDD promoters are not depleted against
their own genome baseline, and it also carries the highest genome-wide promoter
Alu density in the panel. One explanation offered for both is a relatively
recent or unusually active Alu expansion in that lineage: elements accumulating
faster than they are lost would raise density and compress differences between
gene classes. That is testable rather than merely assertable, because
RepeatMasker reports the divergence of each copy from its subfamily consensus
(milliDiv), which is a proxy for age: low values mean recently inserted.

This compares the age distribution of Alu elements across the species for which
a RepeatMasker table is held locally. If the marmoset explanation holds, its Alu
should be systematically younger than those of the other primates, and in
particular younger than those of the other Platyrrhine, the squirrel monkey.

Output: results/mechanism/alu_age_by_species.csv
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "results" / "mechanism"
OUT.mkdir(parents=True, exist_ok=True)

SPECIES = [
    ("hg38", "Human"),
    ("ponAbe3", "Orangutan"),
    ("nomLeu3", "Gibbon"),
    ("rheMac10", "Macaque"),
    ("calJac4", "Marmoset"),
    ("saiBol1", "Squirrel monkey"),
]
YOUNG = 150  # milliDiv below this counts as a recent insertion


def main():
    rows = []
    for key, label in SPECIES:
        path = DATA / key / "rmsk" / "rmsk.txt.gz"
        if not path.exists():
            print(f"[skip] {label}: no local RepeatMasker table")
            continue
        div = []
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 12 and f[11] == "SINE" and f[12] == "Alu":
                    try:
                        div.append(int(f[2]))
                    except ValueError:
                        continue
        d = np.array(div, dtype=float)
        rows.append(dict(species=label, n_alu=len(d),
                         median_milliDiv=float(np.median(d)),
                         mean_milliDiv=round(float(d.mean()), 1),
                         q25=float(np.percentile(d, 25)),
                         pct_young=round(100 * float((d < YOUNG).mean()), 1)))
        print(f"{label:16s} n={len(d):8d}  median milliDiv={rows[-1]['median_milliDiv']:.0f}  "
              f"<{YOUNG}: {rows[-1]['pct_young']:.1f}%")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "alu_age_by_species.csv", index=False)
    pd.set_option("display.width", 200)
    print("\n=== Alu age distribution by species (lower milliDiv = younger) ===")
    print(out.to_string(index=False))
    print(f"\nWrote {OUT / 'alu_age_by_species.csv'}")


if __name__ == "__main__":
    main()
