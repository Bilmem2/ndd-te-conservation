#!/usr/bin/env python3
"""
21_lemur_ortholog.py — Ortholog-validated replication of the mouse-lemur Alu
depletion, matching the manuscript's Ensembl BioMart one2one validation used
for the other non-human species.

Fetches human<->Microcebus murinus orthologs from Ensembl BioMart, restricts
the NDD/HK gene sets to genes with a 1:1 ortholog, and recomputes the Alu
depletion effect size. Concordance with the symbol-based result (r=-0.262)
indicates the finding is robust to ortholog definition.
"""
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
MM = REPO / "results" / "mmur3"
OD = REPO / "data" / "orthologs"
OD.mkdir(parents=True, exist_ok=True)

BIOMART = "https://www.ensembl.org/biomart/martservice"
XML = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
<Dataset name="hsapiens_gene_ensembl" interface="default">
<Filter name="with_mmurinus_homolog" excluded="0"/>
<Attribute name="external_gene_name"/>
<Attribute name="mmurinus_homolog_associated_gene_name"/>
<Attribute name="mmurinus_homolog_orthology_type"/>
</Dataset>
</Query>'''


def fetch_orthologs():
    cache = OD / "mmur3_raw.tsv"
    if cache.exists() and cache.stat().st_size > 1000:
        return pd.read_csv(cache, sep="\t", header=None,
                           names=["human", "lemur", "type"], dtype=str).dropna(subset=["human"])
    data = urllib.parse.urlencode({"query": XML}).encode()
    req = urllib.request.Request(BIOMART, data=data)
    with urllib.request.urlopen(req, timeout=180) as r:
        txt = r.read().decode()
    with open(cache, "w", encoding="utf-8") as f:
        f.write(txt)
    rows = [l.split("\t") for l in txt.strip().splitlines() if l and not l.startswith("Query")]
    df = pd.DataFrame([r for r in rows if len(r) == 3], columns=["human", "lemur", "type"])
    return df


def load_lemur():
    frames = []
    for cat in ("HighConfNDD", "Housekeeping"):
        d = pd.read_csv(MM / f"{cat}_Alu.bed", sep="\t", header=None,
                        names=["chrom", "start", "end", "gene", "category", "alu"])
        d["density"] = d["alu"] / ((d["end"] - d["start"]) / 1000)
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def rstat(hk, nd):
    U, p = stats.mannwhitneyu(hk, nd, alternative="greater")
    u2, _ = stats.mannwhitneyu(hk, nd, alternative="two-sided")
    return p, round(1 - (2 * u2) / (len(hk) * len(nd)), 3)


def main():
    orth = fetch_orthologs()
    print(f"BioMart human-lemur homolog rows: {len(orth)}")
    one2one = set(orth[orth["type"] == "ortholog_one2one"]["human"])
    print(f"  one2one orthologs: {len(one2one)}")

    lem = load_lemur()
    val = lem[lem["gene"].isin(one2one)]
    out = []
    for label, df in [("symbol-based (all mapped)", lem), ("ortholog one2one", val)]:
        hk = df[df["category"] == "Housekeeping"]["density"]
        nd = df[df["category"] == "HighConfNDD"]["density"]
        p, r = rstat(hk, nd)
        out.append(dict(method=label, n_HK=len(hk), n_NDD=len(nd), p_value=p, r=r))
    res = pd.DataFrame(out)
    res.to_csv(MM / "lemur_ortholog_validation.csv", index=False)
    pd.set_option("display.width", 150)
    print("\n=== Mouse lemur Alu depletion: symbol vs ortholog-validated ===")
    print(res.to_string(index=False))
    dr = abs(out[0]["r"] - out[1]["r"])
    print(f"\n|delta r| = {dr:.3f}  ({'concordant' if dr < 0.03 else 'check'})")


if __name__ == "__main__":
    main()
