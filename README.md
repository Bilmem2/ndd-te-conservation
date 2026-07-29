# Conserved SINE and Lineage-Variable LINE-1 Depletion at Neurodevelopmental Disorder Promoters

**Can Sevilmiş** | Department of Molecular Biology and Genetics, Bahçeşehir University, Istanbul, Turkey

**Repository:** https://github.com/Bilmem2/ndd-te-conservation  
**Related work:** Sevilmiş, C. (2026). LINE-1 Depletion at Promoters of Neurodevelopmental Disorder Genes: A Genome-Wide Analysis. *Preprints*. https://doi.org/10.20944/preprints202604.0715.v1

---

## Overview

This repository contains all analysis scripts, processed gene lists, statistical results, and figures for the above manuscript. The analysis spans **nine mammalian genomes** — seven primates (human, orangutan, gibbon, macaque, two New World monkeys — marmoset and **squirrel monkey, _Saimiri boliviensis_** — and the strepsirrhine **gray mouse lemur, _Microcebus murinus_**) plus mouse and dog — and adds **genomic-context confounder controls** (local gene density, recombination rate, joint matched control), an **ENCODE cCRE functional overlay**, and a **fetal-brain regulatory-specificity test** (ENCODE fetal DNase-seq, three brain donors vs three non-neural fetal tissues). Raw genome assemblies, RepeatMasker annotations, GTF files, and the hg38 reference are not included due to size constraints; download instructions are provided below.

---

## Repository Structure

```
.
├── scripts/
│   │   # Numbered 01-31 in pipeline order. Scripts tagged (exploratory) are retained
│   │   # for provenance and are NOT used as evidence in the manuscript.
│   ├── 01_prepare_gene_lists.py       # Curate NDD and Housekeeping gene sets
│   ├── 02_rmsk_to_bed.sh              # Extract Alu/LINE-1/B1B2 BED from RepeatMasker
│   ├── 03_get_promoters.sh            # Extract TSS ± 2 kb promoter windows from GTF
│   ├── 04_split_promoters.sh          # Split promoter BEDs by gene category
│   ├── 05_intersect.sh                # BEDTools intersect: TE count per promoter
│   ├── 06_statistics.py               # (exploratory) preliminary statistics
│   ├── 07_pli_correlation.py          # (exploratory) pLI vs Alu density
│   ├── 08_encode_overlap_v2.py        # (exploratory) CTCF/DNase overlap
│   ├── 09_window_sensitivity.py       # Window size sensitivity (± 0.5–3 kb)
│   ├── 10_cross_disease.py            # Cross-disease specificity (ClinVar)
│   ├── 11_stats_updated.py            # Cross-species statistics (7 species)
│   ├── 12_figures_final.py            # Core figures (Fig 1–6)
│   ├── 13_ortholog_analysis.py        # Ortholog-validated replication (Ensembl BioMart)
│   ├── 14_gnomad_mei.py               # (exploratory) gnomAD polymorphic MEI decomposition (inconclusive)
│   ├── 15_context_controls.py         # Gene density + recombination controls (Fig 7)
│   ├── 16_ccre_overlay.py             # ENCODE cCRE functional overlay (Fig 8)
│   ├── 17_functional_consequence.py   # (exploratory) gnomAD constraint / GTEx expression vs Alu (inconclusive)
│   ├── 18_matched_control.py          # Joint GC + density + recombination matched control
│   ├── 19_rebaseline.py               # Honest genome / expression-matched re-baselining
│   ├── 20_lemur.py                    # Mouse lemur (strepsirrhine) Alu depletion
│   ├── 21_lemur_ortholog.py           # Mouse lemur ortholog validation (BioMart)
│   ├── 22_consolidate.py              # Master cross-species table + results backbone
│   ├── 23_brain_overlay.py            # Fetal-brain DNase regulatory-specificity overlay (Fig 9)
│   ├── 24_null_model.py               # Permutation null model, all 5 species-TE combos
│   ├── 25_genome_baseline.py          # Within-species genome baseline, all species (Table 2)
│   ├── 26_cross_disease_recompute.py  # Repaired ClinVar gene sets + genome-referenced comparison
│   ├── 27_constraint_matched.py       # LOEUF + brain-expression + GC matched control
│   ├── 28_mane_tss.py                 # Canonical (MANE Select) TSS sensitivity analysis
│   ├── 29_b1_b2_split.py              # Mouse B1 vs B2 separately + BH q-values for Table 1
│   ├── 30_squirrel_monkey.py          # Second Platyrrhine (Saimiri boliviensis)
│   ├── 31_pseudoreplication.py        # Thinning clustered/paralogous promoters
│   ├── probe_dosage.py                # (exploratory) dosage-sensitivity continuum probe (did not hold)
│   ├── probe_subfamily.py             # (exploratory) Alu subfamily-age (AluJ/S/Y) probe (inconclusive)
│   ├── fig_lemur_update.py            # Regenerate Fig 1 (7 primates) & Fig 3 (9 species)
│   ├── fig_new.py                     # Generate Fig 7 (context) & Fig 8 (cCRE)
│   ├── fig_brain.py                   # Generate Fig 9 (fetal-brain specificity)
│   ├── fig_null_update.py             # Regenerate Fig 4 (5-panel permutation null model)
│   └── fig_phylo.py                   # Generate Fig S5 (effect size vs divergence time)
│
├── data/
│   ├── gene_lists/
│   │   ├── HighConfNDD_genes.txt      # SFARI Tier 1+2 ∪ ClinGen Epilepsy Definitive/Strong (n=1020)
│   │   ├── Housekeeping_genes.txt     # HRT Atlas ∩ brain TPM≥1, NDD-free (n=1679)
│   │   ├── Cardiovascular_genes.txt   # ClinVar P/LP cardiovascular genes
│   │   └── Mendelian_genes.txt        # ClinVar P/LP broad Mendelian genes
│   │
│   └── orthologs/
│       ├── ensembl_to_symbol.tsv      # Ensembl ID ↔ HGNC symbol mapping
│       └── symbol_to_ensembl.csv      # HGNC symbol ↔ Ensembl ID mapping
│
├── results/
│   ├── statistics_final.csv        # Main results: all species, symbol-based matching
│   ├── statistics_ortholog.csv     # Validation: Ensembl BioMart 1:1 ortholog-based matching
│   ├── null_model_full.csv         # Permutation test, all 5 species-TE combos (n=10,000)
│   ├── encode_overlap_v2.csv       # ENCODE overlap results (exploratory)
│   ├── hg38/                       # Human
│   │   ├── HighConfNDD_Alu.bed
│   │   ├── Housekeeping_Alu.bed
│   │   ├── HighConfNDD_LINE1.bed
│   │   └── Housekeeping_LINE1.bed
│   ├── ponAbe3/                    # Orangutan
│   │   ├── HighConfNDD_Alu.bed
│   │   ├── Housekeeping_Alu.bed
│   │   ├── HighConfNDD_LINE1.bed
│   │   └── Housekeeping_LINE1.bed
│   ├── nomLeu3/                    # Gibbon
│   │   ├── HighConfNDD_Alu.bed
│   │   ├── Housekeeping_Alu.bed
│   │   ├── HighConfNDD_LINE1.bed
│   │   └── Housekeeping_LINE1.bed
│   ├── rheMac10/                   # Macaque
│   │   ├── HighConfNDD_Alu.bed
│   │   ├── Housekeeping_Alu.bed
│   │   ├── HighConfNDD_LINE1.bed
│   │   └── Housekeeping_LINE1.bed
│   ├── calJac4/                    # Marmoset
│   │   ├── HighConfNDD_Alu.bed
│   │   ├── Housekeeping_Alu.bed
│   │   ├── HighConfNDD_LINE1.bed
│   │   └── Housekeeping_LINE1.bed
│   ├── saiBol1/                    # Squirrel monkey (2nd Platyrrhine; Alu only)
│   │   ├── HighConfNDD_Alu.bed
│   │   ├── Housekeeping_Alu.bed
│   │   └── squirrel_stats.csv
│   ├── mmur3/                      # Mouse lemur (strepsirrhine, ~70 Mya; Alu only)
│   │   ├── HighConfNDD_Alu.bed
│   │   ├── Housekeeping_Alu.bed
│   │   ├── lemur_stats.csv
│   │   └── lemur_ortholog_validation.csv
│   ├── mm10/                       # Mouse (B1/B2 SINE analogs of primate Alu + LINE-1)
│   │   ├── HighConfNDD_B1B2.bed
│   │   ├── Housekeeping_B1B2.bed
│   │   ├── HighConfNDD_LINE1.bed
│   │   └── Housekeeping_LINE1.bed
│   ├── canFam4/                    # Dog (LINE-1 only; no Alu in dog)
│   │   ├── HighConfNDD_LINE1.bed
│   │   └── Housekeeping_LINE1.bed
│   ├── sensitivity/
│   │   ├── HighConfNDD_w500.bed
│   │   ├── HighConfNDD_w1000.bed
│   │   ├── HighConfNDD_w2000.bed
│   │   ├── HighConfNDD_w3000.bed
│   │   ├── Housekeeping_w500.bed
│   │   ├── Housekeeping_w1000.bed
│   │   ├── Housekeeping_w2000.bed
│   │   ├── Housekeeping_w3000.bed
│   │   └── window_sensitivity.csv
│   ├── cross_disease/
│   │   ├── Cardiovascular_promoters.bed
│   │   ├── Mendelian_promoters.bed
│   │   └── cross_disease_results.csv
│   ├── context/                    # Gene density + recombination controls (Fig 7)
│   ├── ccre/                       # ENCODE cCRE functional overlay (Fig 8)
│   ├── brain/                      # Fetal-brain DNase specificity overlay (Fig 9)
│   ├── matched/                    # Matched control + honest re-baselining
│   ├── consolidated/              # Master cross-species table + RESULTS_BACKBONE.md
│   ├── functional/                 # Exploratory constraint/expression + subfamily probes
│   └── gnomad_mei/                 # Exploratory polymorphic-MEI decomposition
│
└── figures/
    │   # File names are historical; the bracketed label is the number the figure
    │   # carries in the manuscript, where four of them are supplementary.
    ├── Fig1_Alu_Primates.pdf/.png     # [Fig 1]  Alu depletion across 6 primates (incl. mouse lemur)
    ├── Fig3_Heatmap.pdf/.png          # [Fig 2]  Significance heatmap (8 species × TE class)
    ├── Fig7_ContextControls.pdf/.png  # [Fig 3]  Gene density + recombination + matched control
    ├── Fig8_cCRE.pdf/.png             # [Fig 4]  ENCODE cCRE functional overlay
    ├── Fig9_BrainSpecificity.pdf/.png # [Fig 5]  Fetal-brain regulatory specificity (DNase-seq)
    ├── Fig2_LINE1_Mammals.pdf/.png    # [Fig S1] LINE-1 depletion across 7 mammals
    ├── Fig4_NullModel.pdf/.png        # [Fig S2] Permutation null model validation
    ├── Fig5_CpG.pdf/.png              # [Fig S3] CpG island confounder analysis
    ├── Fig6_GC_Analysis.pdf/.png      # [Fig S4] Promoter GC content analysis
    └── FigS5_PhyloEffect.pdf/.png     # [Fig S5] Effect size vs divergence time
```

---

## Data Sources

Raw data files must be downloaded separately. The following sources were used:

### Genome Annotations (GTF)

| Species | Assembly | Source |
|---------|----------|--------|
| Human | hg38, GENCODE v47 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/ |
| Orangutan | ponAbe3, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/pongo_abelii/ |
| Gibbon | nomLeu3, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/nomascus_leucogenys/ |
| Macaque | rheMac10, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/macaca_mulatta/ |
| Marmoset | calJac4, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/callithrix_jacchus/ |
| Squirrel monkey | SaiBol1.0, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/saimiri_boliviensis_boliviensis/ |
| Mouse lemur | Mmur_3.0, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/microcebus_murinus/ |
| Mouse | mm10, GENCODE vM25 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/ |
| Dog | canFam6, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/canis_lupus_familiaris/ |

### RepeatMasker Annotations

Retrieved from the UCSC Genome Browser for all seven species:

```
https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/database/rmsk.txt.gz
```

Replace `{assembly}` with: `hg38`, `ponAbe3`, `nomLeu3`, `rheMac10`, `calJac4`, `mm10`, `canFam6`

The **mouse lemur** RepeatMasker track (`SINE/Alu`) is instead taken from the UCSC
**GenArk** assembly hub (`GCF_000165445.2`, Mmur_3.0); its RefSeq sequence names are
mapped onto the Ensembl GTF names via the hub's `chromAlias.txt`:

```
https://hgdownload.soe.ucsc.edu/hubs/GCF/000/165/445/GCF_000165445.2/
```

> **Note on the dog assembly folder.** The `data/canFam4/` and `results/canFam4/`
> directories carry this name for historical reasons. The dog data they contain is in
> fact **ROS_Cfam_1.0 (UCSC canFam6)**: the GTF is from Ensembl release 112 and the
> RepeatMasker track from UCSC canFam6. The two are coordinate-compatible (they differ
> only by the `chr` prefix, which the promoter script adds), so the folder name is
> cosmetic only and is deliberately kept as-is to preserve every path reference in the
> scripts and committed results.

### Ortholog mapping

The 1:1 ortholog tables in `data/orthologs/` — `ensembl_to_symbol.tsv`,
`symbol_to_ensembl.csv`, and the per-species `<assembly>_raw.tsv` files consumed by
[`13_ortholog_analysis.py`](scripts/13_ortholog_analysis.py) — were generated from
**Ensembl BioMart (release 112)** and are committed directly to the repository. There is
no automated fetch script on purpose: the BioMart service changes between releases, so
shipping the exact tables used here is what makes the ortholog validation reproducible.

### Additional Data

| Dataset | Source |
|---------|--------|
| SFARI Gene 2.0 (Tier 1+2) | https://sfari.org/resource/sfari-gene |
| ClinGen Epilepsy GCEP | https://clinicalgenome.org |
| Human Phenotype Ontology | https://hpo.jax.org |
| HRT Atlas v1.0 | https://www.housekeeping.unicamp.br |
| GTEx brain expression | https://gtexportal.org |
| gnomAD v4.1 constraint | https://gnomad.broadinstitute.org |
| ClinVar variant summary | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ |
| hg38 CpG islands | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz |
| hg38 reference genome / 2bit | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit |
| ENCODE SCREEN cCRE registry (GRCh38) | https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed |
| ENCODE fetal DNase-seq peaks (GRCh38) | https://www.encodeproject.org — brain: ENCFF955AQD, ENCFF631TDE, ENCFF670PXX; non-neural: ENCFF667IEN (liver), ENCFF362PZG (lung), ENCFF016LYI (stomach) |
| Recombination map (GRCh38, deCODE-derived) | https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ |
| GTEx v8 median TPM by tissue | https://gtexportal.org |
| gnomAD v4.1 SV mobile-element insertions | https://gnomad.broadinstitute.org *(exploratory only)* |

---

## Requirements

Python dependencies are pinned in [`environment.yml`](environment.yml) (conda environment name `bio_master`):

```bash
conda env create -f environment.yml
conda activate bio_master
```

This provides Python 3.10 with pandas, numpy, scipy, matplotlib, seaborn, tqdm, and requests. The revision-stage scripts additionally use **twobitreader** (promoter GC from `hg38.2bit`) and **tabulate** (used by `22_consolidate.py`).

**System dependency (install separately):**

- **BEDTools ≥ 2.31** — used by the original scripts 05 and 08–10, 12. The revision-stage scripts (14–23, `fig_*`, `probe_*`) compute all interval overlaps with pure-Python (NumPy) routines and do **not** require BEDTools.

> RepeatMasker itself is **not** required. The pipeline parses pre-computed UCSC
> RepeatMasker tracks (`rmsk.txt.gz`); it never runs RepeatMasker locally.

---

## Reproducing the Analysis

All scripts resolve their paths **relative to the repository root** (via `BASH_SOURCE`
in the shell scripts and `Path(__file__)` in the Python scripts), so the pipeline runs
from a fresh clone on any machine — no path edits required.

### Quick start

```bash
# 1. Clone
git clone https://github.com/Bilmem2/ndd-te-conservation.git
cd ndd-te-conservation

# 2. Create and activate the conda environment
conda env create -f environment.yml
conda activate bio_master

# 3. Download raw genomes, GTF, and RepeatMasker tracks (not in the repo)
bash scripts/00_download_data.sh

# 4. Run the analysis in the order below.
```

The committed `results/`, `data/gene_lists/`, and `data/orthologs/` mean you can
regenerate **all statistics and figures** (steps 6–10 below) without the multi-GB
downloads. `00_download_data.sh` and steps 2–5 are only needed to rebuild the
intermediate promoter/TE BED files from scratch.

### Pipeline

```bash
# 1. Prepare gene lists
#    (SFARI, ClinGen, HPO, HRT Atlas, GTEx) placed in data/sources/. The outputs
#    are already committed under data/gene_lists/, so this step can be skipped.
python scripts/01_prepare_gene_lists.py

# 2. Extract Alu, LINE-1, and mouse B1/B2 BED files from RepeatMasker annotations
bash scripts/02_rmsk_to_bed.sh

# 3. Extract TSS ± 2 kb promoter windows from GTF files
bash scripts/03_get_promoters.sh

# 4. Split promoter BEDs by gene category (HighConfNDD, Housekeeping)
bash scripts/04_split_promoters.sh

# 5. Count TE overlaps per promoter window (BEDTools intersect)
bash scripts/05_intersect.sh

# 6. Window sensitivity analysis (± 0.5, 1, 2, 3 kb)
python scripts/09_window_sensitivity.py

# 7. Cross-disease specificity analysis (ClinVar).
#    Requires data/clinvar_variants.txt.gz. If ClinVar is not available locally,
#    26_cross_disease_recompute.py reproduces the reported numbers from the
#    committed gene lists instead.
python scripts/10_cross_disease.py

# 8. Final statistics across all 7 species (including mouse B1/B2)
python scripts/11_stats_updated.py

# 9. Generate core figures (Fig 1–6)
python scripts/12_figures_final.py

# 10. Ortholog-validated replication (Ensembl BioMart 1:1)
python scripts/13_ortholog_analysis.py

# 11. Genomic-context controls: gene density + recombination (Fig 7 data)
python scripts/15_context_controls.py

# 12. ENCODE cCRE functional overlay (Fig 8 data)
python scripts/16_ccre_overlay.py

# 12b. Fetal-brain DNase regulatory-specificity overlay (Fig 9 data)
python scripts/23_brain_overlay.py

# 13. Joint matched control + honest re-baselining
python scripts/18_matched_control.py
python scripts/19_rebaseline.py

# 14. Mouse lemur (strepsirrhine) Alu depletion + ortholog validation
python scripts/20_lemur.py
python scripts/21_lemur_ortholog.py

# 15. Permutation null model across all five species-TE combinations (Fig 4 data)
python scripts/24_null_model.py

# 15b. Within-species genome baseline for every species (Table 2).
#      Needs the raw GTF + RepeatMasker tracks from 00_download_data.sh.
python scripts/25_genome_baseline.py

# 15c. Specificity and sensitivity checks: constraint/brain-expression matched
#      control, canonical-TSS robustness, mouse B1 vs B2, FDR for Table 1
python scripts/27_constraint_matched.py
python scripts/28_mane_tss.py
python scripts/29_b1_b2_split.py

# 15d. Second Platyrrhine, resolving whether the marmoset result is species-specific
python scripts/30_squirrel_monkey.py

# 15e. Non-independence check: genomic thinning + one-per-family
python scripts/31_pseudoreplication.py

# 16. Regenerate Fig 1 & 3 (with lemur) and Fig 4 (5-panel null), build Fig 7, 8 & 9
python scripts/fig_lemur_update.py
python scripts/fig_null_update.py
python scripts/fig_new.py
python scripts/fig_brain.py

# 17. Consolidate all results into the master cross-species table + backbone
python scripts/22_consolidate.py
```

The remaining scripts — `06_statistics.py` (preliminary), `07_pli_correlation.py`,
`08_encode_overlap_v2.py`, `14_gnomad_mei.py`, `17_functional_consequence.py`,
`probe_dosage.py`, and `probe_subfamily.py` — contain **exploratory analyses that are
not used as evidence in the manuscript** (the polymorphic-MEI targeting-vs-selection
decomposition and the Alu subfamily-age probe are mentioned in the Discussion only as
inconclusive). They are retained for provenance and full transparency.

---
