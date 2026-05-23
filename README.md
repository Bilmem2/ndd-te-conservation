# Evolutionary Conservation of SINE and LINE-1 Depletion at NDD Gene Promoters Across Primates and Mammals

**Can Sevilmiş** | Department of Molecular Biology and Genetics, Bahçeşehir University, Istanbul, Turkey

**Repository:** https://github.com/Bilmem2/ndd-te-conservation  
**Related work:** Sevilmiş, C. (2026). LINE-1 Depletion at Promoters of Neurodevelopmental Disorder Genes: A Genome-Wide Analysis. *Preprints*. https://doi.org/10.20944/preprints202604.0715.v1

---

## Overview

This repository contains all analysis scripts, processed gene lists, statistical results, and figures for the above manuscript. Raw genome assemblies, RepeatMasker annotations, GTF files, and the hg38 reference FASTA are not included due to size constraints; download instructions are provided below.

---

## Repository Structure

```
.
├── scripts/
│   ├── 01_prepare_gene_lists.py    # Curate NDD, BroadNDD, and Housekeeping gene sets
│   ├── 02_rmsk_to_bed.sh           # Extract Alu, LINE-1, and mouse B1/B2 BED files from RepeatMasker
│   ├── 03_get_promoters.sh         # Extract TSS ± 2 kb promoter windows from GTF (multi-species)
│   ├── 04_get_promoters.sh         # Alternative promoter extraction script (retained for reference)
│   ├── 04_split_promoters.sh       # Split promoter BEDs by gene category
│   ├── 05_intersect.sh             # BEDTools intersect: TE count per promoter window
│   ├── 06_statistics.py            # Preliminary statistics (exploratory)
│   ├── 07_pli_correlation.py       # pLI vs Alu frequency correlation (exploratory)
│   ├── 08_encode_overlap_v2.py     # ENCODE CTCF/DNase overlap analysis (exploratory)
│   ├── 09_window_sensitivity.py    # Window size sensitivity analysis (± 0.5–3 kb)
│   ├── 10_cross_disease.py         # Cross-disease specificity analysis (ClinVar)
│   ├── 11_stats_updated.py         # Final statistics across all 7 species (incl. mouse B1/B2)
│   ├── 12_figures_final.py         # Generate all manuscript figures (Fig 1–6)
│   └── 13_ortholog_analysis.py     # Ortholog-validated replication (Ensembl BioMart 1:1)
│
├── data/
│   ├── gene_lists/
│   │   ├── HighConfNDD_genes.txt   # SFARI Tier 1+2 ∪ ClinGen Epilepsy Definitive/Strong (n=1020)
│   │   ├── BroadNDD_genes.txt      # HPO: autism, seizure, ADHD, ID (n=3360)
│   │   ├── Housekeeping_genes.txt  # HRT Atlas ∩ brain TPM≥1, NDD-free (n=1679)
│   │   ├── Cardiovascular_genes.txt
│   │   └── Mendelian_genes.txt
│   └── orthologs/
│       ├── ensembl_to_symbol.tsv   # Ensembl ID ↔ HGNC symbol mapping
│       └── symbol_to_ensembl.csv   # HGNC symbol ↔ Ensembl ID mapping
│
├── results/
│   ├── statistics_final.csv        # Main results: all species, symbol-based matching
│   ├── statistics_ortholog.csv     # Validation: Ensembl BioMart 1:1 ortholog-based matching
│   ├── null_model_new.csv          # Permutation test results (n=10,000)
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
│   └── cross_disease/
│       ├── Cardiovascular_promoters.bed
│       ├── Mendelian_promoters.bed
│       └── cross_disease_results.csv
│
└── figures/
    ├── Fig1_Alu_Primates.pdf/.png   # Alu/SINE depletion across 5 primates
    ├── Fig2_LINE1_Mammals.pdf/.png  # LINE-1 depletion across 7 mammals
    ├── Fig3_Heatmap.pdf/.png        # Significance heatmap (species × TE class, incl. mouse B1/B2)
    ├── Fig4_NullModel.pdf/.png      # Permutation null model (human/orangutan/macaque Alu; mouse B1/B2)
    ├── Fig5_CpG.pdf/.png            # CpG island confounder analysis (hg38)
    └── Fig6_GC_Analysis.pdf/.png    # Promoter GC content analysis (hg38)
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
| Mouse | mm10, GENCODE vM25 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/ |
| Dog | canFam6, Ensembl 112 | https://ftp.ensembl.org/pub/release-112/gtf/canis_lupus_familiaris/ |

### RepeatMasker Annotations

Retrieved from the UCSC Genome Browser for all seven species:

```
https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/database/rmsk.txt.gz
```

Replace `{assembly}` with: `hg38`, `ponAbe3`, `nomLeu3`, `rheMac10`, `calJac4`, `mm10`, `canFam6`

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
| hg38 reference genome | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz |

---

## Requirements

```
Python >= 3.10
pandas
scipy
matplotlib
seaborn
tqdm
requests
bedtools >= 2.31
```

---

## Reproducing the Analysis

Scripts are provided for transparency and documentation purposes. Full reproduction requires downloading the raw genome assemblies and gene annotation files listed under Data Sources, and placing them in the directory structure described above. Absolute paths in scripts may require adjustment for local environments.

The analysis follows this sequence:

```bash
# 1. Prepare gene lists (SFARI, ClinGen, HPO, HRT Atlas)
python scripts/01_prepare_gene_lists.py

# 2. Extract Alu, LINE-1, and mouse B1/B2 BED files from RepeatMasker annotations
bash scripts/02_rmsk_to_bed.sh

# 3. Extract TSS ± 2 kb promoter windows from GTF files
bash scripts/03_get_promoters.sh

# 4. Split promoter BEDs by gene category (HighConfNDD, BroadNDD, Housekeeping)
bash scripts/04_split_promoters.sh

# 5. Count TE overlaps per promoter window (BEDTools intersect)
bash scripts/05_intersect.sh

# 6. Window sensitivity analysis (± 0.5, 1, 2, 3 kb)
python scripts/09_window_sensitivity.py

# 7. Cross-disease specificity analysis (ClinVar)
python scripts/10_cross_disease.py

# 8. Final statistics across all 7 species (including mouse B1/B2)
python scripts/11_stats_updated.py

# 9. Generate manuscript figures (Fig 1–6)
python scripts/12_figures_final.py

# 10. Ortholog-validated replication (Ensembl BioMart 1:1)
python scripts/13_ortholog_analysis.py
```

Scripts 06–08 (`06_statistics.py`, `07_pli_correlation.py`, `08_encode_overlap_v2.py`) contain exploratory analyses not included in the final manuscript but retained for completeness.
