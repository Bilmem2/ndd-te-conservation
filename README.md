# Conserved SINE and Lineage-Variable LINE-1 Depletion at Neurodevelopmental Disorder Promoters

**Can Sevilmiş** | Department of Molecular Biology and Genetics, Bahçeşehir University, Istanbul, Turkey

**Related work:** Sevilmiş, C. (2026). LINE-1 Depletion at Promoters of Neurodevelopmental Disorder Genes: A Genome-Wide Analysis. *Preprints*. https://doi.org/10.20944/preprints202604.0715.v1

---

## Overview

Analysis code, gene lists, statistical results and figures for the manuscript above.

The study compares transposable element density at the promoters (TSS ± 2 kb) of
high-confidence neurodevelopmental disorder (NDD) genes against brain-expressed
housekeeping genes across **nine mammalian genomes**: seven primates — human, orangutan,
gibbon, macaque, two New World monkeys (marmoset, squirrel monkey) and the strepsirrhine
gray mouse lemur — plus mouse and dog. Beyond the primary contrast it adds within-species
genome baselines, a constraint- and expression-matched control, functional overlays
(ENCODE cCREs, fetal-brain DNase-seq) and a set of robustness and mechanism checks.

Raw genomes, GTFs and RepeatMasker tracks are too large to commit; `00_download_data.sh`
fetches them. Everything needed to regenerate the **statistics and figures** is committed.

---

## Repository layout

```
scripts/     01–41, in pipeline order (see "Reproducing the analysis")
data/
  gene_lists/    HighConfNDD (n=1020), Housekeeping (n=1679), Cardiovascular, Mendelian
  orthologs/     Ensembl BioMart 1:1 ortholog tables (committed; see note below)
results/
  <assembly>/    per-species promoter BEDs with TE counts, one file per gene set and TE class
                 (hg38, ponAbe3, nomLeu3, rheMac10, calJac4, saiBol1, mmur3, mm10, canFam4)
  consolidated/  master cross-species table, BH q-values, results backbone
  context/       gene density + recombination controls
  matched/       matched controls, genome baselines, matching sensitivity
  mechanism/     insertion opportunity, orientation, subfamily age
  sensitivity/   promoter window, canonical TSS, promoter non-independence
  ccre/ brain/   ENCODE cCRE and fetal-brain DNase overlays
  cross_disease/ ClinVar comparison sets
  functional/ gnomad_mei/   exploratory analyses, not used as evidence
figures/     see mapping below
```

Figure files keep their original names; the bracketed label is the number the figure
carries in the manuscript, where four of the nine are supplementary.

| File | Manuscript | Content |
|------|-----------|---------|
| `Fig1_Alu_Primates` | Fig 1 | Alu depletion across seven primates |
| `Fig3_Heatmap` | Fig 2 | Significance heatmap, nine species × TE class |
| `Fig7_ContextControls` | Fig 3 | Gene density, recombination, matched control |
| `Fig8_cCRE` | Fig 4 | ENCODE cCRE overlay |
| `Fig9_BrainSpecificity` | Fig 5 | Fetal-brain regulatory specificity |
| `Fig2_LINE1_Mammals` | Fig S1 | LINE-1 across nine mammals |
| `Fig4_NullModel` | Fig S2 | Permutation null model |
| `Fig5_CpG` | Fig S3 | CpG island stratification |
| `Fig6_GC_Analysis` | Fig S4 | Promoter GC content |
| `FigS5_PhyloEffect` | Fig S5 | Effect size vs divergence time |

---

## Data sources

### Gene annotations

| Species | Assembly | Source |
|---------|----------|--------|
| Human | hg38, GENCODE v47 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/ |
| Mouse | mm10, GENCODE vM25 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/ |
| Orangutan, gibbon, macaque, marmoset, squirrel monkey, mouse lemur, dog | ponAbe3, nomLeu3, rheMac10, calJac4, SaiBol1.0, Mmur\_3.0, canFam6 | Ensembl release 112, `https://ftp.ensembl.org/pub/release-112/gtf/<species>/` |

### RepeatMasker

`https://hgdownload.soe.ucsc.edu/goldenPath/{assembly}/database/rmsk.txt.gz` for
`hg38`, `ponAbe3`, `nomLeu3`, `rheMac10`, `calJac4`, `saiBol1`, `mm10`, `canFam6`.

The **mouse lemur** track comes instead from the UCSC GenArk hub
(`GCF_000165445.2`); its RefSeq sequence names are mapped onto the Ensembl GTF via the
hub's `chromAlias.txt`. Alu was taken from that hub's BED; LINE-1 is parsed from the
hub's `GCF_000165445.2.repeatMasker.out.gz` by `37_lemur_line1.py`, which writes
`data/mmur3/line1_refseq.bed` on first run. For the **squirrel monkey**, Ensembl uses
versioned INSDC accessions (`JH378105.1`) where UCSC uses the unversioned form, so the
suffix is stripped before matching.

Two further downloads are fetched on demand rather than by `00_download_data.sh`:
the hg38 CpG-island track (`cpgIslandExt`, used by `12_figures_final.py` for the CpG
stratification figure) and the mouse lemur RepeatMasker output above. Both are excluded
from version control by size.

> **Dog folder name.** `data/canFam4/` and `results/canFam4/` are named for historical
> reasons; the data they hold is **ROS_Cfam_1.0 (UCSC canFam6)**. The two are
> coordinate-compatible apart from the `chr` prefix that the promoter script adds, so the
> name is cosmetic and is kept to preserve every path reference in the committed results.

### Ortholog tables

The 1:1 ortholog tables in `data/orthologs/` were generated from **Ensembl BioMart
(release 112)** and are committed directly. There is no fetch script on purpose: BioMart
changes between releases, and shipping the exact tables is what makes the ortholog
validation reproducible.

### Other datasets

| Dataset | Source |
|---------|--------|
| SFARI Gene 2.0 (Tier 1+2) | https://sfari.org/resource/sfari-gene |
| ClinGen Epilepsy GCEP | https://clinicalgenome.org |
| HRT Atlas v1.0 (housekeeping) | https://www.housekeeping.unicamp.br |
| GTEx v8 median TPM by tissue | https://gtexportal.org |
| gnomAD v4.1 constraint (LOEUF, pLI) | https://gnomad.broadinstitute.org |
| ClinVar variant summary | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ |
| hg38 CpG islands / 2bit | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/ |
| ENCODE SCREEN cCRE registry (GRCh38) | https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed |
| ENCODE fetal DNase-seq | brain ENCFF955AQD, ENCFF631TDE, ENCFF670PXX; non-neural ENCFF667IEN, ENCFF362PZG, ENCFF016LYI |
| Recombination map (GRCh38, deCODE-derived) | https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ |
| gnomAD v4.1 SV mobile-element insertions | https://gnomad.broadinstitute.org *(exploratory only)* |

---

## Requirements

```bash
conda env create -f environment.yml
conda activate bio_master
```

Python 3.10 with pandas, numpy, scipy, matplotlib, seaborn, plus **twobitreader**
(promoter GC from `hg38.2bit`) and **tabulate**.

**BEDTools ≥ 2.31** must be installed separately; it is used by scripts 05 and 08–11.
Script 12 and everything after it compute interval overlaps in NumPy and do not need it. RepeatMasker itself
is never run — the pipeline parses pre-computed UCSC tracks.

---

## Reproducing the analysis

All scripts resolve paths relative to the repository root, so the pipeline runs from a
fresh clone with no path edits. Because `results/`, `data/gene_lists/` and
`data/orthologs/` are committed, **steps 6 onwards run without the multi-GB downloads**;
steps 1–5 are only needed to rebuild the promoter and TE BED files from scratch.

```bash
git clone https://github.com/Bilmem2/ndd-te-conservation.git
cd ndd-te-conservation
conda env create -f environment.yml && conda activate bio_master
bash scripts/00_download_data.sh          # raw genomes, GTF, RepeatMasker

# 1–5  build gene lists, TE BEDs, promoter windows, TE counts per promoter
python scripts/01_prepare_gene_lists.py
bash   scripts/02_rmsk_to_bed.sh
bash   scripts/03_get_promoters.sh
bash   scripts/04_split_promoters.sh
bash   scripts/05_intersect.sh

# 6    core cross-species statistics and figures
python scripts/11_stats_updated.py
python scripts/12_figures_final.py
python scripts/13_ortholog_analysis.py    # 1:1 ortholog validation

# 7    promoter window sizes and cross-disease sets
python scripts/09_window_sensitivity.py
python scripts/10_cross_disease.py        # needs ClinVar; 26_ reproduces it without
python scripts/26_cross_disease_recompute.py

# 8    genomic-context controls and honest baselines
python scripts/15_context_controls.py
python scripts/18_matched_control.py
python scripts/19_rebaseline.py
python scripts/25_genome_baseline.py

# 9    added species
python scripts/20_lemur.py
python scripts/21_lemur_ortholog.py
python scripts/30_squirrel_monkey.py

# 10   specificity, robustness, mechanism
python scripts/24_null_model.py           # permutation, all species-TE combinations
python scripts/27_constraint_matched.py   # LOEUF + brain expression + GC matched control
python scripts/28_mane_tss.py             # canonical TSS
python scripts/29_b1_b2_split.py          # mouse B1 vs B2, and BH q-values
python scripts/31_pseudoreplication.py    # clustered and paralogous promoters
python scripts/32_matching_sensitivity.py # matching order, metric, covariates
python scripts/33_insertion_opportunity.py# L1 endonuclease site density
python scripts/34_orientation_bias.py     # orientation and subfamily age
python scripts/35_alu_age_by_species.py   # Alu divergence per species
python scripts/36_loeuf_gradient.py       # Alu density across LOEUF deciles
python scripts/37_lemur_line1.py          # mouse lemur LINE-1
python scripts/38_squirrel_ortholog.py    # squirrel monkey 1:1 ortholog control
python scripts/39_squirrel_line1.py       # squirrel monkey LINE-1 (completes the panel)
python scripts/40_gc_analysis.py          # promoter GC + GC-stratified depletion (Fig S4)
python scripts/41_dog_cansine.py          # dog Can-SINE boundary test

# 11   functional overlays
python scripts/16_ccre_overlay.py
python scripts/23_brain_overlay.py

# 12   figures and consolidation
python scripts/fig_lemur_update.py        # Fig 1, Fig 2
python scripts/fig_null_update.py         # Fig S2
python scripts/fig_new.py                 # Fig 3, Fig 4
python scripts/fig_brain.py               # Fig 5
python scripts/fig_phylo.py               # Fig S5
python scripts/22_consolidate.py
```

`06_statistics.py`, `07_pli_correlation.py`, `08_encode_overlap_v2.py`,
`14_gnomad_mei.py`, `17_functional_consequence.py`, `probe_dosage.py` and
`probe_subfamily.py` hold **exploratory analyses that are not used as evidence in the
manuscript**. They are kept for provenance; the polymorphic-MEI decomposition is
mentioned in the Discussion only as inconclusive.

---

## Citation

If you use this code or the derived gene sets, please cite the manuscript. Raw
annotations remain subject to the terms of their original providers.
