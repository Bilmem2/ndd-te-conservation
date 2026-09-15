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
scripts/     00–49 in pipeline order (see "Reproducing the analysis"), plus
             fig_*.py, which draw the figures, and probe_*.py, exploratory
data/
  gene_lists/    HighConfNDD (n=1020), Housekeeping (n=1679), and the ClinVar
                 sets in two versions: Cardiovascular/Mendelian_genes.txt from a
                 substring match, and the _strict.txt files the manuscript uses,
                 built by the explicit P/LP rule in 48_cross_disease_strict.py
  orthologs/     Ensembl BioMart 1:1 ortholog tables (committed; see note below)
results/
  hg38/ ponAbe3/ nomLeu3/ rheMac10/ calJac4/ saiBol1/ mmur3/ mm10/ canFam6/
                 one directory per assembly: promoter BEDs with TE counts, one
                 file per gene set and TE class
  consolidated/  master cross-species table, BH q-values, results backbone
  context/       gene density + recombination controls
  matched/       matched controls, genome baselines, matching sensitivity, paired tests
  mechanism/     insertion opportunity, orientation, subfamily age
  sensitivity/   promoter window, canonical TSS, promoter non-independence
  synteny/       promoter windows transferred between genomes by liftOver
  ccre/ brain/   ENCODE cCRE and fetal-brain DNase overlays
  cross_disease/ ClinVar comparison sets
  functional/ gnomad_mei/   exploratory analyses, not used as evidence
figures/     see mapping below
```

Each file is named for the number it carries in the submission: `Fig<n>` for the
six figures inside the article, `ESM<n>` for those submitted separately as Online
Resources. Online Resource 6 is supplementary text and has no figure.

| File | In the submission | Produced by | Content |
|------|------------------|-------------|---------|
| `Fig1_Alu_Primates` | Figure 1 | `fig_alu_primates.py` | Alu depletion across seven primates |
| `Fig2_Heatmap` | Figure 2 | `fig_heatmap_effect.py` | Effect-size heatmap, nine genomes × TE class |
| `Fig3_ContextControls` | Figure 3 | `fig_new.py` | Gene density, recombination, matched control |
| `Fig4_cCRE` | Figure 4 | `fig_new.py` | ENCODE cCRE overlay |
| `Fig5_BrainSpecificity` | Figure 5 | `fig_brain.py` | Fetal-brain regulatory specificity |
| `Fig6_Flanking` | Figure 6 | `fig_flanking.py` | Promoter versus regional depletion |
| `ESM1_LINE1_Mammals` | Online Resource 1 | `fig_line1_esm1.py` | SINE and LINE-1 density across nine genomes |
| `ESM2_NullModel` | Online Resource 2 | `fig_null_update.py` | Permutation null model |
| `ESM3_CpG` | Online Resource 3 | `12_figures_final.py` | CpG island stratification |
| `ESM4_GC_Analysis` | Online Resource 4 | `40_gc_analysis.py` | Promoter GC content |
| `ESM5_PhyloEffect` | Online Resource 5 | `fig_phylo.py` | Effect size vs divergence time |
| `ESM7_Alu_Boxplots` | Online Resource 7 | (kept from the original submission) | Per-species Alu distributions |

Exactly one script writes each figure. `12_figures_final.py` is the original
monolithic figure script and once wrote four of these as well; those `savefig`
calls are commented out, because running it in pipeline order otherwise reverted
Figure 1 to a five-species boxplot from before the panel was extended and the
heatmap to the colour scale the review objected to.

---

## Data sources

### Gene annotations

| Species | Assembly | Source |
|---------|----------|--------|
| Human | hg38, GENCODE v47 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/ |
| Mouse | mm10, GENCODE vM25 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/ |
| Orangutan, gibbon, macaque, marmoset, squirrel monkey, mouse lemur, dog | ponAbe3, nomLeu3, rheMac10, calJac4, SaiBol1.0, Mmur\_3.0, canFam6 | Ensembl release 112, `https://ftp.ensembl.org/pub/release-112/gtf/<species>/` |

### RepeatMasker

`https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/database/rmsk.txt.gz` for
`hg38`, `ponAbe3`, `nomLeu3`, `rheMac10`, `calJac4`, `saiBol1`, `mm10`, `canFam6`.

The **mouse lemur** track comes instead from the UCSC GenArk hub
(`GCF_000165445.2`); its RefSeq sequence names are mapped onto the Ensembl GTF via the
hub's `chromAlias.txt`. Alu was taken from that hub's BED; LINE-1 is parsed from the
hub's `GCF_000165445.2.repeatMasker.out.gz` by `37_lemur_line1.py`, which writes
`data/mmur3/line1_refseq.bed` on first run. For the **squirrel monkey**, Ensembl uses
versioned INSDC accessions (`JH378105.1`) where UCSC uses the unversioned form, so the
suffix is stripped before matching.

`00_download_data.sh` covers all nine assemblies, including the mouse lemur GenArk
files above and `hg38.2bit`, from which promoter GC content is read directly. The
hg38 CpG-island track (`cpgIslandExt`, used by `12_figures_final.py` for the CpG
stratification figure) is fetched on demand instead, as are the human-only
auxiliary sources in the table below. All are excluded from version control by
size.

> **Dog.** Ensembl's `ROS_Cfam_1.0` and UCSC's `canFam6` are the same assembly.
> The repository uses `canFam6` throughout.

### Ortholog tables

The 1:1 ortholog tables in `data/orthologs/` were generated from **Ensembl BioMart**
and are committed directly. There is no fetch script on purpose: BioMart changes between
releases, and shipping the exact tables is what makes the ortholog validation
reproducible. Each was pulled from the `hsapiens_gene_ensembl` dataset with the
`<species>_homolog_associated_gene_name` and `<species>_homolog_orthology_type`
attributes, filtered to `ortholog_one2one`; the squirrel monkey prefix is
`sbboliviensis`.

### Other datasets

| Dataset | Source |
|---------|--------|
| SFARI Gene 2.0 (Tier 1+2) | https://sfari.org/resource/sfari-gene |
| ClinGen Epilepsy GCEP | https://clinicalgenome.org |
| HRT Atlas v1.0 (housekeeping) | https://www.housekeeping.unicamp.br |
| GTEx v8 median TPM by tissue | https://gtexportal.org |
| gnomAD v4.1 constraint (LOEUF, pLI) | https://gnomad.broadinstitute.org |
| ClinVar variant summary | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ |
| hg38 CpG islands (`cpgIslandExt`) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ |
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

**BEDTools ≥ 2.31** must be installed separately. It is called by `05_intersect.sh`,
`08_encode_overlap_v2.py`, `09_window_sensitivity.py`, `10_cross_disease.py`,
`14_gnomad_mei.py`, `15_context_controls.py` and `16_ccre_overlay.py`. Every other
script computes interval overlaps in NumPy, including `12_figures_final.py`, whose
NumPy counter was checked against `bedtools intersect -c` and returns identical counts.
RepeatMasker itself is never run — the pipeline parses pre-computed UCSC tracks.

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
python scripts/33_insertion_opportunity.py # L1 endonuclease site density
python scripts/34_orientation_bias.py     # orientation and subfamily age
python scripts/35_alu_age_by_species.py   # Alu divergence per species
python scripts/36_loeuf_gradient.py       # Alu density across LOEUF deciles
python scripts/37_lemur_line1.py          # mouse lemur LINE-1
python scripts/38_squirrel_ortholog.py    # squirrel monkey 1:1 ortholog control
python scripts/39_squirrel_line1.py       # squirrel monkey LINE-1 (completes the panel)
python scripts/40_gc_analysis.py          # promoter GC + GC-stratified depletion (ESM 4)
python scripts/41_dog_cansine.py          # dog Can-SINE boundary test
python scripts/42_coverage_robustness.py  # Alu as merged bp coverage, not record counts
python scripts/43_flanking_control.py     # promoter vs flanking windows out to 250 kb

# 11   functional overlays
python scripts/16_ccre_overlay.py
python scripts/23_brain_overlay.py

# 12   analyses added for the revision
python scripts/45_promoter_vs_local.py    # promoter against its own regional background
python scripts/46_paired_matched_tests.py # Wilcoxon, sign-flip permutation, Kerby r
python scripts/47_line1_floor.py          # is an Alu-sized deficit detectable at LINE-1 density
python scripts/48_cross_disease_strict.py # explicit P/LP rule for the ClinVar sets
python scripts/49_syntenic_promoters.py   # liftOver transfer; needs tools/liftOver + chain files

# 13   figures and consolidation
python scripts/fig_alu_primates.py        # Fig 1
python scripts/fig_heatmap_effect.py      # Fig 2
python scripts/fig_new.py                 # Fig 3, Fig 4
python scripts/fig_brain.py               # Fig 5
python scripts/fig_flanking.py            # Fig 6
python scripts/fig_line1_esm1.py          # ESM 1
python scripts/fig_null_update.py         # ESM 2
python scripts/12_figures_final.py        # ESM 3
python scripts/fig_phylo.py               # ESM 5
python scripts/22_consolidate.py
```

Exactly one script writes each figure, and `12_figures_final.py` must run before
`fig_alu_primates.py`, `fig_heatmap_effect.py`, `fig_line1_esm1.py` and
`fig_null_update.py` if it is run at all: it once wrote those four figures too,
and although those writes are now commented out, keeping the order removes the
question. `40_gc_analysis.py` in step 10 draws ESM 4.

Not in the recipe: `06`, `07`, `08` and `17`, superseded by `11`, `16` and `23`;
`14` and `44`, which produce the exploratory gnomAD mobile-element analysis that
the manuscript reports but does not rely on; and `probe_*.py`, which are scratch
checks.

`06_statistics.py`, `07_pli_correlation.py`, `08_encode_overlap_v2.py`,
`14_gnomad_mei.py`, `17_functional_consequence.py`, `probe_dosage.py` and
`probe_subfamily.py` hold **exploratory analyses that are not used as evidence in the
manuscript**. They are kept for provenance; the polymorphic-MEI decomposition is
mentioned in the Discussion only as inconclusive.

---

## Citation

If you use this code or the derived gene sets, please cite the manuscript. Raw
annotations remain subject to the terms of their original providers.
