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

| Path | Contents |
|------|----------|
| `scripts/` | Numbered `00`–`49` in the order the analysis was built; the recipe below marks which of them reproduce the published results. `fig_*.py` draw the figures; `probe_*.py` are exploratory probes |
| `data/gene_lists/` | HighConfNDD (n = 1020) and Housekeeping (n = 1679), plus the ClinVar disease sets in a loose and a strict version — the manuscript uses the strict ones |
| `data/orthologs/` | Ensembl BioMart 1:1 ortholog tables, committed on purpose (see below) |
| `results/` — one per assembly | `hg38`, `ponAbe3`, `nomLeu3`, `rheMac10`, `calJac4`, `saiBol1`, `mmur3`, `mm10`, `canFam6`: promoter BEDs with TE counts, one file per gene set and element class |
| `results/consolidated/` | Master cross-species table, BH *q*-values |
| `results/context/` | Gene density and recombination controls |
| `results/matched/` | Matched controls, genome baselines, matching sensitivity, paired tests |
| `results/mechanism/` | Insertion opportunity, orientation, subfamily age |
| `results/sensitivity/` | Promoter window size, canonical TSS, promoter non-independence |
| `results/synteny/` | Promoter windows transferred between genomes by liftOver |
| `results/ccre/`, `results/brain/` | ENCODE cCRE and fetal-brain DNase overlays |
| `results/cross_disease/` | ClinVar comparison sets |
| `results/functional/`, `results/gnomad_mei/` | Exploratory analyses, not used as evidence |
| `figures/` | One file per figure, named for the number it carries in the submission |

## Figures

Each figure has a single producing script; `fig_new.py` draws two of them. This
table maps figure *files*, which are in `figures/`. The Online Resource
documents themselves are part of the manuscript and are not in this repository.

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
| *(none)* | Online Resource 6 | — | Supplementary **text**, not a figure: the permutation null, CpG stratification and GC analysis reported in full |
| `ESM7_Alu_Boxplots` | Online Resource 7 | (kept from the original submission) | Per-species Alu distributions |

---

## Data sources

### Gene annotations

| Species | Assembly | Source |
|---------|----------|--------|
| Human | hg38, GENCODE v47 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/ |
| Mouse | mm10, GENCODE vM25 | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/ |
| Orangutan, gibbon, macaque, marmoset, squirrel monkey, mouse lemur, dog | ponAbe3, nomLeu3, rheMac10, calJac4, SaiBol1.0, Mmur\_3.0, canFam6 | Ensembl release 112, `https://ftp.ensembl.org/pub/release-112/gtf/` — `00_download_data.sh` maps each assembly to its Ensembl species directory |

Ensembl's `ROS_Cfam_1.0` and UCSC's `canFam6` are the same dog assembly; the
repository uses `canFam6` throughout.

### RepeatMasker

`https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/database/rmsk.txt.gz` for
`hg38`, `ponAbe3`, `nomLeu3`, `rheMac10`, `calJac4`, `saiBol1`, `mm10`, `canFam6`.
RepeatMasker itself is never run; the pipeline parses these pre-computed tracks.

Two assemblies need extra handling. The **mouse lemur** track comes from the UCSC
GenArk hub (`GCF_000165445.2`) instead, and its RefSeq sequence names are mapped
onto the Ensembl GTF through the hub's `chromAlias.txt`; Alu is taken from the
hub's BED and LINE-1 is parsed from `GCF_000165445.2.repeatMasker.out.gz` by
`37_lemur_line1.py`. For the **squirrel monkey**, Ensembl uses versioned INSDC
accessions (`JH378105.1`) where UCSC uses the unversioned form, so the suffix is
stripped before matching.

### Ortholog tables

BioMart content changes between releases, so the tables that were used are
cached in `data/orthologs/` rather than refetched. Each was pulled from the
`hsapiens_gene_ensembl` dataset with the `<species>_homolog_*` attributes and
filtered to `ortholog_one2one`; the squirrel monkey prefix is `sbboliviensis`.

Two of them are committed, `mmur3_raw.tsv` and `saiBol1_raw.tsv`, so
`21_lemur_ortholog.py` and `38_squirrel_ortholog.py` run as they stand. The six
tables that `13_ortholog_analysis.py` needs — orangutan, gibbon, macaque,
marmoset, mouse and dog — are **not** in the repository. That script tries to
download them and, failing that, prints the exact BioMart query and file path
for each. At the last check (22 September 2026) the `martservice` endpoint was
answering "Service unavailable", so expect to build those six by hand from
[BioMart](https://www.ensembl.org/biomart/martview) if you want to rerun this
step. Its published output, `results/statistics_ortholog.csv`, is committed, so
the numbers quoted in Methods 2.6 can be checked without it.

### Other datasets

| Dataset | Source |
|---------|--------|
| SFARI Gene 2.0 (Tier 1+2) | https://sfari.org/resource/sfari-gene |
| ClinGen Epilepsy GCEP | https://search.clinicalgenome.org/kb/affiliate/40005 |
| HPO term gene lists (HP:0000729, HP:0001249, HP:0001250, HP:0007018) | https://hpo.jax.org |
| HRT Atlas v1.0 (housekeeping) | https://housekeeping.unicamp.br |
| GTEx v8 brain-region TPM (`rna_brain_gtex.tsv`) | https://www.proteinatlas.org/about/download |
| gnomAD v4.1 constraint (LOEUF, pLI) | https://gnomad.broadinstitute.org |
| ClinVar variant summary | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ |
| hg38 CpG islands (`cpgIslandExt`) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ |
| ENCODE SCREEN cCRE registry (GRCh38) | https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed |
| ENCODE fetal DNase-seq | brain ENCFF955AQD, ENCFF631TDE, ENCFF670PXX; non-neural ENCFF667IEN, ENCFF362PZG, ENCFF016LYI |
| Recombination map (GRCh38, deCODE-derived) | https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ |
| gnomAD v4.1 SV mobile-element insertions | https://gnomad.broadinstitute.org *(exploratory only)* |

Two inputs are derived from the downloads above rather than fetched directly,
and `00_download_data.sh` does not build them, so the scripts that read them
stop unless you make them first:

| File | Built from | Read by |
|------|-----------|---------|
| `data/hg38/gnomad_mei.tsv` | the `INS:ME:*` rows of `gnomad.v4.1.sv.sites.bed.gz`, keeping chrom, start, end, name, svtype, AN, AC, AF | `14_gnomad_mei.py`, `44_mei_matched.py` |
| `data/hg38/alu_detailed.bed` | the hg38 `rmsk` track, keeping chrom, start, end, repName, milliDiv for `SINE/Alu` | `probe_subfamily.py` |

Both belong to analyses outside the main recipe. The polymorphic-insertion
comparison that `44_mei_matched.py` produces is quoted in the manuscript
(Section 4.7) and is reported there as inconclusive; its output,
`results/gnomad_mei/matched_mei.csv`, is committed, so that number can be
checked without rebuilding the input.

The first five are downloaded by hand rather than by `00_download_data.sh`,
because each portal exports through its own button. None requires an account.
Put them in `data/sources/` under the names `01_prepare_gene_lists.py` expects:

```
SFARI-Gene_genes_<release>_<export>.csv
clingen_epilepsy.csv
genes_for_HP_0000729.txt   genes_for_HP_0001249.txt
genes_for_HP_0001250.txt   genes_for_HP_0007018.txt
Housekeeping_GenesHuman.csv
rna_brain_gtex.tsv
```

`data/sources/` is gitignored, so these are the one set of inputs a fresh clone
has to fetch before step 1. Everything else comes down with
`00_download_data.sh`. The gene lists these produce are committed under
`data/gene_lists/`, so steps 2 onwards run without them.

SFARI names its export after the release date, and its curation moves. Step 1
uses the 01 May 2026 release the published lists were built from when that file
is present. With any other export it prints which one it found, compares the
sets it would produce against the committed ones, writes the gene-by-gene
difference to `results/gene_set_source_diff.tsv`, and **leaves
`data/gene_lists/` untouched** — so a newer SFARI download cannot silently put
the pipeline out of step with the paper. Pass `--overwrite` to rebuild the lists
from whatever sources are present; every downstream result then belongs to that
release, not to the one reported.

---

## Requirements

```bash
conda env create -f environment.yml
conda activate bio_master
```

Python 3.10 with pandas, numpy, scipy, matplotlib, seaborn, **twobitreader** and
**tabulate**. **BEDTools ≥ 2.31** must be installed separately; the scripts that
do not call it compute interval overlaps in NumPy instead, and that routine was
checked against `bedtools intersect -c` on the hg38 CpG-island track. The
syntenic transfer also needs the UCSC `liftOver` binary in `tools/`.

---

## Reproducing the analysis

All scripts resolve paths relative to the repository root, so the pipeline runs
from a fresh clone with no path edits. Because `results/`, `data/gene_lists/` and
`data/orthologs/` are committed, **steps 1–5 are only needed to rebuild the
promoter and TE BED files from scratch.** Later steps still read individual
downloaded files — `hg38.2bit` for GC content, the gnomAD constraint table,
ClinVar, the cCRE registry, the fetal DNase peaks, the recombination map, or the
chain files. `00_download_data.sh` fetches all of them.

```bash
git clone https://github.com/Bilmem2/ndd-te-conservation.git
cd ndd-te-conservation
conda env create -f environment.yml && conda activate bio_master
bash scripts/00_download_data.sh          # raw genomes, GTF, RepeatMasker

# 1–5  build gene lists, TE BEDs, promoter windows, TE counts per promoter
python scripts/01_prepare_gene_lists.py   # see note below before running
bash   scripts/02_rmsk_to_bed.sh
bash   scripts/03_get_promoters.sh
bash   scripts/04_split_promoters.sh
bash   scripts/05_intersect.sh

# 6    core cross-species statistics
python scripts/11_stats_updated.py
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
python scripts/40_gc_analysis.py          # promoter GC + GC-stratified depletion; draws ESM 4
python scripts/41_dog_cansine.py          # dog Can-SINE boundary test
python scripts/50_genome_baseline_by_class.py  # genome baseline per element class, all
                                          # nine genomes (Table 4); needs 20, 37 and 41
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
python scripts/49_syntenic_promoters.py   # liftOver transfer; needs chain files

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

Two Online Resources are absent from this step. ESM 4 is drawn by
`40_gc_analysis.py` in step 10, which performs the GC analysis and its figure
together; ESM 6 is supplementary text and has no figure.

The recipe calls 53 of the 61 scripts. The other eight are kept for provenance:
`06_statistics.py` and `08_encode_overlap_v2.py`, superseded by
`11_stats_updated.py` and `16_ccre_overlay.py`; `07_pli_correlation.py`,
`14_gnomad_mei.py`, `17_functional_consequence.py` and `44_mei_matched.py`,
exploratory and not used as evidence — the polymorphic mobile-element
decomposition they produce is reported in the Discussion only as inconclusive;
and `probe_dosage.py` and `probe_subfamily.py`, feasibility probes written to
find out whether a question was worth pursuing before it was analysed properly.
The subfamily-age question `probe_subfamily.py` scoped is reported in the
Results, computed there with confidence intervals by `34_orientation_bias.py`.

---

## Citation

If you use this code or the derived gene sets, please cite:

> Sevilmiş, C. Conserved SINE and Lineage-Variable LINE-1 Depletion at
> Neurodevelopmental Disorder Promoters. Manuscript under review (2026).
