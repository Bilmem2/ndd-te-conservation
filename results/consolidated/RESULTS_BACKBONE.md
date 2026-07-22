# Results Backbone — Conserved SINE Depletion at NDD Promoters (revised)

_Auto-generated consolidation of all final analyses. Honest, re-baselined numbers._

## 1. Cross-species SINE-class depletion (headline, comparative genomics)

SINE = Alu in primates, B1/B2 in mouse. One-sided Mann-Whitney U (HK>NDD), rank-biserial r. **Mouse lemur (strepsirrhine, ~70 My) is new; ortholog-validated |delta r|=0.001.**

| species    |   mya | TE   |   n_HK |   n_NDD |   median_HK |   median_NDD |     p_value |      r | sig   |
|:-----------|------:|:-----|-------:|--------:|------------:|-------------:|------------:|-------:|:------|
| Human      |     0 | Alu  |   1714 |    1025 |        0.5  |         0.25 | 1.26854e-53 | -0.345 | ***   |
| Orangutan  |    16 | Alu  |   1442 |     910 |        0.5  |         0.25 | 5.11472e-53 | -0.367 | ***   |
| Gibbon     |    20 | Alu  |   1280 |     823 |        0.5  |         0.25 | 7.45258e-31 | -0.291 | ***   |
| Macaque    |    25 | Alu  |   1533 |     949 |        0.5  |         0.25 | 1.08478e-60 | -0.384 | ***   |
| Marmoset   |    40 | Alu  |   1340 |     908 |        0.5  |         0.5  | 9.01879e-10 | -0.147 | ***   |
| MouseLemur |    70 | Alu  |   1324 |     899 |        0.25 |         0    | 2.70443e-29 | -0.262 | ***   |
| Mouse      |    90 | B1B2 |   1605 |     960 |        0.75 |         0.25 | 1.00897e-70 | -0.413 | ***   |

## 2. LINE-1 (internal contrast: weak, lineage-variable)

| species   |   mya | TE    |   n_HK |   n_NDD |   median_HK |   median_NDD |    p_value |      r | sig   |
|:----------|------:|:------|-------:|--------:|------------:|-------------:|-----------:|-------:|:------|
| Human     |     0 | LINE1 |   1714 |    1025 |        0    |         0    | 0.00371898 | -0.049 | **    |
| Orangutan |    16 | LINE1 |   1442 |     910 |        0    |         0    | 0.00965963 | -0.044 | **    |
| Gibbon    |    20 | LINE1 |   1280 |     823 |        0    |         0    | 0.272396   | -0.013 | ns    |
| Macaque   |    25 | LINE1 |   1533 |     949 |        0    |         0    | 0.00399077 | -0.048 | **    |
| Marmoset  |    40 | LINE1 |   1340 |     908 |        0.25 |         0.25 | 0.736254   |  0.015 | ns    |
| Mouse     |    90 | LINE1 |   1605 |     960 |        0    |         0    | 0.0591572  | -0.025 | ns    |
| Dog       |    95 | LINE1 |   1497 |     880 |        0.25 |         0.25 | 0.853864   |  0.025 | ns    |

## 3. Supporting rigor (all human hg38)

**3a. Multivariate matched control (GC + gene density + recombination):** depletion persists after joint matching.

| comparison                  |   n_HK |   n_NDD |   median_alu_HK |   median_alu_NDD |   mean_alu_HK |   mean_alu_NDD |     p_value |      r |
|:----------------------------|-------:|--------:|----------------:|-----------------:|--------------:|---------------:|------------:|-------:|
| NDD vs ALL housekeeping     |   1712 |    1021 |             0.5 |             0.25 |         0.665 |          0.381 | 1.24516e-53 | -0.345 |
| NDD vs MATCHED housekeeping |   1021 |    1021 |             0.5 |             0.25 |         0.581 |          0.381 | 1.96167e-27 | -0.269 |

**3b. Honest re-baselining (vs genome & expression-breadth+GC-matched control):** the dramatic vs-housekeeping ratio was partly housekeeping Alu-enrichment; a modest (~18%) NDD-specific fixed-Alu deficit survives; polymorphic shows no deficit.

| reference                   |   fixed_NDD |   fixed_ratio | fixed_CI                               |   poly_NDD |   poly_ratio | poly_CI                                |
|:----------------------------|------------:|--------------:|:---------------------------------------|-----------:|-------------:|:---------------------------------------|
| GENOME (all PC genes)       |       0.377 |         0.773 | ('—', '—')                             |     0.0488 |        1.007 | ('—', '—')                             |
| MATCHED non-disease control |       0.377 |         0.817 | (np.float64(0.738), np.float64(0.904)) |     0.0488 |        1.089 | (np.float64(0.882), np.float64(1.333)) |

**3c. Context (NDD vs HK differ, but depletion survives strata):**

| metric       |   n_NDD |   n_HK |   median_NDD |   median_HK |   mean_NDD |   mean_HK |     p_value |      r |
|:-------------|--------:|-------:|-------------:|------------:|-----------:|----------:|------------:|-------:|
| gene_density |    1025 |   1714 |        8     |      15     |     13.155 |    18.944 | 7.46006e-47 |  0.328 |
| recomb_cM_Mb |    1021 |   1712 |        0.729 |       0.561 |      1.273 |     0.964 | 1.13031e-07 | -0.121 |

Alu depletion within gene-density quartiles:

| stratum   | range       |   n_HK |   n_NDD |   median_HK |   median_NDD |     p_value |      r |
|:----------|:------------|-------:|--------:|------------:|-------------:|------------:|-------:|
| Q1        | 0.00-6.00   |    239 |     446 |        0.25 |         0    | 1.3528e-14  | -0.329 |
| Q2        | 6.00-13.00  |    480 |     205 |        0.5  |         0.25 | 2.00413e-08 | -0.26  |
| Q3        | 13.00-25.00 |    490 |     194 |        0.75 |         0.5  | 5.8261e-06  | -0.212 |
| Q4        | 25.00-69.00 |    505 |     180 |        0.75 |         0.5  | 3.11055e-06 | -0.224 |

Alu depletion within recombination quartiles:

| stratum   | range      |   n_HK |   n_NDD |   median_HK |   median_NDD |     p_value |      r |
|:----------|:-----------|-------:|--------:|------------:|-------------:|------------:|-------:|
| Q1        | 0.00-0.22  |    457 |     227 |         0.5 |         0.25 | 1.75404e-09 | -0.273 |
| Q2        | 0.22-0.62  |    445 |     238 |         0.5 |         0.25 | 5.87143e-09 | -0.26  |
| Q3        | 0.62-1.51  |    447 |     236 |         0.5 |         0.25 | 2.17688e-16 | -0.371 |
| Q4        | 1.51-18.86 |    363 |     320 |         0.5 |         0    | 2.8647e-24  | -0.436 |

**3d. Functional overlay (ENCODE cCRE):** NDD promoters are regulatory-active; Alu-free promoters are more cCRE-dense.

| cCRE_group   |   median_NDD |   median_HK |   mean_NDD |   mean_HK |     p_value |      r | direction   |
|:-------------|-------------:|------------:|-----------:|----------:|------------:|-------:|:------------|
| PLS          |         0.5  |         0.5 |      0.466 |     0.467 | 0.478429    |  0.015 | NDD<HK      |
| ELS          |         1.75 |         2   |      1.757 |     1.857 | 0.00435627  |  0.064 | NDD<HK      |
| CTCF         |         0    |         0   |      0.009 |     0.003 | 1.34991e-05 | -0.024 | NDD>HK      |
| active_all   |         2.5  |         2.5 |      2.303 |     2.387 | 0.258861    |  0.026 | NDD<HK      |

Alu-free vs Alu-positive regulatory density:

| cCRE_group   |   n_alu_free |   n_alu_pos |   median_free |   median_pos |   p_free_gt_pos |      r |
|:-------------|-------------:|------------:|--------------:|-------------:|----------------:|-------:|
| PLS          |          704 |        2035 |             2 |            2 |     1.2189e-06  | -0.113 |
| ELS          |          704 |        2035 |             8 |            7 |     1.59768e-11 | -0.166 |
| CTCF         |          704 |        2035 |             0 |            0 |     0.0833553   | -0.008 |
| active_all   |          704 |        2035 |            11 |           10 |     1.95958e-20 | -0.23  |

## 4. Mechanism probes (reported honestly as inconclusive)

**4a. gnomAD MEI targeting-vs-selection decomposition** (deflated after honest re-baselining; AFS not significant vs matched control):

| insertion_class             |   dens_HK |   dens_NDD |   ratio_NDD_HK | ci95           |
|:----------------------------|----------:|-----------:|---------------:|:---------------|
| FIXED Alu (RepeatMasker)    |    0.666  |     0.383  |          0.575 | [0.526, 0.627] |
| POLYMORPHIC Alu (gnomAD SV) |    0.0519 |     0.0483 |          0.93  | [0.776, 1.110] |

**4b. Alu subfamily-age** (depletion ~uniform across AluJ/S/Y => no distinctive age-dependent selection signature; stably refractory environment):

| subfamily   |   mean_age_milliDiv |   dens_HK |   dens_NDD |   ratio_NDD_HK |   p_depletion |      r |
|:------------|--------------------:|----------:|-----------:|---------------:|--------------:|-------:|
| AluY        |                68.6 |     0.066 |      0.037 |          0.564 |   5.54677e-10 | -0.094 |
| AluS        |               111.9 |     0.398 |      0.231 |          0.579 |   8.84602e-43 | -0.3   |
| AluJ        |               165   |     0.164 |      0.097 |          0.592 |   3.10976e-19 | -0.174 |

## 5. Honest framing / limitations

- Phenomenon of Alu-depletion at developmental promoters + Alu/B1 convergence is **prior art** (Polak & Domany 2006; Tsirigos & Rigoutsos 2009; Simons/Mattick TFRs). Our contribution = systematic cross-species quantification incl. a strepsirrhine, modern clinical NDD sets, and comprehensive confounder control. Do NOT claim first discovery of the phenomenon.

- Mechanism (targeting vs post-insertion selection) **remains open**: population (gnomAD MEI) and subfamily-age probes are inconclusive once properly controlled.

- Magnitude is modest vs an honest (genome / breadth-matched) baseline; the large vs-housekeeping effect was partly housekeeping Alu-enrichment (Eller et al. 2007).

- Functional validation (wet-lab) is future work / capstone.
