# Plotting Scripts

This directory contains all final visualization scripts used in the DNA repair and damage analysis pipeline.

---

## 📊 R Scripts (`plots/R/`)

| Script | Description |
|--------|-------------|
| `plot_tssquartile_damage.R` | TSS metaprofile of Damage-seq by expression quartile and strand |
| `plot_tssquartile_xr.R` | Same as above for XR-seq (adjusts y-axis per strain) |
| `TSNTSviolin_expression_quartiles.R` | Violin plots comparing TS vs NTS strand RPKM across quartiles |
| `timecoursetstototalreadsratio_xr_expressionquartiles.R` | XR-seq TS/(TS+NTS) ratio across timepoints and expression quartiles |
| `spearmancorr2kbwindowhclust.R` | Spearman correlation + hierarchical clustering heatmap from 2kb RPKM/RPM matrices |
| `boxplot_epigenome_XR.R` | XR-seq signals stratified by histone-mark tertiles (2 kb bins; simulation-normalized; all strains) |
| `boxplot_epigenome_Damage.R` | Damage-seq signals stratified by histone-mark tertiles (2 kb bins; simulation-normalized; wt only) |

---

## 📊 Python Scripts — ATAC-seq (`plots/atacseq/`)

| Script | Description |
|--------|-------------|
| `atacseq_damage_allpeaks.py` | Damage-seq over all ATAC peaks (log₂ sim-normalized) |
| `atacseq_xr_allpeaks.py` | XR-seq over all ATAC peaks (RPM) |
| `*_genicintergenic.py` | Stratified by genic vs intergenic peaks |
| `*_quartiles.py` | Stratified by expression quartiles |

---

## 📊 Python Scripts — Epigenome (`plots/epigenome/`)

| Script | Description |
|--------|-------------|
| `epigenome_damage_timecourse.py` | Damage-seq over ChIP-seq marks (log₂ sim-normalized; stranded) |
| `epigenome_xr_timecourse.py` | XR-seq over ChIP-seq marks (same logic; all strains) |

---

## 📊 Python Scripts — TSS/TES Meta-Profiles (`plots/TSS/`)

| Script | Description |
|--------|-------------|
| `timecourseTSS.py` | Plots time course TSS metaprofiles by expression quartile (upper/lower) for XR |
| `TSSTESgenomewideTCR.py` | Plots side-by-side TS vs NTS at TSS and TES for each sample |

---

## 📦 Requirements

### R
- `ggplot2`
- `dplyr`
- `ggpubr`
- `pheatmap`
- `purrr`

### Python
- `pandas`, `numpy`
- `matplotlib`
- (Optional: `seaborn`, `scipy`)

---

## 📁 Outputs

All plots are saved as `.pdf` files in their respective `../output/.../` folders with high-resolution, color-safe formatting. Confidence intervals (CI95) are shaded. TS/NTS strand orientation is consistently handled.

---

## Author

Cansu Kose, 2025 – UNC Chapel Hill | Sancar Lab
