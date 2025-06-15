# Intersect Scripts

This folder contains a collection of SLURM-compatible `bedtools intersect` scripts designed to compare Damage-seq and XR-seq datasets with various genomic features including ATAC-seq peaks, histone marks, gene quartiles, and genome-wide bins.

Each script performs parallelized intersection analysis and outputs bin-level counts or filtered overlaps.

---

## 📜 Script Overview

| Script | Description |
|--------|-------------|
| `atacseq_intersect.sh` | Intersects Damage-seq and XR-seq `.bed` files with ATAC-seq peaks (all peaks, genic/intergenic, expression quartiles) |
| `epigenome_intersect.sh` | Intersects Damage/XR `.bed` files (and their simulations) with 5kb binned ChIP-seq histone marks (H3K4me3, H3K4me1, H3K27me3, H3K36me3) |
| `genebodyexpressionquartiles_intersect_fordamageseq.sh` | Intersects Damage-seq samples with gene body quartiles in a strand-specific way |
| `genebodyexpressionquartiles_intersect_forXRseq.sh` | Same as above, but for XR-seq samples |
| `genomewide2kbwindows_intersect.sh` | Intersects XR-seq samples with 2kb genome-wide windows and computes RPKM |
| `TSSandTES_binnedfiles_intersect.sh` | Intersects Damage-seq and XR-seq samples with 100-bin TSS and TES windows (strand-specific) |
| `TSSquartiles_binnedfiles_intersect_damageseq.sh` | Intersects Damage-seq with TSS-based expression quartiles (150 bins) |
| `TSSquartiles_binnedfiles_intersect_XRseq.sh` | Intersects XR-seq with TSS-based expression quartiles (150 bins) |

---

## 🧰 Requirements

- SLURM job scheduler
- `bedtools`
- Properly formatted `.bed` files from Damage-seq and XR-seq pipelines
- Pre-processed annotation BEDs for:
  - ATAC-seq (peak bins, quartiles)
  - ChIP-seq 5kb bins
  - TSS/TES regions
  - Gene body quartiles
  - 2kb windows

---

## 📁 Output

Each script writes to a dedicated `../output/` subdirectory. Files typically include:

- `*_readCount.txt`: total lines in the input `.bed` (used for normalization)
- `*_intersect.bed`: output from `bedtools intersect -c`
- SLURM logs: `slurm-*.out`

---

## 🧪 Notes

- Scripts use `-F 0.5` or `-F 0.2` overlap fraction and may use `-s`/`-S` for strand specificity.
- Simulated data intersections are supported in `epigenome_intersect.sh`.

## Author

Cansu Kose, 2025 – UNC Chapel Hill | Sancar Lab

