# C. elegans DNA Repair and Damage Analysis Pipeline

This repository provides a complete and modular pipeline for analyzing DNA damage and repair in *Caenorhabditis elegans*, including Damage-seq, XR-seq, ATAC-seq, ChIP-seq, RNA-seq, and annotation processing. Each module is stored in a separate folder with its own `README.md` explaining scripts, dependencies, and outputs.

---

## 📁 Repository Structure

| Folder             | Description | Documentation |
|--------------------|-------------|----------------|
| `wormbaseWS295/`   | Processes WormBase WS295 GFF2 annotations into BED format | [README_annotationprocess.md](wormbaseWS295/README_annotationprocess.md) |
| `xrseqdata/`       | Complete XR-seq processing and normalization pipeline | [README_XRseq.md](xrseqdata/README_XRseq.md) |
| `damageseqdata/`   | Full Damage-seq processing including adapter trimming, filtering, and QC | [README_damage_seq_pipeline.md](damageseqdata/README_damage_seq_pipeline.md) |
| `atacseqdata/`     | ATAC-seq peak reprocessing and binning by genic/intergenic regions or expression | [README_atacseq.md](atacseqdata/README_atacseq.md) |
| `epigenome/`       | ChIP-seq raw processing, peak calling, and binned region generation | [README_epigenome.md](epigenome/README_epigenome.md) |
| `RNAseqdata/`      | RNA-seq processing, TPM filtering, expression quartile assignment | [README_rnaseq.md](RNAseqdata/README_rnaseq.md) |
| `intersect/`       | `bedtools intersect` scripts for combining Damage/XR-seq data with genomic features | [README_intersect.md](intersect/README_intersect.md) |
| `plots/`           | Final visualization scripts (R and Python) for strand bias, TSS metaprofiles, clustering | [README_plots.md](plots/README_plots.md) |

---

## ⚙️ Dependencies

- **Core tools:** `bedtools`, `samtools`, `bowtie2`, `cutadapt`, `macs2`, `STAR`, `deeptools`, `sra-tools`
- **R packages:** `ggplot2`, `dplyr`, `ggpubr`, `pheatmap`, `tidyr`, `wesanderson`, `ggthemes`, `purrr`
- **Python packages:** `pandas`, `numpy`, `matplotlib`, `scipy`, `seaborn` (optional)

Tested on a SLURM-based HPC environment using the `WBcel235` (*C. elegans* ce11) genome.

---

## 📤 Upload Instructions

```bash
cd Celegans_DNARepair_Pipeline/
git init
git add .
git commit -m "Initial commit: modular DNA repair and damage pipeline"
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME.git
git branch -M main
git push -u origin main
```

---

## 👩‍🔬 Author

Cansu Kose, 2025 – UNC Chapel Hill | Sancar Lab

If you use or adapt this pipeline, please cite appropriately or contact me via GitHub.
