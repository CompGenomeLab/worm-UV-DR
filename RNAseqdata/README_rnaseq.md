## RNAseqdata

This folder contains scripts to process *C. elegans* RNA-seq data from raw FASTQ files through gene-level filtering and quartile-based categorization.

### 📜 Pipeline Overview

1. **Alignment & Counting**
   - `RNA1_rnaseqbioinformatics.sh`: SLURM pipeline for STAR alignment, SAM filtering, and gene-level counting via FeatureCounts.
   - `filter_sam.pl`: Perl script used to filter SAM files before converting to BAM.

2. **Expression-Based Filtering**
   - `RNA2_filternotexpressed.R`: Computes TPM and filters out lowly expressed genes.

3. **Gene Filtering by Length/Distance**
   - `RNA3_filtergenesinquartiles.sh`: Removes overlapping/short genes, retains those ≥2kb and ≥500bp apart.

4. **Quartile Assignment**
   - `quartile_split.R`: Splits retained genes into 4 expression quartiles.

5. **TSS/TES Binning**
   - `RNA4_process_quartilebed.sh`: Adds slop, generates 150-bin TSS and TES windows per quartile.

### 📁 Inputs
- Paired-end FASTQ files (`*_1.fq`, `*_2.fq`)
- GTF: `c_elegans.PRJNA13758.WS295.canonical_geneset.gtf`
- BED: WS295 gene annotation in BED format

### 💻 Requirements
- `STAR`, `samtools`, `bedtools`, `featureCounts`, `R`, `data.table`, `perl`
- Executed in SLURM-managed HPC environment

### ⚠️ Notes
- Update paths in each script (e.g., chromosome size files, working directory)
- Run `RNA1_rnaseqbioinformatics.sh` first before subsequent steps


## Author

Cansu Kose, 2025 – UNC Chapel Hill | Sancar Lab


