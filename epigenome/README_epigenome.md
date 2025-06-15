## epigenome

This folder contains a full pipeline to process ChIP-seq data for four histone marks (H3K4me3, H3K4me1, H3K27me3, H3K36me3) in *C. elegans* from raw SRA files to peak binning.

---

### 🧾 Pipeline Steps

#### 1. Download Data
- **`downloadrawdata.sh`**: Downloads raw SRA files using `prefetch` and converts them to FASTQ using `fasterq-dump`.

#### 2. Preprocess & Align
- **`processrawdata.sh`**:
  - Aligns reads with Bowtie2
  - Sorts, removes duplicates, merges replicates
  - Calls peaks with MACS2 (broad or narrow depending on histone)
  - Generates normalized BigWig files

#### 3. Generate Binned Peak Regions
- **`binnedepigenome.sh`**:
  - Takes called peak files
  - Centers each peak and expands ±2.5kb (5kb region)
  - Bins each into 500 × 10bp windows

---

### ⚙️ Requirements

- SLURM cluster with:
  - `bowtie2`, `samtools`, `macs2`, `bedtools`, `deeptools`
  - `sra-tools` for download
- Genome:
  - `ce11.chrom.sizes`
  - Bowtie2 index (e.g., `WBcel235`)

---

### 📁 Expected Folder Outputs

- `fastq/`: Paired-end FASTQ files
- `sra/`: Raw `.sra` files
- `peaks_chipseq/`: Sorted BAMs, MACS2 peaks, BigWigs, and 5kb binned BEDs

> Note: Adjust paths (`/pathto/`) inside scripts before running.

---

### 📊 Outputs

- Merged BAMs
- Narrow/broad MACS2 peak BEDs
- RPKM-normalized BigWigs
- `*_5kb_500bins.bed`: 10bp-binned peaks for metaprofiles

## Author

Cansu Kose, 2025 – UNC Chapel Hill | Sancar Lab


