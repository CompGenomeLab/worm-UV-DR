## atacseqrawdata

This folder contains a modular pipeline for processing raw ATAC-seq data from FASTQ files to peak calling and nucleosome positioning, using the C. elegans ce11 genome.

### 📜 Pipeline Summary

1. **Prepare Reference Genome (`0_prepare_genome_ce11.sh`)**
   - Downloads the ce11 reference genome (FASTA)
   - Builds Bowtie2 index
   - Indexes the FASTA with samtools

2. **Adaptor Trimming & Deduplication (`1_trim_dedup.sh`)**
   - Uses `fastp` to trim adaptors and remove PCR duplicates from paired-end FASTQ files

3. **Alignment (`2_alignment.sh`)**
   - Aligns trimmed reads to ce11 genome with Bowtie2
   - Sorts and indexes BAM files

4. **Post-Alignment Processing (`3_post_alignment.sh`)**
   - Generates alignment statistics (idxstats, flagstat)
   - Removes mitochondrial reads
   - Adds read group info
   - Marks and quantifies duplicates (Picard)
   - Filters for high-quality, properly paired reads

5. **Blacklist Filtering (`4_remove_blacklisted.sh`)**
   - Removes reads overlapping blacklisted genomic regions (using BED file)
   - Sorts and indexes the filtered BAM

6. **Peak Calling (`5_macs2.sh`)**
   - Calls ATAC-seq peaks using MACS2 on blacklist-filtered BAM files

7. **Nucleosome Positioning (`6_nucleoatac.sh`)**
   - Runs NucleoATAC to infer nucleosome positions and generate nucleosome/ATAC signal tracks

### ⚙️ Requirements

- `bowtie2`
- `samtools`
- `fastp`
- `picard`
- `bedtools`
- `macs2`
- `nucleoatac`
- `awk`, `wget`
- ce11 blacklist BED file

### 📦 Outputs

- `alignment/reference_ce11/ce11.fa`, `.genome.build`: Reference genome and index
- `adaptor_trimming/`: Trimmed and deduplicated FASTQ files
- `alignment/*.bam`, `*.bai`: Aligned, filtered, and indexed BAM files
- `alignment/*.idxstats`, `*.flagstat`, `*.dup_ce11.metrics`: Alignment and duplication stats
- `macs2/*.broadPeak`, `*.summits.bed`: MACS2 peak calls
- `nucleoatac/*_nucpos_org.bed`, `*_atac_org.bed`: Nucleosome and ATAC signal tracks


## Author

Cem Azgari, 2025 – Sabancı University | Adebali Lab 