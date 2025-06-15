#!/bin/bash
# === SETUP ===
# Load SRA Toolkit module if on cluster (optional)
# module load sra-tools

# Create directories
mkdir -p sra fastq

# Move into sra directory
cd sra

# List of your samples
samples=(
SRR7164190
SRR7164191
SRR7164179
SRR7164178
SRR7164167
SRR7164166
SRR7164155
SRR7164154
)

# Loop through each sample and download
for sample in "${samples[@]}"; do
    echo "[INFO] Downloading $sample"
    prefetch "$sample"
    echo "[INFO] Converting $sample to FASTQ"
    fasterq-dump "$sample" --split-files --outdir ../fastq
done

echo "[DONE] All FASTQs are saved in 'fastq/' directory."
