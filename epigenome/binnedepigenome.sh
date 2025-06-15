#!/bin/bash
# binnedepigenome.sh
# Script to generate 5kb windows around centered ChIP-seq peaks, binned into 10 bp intervals

# === SETTINGS ===
INPUT_DIR="./peaks_chipseq"
OUTPUT_DIR="./peaks_chipseq"
GENOME_FILE="/pathto/ce11.chrom.sizes"  # Replace with actual genome file
mkdir -p "$OUTPUT_DIR"

# === Define peak files to process ===
PEAK_FILES=(
  H3K4me3_peaks.bed
  H3K4me1_peaks.bed
  H3K27me3_peaks.bed
  H3K36me3_peaks.bed
)

# === Load bedtools module ===
module load bedtools

# === Process each peak file ===
for peak_file in "${PEAK_FILES[@]}"; do
  basename="${peak_file%.*}"
  echo "Processing $peak_file..."

  sort -k1,1 -k2,2n "${INPUT_DIR}/${peak_file}" > "${OUTPUT_DIR}/${basename}_sorted.bed"

  awk 'BEGIN{OFS="\t"} {
    center = int(($2+$3)/2);
    print $1, center, center+1, ($4 != "" ? $4 : "peak" NR)
  }' "${OUTPUT_DIR}/${basename}_sorted.bed" > "${OUTPUT_DIR}/${basename}_centered.bed"

  bedtools slop -i "${OUTPUT_DIR}/${basename}_centered.bed" -g "$GENOME_FILE" -b 2500 > "${OUTPUT_DIR}/${basename}_centered_5kb.bed"

  awk '($3 - $2) >= 4999' "${OUTPUT_DIR}/${basename}_centered_5kb.bed" > "${OUTPUT_DIR}/${basename}_centered_5kb_filtered.bed"

  awk 'BEGIN{OFS="\t"} {
    peak_id=$4;
    for(i=0; i<500; i++) {
      start=$2+(i*10);
      end=start+10;
      print $1, start, end, peak_id, i+1;
    }
  }' "${OUTPUT_DIR}/${basename}_centered_5kb_filtered.bed" > "${OUTPUT_DIR}/${basename}_5kb_500bins.bed"

  echo "[✅ DONE] ${basename}_5kb_500bins.bed"
done
