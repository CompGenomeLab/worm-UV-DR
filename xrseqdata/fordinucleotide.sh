#!/bin/bash

# Listing samples in an array

SAMPLES=()
for i in "${SAMPLE_DIR}"./*_20senssim.bed; do
  i=${i##*/}
  SAMPLES+=("${i%.bed}")
done

echo Your samples are:
for SAMPLE in "${SAMPLES[@]}"; do
  echo "${SAMPLE}"
done

read -r -s -p $'Press enter to continue.\n'


# Modules
echo Loading bedtools
module load bedtools


# Filter by length

for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --mem=32g --time=24:00:00 --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-length_filter.out" --wrap="awk '{if(\$3-\$2>=24 && \$3-\$2<=24){print}}' ${SAMPLE}.bed > dinucleotide/${SAMPLE}_24.bed"
done

# Get FASTA from BED

for SAMPLE in "${SAMPLES[@]}"; do
  sbatch --dependency=singleton  --mem=32g --time=24:00:00 --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-getfasta.out" --wrap="bedtools getfasta -s -fi /users/c/a/cansuk/seq/ce/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa -bed dinucleotide/${SAMPLE}_24.bed -fo dinucleotide/${SAMPLE}_24_bed.fa"
done
