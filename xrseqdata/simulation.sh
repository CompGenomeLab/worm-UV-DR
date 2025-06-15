#!/bin/bash

# Listing samples in an array

SAMPLES=()
for i in "${SAMPLE_DIR}"/merged_bed/*merged.bed; do
  i=${i##*/}
  SAMPLES+=("${i%.bed}")
done

echo Your samples are:
for SAMPLE in "${SAMPLES[@]}"; do
  echo "${SAMPLE}"
done

read -r -s -p $'Press enter to continue.\n'


# Filter by length

for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --mem=32g --time=24:00:00 --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-sim.out" --wrap="boquila --fasta ${SAMPLE}.fa --bed simulation/${SAMPLE}_simulation.bed --ref /users/c/a/cansuk/seq/ce/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa --seed 1 --sens 20 --kmer 2 --regions regions.ron > simulation/${SAMPLE}_simulation.fa"
done

