#!/bin/bash

#SBATCH --job-name=generate_bw
#SBATCH --output=logs/master_bigwig_%j.out
#SBATCH --error=logs/master_bigwig_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --partition=general

# Load required modules
module load ucsctools
module load bedtools

# Paths
chrom_sizes="/users/c/a/cansuk/seq/ce/ce11.chrom.sizes"
fai_file="/users/c/a/cansuk/seq/ce/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai"
merged_dirs=(
  "/work/users/c/a/cansuk/2025JanuaryCeleganspaperworkflow/damageseqdata/merged_bed/simulation"
  "/work/users/c/a/cansuk/2025JanuaryCeleganspaperworkflow/xrseqdata/merged_bed/simulation"
)
output_dir="/work/users/c/a/cansuk/2025JanuaryCeleganspaperworkflow/IGV_bigwig/simulation"
mkdir -p "$output_dir" logs

# Loop over files
for dir in "${merged_dirs[@]}"; do
  for bed in "$dir"/*simulation.bed; do
    [[ -f "$bed" ]] || continue
    base=$(basename "$bed" .bed)
    unsorted_bed="${output_dir}/${base}_noMtDNA.bed"
    sorted_bed="${output_dir}/${base}_noMtDNA_sorted.bed"

    # Skip if unsorted bed is missing
    [[ -f "$unsorted_bed" ]] || { echo "Missing $unsorted_bed"; continue; }

    # 1. Submit sorting job
    sort_jobid=$(sbatch --parsable --mem=64g --time=12:00:00 --output=logs/${base}_sort.out --wrap="
      bedtools sort -i $unsorted_bed > $sorted_bed
    ")

    # 2A. Plus strand: genomecov from sorted bed
    plus_raw_jobid=$(sbatch --parsable --mem=64g --time=12:00:00 --dependency=afterok:$sort_jobid --output=logs/${base}_plus_raw.out --wrap="
      total=\$(wc -l < $sorted_bed)
      scale=\$(awk -v total=\$total 'BEGIN{printf \"%.8f\", 1000000/total}')
      bedtools genomecov -i $sorted_bed -g $fai_file -bg -strand + -scale \$scale \
        > ${output_dir}/${base}_plus_unsorted.bedGraph
    ")

    # 2B. Plus strand: sort + bigwig
    sbatch --mem=64g --time=12:00:00 --dependency=afterok:$plus_raw_jobid --output=logs/${base}_plus_bw.out --wrap="
      bedtools sort -faidx $chrom_sizes -i ${output_dir}/${base}_plus_unsorted.bedGraph > ${output_dir}/${base}_plus_sorted.bedGraph &&
      bedGraphToBigWig ${output_dir}/${base}_plus_sorted.bedGraph $chrom_sizes ${output_dir}/${base}_plus.bw
    "

    # 3A. Minus strand: genomecov from sorted bed
    minus_raw_jobid=$(sbatch --parsable --mem=64g --time=12:00:00 --dependency=afterok:$sort_jobid --output=logs/${base}_minus_raw.out --wrap="
      total=\$(wc -l < $sorted_bed)
      scale=\$(awk -v total=\$total 'BEGIN{printf \"%.8f\", 1000000/total}')
      bedtools genomecov -i $sorted_bed -g $fai_file -bg -strand - -scale \$scale \
        | awk '{OFS=\"\t\"; \$4 = -\$4; print}' > ${output_dir}/${base}_minus_unsorted.bedGraph
    ")

    # 3B. Minus strand: sort + bigwig
    sbatch --mem=64g --time=12:00:00 --dependency=afterok:$minus_raw_jobid --output=logs/${base}_minus_bw.out --wrap="
      bedtools sort -faidx $chrom_sizes -i ${output_dir}/${base}_minus_unsorted.bedGraph > ${output_dir}/${base}_minus_sorted.bedGraph &&
      bedGraphToBigWig ${output_dir}/${base}_minus_sorted.bedGraph $chrom_sizes ${output_dir}/${base}_minus.bw
    "

  done
done
