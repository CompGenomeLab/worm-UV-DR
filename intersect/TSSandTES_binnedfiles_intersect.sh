#!/bin/bash

# Set the base directories
XRSEQ_DIR="../xrseqdata/merged_bed"
DAMAGESEQ_DIR="../damageseqdata/merged_bed"
OUTPUT_DIR="../output/genomewidetcr"  # Change this to your desired output directory
GENE_LIST_DIR="../wormbaseWS295"
QUARTILEGENELIST_DIR="../RNAseqdata"


# Find all .bed files ending in simulation.bed under XR-seq and Damage-seq directories
XR_SAMPLES=($(find "${XRSEQ_DIR}" -name "*merged.bed" | sed 's/\.bed$//'))
DAMAGE_SAMPLES=($(find "${DAMAGESEQ_DIR}" -name "*merged.bed" | sed 's/\.bed$//'))

# Combine both sample lists
SAMPLES=("${XR_SAMPLES[@]}" "${DAMAGE_SAMPLES[@]}")

# Ensure output directory exists
mkdir -p "${OUTPUT_DIR}"

module load bedtools

# Loop through each sample and submit jobs to Slurm
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_NAME=$(basename "${SAMPLE}")
    
    # Submit intersect jobs for TES
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-div_intersect_tesNTS.out" \
        --wrap="bedtools intersect -c -a ${GENE_LIST_DIR}/TESce11_processed_genes_final_100bins.bed -b ${SAMPLE}.bed -wa -s -F 0.2 > ${OUTPUT_DIR}/${SAMPLE_NAME}_div_tesNTS.bed"
    sbatch --mem=32g -t 120 --dependency=singleton --job-name="${SAMPLE_NAME}" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-div_intersect_tesTS.out" \
        --wrap="bedtools intersect -c -a ${GENE_LIST_DIR}/TESce11_processed_genes_final_100bins.bed -b ${SAMPLE}.bed -wa -S -F 0.2 > ${OUTPUT_DIR}/${SAMPLE_NAME}_div_tesTS.bed"

    # Submit intersect jobs for TSS
    sbatch --mem=32g -t 120 --dependency=singleton --job-name="${SAMPLE_NAME}" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-div_intersect_tssNTS.out" \
        --wrap="bedtools intersect -c -a ${GENE_LIST_DIR}/TSSce11_processed_genes_final_100bins.bed -b ${SAMPLE}.bed -wa -s -F 0.2 > ${OUTPUT_DIR}/${SAMPLE_NAME}_div_tssNTS.bed"
    sbatch --mem=32g -t 120 --dependency=singleton --job-name="${SAMPLE_NAME}" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-div_intersect_tssTS.out" \
        --wrap="bedtools intersect -c -a ${GENE_LIST_DIR}/TSSce11_processed_genes_final_100bins.bed -b ${SAMPLE}.bed -wa -S -F 0.2 > ${OUTPUT_DIR}/${SAMPLE_NAME}_div_tssTS.bed"
   
    # Submit job for read count
    sbatch --mem=32g -t 120 --dependency=singleton --job-name="${SAMPLE_NAME}" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
        --wrap="grep -c \"^\" ${SAMPLE}.bed > ${OUTPUT_DIR}/${SAMPLE_NAME}_readCount.txt"
done
