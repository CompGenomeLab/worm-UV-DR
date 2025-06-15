#!/bin/bash

# Set the base directories
XRSEQ_DIR="../damageseqdata/merged_bed"
OUTPUT_DIR="../output/genebodyintersect/damage"  # Change this to your desired output directory
GENELIST_DIR="../RNAseqdata"

# Ensure output directory exists
mkdir -p "${OUTPUT_DIR}"

# Load the required module
module load bedtools

# Find all .bed files in the XRSEQ_DIR
SAMPLES=($(find "${XRSEQ_DIR}" -name "*merged.bed"))

# Loop through each sample and submit jobs to Slurm
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_NAME=$(basename "${SAMPLE}" .bed)  # Extract sample name without extension

    # Extract the strain name (wt, csb, or xpc)
    STRAIN=$(echo "${SAMPLE_NAME}" | cut -d '_' -f 1)

    echo "Processing sample: ${SAMPLE_NAME} (Strain: ${STRAIN})"
    
    # Submit job for read count (unique for each sample)
    JOB_READCOUNT=$(sbatch --parsable --mem=32g -t 120 --job-name="${SAMPLE_NAME}_readcount" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
        --wrap="grep -c \"^\" ${SAMPLE} > ${OUTPUT_DIR}/${SAMPLE_NAME}_readCount.txt")

    # Define quartile files for the strain
    QUARTILE_FILES=(
        "${GENELIST_DIR}/wt_0_25_quartile.bed"
        "${GENELIST_DIR}/wt_25_50_quartile.bed"
        "${GENELIST_DIR}/wt_50_75_quartile.bed"
        "${GENELIST_DIR}/wt_75_100_quartile.bed"
    )

    # Loop through quartile files and submit jobs
    for QUARTILE_FILE in "${QUARTILE_FILES[@]}"; do
        if [[ ! -f "${QUARTILE_FILE}" ]]; then
            echo "Quartile file not found: ${QUARTILE_FILE}"
            continue
        fi

        QUARTILE_NAME=$(basename "${QUARTILE_FILE}" .bed)  # Extract the quartile file name without extension

        # Submit intersect jobs for TES (NTS and TS)
        JOB_NTS=$(sbatch --parsable --mem=32g -t 120 --dependency=afterok:${JOB_READCOUNT} --job-name="${SAMPLE_NAME}_${QUARTILE_NAME}_NTS" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-${QUARTILE_NAME}-genesNTS.out" \
            --wrap="bedtools intersect -c -a ${QUARTILE_FILE} -b ${SAMPLE} -wa -s -F 0.5 > ${OUTPUT_DIR}/${SAMPLE_NAME}_${QUARTILE_NAME}_NTS.bed")
        
        JOB_TS=$(sbatch --parsable --mem=32g -t 120 --dependency=afterok:${JOB_READCOUNT} --job-name="${SAMPLE_NAME}_${QUARTILE_NAME}_TS" --output="${OUTPUT_DIR}/slurm-%j-${SAMPLE_NAME}-${QUARTILE_NAME}-genesTS.out" \
            --wrap="bedtools intersect -c -a ${QUARTILE_FILE} -b ${SAMPLE} -wa -S -F 0.5 > ${OUTPUT_DIR}/${SAMPLE_NAME}_${QUARTILE_NAME}_TS.bed")
    done
done