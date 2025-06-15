#!/bin/bash

########## 1) SETUP DIRECTORIES ##########
# Input BED file directories
XRSEQ_DIR="../xrseqdata/merged_bed"
DAMAGESEQ_DIR="../damageseqdata/merged_bed"
XRSEQ_SIMULATION_DIR="../xrseqdata/merged_bed/simulation"
DAMAGESEQ_SIMULATION_DIR="../damageseqdata/merged_bed/simulation"
# Output directories
XRSEQ_OUT_DIR="../output/epigenome/intersect/xrseq/5kb/stranded"
DAMAGESEQ_OUT_DIR="../output/epigenome/intersect/damageseq/5kb/stranded"
XRSEQ_SIMULAITON_OUT_DIR="../output/epigenome/intersect/xrseq/5kb/stranded/simulation"
DAMAGESEQ_SIMULATION_OUT_DIR="../output/epigenome/intersect/damageseq/5kb/stranded/simulation"

mkdir -p "${XRSEQ_OUT_DIR}" "${DAMAGESEQ_OUT_DIR}" "XRSEQ_SIMULAITON_OUT_DIR" "DAMAGESEQ_SIMULATION_OUT_DIR"

# ATAC/chIP-seq epigenome bins
EPIGENOME_LIST_DIR="/work/users/c/a/cansuk/2025JanuaryCeleganspaperworkflow/epigenome/sra/peaks_chipseq"

# Histone marks to loop over
MARKS=("H3K4me3" "H3K4me1" "H3K27me3" "H3K36me3")

########## 2) LOAD MODULES ##########
module load bedtools



########## 4) LOOP OVER XR-SEQ FILES ##########
find "${XRSEQ_DIR}" -maxdepth 1 -type f -name "*us.bed" | while read -r SAMPLE_PATH; do
    [[ -e "$SAMPLE_PATH" ]] || continue
    SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .bed)

    # Read count
    sbatch --parsable --mem=32g -t 120 --job-name="${SAMPLE_NAME}_readcount" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
      --wrap="grep -c '^' ${SAMPLE_PATH} > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_readCount.txt"

    # Intersect with each histone mark
    for MARK in "${MARKS[@]}"; do
        sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}_${MARK}" \
          --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-${MARK}.out" \
          --wrap="bedtools intersect -c \
            -a ${EPIGENOME_LIST_DIR}/${MARK}_sorted_peaks_5kb_500bins.bed \
            -b ${SAMPLE_PATH} -wa \
            > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_${MARK}_5kb.bed"
    done

done
########## 4) LOOP OVER DAMAGE-SEQ FILES ##########
find "${DAMAGESEQ_DIR}" -maxdepth 1 -type f -name "*us.bed" | while read -r SAMPLE_PATH; do
    [[ -e "$SAMPLE_PATH" ]] || continue
    SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .bed)

    # Read count
    sbatch --parsable --mem=32g -t 120 --job-name="${SAMPLE_NAME}_readcount" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
      --wrap="grep -c '^' ${SAMPLE_PATH} > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_readCount.txt"

    # Intersect with each histone mark
    for MARK in "${MARKS[@]}"; do
        sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}_${MARK}" \
          --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-${MARK}.out" \
          --wrap="bedtools intersect -c \
            -a ${EPIGENOME_LIST_DIR}/${MARK}_sorted_peaks_5kb_500bins.bed \
            -b ${SAMPLE_PATH} -wa \
            > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_${MARK}_5kb.bed"
    done

done
########## 4) LOOP OVER XRSEQ_SIMULATION FILES ##########
find "${XRSEQ_SIMULAITON_DIR}" -maxdepth 1 -type f -name "*us.bed" | while read -r SAMPLE_PATH; do
    [[ -e "$SAMPLE_PATH" ]] || continue
    SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .bed)

    # Read count
    sbatch --parsable --mem=32g -t 120 --job-name="${SAMPLE_NAME}_readcount" \
      --output="${XRSEQ_SIMULAITON_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
      --wrap="grep -c '^' ${SAMPLE_PATH} > ${XRSEQ_SIMULAITON_OUT_DIR}/${SAMPLE_NAME}_readCount.txt"

    # Intersect with each histone mark
    for MARK in "${MARKS[@]}"; do
        sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}_${MARK}" \
          --output="${XRSEQ_SIMULAITON_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-${MARK}.out" \
          --wrap="bedtools intersect -c \
            -a ${EPIGENOME_LIST_DIR}/${MARK}_sorted_peaks_5kb_500bins.bed \
            -b ${SAMPLE_PATH} -wa \
            > ${XRSEQ__SIMULAITONOUT_DIR}/${SAMPLE_NAME}_${MARK}_5kb.bed"
    done

done
########## 4) LOOP OVER DAMAGE-SEQ_SIMULATION FILES ##########
find "${DAMAGESEQ_SIMULAITON_DIR}" -maxdepth 1 -type f -name "*us.bed" | while read -r SAMPLE_PATH; do
    [[ -e "$SAMPLE_PATH" ]] || continue
    SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .bed)

    # Read count
    sbatch --parsable --mem=32g -t 120 --job-name="${SAMPLE_NAME}_readcount" \
      --output="${DAMAGESEQ_SIMULAITON_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
      --wrap="grep -c '^' ${SAMPLE_PATH} > ${DAMAGESEQ_SIMULAITON_OUT_DIR}/${SAMPLE_NAME}_readCount.txt"

    # Intersect with each histone mark
    for MARK in "${MARKS[@]}"; do
        sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}_${MARK}" \
          --output="${DAMAGESEQ_SIMULAITON_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-${MARK}.out" \
          --wrap="bedtools intersect -c \
            -a ${EPIGENOME_LIST_DIR}/${MARK}_sorted_peaks_5kb_500bins.bed \
            -b ${SAMPLE_PATH} -wa \
            > ${DAMAGESEQ_SIMULAITON_OUT_DIR}/${SAMPLE_NAME}_${MARK}_5kb.bed"
    done

done
