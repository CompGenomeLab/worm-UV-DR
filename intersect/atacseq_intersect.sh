#!/bin/bash

########## 1) SETUP DIRECTORIES ##########
# Directories containing .bed files
XRSEQ_DIR="../xrseqdata/merged_bed"
DAMAGESEQ_DIR="../damageseqdata/merged_bed"

# Output directories (created if they do not exist)
XRSEQ_OUT_DIR="../output/atacseqintersects/xrseq"
DAMAGESEQ_OUT_DIR="../output/atacseqintersects/damageseq"

mkdir -p "${XRSEQ_OUT_DIR}"
mkdir -p "${DAMAGESEQ_OUT_DIR}"

# Directory containing ATAC-seq peak lists
ATAC_LIST_DIR="../atacseqdata"

########## 2) LOAD MODULES (IF NEEDED) ##########
module load bedtools

########## 3) LOOP OVER XR-SEQ .BED FILES ##########
for SAMPLE_PATH in "${XRSEQ_DIR}"/*merged.bed; do
    
    # If the glob doesn't match anything, skip
    [[ -e "$SAMPLE_PATH" ]] || continue
    
    SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .bed)
    
    # Submit intersect jobs for various ATAC-seq peak categories
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-allpeaks.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_allpeaks.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-promgenic.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_genic_2kb_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_promotergenic.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-intergenic.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_intergenic_2kb_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_intergenic.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-0_25.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_0_25_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_0_25.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-25_50.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_25_50_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_25_50.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-50_75.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_50_75_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_50_75.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-75_100.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_75_100_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_75_100.bed"
    
    # Submit job to get read count (total lines in the .bed file)
    sbatch --mem=32g -t 120 --dependency=singleton --job-name="${SAMPLE_NAME}" \
      --output="${XRSEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
      --wrap="grep -c '^' ${SAMPLE_PATH} > ${XRSEQ_OUT_DIR}/${SAMPLE_NAME}_readCount.txt"
done

########## 4) LOOP OVER DAMAGE-SEQ .BED FILES ##########
for SAMPLE_PATH in "${DAMAGESEQ_DIR}"/*merged.bed; do
    
    # If the glob doesn't match anything, skip
    [[ -e "$SAMPLE_PATH" ]] || continue
    
    SAMPLE_NAME=$(basename "${SAMPLE_PATH}" .bed)
    
    # Submit intersect jobs for various ATAC-seq peak categories
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-allpeaks.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_allpeaks.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-promgenic.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_genic_2kb_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_promotergenic.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-intergenic.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_intergenic_2kb_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_intergenic.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-0_25.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_0_25_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_0_25.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-25_50.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_25_50_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_25_50.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-50_75.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_50_75_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_50_75.bed"
    
    sbatch --mem=32g -t 120 --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-75_100.out" \
      --wrap="bedtools intersect -c \
         -a ${ATAC_LIST_DIR}/ce11_L1atacseqpeaks_2kb_75_100_200bins.bed \
         -b ${SAMPLE_PATH} -wa \
         > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_75_100.bed"
    
    # Submit job to get read count
    sbatch --mem=32g -t 120 --dependency=singleton --job-name="${SAMPLE_NAME}" \
      --output="${DAMAGESEQ_OUT_DIR}/slurm-%j-${SAMPLE_NAME}-read_count.out" \
      --wrap="grep -c '^' ${SAMPLE_PATH} > ${DAMAGESEQ_OUT_DIR}/${SAMPLE_NAME}_readCount.txt"
done
