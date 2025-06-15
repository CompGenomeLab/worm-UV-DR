
#!/bin/bash

# this scirpt should be run by this command: bash RNAseq_bioinformatics.sh -d ./
# STAR index, annotation file (.gtf), and filter_sam.pl should be in the same directory.

# Command log
echo bash "$0" "$@" > command_log.txt

while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
      -d|--dir)
          SAMPLE_DIR="$2"
          shift 2;;
done

# Listing samples in an array

SAMPLES=()
for i in "${SAMPLE_DIR}"/*_1.fq; do
  i=${i##*/}
  SAMPLES+=("${i%_1.fq}")
done

echo Your samples are:
for SAMPLE in "${SAMPLES[@]}"; do
  echo "${SAMPLE}"
done

read -r -s -p $'Press enter to continue.\n'

# Modules

echo Loading modules

module load samtools
module load gcc
module load star
module load perl
module load featurecounts
module load bedtools

#star alignment
for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-star.out" --wrap="STAR --genomeDir CE_star --readFilesIn ${SAMPLE}_1.fq ${SAMPLE}_2.fq --outFilterIntronMotifs RemoveNoncanonicalUnannotated –runThreadN 2 --outFileNamePrefix ${SAMPLE}_"
done
#filtering
for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-filter.out" --wrap="perl filter_sam.pl ${SAMPLE}_Aligned.out.sam ${SAMPLE}_filtered.sam"
done

#remove unmappedreads and convert sam to bam
for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-samtobam.out" --wrap="samtools view -bS -F 4 ${SAMPLE}_filtered.sam > ${SAMPLE}.bam"
done

#sort bam file
for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-sortbam.out" --wrap="samtools sort -O bam -o ${SAMPLE}_sorted.bam ${SAMPLE}_final.bam"
done

#get the read counts by featurecounts
for SAMPLE in "${SAMPLES[@]}"; do
    sbatch --dependency=singleton --job-name="${SAMPLE}" --output="slurm-%j-${SAMPLE}-featurecounts.out" --wrap="featureCounts -p -t exon -g gene_id -a ../wormbaseWS295/c_elegans.PRJNA13758.WS295.canonical_geneset.gtf -o ${SAMPLE}_count.txt ${SAMPLE}_sorted.bam"
done
