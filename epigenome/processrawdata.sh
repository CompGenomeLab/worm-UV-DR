#!/bin/bash

# === SETTINGS ===
FASTQ_DIR="./fastq"          # Directory where fastqs are
REF_GENOME_INDEX="/pathto/Bowtie2Index/WBcel235"  # Bowtie2 index prefix
RESULT_DIR="./peaks_chipseq"
mkdir -p "$RESULT_DIR"

# === MODULES ===
echo Loading modules...
module load gcc/11.2.0
module load bowtie2
module load samtools
module load macs/2.2.9.1
module load bedtools
module load deeptools


# === SAMPLE-HISTONE MAP ===
declare -A histone_map
histone_map["SRR7164154"]="H3K4me3"
histone_map["SRR7164155"]="H3K4me3"
histone_map["SRR7164166"]="H3K4me1"
histone_map["SRR7164167"]="H3K4me1"
histone_map["SRR7164178"]="H3K36me3"
histone_map["SRR7164179"]="H3K36me3"
histone_map["SRR7164190"]="H3K27me3"
histone_map["SRR7164191"]="H3K27me3"


# === GROUP SAMPLES BY HISTONE ===
declare -A samples_by_histone
for sample in "${!histone_map[@]}"; do
    histone=${histone_map[$sample]}
    samples_by_histone[$histone]+="${sample} "
done

# === ALIGN EACH SAMPLE ===
for sample in "${!histone_map[@]}"; do
  histone=${histone_map[$sample]}
  sbatch --mem=64g -n 4 -t 1-00:00:00 --job-name="${sample}_bowtie2_align" --output="slurm-%j-${sample}_align.out" --wrap="
    bowtie2 -x ${REF_GENOME_INDEX} -U ${FASTQ_DIR}/${sample}_1.fastq -S ${sample}.sam --very-sensitive;
    samtools view -bS ${sample}.sam | samtools sort -o ${sample}_sorted.bam;
    samtools index ${sample}_sorted.bam;
    rm ${sample}.sam
 "
done

# === REMOVE DUPLICATES ===
for sample in "${!histone_map[@]}"; do
  sbatch --dependency=singleton --mem=64g -t 12:00:00 --job-name="${sample}_rmdup" --output="slurm-%j-${sample}_rmdup.out" --wrap="
    samtools markdup -r ${sample}_sorted.bam ${sample}_sorted_rmdup.bam;
    samtools index ${sample}_sorted_rmdup.bam
  "
done

# === MERGE REPLICATES BY HISTONE MARK ===
for histone in "${!samples_by_histone[@]}"; do
  samples="${samples_by_histone[$histone]}"
  merged_name="${histone}_merged"
  
  sbatch --dependency=singleton --mem=64g -t 12:00:00 --job-name="merge_${histone}" --output="slurm-%j-merge_${histone}.out" --wrap="
    samtools merge ${merged_name}.bam ${samples// /_sorted_rmdup.bam };
    samtools sort -o ${merged_name}_sorted.bam ${merged_name}.bam;
    samtools index ${merged_name}_sorted.bam
  "
done

# === PEAK CALLING ===
for histone in "${!samples_by_histone[@]}"; do
  merged_name="${histone}_merged"

  if [[ "$histone" == "H3K27me3" || "$histone" == "H3K36me3" ]]; then
    # Broad peaks
    sbatch --dependency=singleton --mem=32g -t 6:00:00 --job-name="macs2_${histone}_broad" --output="slurm-%j-${histone}_macs2_broad.out" --wrap="
      macs2 callpeak -t ${merged_name}_sorted.bam -f BAM -g 9e7 -n ${histone} --broad --broad-cutoff 0.1 --outdir ${RESULT_DIR}
    "
  elif [[ "$histone" == "H3K4me1" || "$histone" == "H3K4me3" ]]; then
    # Narrow peaks with manual fragment size
    sbatch --dependency=singleton --mem=32g -t 6:00:00 --job-name="macs2_${histone}_narrow" --output="slurm-%j-${histone}_macs2_narrow.out" --wrap="
      macs2 callpeak -t ${merged_name}_sorted.bam -f BAM -g 9e7 -n ${histone} --nomodel --extsize 147 --outdir ${RESULT_DIR}
    "
  else
    # Fallback (should not occur)
    echo "⚠️ Unknown histone mark: $histone"
  fi
done


# === CREATE BIGWIGS ===
for histone in "${!samples_by_histone[@]}"; do
  merged_name="${histone}_merged"
  
  sbatch --dependency=singleton --mem=64g -t 6:00:00 --job-name="bamCov_${histone}" --output="slurm-%j-${histone}_bamCoverage.out" --wrap="
    bamCoverage -b ${merged_name}_sorted.bam -o ${RESULT_DIR}/${histone}_merged.bw --normalizeUsing RPKM --binSize 10
  "
done

echo "[DONE] All jobs submitted!"
