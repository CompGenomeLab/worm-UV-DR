#!/bin/bash

# module load samtools
# module load bowtie2

WORKDIR="./"

sample_path="${WORKDIR}/adaptor_trimming"
sample="" # Basename of the sample to be processed without the suffix

output_path="${WORKDIR}/alignment"
bt2idx="${output_path}/reference_ce11/ce11.genome.build"

echo "Processing sample $sample , alignment starts"
bowtie2 \
    -p 16 \
    --local \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    -I 25 \
    -X 700 \
    -x $bt2idx \
    -1 $sample_path/${sample}_fastp_dedup_1.fastq.gz \
    -2 $sample_path/${sample}_fastp_dedup_2.fastq.gz | \
    samtools view -bS - \
    > $output_path/${sample}_ce11.bam

echo "Alignment done, sorting and indexing starts"
samtools sort $output_path/${sample}_ce11.bam -o $output_path/${sample}_sorted_ce11.bam 

echo "Sorting done, indexing starts"
samtools index $output_path/${sample}_sorted_ce11.bam

echo "Done"
