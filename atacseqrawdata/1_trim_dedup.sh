#!/bin/bash

# module load fastp

WORKDIR="./"

sample_path="" # Path to the directory containing the sample files

# Sample names are expected to be in the format: sample_1.fastq.gz and sample_2.fastq.gz
sample="" # Basename of the sample to be processed without the suffix _1/2.fastq.gz

output_path="${WORKDIR}/adaptor_trimming"
mkdir -p $output_path

echo "Processing sample $sample"    

fastp \
    --thread 16 \
    --in1 $sample_path/${sample}_1.fastq.gz \
    --in2 $sample_path/${sample}_2.fastq.gz \
    --out1 $output_path/${sample}_fastp_dedup_1.fastq.gz \
    --out2 $output_path/${sample}_fastp_dedup_2.fastq.gz \
    --dedup --detect_adapter_for_pe \
    -h $output_path/${sample}_fastp_dedup_trimming.html \
    -R $output_path/${sample}_fastp_dedup_trimming.out

echo "Trimming $sample with removing duplicates is done"
