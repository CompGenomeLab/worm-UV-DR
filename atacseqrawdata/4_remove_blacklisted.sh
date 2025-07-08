#!/bin/bash

# module load samtools
# module load bedtools

WORKDIR="./"

sample_path="${WORKDIR}/alignment"
sample="" # Basename of the sample to be processed without the suffix

blacklist="" # blacklist regions file 

echo "Processing sample $sample, removing blacklist regions"

#Remove reads within the blacklist regions
bedtools intersect \
    -nonamecheck \
    -v \
    -abam $sample_path/${sample}.filtered_ce11.bam \
    -b $blacklist \
    > $sample_path/${sample}.tmp.bam

#Sort and index the bam file
samtools sort \
    -O bam \
    -o $sample_path/${sample}.blacklist-filtered_ce11.bam \
    $sample_path/${sample}.tmp.bam

samtools index $sample_path/${sample}.blacklist-filtered_ce11.bam

rm $sample_path/${sample}.tmp.bam

echo "done"