#!/bin/bash

# module load samtools
# module load picard

WORKDIR="./"

sample_path="${WORKDIR}/alignment"
sample="" # Basename of the sample to be processed without the suffix

echo "Processing sample $sample, post alignment processing starts"

#Generate the idxstats report
samtools idxstats $sample_path/${sample}_sorted_ce11.bam \
    > $sample_path/${sample}_sorted_ce11.idxstats

#Generate the flagstat report
samtools flagstat $sample_path/${sample}_sorted_ce11.bam \
    > $sample_path/${sample}_sorted_ce11.flagstat

#Remove reads aligned to the mitochondria
samtools view -h $sample_path/${sample}_sorted_ce11.bam | \
    grep -v chrM | \
    samtools sort -O bam -o $sample_path/${sample}_sorted_ce11.rmChrM.bam -T .

# echo "reheading the bam file"
samtools addreplacerg \
    -@ 16 \
    -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" \
    -o $sample_path/${sample}_sorted_ce11.rmChrM.reheaded.bam \
    $sample_path/${sample}_sorted_ce11.rmChrM.bam

echo "marking duplicates"
picard MarkDuplicates \
    QUIET=true \
    INPUT=$sample_path/${sample}_sorted_ce11.rmChrM.reheaded.bam \
    OUTPUT=$sample_path/${sample}.marked_ce11.bam \
    METRICS_FILE=$sample_path/${sample}.dup_ce11.metrics \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR=.

# View the % of duplicates
head -n 8 $sample_path/${sample}.dup_ce11.metrics | \
    cut -f 7,9 | \
    grep -v ^# | \
    tail -n 2 \
    > $sample_path/${sample}.dup_ce11.txt

echo "removing duplicates and filtering"
samtools view -h -b -f 2 -F 1548 -q 30 $sample_path/${sample}.marked_ce11.bam | \
samtools sort -o $sample_path/${sample}.filtered_ce11.bam

samtools index $sample_path/${sample}.filtered_ce11.bam

echo "removing unnecessary bamfiles"

rm -r $sample_path/${sample}_sorted_ce11.rmChrM.bam
rm -r $sample_path/${sample}.marked_ce11.bam
rm -r $sample_path/${sample}.marked_ce11.bam.bai
rm -r $sample_path/${sample}_ce11.bam
rm -r $sample_path/${sample}_ce11.bam.bai
rm -r $sample_path/${sample}_sorted_ce11.rmChrM.reheaded.bam

echo "done"