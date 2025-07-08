#!/bin/bash

# module load nucleoatac

WORKDIR="./"

sample_path="${WORKDIR}/macs2"
sample="" # Basename of the sample to be processed without the suffix

nucleoatacpath="${WORKDIR}/nucleoatac"
mkdir -p $nucleoatacpath

echo "Processing sample $sample, finding atac-seq peaks"
nucleoatac run \
    --bed $sample_path/${sample}_peaks.broadPeak \
    --bam $sample_path/${sample}.blacklist-filtered_ce11.bam \
    --fasta ${WORKDIR}/alignment/reference_ce11/ce11.fa \
    --out $nucleoatacpath/${sample}_atac \
    --cores 16 &> $nucleoatacpath/${sample}_nucleoatac.log

echo "done"

gunzip -c $nucleoatacpath/${sample}_atac.nucpos.bed.gz | \
    awk '{{print $1"\\t"$2"\\t"$3"\\t"".""\\t"$5}}' \
    > $nucleoatacpath/${sample}_atac.nucpos_org.bed

awk '{{print $1"\\t"$2"\\t"$3"\\t"".""\\t"$5}}' \
    $sample_path/${sample}_summits.bed \
    > $nucleoatacpath/${sample}_atac.atac_org.bed
