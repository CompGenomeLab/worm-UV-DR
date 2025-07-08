#!/bin/bash

# module load macs2

WORKDIR="./"

sample_path="${WORKDIR}/alignment"
sample="" # Basename of the sample to be processed without the suffix

macs2path="${WORKDIR}/macs2"
mkdir -p $macs2path

echo "Processing sample $sample, finding atac-seq peaks"
macs2 callpeak \
    --format BAM \
    --bdg \
    --SPMR \
    --nomodel \
    --shift -37 \
    --extsize 73 \
    -g ce \
    --keep-dup all \
    -q 0.05 \
    -n ${sample} \
    -t $sample_path/${sample}.blacklist-filtered_ce11.bam \
    --outdir $macs2path/ \
    2> $macs2path/${sample}_macs2.log

echo "done"
