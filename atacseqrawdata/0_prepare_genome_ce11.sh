#!/bin/bash

# module load samtools
# module load bowtie2

WORKDIR="./"

reference_path="${WORKDIR}/alignment/reference_ce11"
mkdir -p ${reference_path}

wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz -P ${reference_path}/.
gunzip ${reference_path}/ce11.fa.gz

bowtie2-build --threads 32 ${reference_path}/ce11.fa ${reference_path}/ce11.genome.build

samtools faidx ${reference_path}/ce11.fa
