#!/usr/bin/env python3

# Align paired-end reads to the C. elegans ce11 genome using Bowtie2

import os
from glob import glob

samples = []
for file in glob("*cut.fastq.gz"):
    if "_R1_" in file:
        head = file.split("_R1_")[0]
    elif "_R2_" in file:
        head = file.split("_R2_")[0]
    else:
        continue
    samples.append(head)

samples = list(set(samples))

for sp in samples:
    input1 = f"{sp}_R1_cut.fastq.gz"
    input2 = f"{sp}_R2_cut.fastq.gz"
    output = f"{sp}.sam"
    job_name = f"Map_{sp}"
    log_file = f"Map_{sp}.out"

    cmd = (
        f"sbatch -t 10-00:00:00 --mem=160000 -n 4 -J {job_name} -o {log_file} "
        f"--wrap=\"bowtie2 -q --phred33 --local -p 4 --seed 123 --no-mixed "
        f"-x ce11_bowtie2_index -1 {input1} -2 {input2} -S {output}\""
    )
    print(f"Submitting Bowtie2 mapping job for: {sp}")
    os.system(cmd)
