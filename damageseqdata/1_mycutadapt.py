#!/usr/bin/env python3

import os
from glob import glob

# Collect sample prefixes by matching paired-end reads
samples = []
for filepath in glob("./*.fastq.gz"):
    if "_R1.fastq.gz" in filepath:
        prefix = filepath.replace("_R1.fastq.gz", "")
        samples.append(prefix)
    elif "_R2.fastq.gz" in filepath:
        prefix = filepath.replace("_R2.fastq.gz", "")
        samples.append(prefix)

samples = list(set(samples))

# Submit cutadapt job for each sample
for sp in samples:
    input1 = f"{sp}_R1.fastq.gz"
    input2 = f"{sp}_R2.fastq.gz"
    output1 = f"{sp}_R1_cut.fastq.gz"
    output2 = f"{sp}_R2_cut.fastq.gz"
    job_name = f"cut_{sp}"
    log_file = f"cut_{sp}.out"

    if os.path.exists(input1) and os.path.exists(input2):
        cmd = (
            f"sbatch -t 1-00:00:00 --mem=30000 -n 8 -J {job_name} -o {log_file} "
            f"--wrap=\"cutadapt --discard-trimmed -g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT "
            f"-o {output1} -p {output2} {input1} {input2}\""
        )
        print(f"Submitting cutadapt job for sample: {sp}")
        os.system(cmd)
    else:
        print(f"Missing input files for sample: {sp}. Skipping.")
