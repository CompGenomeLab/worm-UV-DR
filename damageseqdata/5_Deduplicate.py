#!/usr/bin/python

# Required modules: picard, SLURM

import os
from glob import glob

# Loop through all sorted BAM files
for file in glob("./*_samSorted.bam"):
    file = file.strip()
    sample_name = file.split("/")[-1].replace("_samSorted.bam", "")
    
    # Define file paths
    dedup_bam = f"{sample_name}_dedup.bam"
    metrics_file = f"{sample_name}_dedup.metrics.txt"
    
    # SLURM job names and logs
    dedup_job = f"dedup_{sample_name}"
    dedup_out = f"{dedup_job}.out"
    
    # MarkDuplicates command
    dedup_command = (
        f"sbatch -t 2-00:00:00 --mem=80000 --cpus-per-task=4 -J {dedup_job} -o {dedup_out} "
        f"--wrap=\"(picard MarkDuplicates INPUT={file} OUTPUT={dedup_bam} METRICS_FILE={metrics_file} TMP_DIR=./tmp && "
        f"echo 'Duplicate removal success' || (echo 'Duplicate removal failed'; exit 1))\""
    )
    
    # Submit the job
    print(f"Submitting deduplication job for {sample_name}...")
    os.system(dedup_command)
