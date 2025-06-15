#!/usr/bin/python

# Required modules: samtools, SLURM

import os
from glob import glob

# Loop through all BAM files
for file in glob("./*.bam"):
    file = file.strip()
    sample_name = file.split("/")[-1].replace(".bam", "")
    
    # Define file paths
    header_file = f"{sample_name}_header.txt"
    sorted_bam = f"{sample_name}_samSorted.bam"
    
    # SLURM job names and logs
    parse_sort_job = f"parse_sort_{sample_name}"
    parse_sort_out = f"{parse_sort_job}.out"
    
    # Parse and sort command
    parse_sort_command = (
        f"sbatch -t 1-00:00:00 --mem=80000 --cpus-per-task=8 -J {parse_sort_job} -o {parse_sort_out} "
        f"--wrap=\"(samtools view -H {file} | grep -v '^@PG' > {header_file} && echo 'Header parsing success' || "
        f"(echo 'Header parsing failed'; exit 1)) && "
        f"(samtools reheader {header_file} {file} | samtools sort -o {sorted_bam} -@ 8 && echo 'Reheader and sorting success' || "
        f"(echo 'Reheader and sorting failed'; exit 1))\""
    )
    
    # Submit the job
    print(f"Submitting parse and sort job for {sample_name}...")
    os.system(parse_sort_command)
