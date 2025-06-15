#!/usr/bin/env python3

# Convert paired SAM files to BAM using samtools (filter for properly paired reads)

import os
from glob import glob

for file in glob("*.sam"):
    sample = os.path.basename(file).split(".")[0]
    bam_output = f"{sample}.bam"
    job_name = f"sambam_{sample}"
    log_file = f"sambam_{sample}.out"

    cmd = (
        f"sbatch -t 1-00:00:00 --mem=80000 -J {job_name} -o {log_file} "
        f"--wrap=\"samtools view -bf 0x2 -o {bam_output} {file}\""
    )
    print(f"Submitting SAM to BAM conversion job for: {sample}")
    os.system(cmd)
