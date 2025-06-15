#!/bin/bash

# Path to the genome index file
GENOME_INDEX="/users/c/a/cansuk/seq/human/genome.fa.fai"  # Replace with your genome index file

# Loop through all _formatted.bed files
for file in *_formatted.bed; do
    # Extract the base name without the suffix
    sample_name=$(basename "$file" _formatted.bed)
    
    # Define output file name
    final_output="${sample_name}_10.bed"
    log_file="process_${sample_name}.out"

    echo "Submitting SLURM job for sample: $sample_name"

    # Submit the job to SLURM
    job_id=$(sbatch --parsable -t 1-00:00:00 --mem=8000 -n 1 -J process_${sample_name} -o ${log_file} \
        --wrap="bedtools flank -i $file -g $GENOME_INDEX -l 6 -r 0 -s | \
        bedtools slop -g $GENOME_INDEX -l 0 -r 4 -s | \
        awk '{ if (\$3 - \$2 == 10) { print } }' > $final_output")

    if [[ $? -eq 0 ]]; then
        echo "Submitted job for $sample_name with Job ID: $job_id"
    else
        echo "Failed to submit job for $sample_name. Check SLURM configuration and logs."
    fi
done

echo "All SLURM jobs submitted."
