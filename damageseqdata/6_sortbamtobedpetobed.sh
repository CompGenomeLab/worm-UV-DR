#!/bin/bash

# Enable debugging
set -e  # Exit on errors
set -x  # Print commands for debugging

# Collect all deduplicated BAM files
bam_files=(*_dedup.bam)

# Check if any files exist
if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "No *_dedup.bam files found in the directory."
    exit 1
fi

for bam_file in "${bam_files[@]}"; do
    # Extract the sample name (remove "_dedup.bam" suffix)
    sample_name=${bam_file%_dedup.bam}
    
    # Define output files
    sorted_bam="${sample_name}_pe.bam"
    bedpe_file="${sample_name}.bedpe"
    formatted_bedpe_file="${sample_name}_formatted.bed"
    sort_log="sort_${sample_name}.out"
    bedpe_log="bedpe_${sample_name}.out"
    format_log="format_${sample_name}.out"

    # SLURM job for sorting the BAM file by name
    sort_job=$(sbatch --parsable -t 1-00:00:00 --mem=80000 -n 4 -J sort_${sample_name} -o ${sort_log} \
        --wrap="samtools sort -n -o ${sorted_bam} ${bam_file}")
    if [[ $? -eq 0 ]]; then
        echo "Submitted sorting job for ${bam_file} with Job ID: ${sort_job}"
    else
        echo "Failed to submit sorting job for ${bam_file}. Check SLURM setup and logs."
        continue
    fi

    # SLURM job for converting to BEDPE
    bedpe_job=$(sbatch --parsable -t 1-00:00:00 --mem=80000 -n 4 -J bedpe_${sample_name} -o ${bedpe_log} \
        --dependency=afterok:${sort_job} \
        --wrap="bedtools bamtobed -bedpe -mate1 -i ${sorted_bam} > ${bedpe_file}")
    if [[ $? -eq 0 ]]; then
        echo "Submitted BEDPE conversion job for ${sorted_bam} with Job ID: ${bedpe_job}"
    else
        echo "Failed to submit BEDPE conversion job for ${sorted_bam}. Check SLURM setup and logs."
        continue
    fi

    # SLURM job for formatting the BEDPE file
    format_job=$(sbatch --parsable -t 1-00:00:00 --mem=8000 -n 1 -J format_${sample_name} -o ${format_log} \
        --dependency=afterok:${bedpe_job} \
        --wrap="awk '{if (\$9 == \"+\") print \$1 \"\\t\" \$2 \"\\t\" \$6 \"\\t\" \$7 \"\\t\" \$8 \"\\t\" \$9; else if (\$9 == \"-\") print \$1 \"\\t\" \$5 \"\\t\" \$3 \"\\t\" \$7 \"\\t\" \$8 \"\\t\" \$9}' ${bedpe_file} > ${formatted_bedpe_file}")
    if [[ $? -eq 0 ]]; then
        echo "Submitted BEDPE formatting job for ${bedpe_file} with Job ID: ${format_job}"
    else
        echo "Failed to submit BEDPE formatting job for ${bedpe_file}. Check SLURM setup and logs."
        continue
    fi
done

echo "All SLURM jobs submitted."
