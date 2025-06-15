#!/bin/bash

module load bedtools
module load r

# Define the input and output directory
DATA_DIR="/pathtoworkingdirectory"

# Change to the data directory to process files
cd "$DATA_DIR" || exit 1

# Loop through each file ending with _filtered_gene_data.bed
for INPUT_FILE in *filtered_gene_data.bed; do
    echo "Processing $INPUT_FILE..."

    # Correct start and end coordinates for negative strand
    awk 'BEGIN {OFS="\t"} {if ($4 == "-") {temp=$2; $2=$3; $3=temp} print}' "$INPUT_FILE" > "${INPUT_FILE%.bed}_corrected.bed"

    # Sort the corrected input file
    sort -k1,1 -k2,2n "${INPUT_FILE%.bed}_corrected.bed" > "${INPUT_FILE%.bed}_sorted.bed"

    # Find overlapping regions excluding self-overlaps
    bedtools intersect -a "${INPUT_FILE%.bed}_sorted.bed" -b "${INPUT_FILE%.bed}_sorted.bed" -wo | awk '$2 != $9 || $3 != $10 {print $0}' > "${INPUT_FILE%.bed}_overlaps_excluding_self.bed"

    # Extract gene names from the overlaps and find unique gene names
    awk '{print $5}' "${INPUT_FILE%.bed}_overlaps_excluding_self.bed" | sort | uniq > "${INPUT_FILE%.bed}_unique_overlapping_genes.txt"

    # Exclude overlapping genes from the all genes list
    grep -v -f "${INPUT_FILE%.bed}_unique_overlapping_genes.txt" "${INPUT_FILE%.bed}_sorted.bed" > "${INPUT_FILE%.bed}_nooverlap.bed"

    # Filter for genes longer than 2kb
    awk '$3 - $2 > 2000' "${INPUT_FILE%.bed}_nooverlap.bed" > "${INPUT_FILE%.bed}_longerthan2kb.bed"
    
    # Select relevant columns
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $6, $7, $4}' "${INPUT_FILE%.bed}_longerthan2kb.bed" > "${INPUT_FILE%.bed}_longerthan2kb_selected_columns.bed"
    
    # Filter genes at least 500 nucleotides away from each other
    bedtools closest -a "${INPUT_FILE%.bed}_longerthan2kb_selected_columns.bed" -b "${INPUT_FILE%.bed}_longerthan2kb_selected_columns.bed" -d -N -t first | awk -v dist="500" '$13 >= dist' | cut -f 1-6 > "${INPUT_FILE%.bed}_500bpaway.bed"

    # Rename final output files for clarity
    BASENAME=$(basename "$INPUT_FILE" _filtered_gene_data.bed)
    mv "${INPUT_FILE%.bed}_500bpaway.bed" "${BASENAME}_500bpaway.bed"

    # Remove intermediate files
    rm "${INPUT_FILE%.bed}_unique_overlapping_genes.txt"
    rm "${INPUT_FILE%.bed}_overlaps_excluding_self.bed"
    rm "${INPUT_FILE%.bed}_longerthan2kb.bed"
    rm "${INPUT_FILE%.bed}_longerthan2kb_selected_columns.bed"
    rm "${INPUT_FILE%.bed}_corrected.bed"
    rm "${INPUT_FILE%.bed}_nooverlap.bed"
    
    echo "Processed $INPUT_FILE"
done

# Pass filtered file to R script and confirm execution
echo "Running R script for quartile analysis..."
Rscript quartile_split.R "$DATA_DIR"
if [ $? -eq 0 ]; then
    echo "R script completed successfully."
else
    echo "R script failed to run. Check for errors."
    exit 1
fi
