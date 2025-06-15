#!/bin/bash

# Load bedtools module
module load bedtools

# Input BED file
INPUT_BED="c_elegansPRJNA13758_WS295_annotations.bed"
CHROM_SIZES="/pathto/ce11.chrom.sizes"

# Declare thresholds
declare -i GENE_LONGER_THAN=2000
declare -i DISTANCE_GREATER_THAN=500

# Intermediate and final file names
TRANSFORMED_BED="transformed_genes.bed"
CLOSEST_BED="closest_genes.bed"
FILTERED_BED="filtered_genes.bed"
FINAL_BED="processed_genes_dis500longer2kbselected.bed"
SLOP_BED="${FINAL_BED}_slop500.bed"
TSS_PREFIX="TSSce11_"
TES_PREFIX="TESce11_"

# Step 1: Transform the BED file and exclude MtDNA
echo "Step 1: Transforming the BED file..."
awk 'BEGIN {OFS="\t"} $1 != "MtDNA" {print $1, $2, $3, $6, 0, $4, $5}' "$INPUT_BED" > "$TRANSFORMED_BED"

# Debug: Check the transformed file
echo "Debug: Transformed BED file:"
head "$TRANSFORMED_BED"

# Step 2: Find closest genes using bedtools
echo "Step 2: Finding closest genes..."
bedtools closest -a "$TRANSFORMED_BED" -b "$TRANSFORMED_BED" -d -N -t first > "$CLOSEST_BED"

# Debug: Check the closest file
echo "Debug: Closest BED file:"
head "$CLOSEST_BED"

# Step 3: Filter by distance
echo "Step 3: Filtering by distance..."
awk -v dist="$DISTANCE_GREATER_THAN" '$15 >= dist' "$CLOSEST_BED" > "$FILTERED_BED"

# Debug: Check the filtered file
echo "Debug: Filtered BED file:"
head "$FILTERED_BED"

# Step 4: Filter by gene length
echo "Step 4: Filtering by gene length..."
awk -v gene_length="$GENE_LONGER_THAN" '$3 - $2 > gene_length' "$FILTERED_BED" | cut -f 1-6 > "$FINAL_BED"

# Debug: Check the final file
echo "Debug: Final BED file:"
head "$FINAL_BED"

# Clean up intermediate files (optional)
echo "Cleaning up intermediate files..."
rm -f "$TRANSFORMED_BED" "$CLOSEST_BED" "$FILTERED_BED"

echo "Final BED file created: $FINAL_BED"

# Step 5: Add slop
bedtools slop -i "$FINAL_BED" -g "$CHROM_SIZES" -b 500 > "$SLOP_BED"

echo "500 bp slop added to final BED file: $SLOP_BED"

# Step 6: Generate TSS and TES files
awk -v OFS='\t' -v tss_prefix="$TSS_PREFIX" -v tes_prefix="$TES_PREFIX" '
{
    if ($6 == "+") {
        tss_start = $2
        tss_end = $2 + 1500
        tes_start = $3 - 1500
        tes_end = $3
    } else if ($6 == "-") {
        tss_start = $3 - 1500
        tss_end = $3
        tes_start = $2
        tes_end = $2 + 1500
    }
    print $1, tss_start, tss_end, $4, $5, $6 > tss_prefix"processed_genes.bed"
    print $1, tes_start, tes_end, $4, $5, $6 > tes_prefix"processed_genes.bed"
}' "$SLOP_BED"

# Step 7: Generate 100 windows for TSS and TES
for file in {${TSS_PREFIX},${TES_PREFIX}}*.bed; do
    echo "Processing $file..."
    base_name=$(basename "$file" ".bed")
    output_bed="${base_name}_100bins.bed"

    awk -v OFS='\t' '{for(i=1; i<=100; i++) print $4, $5, $6}' "$file" > temp_columns.txt

    # Step 2: Use bedtools makewindows
    bedtools makewindows -b "$file" -n 100 > temp_windows.bed

    # Step 3: Paste the columns from step 1 to the output of makewindows
    paste temp_windows.bed temp_columns.txt > "$output_bed"

    # Clean up temporary files
    rm temp_columns.txt temp_windows.bed

    echo "Processed 100-bin BED file saved as $output_bed"
done

# Step 8: Reverse windows for negative strand genes and add bin numbers
for input_file in *_100bins.bed; do
    output_file_name="${input_file/_100bins.bed/_final_100bins.bed}"
    echo "Processing $input_file -> $output_file_name..."

    current_gene=""
    current_bins=()
    reverse_order=false

    write_bins() {
        bin_number=1  # Reset bin numbering for each gene
        if [ "$reverse_order" = true ]; then
            # Reverse the order of bins for negative strand genes
            for ((i=${#current_bins[@]}-1; i>=0; i--)); do
                echo -e "${current_bins[i]}\t$bin_number" >> "$output_file_name"
                ((bin_number++))
            done
        else
            for bin in "${current_bins[@]}"; do
                echo -e "$bin\t$bin_number" >> "$output_file_name"
                ((bin_number++))
            done
        fi
    }

    while IFS=$'\t' read -r chrom start end col4 col5 strand; do
        gene="$col4"
        if [[ "$gene" != "$current_gene" ]]; then
            if [[ -n "$current_gene" ]]; then write_bins; fi
            current_bins=()
            reverse_order=false
            [[ "$strand" == "-" ]] && reverse_order=true
            current_gene="$gene"
        fi
        current_bins+=("$chrom\t$start\t$end\t$col4\t$col5\t$strand")
    done < "$input_file"

    if [[ -n "$current_gene" ]]; then write_bins; fi

done

echo "Bin numbering completed."
