#!/bin/bash

# Load bedtools module
module load bedtools

# Chromosome sizes file
CHROM_SIZES="/pathtochromsizefile"

# TSS and TES file prefixes
TSS_PREFIX="TSS_"
TES_PREFIX="TES_"

# Number of bins for TSS and TES
NUM_BINS=150

# Step 1: Add slop to quartile BED files
for file in *quartile.bed; do
    base_name=$(basename "$file" ".bed")
    slop_file="${base_name}_slop500.bed"

    echo "Adding 500 bp slop to $file -> $slop_file"
    bedtools slop -i "$file" -g "$CHROM_SIZES" -b 500 > "$slop_file"

    # Step 2: Generate TSS and TES files
    tss_file="${TSS_PREFIX}${base_name}.bed"
    tes_file="${TES_PREFIX}${base_name}.bed"

    echo "Generating TSS and TES files for $slop_file -> $tss_file and $tes_file"
    awk -v OFS='\t' -v tss_file="$tss_file" -v tes_file="$tes_file" '
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
      print $1, tss_start, tss_end, $4, $5, $6 > tss_file
      print $1, tes_start, tes_end, $4, $5, $6 > tes_file
    }' "$slop_file"
done

# Step 7: Generate 150 windows for TSS and TES
for file in ${TSS_PREFIX}*.bed ${TES_PREFIX}*.bed; do
    echo "Processing $file..."
    base_name=$(basename "$file" ".bed")
    output_bed="${base_name}_150bins.bed"

    # Step 1: Use bedtools makewindows
    bedtools makewindows -b "$file" -n "$NUM_BINS" > temp_windows.bed

    # Step 2: Repeat metadata columns (4, 5, 6) for each bin
    awk -v OFS='\t' '{for(i=1; i<=150; i++) print $4, $5, $6}' "$file" > temp_columns.txt

    # Step 3: Paste the columns from step 1 to the output of makewindows
    paste temp_windows.bed temp_columns.txt > "$output_bed"

    # Clean up temporary files
    rm temp_columns.txt temp_windows.bed

    echo "Processed 150-bin BED file saved as $output_bed"
done

# Step 8: Reverse windows for negative strand genes and add bin numbers
for input_file in *_150bins.bed; do
    output_file_name="${input_file/_150bins.bed/_final_150bins.bed}"
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
