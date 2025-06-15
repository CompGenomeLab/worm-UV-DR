#!/bin/bash

# Set working directory
cd ./merged_bed || { echo "Directory ./merged_bed not found"; exit 1; }

# Loop over all *_merged.bed files
for FILE in *_merged.bed; do
    BASENAME="${FILE%.bed}"

    # Separate by strand
    awk '$6 == "+" { print }' "$FILE" > "${BASENAME}_plus.bed"
    awk '$6 == "-" { print }' "$FILE" > "${BASENAME}_minus.bed"

    echo "Processed $FILE → ${BASENAME}_plus.bed and ${BASENAME}_minus.bed"
done
