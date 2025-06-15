#!/usr/bin/env bash

# ==============================================================================
# Unified Script: process_atac_data.sh
#
# Major Steps:
#   1) Download eLife data -> create ce10_L1atacpeaks.bed (cols 1,2,3,5).
#   2) LiftOver -> ce11_L1atacseqpeaks.bed, strip 'chr' prefix.
#   3) Slop by 925 -> ce11_L1atacseqpeaks_2kb.bed, sort by col4, split quartiles,
#      coordinate-sort quartiles, bin into 200 segments.
#   4) Reformat annotation to 6-column BED (strand in col6), excluding MtDNA.
#   5) Create genic+promoter(1) vs intergenic(2) annotation using 5′ slop (200 nt).
#   6) Split ce11_L1atacseqpeaks.bed into genic vs intergenic if ≥60% overlap
#      => ce11_L1atacseqpeaks_genic.bed, ce11_L1atacseqpeaks_intergenic.bed
#   7) Slop each genic/intergenic set by 925, coordinate-sort, then 200-bin them.
# ------------------------------------------------------------------------------
# Requirements:
#   - UCSC liftOver (module load ucsctools)
#   - bedtools
#   - AWK
#   - ce11.chrom.sizes in /users/c/a/cansuk/seq/ce/ce11.chrom.sizes
#   - A gene annotation with columns: (1=chrom, 2=start, 3=end, 4=strand, 5=geneID, 6=geneName, 7=geneType) or similar.
#     We'll rearrange so that col6 is the strand for bedtools -s, and we'll remove any 'MtDNA' entries.
# ==============================================================================

set -euo pipefail  # Bash strict mode (optional but recommended)

# ----------------------------- CONFIGURE PATHS --------------------------------
WORKDIR="./"

ELIFE_URL="https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzczNDQvZWxpZmUtMzczNDQtZmlnMS1kYXRhMS12Mi50eHQ-/elife-37344-fig1-data1-v2.txt?_hash=S6tOfhfttAOIosQvxZIyFwFKgJ5yfYfpnXuvXIwyv%2BQ%3D"
CHAIN_URL="ftp://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"

# Genome sizes (ce11) & original annotation
GENOME_SIZES="/pathto/ce11.chrom.sizes"
ANNOT_BED_ORIG="../wormbaseWS295/c_elegansPRJNA13758_WS295_annotations.bed"

# We will produce a 6-column annotation for strand-aware bedtools calls:
ANNOT_BED_6COL="c_elegansPRJNA13758_WS295_annotations_6col.bed"

# Output file names (main set)
BED_CE10="ce10_L1atacpeaks.bed"
BED_CE11="ce11_L1atacseqpeaks.bed"
BED_CE11_UNMAPPED="unMapped.bed"
BED_CE11_SLOP="ce11_L1atacseqpeaks_2kb.bed"  # after 925-bp slop

# Quartile files
Q1="ce11_L1atacseqpeaks_2kb_0_25.bed"
Q2="ce11_L1atacseqpeaks_2kb_25_50.bed"
Q3="ce11_L1atacseqpeaks_2kb_50_75.bed"
Q4="ce11_L1atacseqpeaks_2kb_75_100.bed"

# Annotation labeling
GENIC_INTERGENIC="c_elegansPRJNA13758_WS295_genicpromoter1_intergenic2.bed"

# Split peaks
PEAKS_CE11="ce11_L1atacseqpeaks.bed"  # reference to the post-liftOver file
OUT_GENIC="ce11_L1atacseqpeaks_genic.bed"
OUT_INTERG="ce11_L1atacseqpeaks_intergenic.bed"

# Genic/Intergenic slop + bin
GENIC_SLOP="ce11_L1atacseqpeaks_genic_2kb.bed"
INTERG_SLOP="ce11_L1atacseqpeaks_intergenic_2kb.bed"
GENIC_200="ce11_L1atacseqpeaks_genic_2kb_200bins.bed"
INTERG_200="ce11_L1atacseqpeaks_intergenic_2kb_200bins.bed"


# ------------------------------------------------------------------------------
# Function: make_200bins
#   Divides each region in a BED into 200 sub-intervals.
#   Columns in the output: chrom, start, end, peakID, bin#
# ------------------------------------------------------------------------------
make_200bins() {
  local infile="$1"
  local outfile="$2"
  echo "  Generating 200 bins for ${infile} -> ${outfile}"

  awk '
  BEGIN { peakCount=0 }
  {
    peakCount++
    chrom = $1
    start = $2
    end   = $3
    # $4 might be enrichment; ignore for binning
    regLen = end - start
    if (regLen <= 0) { next }

    binSize   = int(regLen / 200)
    remainder = regLen % 200
    currentStart = start

    for (i = 1; i <= 200; i++) {
      thisSize = binSize
      if (i <= remainder) {
        thisSize++
      }

      binStart = currentStart
      binEnd   = binStart + thisSize

      # Output: chrom, binStart, binEnd, peakID, binNumber
      print chrom, binStart, binEnd, "peak" peakCount, i

      currentStart = binEnd
    }
  }' OFS="\t" "${infile}" > "${outfile}"
}


# ------------------------------------------------------------------------------
# PART 1: Prepare WORKDIR & Download Data
# ------------------------------------------------------------------------------
mkdir -p "${WORKDIR}"
cd "${WORKDIR}" || { echo "Cannot cd into ${WORKDIR}"; exit 1; }

echo "Downloading eLife text file..."
wget -O elifedata.txt "${ELIFE_URL}"

echo "Creating ${BED_CE10} by extracting columns 1,2,3,5 (excluding col5=0)..."
awk 'NR>1 && $5 != 0 { OFS="\t"; print $1, $2, $3, $5 }' elifedata.txt > "${BED_CE10}"


# ------------------------------------------------------------------------------
# PART 2: Download chain & run liftOver -> ce11
# ------------------------------------------------------------------------------
echo "Downloading chain file..."
wget -O ce10ToCe11.over.chain.gz "${CHAIN_URL}"

echo "Loading UCSC Tools (liftOver)..."
module load ucsctools 2>/dev/null || {
  echo "WARNING: 'module load ucsctools' failed. Ensure 'liftOver' is in PATH."
}

echo "Running liftOver -> ${BED_CE11}..."
liftOver "${BED_CE10}" \
         ce10ToCe11.over.chain.gz \
         "${BED_CE11}" \
         "${BED_CE11_UNMAPPED}"


# ------------------------------------------------------------------------------
# PART 3: Remove 'chr' prefix from the first column
# ------------------------------------------------------------------------------
echo "Removing 'chr' prefix from first column in ${BED_CE11}..."
awk '{ sub(/^chr/, "", $1); print $0 }' OFS="\t" "${BED_CE11}" > tmp_nochr.bed
mv tmp_nochr.bed "${BED_CE11}"


# ------------------------------------------------------------------------------
# PART 4: Slop by 925 -> ce11_L1atacseqpeaks_2kb.bed, then sort by col4 -> quartiles
# ------------------------------------------------------------------------------
echo "Slopping each interval by 925 -> ${BED_CE11_SLOP}..."
bedtools slop \
  -i "${BED_CE11}" \
  -g "${GENOME_SIZES}" \
  -b 925 \
  > "${BED_CE11_SLOP}"

# Sort by enrichment (col4) to split into quartiles
echo "Sorting ${BED_CE11_SLOP} by enrichment (col4 ascending)..."
sort -k4,4n "${BED_CE11_SLOP}" > tmp_sorted.bed

echo "Overwriting slop file with col4-sorted version..."
mv tmp_sorted.bed "${BED_CE11_SLOP}"

# Split into quartiles
TOTAL=$(wc -l < "${BED_CE11_SLOP}")
Q=$((TOTAL / 4))

head -n "${Q}"                "${BED_CE11_SLOP}" > "${Q1}"
sed -n "$((Q+1)),$((2*Q))p"   "${BED_CE11_SLOP}" > "${Q2}"
sed -n "$((2*Q+1)),$((3*Q))p" "${BED_CE11_SLOP}" > "${Q3}"
sed -n "$((3*Q+1)),$((TOTAL))p" "${BED_CE11_SLOP}" > "${Q4}"

# ------------------------------------------------------------------------------
# Coordinate-sort each quartile by chrom/start (rather than by col4)
# ------------------------------------------------------------------------------
echo "Coordinate-sorting each quartile file..."
for QFILE in "${Q1}" "${Q2}" "${Q3}" "${Q4}"
do
  echo "  Sorting $QFILE by coordinate..."
  sort -k1,1 -k2,2n "$QFILE" > tmp_q.bed
  mv tmp_q.bed "$QFILE"
done

# ------------------------------------------------------------------------------
# (Optional) Coordinate-sort the entire file for binning
# ------------------------------------------------------------------------------
echo "Coordinate-sorting the entire ${BED_CE11_SLOP} file..."
sort -k1,1 -k2,2n "${BED_CE11_SLOP}" > tmp_coord.bed
mv tmp_coord.bed "${BED_CE11_SLOP}"

# ------------------------------------------------------------------------------
# Generate 200 bins for the main 2kb file + each quartile
# ------------------------------------------------------------------------------
echo "Creating 200-bin files from coordinate-sorted intervals..."

make_200bins "${BED_CE11_SLOP}" "ce11_L1atacseqpeaks_2kb_200bins.bed"
make_200bins "${Q1}" "ce11_L1atacseqpeaks_2kb_0_25_200bins.bed"
make_200bins "${Q2}" "ce11_L1atacseqpeaks_2kb_25_50_200bins.bed"
make_200bins "${Q3}" "ce11_L1atacseqpeaks_2kb_50_75_200bins.bed"
make_200bins "${Q4}" "ce11_L1atacseqpeaks_2kb_75_100_200bins.bed"



# ------------------------------------------------------------------------------
# PART 4b: Reformat original annotation to 6 columns (strand in col6), remove MtDNA
# ------------------------------------------------------------------------------
echo "Reformatting original annotation to 6-col BED with strand in col6, excluding MtDNA..."

# Original (example) columns: (1=chrom, 2=start, 3=end, 4=strand, 5=geneID, 6=geneName, 7=geneType)
# We want (1=chrom, 2=start, 3=end, 4=name, 5=score=0, 6=strand),
# and also skip lines where $1 == "MtDNA".

awk '{OFS="\t"; print $1, $2, $3, $5, 0, $4}' "${ANNOT_BED_ORIG}" \
  | awk '$1 != "MtDNA"' \
  > "${ANNOT_BED_6COL}"


# ------------------------------------------------------------------------------
# PART 5: Create genic+promoter(1) vs intergenic(2) annotation
# ------------------------------------------------------------------------------
echo "Creating genic+promoter (1) and intergenic (2) annotation using 5′ slop..."

echo "Slopping each gene 200 nt on the 5′ side..."
bedtools slop \
  -i "${ANNOT_BED_6COL}" \
  -g "${GENOME_SIZES}" \
  -l 200 -r 0 -s \
  > genes_promoter_5prime_200.bed

echo "Sorting slopped file for complement..."
sort -k1,1 -k2,2n genes_promoter_5prime_200.bed > tmp_sorted.bed
mv tmp_sorted.bed genes_promoter_5prime_200.bed

echo "Computing complement (intergenic) regions..."
bedtools complement \
  -i genes_promoter_5prime_200.bed \
  -g "${GENOME_SIZES}" \
  > intergenic.bed

# Step C: Label each set
awk '{OFS="\t"; print $1, $2, $3, 1}' genes_promoter_5prime_200.bed \
  > genic_1.bed
awk '{OFS="\t"; print $1, $2, $3, 2}' intergenic.bed \
  > intergenic_2.bed

# Combine + sort => final annotation
cat genic_1.bed intergenic_2.bed \
  | sort -k1,1 -k2,2n \
  > "${GENIC_INTERGENIC}"

# Cleanup intermediate files
rm -f genic_1.bed intergenic_2.bed genes_promoter_5prime_200.bed intergenic.bed

echo "Annotation file created: ${GENIC_INTERGENIC} (1=genic+promoter, 2=intergenic)"


# ------------------------------------------------------------------------------
# PART 6: Split ce11_L1atacseqpeaks.bed => genic vs intergenic (≥60% overlap)
# ------------------------------------------------------------------------------
echo "Splitting ${PEAKS_CE11} into genic vs intergenic with -f 0.6..."

bedtools intersect \
  -a "${PEAKS_CE11}" \
  -b "${GENIC_INTERGENIC}" \
  -f 0.60 \
  -wa -wb \
  > tmp_peaks_annot.bed

awk '$8 == 1 { OFS="\t"; print $1, $2, $3, $4 }' tmp_peaks_annot.bed \
  > "${OUT_GENIC}"

awk '$8 == 2 { OFS="\t"; print $1, $2, $3, $4 }' tmp_peaks_annot.bed \
  > "${OUT_INTERG}"

rm -f tmp_peaks_annot.bed

echo "Created:"
echo "  - ${OUT_GENIC}      (peaks with >=60% overlap to label=1 => genic+promoter)"
echo "  - ${OUT_INTERG}     (peaks with >=60% overlap to label=2 => intergenic)"


# ------------------------------------------------------------------------------
# PART 7: Slop Genic/Intergenic peaks by 925, sort, bin
# ------------------------------------------------------------------------------
echo "Slopping genic/intergenic peak sets by 925..."

bedtools slop -i "${OUT_GENIC}"  -g "${GENOME_SIZES}" -b 925 > "${GENIC_SLOP}"
bedtools slop -i "${OUT_INTERG}" -g "${GENOME_SIZES}" -b 925 > "${INTERG_SLOP}"

echo "Sorting slopped genic/intergenic..."
sort -k1,1 -k2,2n "${GENIC_SLOP}"   > tmp_sorted.bed && mv tmp_sorted.bed "${GENIC_SLOP}"
sort -k1,1 -k2,2n "${INTERG_SLOP}" > tmp_sorted.bed && mv tmp_sorted.bed "${INTERG_SLOP}"

echo "Creating 200-bin files for genic/intergenic sets..."
make_200bins "${GENIC_SLOP}"  "${GENIC_200}"
make_200bins "${INTERG_SLOP}" "${INTERG_200}"

# ------------------------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------------------------
echo
echo "All steps complete!"
echo "Main peak sets generated:"
echo "  - ${BED_CE10} (ce10 coords)"
echo "  - ${BED_CE11} (liftOver -> ce11, no 'chr')"
echo "  - ${BED_CE11_SLOP} (slopped by 925, sorted by col4 => quartiles)"
echo "Quartiles: ${Q1}, ${Q2}, ${Q3}, ${Q4}"
echo "200-bin files for quartiles & main 2kb: ce11_L1atacseqpeaks_2kb_*"
echo
echo "Annotation reformat => ${ANNOT_BED_6COL} (strand in col6, MtDNA removed)"
echo "Promoter vs intergenic annotation => ${GENIC_INTERGENIC} (1=genic+promoter, 2=intergenic)"
echo "Genic/Intergenic peak split:"
echo "  - ${OUT_GENIC}"
echo "  - ${OUT_INTERG}"
echo "Slop & bin of genic/intergenic:"
echo "  - ${GENIC_SLOP}, ${INTERG_SLOP}"
echo "  - ${GENIC_200}, ${INTERG_200}"
echo "Done!"     in thsi scirpt when it do # Overwrite main slop file with sorted version
mv tmp_sorted.bed "${BED_CE11_SLOP}"  it does not double is not it? it deletes the bed ce11 slop before this