## atacseqdata

This folder contains a unified script to process L1-stage ATAC-seq peaks from eLife and bin them into normalized genomic intervals for downstream analysis.

### 📜 Pipeline Summary (`process_atac_data.sh`)

1. **Download & Format Peaks**
   - Downloads eLife peaks (ce10)
   - Extracts enriched peaks and reformats as 4-column BED

2. **LiftOver to ce11**
   - Converts coordinates to ce11 using UCSC `liftOver`

3. **Slop + Quartile Split**
   - Slops ±925 bp → 2kb peaks
   - Sorts by enrichment → splits into quartiles
   - Sorts by coordinate and bins each into 200 windows

4. **Genic vs Intergenic Annotation**
   - Reformats GTF-derived gene annotations
   - Applies 5′ slop (200bp) to genes
   - Labels genic+promoter (1) vs intergenic (2)
   - Splits peaks based on ≥60% overlap

5. **Final Binning**
   - Bins each genic/intergenic category into 200 segments

### ⚙️ Requirements
- `bedtools`
- `liftOver` from UCSC tools
- `awk`, `wget`, `sort`
- Genome size file: `ce11.chrom.sizes`
- Gene annotation BED (from WS295)

### 📦 Outputs
- `ce11_L1atacseqpeaks_2kb_*.bed`: Slopped ATAC peaks
- `*_200bins.bed`: Binned windows for metaprofile plotting
- `ce11_L1atacseqpeaks_genic.bed`, `intergenic.bed`: Split peak files


## Author

Cansu Kose, 2025 – UNC Chapel Hill | Sancar Lab


