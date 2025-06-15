## wormbaseWS295

This folder contains scripts used to process WormBase WS295 GFF2 gene annotation files for *C. elegans*.

### Scripts

- `1_gff2bed.py`: Converts WormBase WS295 `.gff2` annotation files to BED format with chromosome mapping and gene feature extraction.
- `2_process_bed.sh`: 
  - Filters for long genes (default: >2kb) separated by at least 500 bp
  - Generates slopped BED files
  - Creates strand-aware 100-bin TSS and TES regions

### Requirements
- Python 3
- `bedtools`
- Shell environment with `awk`, `bash`

### Inputs
- `c_elegans.PRJNA13758.WS295.annotations.gff2`
- Chromosome size file (`ce11.chrom.sizes`)

### Output
- Final processed gene BEDs and binned TSS/TES BED files



## Author

Cansu Kose, 2025 – UNC Chapel Hill | Sancar Lab


