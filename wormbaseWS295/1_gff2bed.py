#!/usr/bin/env python3

# Define chromosome mapping from GFF2 to BED format
chromosome_mapping = {
    "CHROMOSOME_I": "I",
    "CHROMOSOME_II": "II",
    "CHROMOSOME_III": "III",
    "CHROMOSOME_IV": "IV",
    "CHROMOSOME_V": "V",
    "CHROMOSOME_X": "X",
    "CHROMOSOME_MtDNA": "MtDNA"
}

# Define input and output file paths
input_gff2_file = 'c_elegans.PRJNA13758.WS295.annotations.gff2'
output_bed_file = 'c_elegansPRJNA13758_WS295_annotations.bed'

# Open input GFF2 file and output BED file
with open(input_gff2_file, 'r') as infile, open(output_bed_file, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):  # Skip header and comment lines
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:  # Ensure there are enough columns
            continue
        attributes = fields[8]
        if not attributes.startswith("Gene "):  # Skip lines where the 9th column doesn't start with 'Gene'
            continue
        
        # Map chromosome name to BED format
        chrom = chromosome_mapping.get(fields[0], fields[0])  # Default to original if no mapping found
        start = int(fields[3]) - 1  # Convert 1-based start position to 0-based (BED format)
        end = fields[4]
        strand = fields[6]
        
        # Extract Gene ID, Locus, Sequence Name, and Biotype
        gene_id = ""
        locus = ""
        sequence_name = ""
        biotype = ""
        for attribute in attributes.split(';'):
            attribute = attribute.strip()
            if attribute.startswith("Gene "):
                gene_id = attribute.split('"')[1]  # Extract value inside quotes
            elif attribute.startswith("Locus "):
                locus = attribute.split('"')[1]  # Extract value inside quotes
            elif attribute.startswith("Sequence_name "):
                sequence_name = attribute.split('"')[1]  # Extract value inside quotes
            elif attribute.startswith("Biotype "):
                biotype = attribute.split('"')[1]  # Extract value inside quotes
        
        # Use Locus if available, otherwise fall back to Sequence Name
        locus_or_sequence = locus if locus else sequence_name

        # Write to BED file
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            chrom, start, end, strand, gene_id, locus_or_sequence, biotype
        ))
