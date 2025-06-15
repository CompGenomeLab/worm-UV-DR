import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py <input_fasta_file> <output_bed_file>")
    sys.exit(1)

# Assign input and output file names from command line arguments
fasta_file = sys.argv[1]
bed_file = sys.argv[2]

with open(fasta_file, 'r') as fasta, open(bed_file, 'w') as bed:
    for line in fasta:
        if line.startswith('>'):  # Header line
            # Example header: >I:0-10(+)
            parts = line[1:].strip().split(':')  # Remove '>' and split by ':'
            chrom = parts[0]
            coords = parts[1].split('(')[0]  # Split by '(' and take the first part
            start, end = coords.split('-')
            strand = parts[1][-2]  # Assuming strand info is always present as either '+' or '-'
            bed.write(f"{chrom}\t{start}\t{end}\t.\t1\t{strand}\n")
