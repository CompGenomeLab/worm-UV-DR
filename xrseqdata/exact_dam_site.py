import os
import glob
import logging
import argparse

from datetime import datetime

parser = argparse.ArgumentParser(description='Map repair and damage signals on given region file.')
parser.add_argument('--input', type=str, help='input folder')
parser.add_argument('--output', type=str, help='output folder')
args = parser.parse_args()

# Configure logging
logging.basicConfig(filename='damage_site_processing.log', level=logging.INFO, 
                    format='%(asctime)s - %(message)s', datefmt='%a, %d %b %Y %H:%M:%S %z')

def process_file(input_file, output_file, strand, seq_type):
    """Process a single file to extract exact damage sites."""
    try:
        logging.info(f"Extracting exact damage sites ({strand} strand) for {seq_type}...")
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                other_data = '\t'.join(parts[3:])

                if strand == 'plus':
                    if seq_type == 'xr':
                        new_start = end - 8
                        new_end = end - 6
                    elif seq_type == 'ds':
                        new_start = end - 6
                        new_end = end - 4
                elif strand == 'minus':
                    if seq_type == 'xr':
                        new_start = start + 6
                        new_end = start + 8
                    elif seq_type == 'ds':
                        new_start = start + 4
                        new_end = start + 6

                outfile.write(f"{chrom}\t{new_start}\t{new_end}\t{other_data}\n")

        logging.info("Success!")
    except Exception as e:
        logging.error(f"Process failed: {e}")
        if os.path.exists(output_file):
            os.remove(output_file)
        raise

def main(input_dir, output_dir):
    """Process all files in the input directory."""
    # Define file patterns
    patterns = {
        'xr': {
            'plus': ['*_xr_merged_plus.bed', '*_xr_merged_sim_plus.bed'],
            'minus': ['*_xr_merged_minus.bed', '*_xr_merged_sim_minus.bed']
        },
        'ds': {
            'plus': ['*_ds_merged_plus.bed', '*_ds_merged_sim_plus.bed'],
            'minus': ['*_ds_merged_minus.bed', '*_ds_merged_sim_minus.bed']
        }
    }

    for seq_type, strands in patterns.items():
        for strand, pattern_list in strands.items():
            # Find all matching files
            for pattern in pattern_list:
                files = glob.glob(os.path.join(input_dir, pattern))
                for file in files:
                    output_file = os.path.join(output_dir, os.path.basename(file))

                    # Process the file
                    process_file(file, output_file, strand, seq_type)

if __name__ == "__main__":
    input_directory = args.input
    output_directory = args.output
    os.makedirs(output_directory, exist_ok=True)
    main(input_directory, output_directory)