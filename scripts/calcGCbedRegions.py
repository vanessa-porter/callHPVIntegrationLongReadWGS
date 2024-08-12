#!/usr/bin/env python
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from concurrent.futures import ThreadPoolExecutor, as_completed

class inputFileError(ValueError):
    '''raise this when there's a mistake in an input file'''

parser = argparse.ArgumentParser(description = "Calculates the GC content in a window size around a bed file")

parser.add_argument("-f", "--fasta", help = "Path to fasta file")
parser.add_argument("-b", "--bed", help = "Path to bed file")
parser.add_argument("-w", "--window", help = "Window size to test (default = 1,000,000 bp)", default=1000000, type=int)
parser.add_argument("-t", "--threads", help = "Number of threads", default=10, type=int)
parser.add_argument("-o", "--out", help = "Output file name")
parser = parser.parse_args()

# Input numbers
num_threads = parser.threads 
window_size = parser.window

# Load the BED file into a pandas DataFrame
bed_file = parser.bed
bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end'])

# Load the reference genome into a dictionary
fasta_file = parser.fasta
reference_genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

print("Files processed")

# Expand each region to 1 Mb around the center
def expand_region(row, total_size=window_size):
    chrom, start, end = row['chrom'], row['start'], row['end']
    region_center = (start + end) // 2
    new_start = max(0, region_center - total_size // 2)
    new_end = new_start + total_size
    return pd.Series([chrom, new_start, new_end])

expanded_bed_df = bed_df.apply(expand_region, axis=1)
expanded_bed_df.columns = ['chrom', 'start', 'end']

print("Bed file expanded to window size")

# Function to get sequence from the reference genome
def get_sequence(chrom, start, end):
    if chrom in reference_genome:
        chrom_length = len(reference_genome[chrom])
        # Ensure the end does not exceed chromosome length
        end = min(end, chrom_length)
        return reference_genome[chrom].seq[start:end]
    else:
        return None

# Function to calculate GC content for a single region
def calculate_gc_content_for_region(row):
    chrom, start, end = row['chrom'], row['start'], row['end']
    sequence = get_sequence(chrom, start, end)
    if sequence:
        gc_content = gc_fraction(sequence) * 100  # convert to percentage
        return gc_content
    else:
        return None

# Function to calculate GC content for all regions in parallel
def calculate_gc_content_for_regions_parallel(bed_df, max_workers=num_threads):
    gc_contents = [None] * len(bed_df)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {executor.submit(calculate_gc_content_for_region, row): idx for idx, row in bed_df.iterrows()}
        
        for future in as_completed(future_to_index):
            idx = future_to_index[future]
            try:
                gc_contents[idx] = future.result()
            except Exception as exc:
                gc_contents[idx] = None
                print(f"Exception occurred for index {idx}: {exc}")

    return gc_contents

print("Functions made. Running function ...")

# Calculate GC content for each expanded region specified in the BED file in parallel
expanded_bed_df['gc_content'] = calculate_gc_content_for_regions_parallel(expanded_bed_df)

# Save the results to a new file
output_file = parser.out
expanded_bed_df.to_csv(output_file, sep='\t', index=False, header=False)

print("GC content calculation for specified regions completed. Results saved to:", output_file)