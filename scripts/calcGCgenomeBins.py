#!/usr/bin/env python

import pandas as pd
import argparse
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import concurrent.futures

class inputFileError(ValueError):
    '''raise this when there's a mistake in an input file'''

parser = argparse.ArgumentParser(description = "Converts a tab-separated file of sample information to a samples.yaml file")

parser.add_argument("-f", "--fasta", help = "Path to fasta file")
parser.add_argument("-w", "--window", help = "Window size to test (default = 1,000,000 bp)", default=1000000, type=int)
parser.add_argument("-s", "--slide", help = "Sliding size", default=100000, type=int)
parser.add_argument("-t", "--threads", help = "Number of threads")
parser.add_argument("-o", "--out", help = "Output file name")
parser = parser.parse_args()

# Load the reference genome into a dictionary
num_threads = int(parser.threads)  # Set this to the desired number of threads
fasta_file = parser.fasta
reference_genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Filter the reference genome to include only autosomes (chromosomes 1 to 22)
autosomes = {f'chr{i}': reference_genome[f'chr{i}'] for i in range(1, 23) if f'chr{i}' in reference_genome}

print("Processed reference genome")

# Define sliding window parameters
window_size = parser.window
slide_size = parser.slide

# Function to generate sliding windows for a given chromosome length
def generate_windows(chrom_length, window_size, slide_size):
    windows = []
    for start in range(0, chrom_length - window_size + 1, slide_size):
        end = start + window_size
        windows.append((start, end))
    return windows

# Function to calculate GC content for a given window
def calculate_gc_content(chrom, start, end):
    sequence = reference_genome[chrom].seq[start:end]
    gc_content = gc_fraction(sequence) * 100  # convert to percentage
    return [chrom, start, end, gc_content]

# Function to calculate GC content for each chromosome in parallel
def calculate_gc_content_sliding_windows_parallel(reference_genome, window_size, slide_size, max_workers=num_threads):
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for chrom in reference_genome:
            chrom_length = len(reference_genome[chrom])
            windows = generate_windows(chrom_length, window_size, slide_size)
            for start, end in windows:
                futures.append(executor.submit(calculate_gc_content, chrom, start, end))
        
        for future in futures:
            results.append(future.result())
    return results

# Calculate GC content in sliding windows across the autosomes in parallel
gc_content_results = calculate_gc_content_sliding_windows_parallel(autosomes, window_size, slide_size)

print("GC content calculated")

# Convert results to a pandas DataFrame
gc_content_df = pd.DataFrame(gc_content_results, columns=['chrom', 'start', 'end', 'gc_content'])

# Save the results to a new file
output_file = parser.out
gc_content_df.to_csv(output_file, sep='\t', index=False, header=False)

print("GC content calculation for sliding windows completed. Results saved to:", output_file)
