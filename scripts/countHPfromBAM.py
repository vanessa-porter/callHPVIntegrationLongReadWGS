#!/usr/bin/env python
import os
import pysam
import argparse
from collections import Counter

parser = argparse.ArgumentParser(description='Count HP tags from a BAM file.')
parser.add_argument('bam', type=str, help='Path to the BAM file')

args = parser.parse_args()
bam = args.bam

if not bam:
    raise ValueError("The BAM file path is not provided. Please provide the path to your BAM file.")
try:
    bamfile = pysam.AlignmentFile(bam, "rb")
except Exception as e:
    raise ValueError(f"Failed to open the BAM file: {e}")

# Open the BAM file
bamfile = pysam.AlignmentFile(bam, "rb")

# Initialize a counter
hp_counter = Counter()

# Iterate over each read in the BAM file
for read in bamfile.fetch():
    if read.has_tag("HP"):  
        hp_tag = read.get_tag("HP")  
        hp_counter[hp_tag] += 1  

# Determine which HP tag has the most reads
most_common_hp = hp_counter.most_common(1)[0]

# Print 
print(f"HP{most_common_hp[0]}")

# Close the BAM file
bamfile.close()