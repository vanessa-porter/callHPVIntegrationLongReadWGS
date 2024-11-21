#!/usr/bin/env python
import pandas as pd
import sys

def calculate_average_methylation(input_file, output_file):
    # Load the data from the input file
    df = pd.read_csv(input_file, sep='\t', header=None)
    df.columns = ['chromosome', 'start', 'end', 'column4', 'column5', 'column6', 'methylation_frequency', 'column8', 'cpg_island_id']
    # Calculate mean
    df['cpg_island_id'] = df['cpg_island_id'].str.strip()
    result = df.groupby('cpg_island_id')['methylation_frequency'].mean().reset_index()
    # Save 
    result.to_csv(output_file, index=False)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        calculate_average_methylation(input_file, output_file)
