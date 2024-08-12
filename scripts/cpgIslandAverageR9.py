#!/usr/bin/env python
import pandas as pd
import sys

def calculate_average_methylation(input_file, output_file):
    # Load the data from the input file
    df = pd.read_csv(input_file, sep='\t', header=None)

    # Assign column names based on the structure of your data, excluding the unnecessary columns
    df.columns = ['chromosome', 'start', 'end', 'column4', 'column5', 'column6', 'methylation_frequency', 'column8', 'cpg_island_id']

    # Remove leading/trailing whitespaces from CpG island IDs
    df['cpg_island_id'] = df['cpg_island_id'].str.strip()

    # Group by the CpG island ID and calculate the average methylation frequency
    result = df.groupby('cpg_island_id')['methylation_frequency'].mean().reset_index()

    # Save the result to the output file
    result.to_csv(output_file, index=False)

    # Print a success message
    print(f"Average methylation frequencies have been saved to {output_file}")
    print(f"Number of unique CpG islands in output: {len(result)}")

if __name__ == "__main__":
    # Ensure the script is called with the correct number of arguments
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        calculate_average_methylation(input_file, output_file)
