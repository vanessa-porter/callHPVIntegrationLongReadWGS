#!/usr/bin/env python
import pandas as pd
import numpy as np
import csv
import argparse
import concurrent.futures
import pybedtools
import random
from pybedtools import BedTool
from scipy.stats import mannwhitneyu
from sklearn.utils import resample

class inputFileError(ValueError):
    '''raise this when there's a mistake in an input file'''

parser = argparse.ArgumentParser(description = "Caculates DMR density at control regions that mimic a specified GC content")

parser.add_argument("-d", "--dmr_paths", help = "Path to file with all the DMR files")
parser.add_argument("-v", "--hpv_paths", help = "Path to file with all the HPV distance event files")
parser.add_argument("-w", "--window", help = "Window size to test (default = 1,000,000 bp)", default=1000000, type=int)
parser.add_argument("-t", "--test", help = "GC content from HPV regions (bedfile)")
parser.add_argument("-c", "--ctrl", help = "GC content from window regions (bedfile)")
parser.add_argument("-p", "--threads", help = "Number of threads")
parser.add_argument("-o", "--out", help = "Output folder path")
parser = parser.parse_args()

##### ---------------------------------------------------------------------
##### GET THE DMR DENSITY IN THE HPV EXPANDED POSITIONS 
##### ---------------------------------------------------------------------
# input objects
genome = '/projects/hpv_nanopore_prj/refs/hg38_chromsizes_autosomes.txt'
window = parser.window
half = int(window/2)

# read in the file locations for all the DMR locations
with open(parser.dmr_paths) as f:
    reader = csv.reader(f, delimiter="\t")
    allDMRFiles = list(reader)

allDmrF = [item for sublist in allDMRFiles for item in sublist] # flatten list
allDmr = {path.split('/')[1]: path for path in allDmrF}

# read in the file locations for all the HPV locations
with open(parser.hpv_paths) as f:
    reader = csv.reader(f, delimiter="\t")
    allHPVFiles = list(reader)

allHPVF = [item for sublist in allHPVFiles for item in sublist] # flatten list
allHPV = {path.split('/')[1]: path for path in allHPVF}

# all the integrated samples
intSamples = list(allHPV.keys())
allSamples = list(allDmr.keys())

#expand bed
def expand_region(feature):
    midpoint = (feature.start + feature.end) // 2
    new_start = max(midpoint - half, 0)  # Ensure start is not negative
    new_end = midpoint + half
    return pybedtools.Interval(feature.chrom, new_start, new_end)

# for each sample, get the DMR density of the HPV region
alldHPVList = []
for s in intSamples:
    hpvFile = allHPV.get(s)
    dmrFile = allDmr.get(s)
    # bedtools objects
    bHPV = BedTool(hpvFile)
    bDMR = BedTool(dmrFile)
    # expand the HPV regions
    expanded_regions = [expand_region(feature) for feature in bHPV]
    expanded_bed = [pybedtools.BedTool([interval]) for interval in expanded_regions]
    #intersect DMRs with the expanded HPV locations
    dHPV = []
    for r in expanded_bed:
        region = r[0]
        size = region.end - region.start
        # GET THE REGION INFO 
        chrom = region.chrom
        start = region.start
        end = region.end
        region_name = f'{chrom}:{start}-{end}'
        # Intersect DMRs with the current expanded region
        sub = bDMR.intersect(r, wa=True)
        # Calculate the sum of lengths of intersected DMRs
        dmr_lengths = sum(i.end - i.start for i in sub)
        # Calculate density
        density = dmr_lengths / size
        d = {'sample': [s], 'region': region_name, 'density': density}
        dHPV.append(d)
    alldHPVList.append(dHPV)

alldHPV = [item for sublist in alldHPVList for item in sublist]

# Initialize dictionaries for aggregating values
hpv_dict = {'sample': [], 'region': [], 'density': []}

# Aggregate values
for item in alldHPV:
    hpv_dict['sample'].extend(item['sample'])
    hpv_dict['region'].append(item['region'])
    hpv_dict['density'].append(item['density'])

# Make pd dataframe
HPVdf = pd.DataFrame(hpv_dict)
# Remove rows where 'region' column contains "chrX"
HPVdf_XRm = HPVdf[~HPVdf['region'].str.contains('chrX')]

##### -----------------------------------------------------------------------------
##### CHOOSE 10,000 RANDOM POSITIONS THAT CREATE THE SAME GC CONTENT DISTRIBUTION
##### -----------------------------------------------------------------------------
# Load the precomputed regions with GC content
precomputed_file = parser.ctrl
precomputed_df = pd.read_csv(precomputed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'gc_content'])

# Load the expanded regions with GC content
expanded_regions_file = parser.test
expanded_regions_df = pd.read_csv(expanded_regions_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'gc_content'])

# add the GC content to the HPV regions
expanded_regions_df['region'] = expanded_regions_df.apply(
    lambda row: f"{row['chrom']}:{row['start']}-{row['end']}", axis=1
)
HPVdf_merge = pd.merge(HPVdf_XRm, expanded_regions_df[['region', 'gc_content']], on='region', how='left')

# Function to create a GC content distribution histogram
def create_gc_content_histogram(gc_content, bins=500):
    hist, bin_edges = np.histogram(gc_content, bins=bins, range=(0, 100), density=True)
    return hist, bin_edges

# Create histogram for expanded regions
expanded_hist, expanded_bin_edges = create_gc_content_histogram(expanded_regions_df['gc_content'])

# Function to calculate similarity between two histograms using the K-S test
def calculate_similarity(hist1, hist2):
    from scipy.stats import ks_2samp
    return ks_2samp(hist1, hist2).statistic

# Function to sample random regions to match GC content distribution
def sample_matching_gc_content_regions(sliding_windows_df, target_hist, bin_edges, sample_size=10000, max_iterations=1000):
    sliding_windows_df = sliding_windows_df[sliding_windows_df['gc_content'] > 0]
    best_sample = None
    best_similarity = float('inf')
    for _ in range(max_iterations):
        sample = resample(sliding_windows_df, n_samples=sample_size)
        sample_hist, _ = np.histogram(sample['gc_content'], bins=bin_edges, density=True)
        similarity = calculate_similarity(target_hist, sample_hist)
        if similarity < best_similarity:
            best_similarity = similarity
            best_sample = sample
    return best_sample

# Sample random regions to match GC content distribution of the expanded regions
matched_ctrl_df = sample_matching_gc_content_regions(precomputed_df, expanded_hist, expanded_bin_edges, sample_size=10000)

##### -----------------------------------------------------------------------------
##### CALCULATE THE DMR DENSITY AT THE CONTROL REGIONS (CHOOSING A RANDOM SAMPLE)
##### -----------------------------------------------------------------------------

# Function to calculate DMR density for a single region
def calculate_dmr_density_for_region(row, all_dmr_files):
    sample = random.choice(allSamples)
    chrom, start, end = row['chrom'], row['start'], row['end']
    region_name = f'{chrom}:{start}-{end}'
    region = pybedtools.BedTool(f"{chrom}\t{start}\t{end}", from_string=True)
    dmr_file = all_dmr_files.get(sample)
    # get density
    dmr_bed = pybedtools.BedTool(dmr_file)
    intersected = dmr_bed.intersect(region, wa=True)
    dmr_lengths = sum(i.end - i.start for i in intersected)
    density = dmr_lengths / (end - start)
    d = {'sample': sample, 'region': region_name, 'density': density}
    return d

# Function to process each region in parallel
def process_regions_in_parallel(regions_df, all_dmr_files, max_workers=10):
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(calculate_dmr_density_for_region, row, all_dmr_files) for _, row in regions_df.iterrows()]
        for future in concurrent.futures.as_completed(futures):
            try:
                results.append(future.result())
            except Exception as exc:
                print(f"Exception occurred: {exc}")
    return results

# Calculate DMR densities for each selected region in parallel
dmr_densities = process_regions_in_parallel(matched_ctrl_df, allDmr)

# Convert the DMR densities to a DataFrame
dmr_densities_df = pd.DataFrame(dmr_densities)

# add the GC content to the control regions
precomputed_df['region'] = precomputed_df.apply(
    lambda row: f"{row['chrom']}:{row['start']}-{row['end']}", axis=1
)
dmr_densities_merge = pd.merge(dmr_densities_df, precomputed_df[['region', 'gc_content']], on='region', how='left')

##### ---------------------------------------------------------------------
##### TEST THE DIFFERENCE IN MEANS AND SAVE DATAFRAMES 
##### ---------------------------------------------------------------------

# Extract density columns
density_test = HPVdf_XRm['density']
density_ctrl = dmr_densities_merge['density']

# Wilcox test
stat, p_value = mannwhitneyu(density_test, density_ctrl, alternative='two-sided')

# save the dataframes
out_hpv = parser.out + '/hpvDensity' + str(window) + '.txt' 
out_ctrl = parser.out + '/shuffleDensity' + str(window) + '.txt'
HPVdf_merge.to_csv(path_or_buf = out_hpv, sep = "\t", header = True, index = False)
dmr_densities_merge.to_csv(path_or_buf = out_ctrl, sep = "\t", header = True, index = False)

# Print the p value
print(p_value)