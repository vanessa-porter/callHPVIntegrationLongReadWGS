#!/usr/bin/env python

import concurrent.futures
import pybedtools
from pybedtools import BedTool
from itertools import combinations
import os
import csv
import random
import pandas as pd
import argparse
from scipy.stats import mannwhitneyu

##### ---------------------------------------------------------------------
##### READ IN THE FILES AND PROCESS
##### ---------------------------------------------------------------------

class inputFileError(ValueError):
    '''raise this when there's a mistake in an input file'''

parser = argparse.ArgumentParser(description = "Converts a tab-separated file of sample information to a samples.yaml file")

parser.add_argument("-d", "--dmr_paths", help = "Path to file with all the DMR files")
parser.add_argument("-v", "--hpv_paths", help = "Path to file with all the HPV distance event files")
parser.add_argument("-w", "--window", help = "Window size to test (default = 1,000,000 bp)", default=1000000, type=int)
parser.add_argument("-t", "--threads", help = "Number of threads")
parser.add_argument("-o", "--out", help = "Output folder path")
parser = parser.parse_args()

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

##### ---------------------------------------------------------------------
##### GET THE DMR DENSITY IN THE HPV POSITIONS (+/- 500 kb from the centre)
##### ---------------------------------------------------------------------

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

##### ---------------------------------------------------------------------
##### GET DENSITY VALUES FOR 10,000 SHUFFLED POSITIONS IN RANDOM SAMPLES
##### ---------------------------------------------------------------------

# generate a random bed file the same size as the testing window (position doesn't matter, will be shuffled)
test_start = 50000000 
test_end = 50000000 + window
test_bed = "chr1" + str(test_start) + str(test_end)
bed_size = f"chr1\t{test_start}\t{test_end}"
bSample = pybedtools.BedTool(bed_size, from_string=True)

def process_iteration(i):
    rsample = random.choice(allSamples)
    sh = bSample.shuffle(g = genome)  # get the shuffled position
    # GET THE REGION INFO 
    interval = sh[0]
    chrom = interval.chrom
    start = interval.start
    end = interval.end
    region_name = f'{chrom}:{start}-{end}'
    # GET THE DENSITY OF THE TEST SAMPLE AT THE SHUFFLED POSITION
    dmrFile = allDmr.get(rsample)
    b = BedTool(dmrFile)
    sub = b.intersect(sh, wa=True)
    dmr_lengths = sum(i.end - i.start for i in sub)
    dT = dmr_lengths / window
    # save the densities
    d = {'sample': [rsample], 'region': region_name, 'density': dT}
    return d

# Initialize results lists and append initial results
allShuffleD = []

# Number of threads
num_threads = int(parser.threads)  # Set this to the desired number of threads

# Using ThreadPoolExecutor to parallelize the for loop
with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = []
    for i in range(10000):
        futures.append(executor.submit(process_iteration, i))
    for future in concurrent.futures.as_completed(futures):
        d = future.result()
        allShuffleD.append(d)

# Initialize dictionaries for aggregating values
sh_dict = {'sample': [], 'region': [], 'density': []}

# Aggregate values
for item in allShuffleD:
    sh_dict['sample'].extend(item['sample'])
    sh_dict['region'].append(item['region'])
    sh_dict['density'].append(item['density'])

# Make pd dataframe
sh_df = pd.DataFrame(sh_dict)

##### ---------------------------------------------------------------------
##### TEST THE DIFFERENCE IN MEANS AND SAVE DATAFRAMES 
##### ---------------------------------------------------------------------

# Extract density columns
density_test = HPVdf_XRm['density']
density_ctrl = sh_df['density']

# Wilcox test
stat, p_value = mannwhitneyu(density_test, density_ctrl, alternative='two-sided')

# save the dataframes
out_hpv = parser.out + '/hpvDensity' + str(window) + '.txt' 
out_ctrl = parser.out + '/shuffleDensity' + str(window) + '.txt'
HPVdf_XRm.to_csv(path_or_buf = out_hpv, sep = "\t", header = True, index = False)
sh_df.to_csv(path_or_buf = out_ctrl, sep = "\t", header = True, index = False)

# Print the p value
print(p_value)
