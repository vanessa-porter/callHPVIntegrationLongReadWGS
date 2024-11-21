#!/usr/bin/env python

import concurrent.futures
from pybedtools import BedTool
from itertools import combinations
import os
import csv
import random
import pandas as pd
import argparse

##### ---------------------------------------------------------------------
##### READ IN THE FILES AND PROCESS
##### ---------------------------------------------------------------------

class inputFileError(ValueError):
    '''raise this when there's a mistake in an input file'''

parser = argparse.ArgumentParser(description = "Permutes a DMR window size and calculated the density fold change in 1000 random positions")

parser.add_argument("-d", "--dmr_paths", help = "Path to test file with all the DMR files")
parser.add_argument("-b", "--region", help = "Path to region bed file to be tested")
parser.add_argument("-s", "--sample", help = "Sample name of the tested sample")
parser.add_argument("-t", "--threads", help = "Number of threads")
parser.add_argument("-o", "--out", help = "Output folder path")
parser = parser.parse_args()

# make the inputs bedtools objects
hs = BedTool(parser.region)
genome = 'tables/hg38_chromsizes_autosomes.txt'

# read in the file locations for all the DMR locations
with open(parser.dmr_paths) as f:
    reader = csv.reader(f, delimiter="\t")
    allDMRFiles = list(reader)

allDmrF = [item for sublist in allDMRFiles for item in sublist] # flatten list

# remove sample from the list
ctrlDMRFiles = list(filter(lambda x:parser.sample not in x, allDmrF))
testDMRFile = list(filter(lambda x:parser.sample in x, allDmrF))

##### ---------------------------------------------------------------------
##### GET THE DMR DENSITY IN THE HPV HOTSPOT IN TEST AND CONTROLS
##### ---------------------------------------------------------------------

# Get the size of the HPV-containing hotspot
hsVal = hs[0]
size = hsVal.stop - hsVal.start

# TEST SAMPLE
b = BedTool(testDMRFile[0])
sub = b.intersect(hs, wa=True)
dmr = [i.stop - i.start for i in sub]
s = sum(dmr)
dT = s/size

## OTHER SAMPLES
dC = []
for s in ctrlDMRFiles:
    b = BedTool(s)
    sub = b.intersect(hs, wa=True)
    dmr = [i.stop - i.start for i in sub]
    s = sum(dmr)
    density = s/size
    dC.append(density)

d = {'test': dT, 'controls': dC, 'region': 'hpv_region'}

##### ---------------------------------------------------------------------
##### GET JACCARD AND DENSITY VALUES FOR 1000 SHUFFLED POSITIONS
##### --------------------------------------------------------------------

def process_iteration(i, genome, testDMRFile, ctrlDMRFiles, size):
    print(i)
    sh = hs.shuffle(g = genome)  # get the shuffled position
    region = 'region' + str(i + 1)
    
    # GET THE DENSITY OF THE TEST SAMPLE AT THE SHUFFLED POSITION
    b = BedTool(testDMRFile[0])
    sub = b.intersect(sh, wa=True)
    dmr = [i.stop - i.start for i in sub]
    s = sum(dmr)
    dT = s / size
    
    # GET THE DENSITY OF THE CONTROL SAMPLES AT THE SHUFFLED POSITION
    dC = []
    for s in ctrlDMRFiles:
        b = BedTool(s)
        sub = b.intersect(sh, wa=True)
        dmr = [i.stop - i.start for i in sub]
        s = sum(dmr)
        density = s / size
        dC.append(density)
    
    # save the densities
    d = {'test': dT, 'controls': dC, 'region': region}
    
    return d

# Initialize results lists and append initial results
allShuffleD = [d]

# Number of threads
num_threads = int(parser.threads)  # Set this to the desired number of threads

# Using ThreadPoolExecutor to parallelize the for loop
with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = []
    for i in range(1000):
        futures.append(executor.submit(process_iteration, i, genome, testDMRFile, ctrlDMRFiles, size))
    
    for future in concurrent.futures.as_completed(futures):
        d = future.result()
        allShuffleD.append(d)

# CONVERT TO PANDAS DATAFRAME
dfD = pd.concat(pd.DataFrame(i) for i in allShuffleD)

# save the dataframes
out1 = parser.out + '/densityPermuteTest.txt'
dfD.to_csv(path_or_buf = out1, sep = "\t", header = True, index = False)
