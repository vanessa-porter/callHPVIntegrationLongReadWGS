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

parser = argparse.ArgumentParser(description = "Permutes a DMR window size and calculated the density and jaccard index fold change in 1000 random positions")

parser.add_argument("-d", "--dmr_paths", help = "Path to test file with all the DMR files")
parser.add_argument("-b", "--region", help = "Path to region bed file to be tested")
parser.add_argument("-s", "--sample", help = "Sample name of the tested sample")
parser.add_argument("-t", "--threads", help = "Number of threads")
parser.add_argument("-o", "--out", help = "Output folder path")
parser = parser.parse_args()

# make the inputs bedtools objects
hs = BedTool(parser.region)
genome = '/projects/hpv_nanopore_prj/refs/hg38_chromsizes_autosomes.txt'

# read in the file locations for all the DMR locations
with open(parser.dmr_paths) as f:
    reader = csv.reader(f, delimiter="\t")
    allDMRFiles = list(reader)

allDmrF = [item for sublist in allDMRFiles for item in sublist] # flatten list

# remove sample from the list
ctrlDMRFiles = list(filter(lambda x:parser.sample not in x, allDmrF))
testDMRFile = list(filter(lambda x:parser.sample in x, allDmrF))

# Get all the pairwise combinations of the test vs the controls
testPairs = [s + "," + testDMRFile[0] for s in ctrlDMRFiles]

# Get all the pairwise combinations of the controls
ctrlP = [",".join(map(str, comb)) for comb in list(combinations(ctrlDMRFiles, 2))]
ctrlPairs = random.sample(ctrlP, len(testPairs)) # pick the same number of combinations as the test 

##### ---------------------------------------------------------------------
##### GET THE JACCARD FOR THE HPV HOTSPOT IN TEST AND CONTROL PAIRS
##### ---------------------------------------------------------------------

## TEST PAIRS
testJ = []
for p in testPairs:
    pairs = p.split(",")
    pairsBT = [BedTool(i) for i in pairs]
    shDMRs = [i.intersect(hs, wa=True) for i in pairsBT]
    smpl1 = shDMRs[0]
    smpl2 = shDMRs[1]
    j = smpl1.jaccard(smpl2)
    jaccard = j.get('jaccard')
    testJ.append(jaccard)

## OTHER PAIRS
ctrlJ = []
for p in ctrlPairs:
    pairs = p.split(",")
    pairsBT = [BedTool(i) for i in pairs]
    shDMRs = [i.intersect(hs, wa=True) for i in pairsBT]
    smpl1 = shDMRs[0]
    smpl2 = shDMRs[1]
    j = smpl1.jaccard(smpl2)
    jaccard = j.get('jaccard')
    ctrlJ.append(jaccard)

j = {'test': testJ, 'controls': ctrlJ, 'region': 'hpv_region'}

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

def process_iteration(i, genome, testPairs, ctrlPairs, testDMRFile, ctrlDMRFiles, size):
    print(i)
    sh = hs.shuffle(g = genome)  # get the shuffled position
    region = 'region' + str(i + 1)
    
    # GET THE JACCARD IN THE SHUFFLED REGION FOR EACH TEST PAIR
    testJ = []
    for p in testPairs:
        pairs = p.split(",")
        pairsBT = [BedTool(i) for i in pairs]
        shDMRs = [i.intersect(sh, wa=True) for i in pairsBT]
        smpl1 = shDMRs[0]
        smpl2 = shDMRs[1]
        j = smpl1.jaccard(smpl2)
        jaccard = j.get('jaccard')
        testJ.append(jaccard)
    
    # GET THE JACCARD IN THE SHUFFLED REGION FOR EACH CONTROL PAIR
    ctrlJ = []
    for p in ctrlPairs:
        pairs = p.split(",")
        pairsBT = [BedTool(i) for i in pairs]
        shDMRs = [i.intersect(sh, wa=True) for i in pairsBT]
        smpl1 = shDMRs[0]
        smpl2 = shDMRs[1]
        j = smpl1.jaccard(smpl2)
        jaccard = j.get('jaccard')
        ctrlJ.append(jaccard)
    
    # save the jaccard indices
    j = {'test': testJ, 'controls': ctrlJ, 'region': region}
    
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
    
    return j, d

# Initialize results lists and append initial results
allShuffleJ = [j]
allShuffleD = [d]

# Number of threads
num_threads = int(parser.threads)  # Set this to the desired number of threads

# Using ThreadPoolExecutor to parallelize the for loop
with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = []
    for i in range(1000):
        futures.append(executor.submit(process_iteration, i, genome, testPairs, ctrlPairs, testDMRFile, ctrlDMRFiles, size))
    
    for future in concurrent.futures.as_completed(futures):
        j, d = future.result()
        allShuffleJ.append(j)
        allShuffleD.append(d)

# CONVERT TO PANDAS DATAFRAME
dfJ = pd.concat(pd.DataFrame(i) for i in allShuffleJ)
dfD = pd.concat(pd.DataFrame(i) for i in allShuffleD)

# save the dataframes
out1 = parser.out + '/densityPermuteTest.txt'
out2 = parser.out + '/jaccardPermuteTest.txt'
dfD.to_csv(path_or_buf = out1, sep = "\t", header = True, index = False)
dfJ.to_csv(path_or_buf = out2, sep = "\t", header = True, index = False)
