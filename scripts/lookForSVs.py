#!/usr/bin/env python

import sys
import os
import pandas as pd
import itertools
import argparse

# Check the number of command line arguments
if not len(sys.argv)==5:
    print("\nError:\tincorrect number of command-line arguments")
    print("Syntax:\tlookForSVs.py [Input VCF] [Read Name Text File] [Output VCF] [Output BED]\n")
    sys.exit()

# File input
fileInput = open(sys.argv[1], "r")

with open(sys.argv[2]) as f:
    reads = f.read().splitlines()

# File output for split multi-alleles
fileOutput = open(sys.argv[3], "w")

### MAKE A VCF FILE CONTAINING ONLY THE ASSOCIATED SVs
bed = []
for strLine in fileInput:
    # Strip the endline character from each input line
    strLine = strLine.rstrip("\n")

    # The '#' character in VCF format indicates that the line is a header. Ignore these and just output to the new file
    if strLine.startswith("#"):
        fileOutput.write(strLine + "\n")
    else:
        # Split the tab-delimited line into an array and then grab the info 
        strArray = [splits for splits in strLine.split("\t")]
        i = strArray[7]
        info = [s for s in i.split(";")]
        info = info[1:]

        # get the read name information
        dInfo = dict(s.split('=') for s in info)
        r = [s for s in dInfo['RNAMES'].split(",")]

        # only keep SVs
        t = strArray[4]
        keep = ['<DEL>', '<DUP>', '<INV>', '<INVDUP>', '<TRA>', '<INS>']

        # if the SV has read names from the text file, output the line of the vcf
        overlap = len(set(r) & set(reads))

        if (overlap > 2) and (t in keep):
            # get for the bed file of the breakpoints
            chr1 = strArray[0]
            pos1 = strArray[1]
            chr2 = dInfo['CHR2']
            pos2 = dInfo['END']

            if (pos1 != pos2):
                fileOutput.write("\t".join(strArray) + "\n")
                # add to the list of chromosomes and positions involved in the SV cluster
                bed.append({chr1: pos1})
                bed.append({chr2: pos2})

### CREATE A BED FILE OF THE REGIONS WITH SV CHANGES
# all the associated chromosomes
all_chr = list(set().union(*(d.keys() for d in bed)))

# make a bed file for each chromosome
b = []
for c in all_chr:
    filt = [key for key in bed if c in key]
    l = []
    for i in filt:
        l.append(i[c])
    nums = [int(x) for x in l]
    nums.sort()

    # pad the SVs by 1 kb for copy number
    if "HPV" in c:
        pos = nums
    elif len(nums) == 1:
        pos = [nums[0] - 1000] + nums + [nums[0] + 1000]
    else:
        pos = [nums[0] - 1000] + nums + [nums[len(nums)-1] + 1000]
        
    # create a bed file 
    chrCol = list(itertools.repeat(c, len(pos)-1))
    startCol = pos[0:(len(pos)-1)]
    endCol = pos[1:(len(pos))]
    bedD = {"chr": chrCol, "start": startCol, "end": endCol}
    bedDF = pd.DataFrame(bedD)
    b.append(bedDF)

outDF = pd.concat(b)


pd.DataFrame.to_csv(outDF, path_or_buf = sys.argv[4], sep = "\t", header = False, index = False)