#!/usr/bin/env python

import sys
import os
import pandas as pd
import numpy as np
import itertools
from collections import Counter

# Check the number of command line arguments
if not len(sys.argv)==5:
    print("\nError:\tincorrect number of command-line arguments")
    print("Syntax:\tlookForSVs.py [Input VCF] [Read Name Text File] [Output VCF] [Output BED]\n")
    sys.exit()


# read in the read names of the event
with open("events/event1.txt") as f:
    reads = f.read().splitlines()

# read in the paf table
paf = pd.read_csv("bam/hpv_reads.paf", sep='\t', lineterminator='\n', names = ["qname","qlen","qstart","qend","strand","tname","tlen","tstart","tend","nmatch","alen","mapq","NM","ms","AS","nn","tp","cm","s1", "s2", "de","rl","cg" ])

# subset the columns
paf = paf.loc[:,'qname':'mapq']
p = paf[paf["mapq"] == 60]

p['tstart2'] = p['strand'].apply(lambda x: p['tstart'] if x == "+" else p['tend'])

# count the reoccuring breaks
bp1 = list(p.tstart)
count1 = Counter(bp1)
c1 = Counter(i for i in count1.elements() if count1[i] >= 2)
starts = dict(c1)
ls = list(starts.keys())

bp2 = list(p.tend)
count2 = Counter(bp2)
c2 = Counter(i for i in count2.elements() if count2[i] >= 2)
ends = dict(c2)
le = list(ends.keys())

# identify all read names
reads = list(set(p.qname))

# loop through each read 
for r in reads:
    psub = p[p["qname"] == r]
    psort = psub.sort_values(by=['qstart'])
    pdic = psub.to_dict('records')


