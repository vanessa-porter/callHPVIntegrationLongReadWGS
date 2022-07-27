import pandas as pd
import os
from collections import Counter

# read in the file
readpath = os.environ.get("READS")
reads = pd.read_csv(readpath, sep='\t', lineterminator='\n', names = ["chr", "start", "end", "read", "MAPQ", "strand"])

# get the filter values
counter = Counter(reads.chr)
count = Counter({k: c for k, c in counter.items() if c > 3})
lst = [item[0] for item in list(count.items())]
chr = list(filter(lambda x:'chr' in x, lst))

# filter the bed file
filt = reads.chr.isin(chr)
reads_filt = reads[filt]

print('\n'.join(reads_filt.to_string(index = False).split('\n')[1:]))
