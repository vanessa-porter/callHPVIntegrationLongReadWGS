import pandas as pd
import os

sumpath = os.environ.get("SUMMARY")
outdir = os.environ.get("OUTDIR")
readbedpath = os.environ.get("READS")
summary = pd.read_csv(sumpath, sep='\t', lineterminator='\n', names = ["chr","pos","HPVchr","HPVpos","HPV.site","event","VAF","read.depth"])
reads = pd.read_csv(readbedpath, sep='\t', lineterminator='\n', names = ["chr","start","end","read", "mapq", "strand"])

# check the number of sites in the event
n = summary.shape[0]
chr = set(summary['chr'].to_list()) 
nchr = len(chr)

# subset reads by mapq value
subreads = reads.loc[reads['mapq'] == 60]

if (n == 2) and (nchr == 1):
    chrom = summary['chr'][0]
    sub = subreads.loc[subreads['chr'] == chrom]
    n1 = min(sub.start)
    n2 = max(sub.end)
    data = [chrom, n1, n2]
    bed = pd.DataFrame.transpose(pd.DataFrame(data))
    pd.DataFrame.to_csv(bed, path_or_buf = outdir + "read_range.bed", sep = "\t", header = False, index = False)