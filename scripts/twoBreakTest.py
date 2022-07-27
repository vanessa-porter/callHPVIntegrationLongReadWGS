import pandas as pd
import os

sumpath = os.environ.get("SUMMARY")
outpath = os.environ.get("OUT")
reppath = os.environ.get("REP")
summary = pd.read_csv(sumpath, sep='\t', lineterminator='\n', names = ["chr","pos","HPVchr","HPVpos","HPV.site","event","VAF","read.depth"])
rep = pd.read_csv(reppath, sep='\t', lineterminator='\n', comment='#', names = ["read", "tool", "search", "start", "end", "score", "strand", "P", "motif"])

## -----------------------------------------------------------------------------------------
## Step 1: Check for repeats
## -----------------------------------------------------------------------------------------

centromere = ["ALR", "HSAT4", "HSAT5", "HSAT6", "SST1"]
telomere = ["TAR1"]

# search the reads for telomere/centromere repeats
any_t = rep[rep.motif.str.contains('|'.join(telomere))]
any_c = rep[rep.motif.str.contains('|'.join(centromere))]

# count n reads with repeat
nt = len(set(any_t['read'].to_list()))
nc = len(set(any_c['read'].to_list()))

# percent of reads
n = len(set(rep['read'].to_list()))
pt = nt / n
pc = nc / n

if (pt > 0.5):
  print("telomere_repeats")
elif (pc > 0.5):
  print("centromere_repeats")
else:
  print("no_telomere_centromere_repeats")

## -----------------------------------------------------------------------------------------
## Step 2: Check number of breakpoints
## -----------------------------------------------------------------------------------------

nsites = summary.shape[0]
chrom = set(summary['chr'].to_list()) 
nchr = len(chrom)

if nsites == 2 and nchr == 1:
  print("two-breakpoint_integration")
  a = min(summary.pos)
  b = max(summary.pos)
  chr = summary.chr[0]
  dist = b - a
  if (dist > 20):
    out = {'chr': [chr, chr], 'start': [(a - 11), (b + 10)], 'end': [(a - 10), (b + 11)]}
    outDF = pd.DataFrame(data=out)
    ins = {'chr': [chr,chr], 'start': [(a + 10), (b - 11)], 'end': [(a + 11), (b - 10)]}
    insDF = pd.DataFrame(data=ins)
    bed = pd.concat([outDF, insDF])
    bedSorted = bed.sort_values('start')
    pd.DataFrame.to_csv(bedSorted, path_or_buf = outpath, sep = "\t", header = False, index = False)
  else:
    middle = a + round(int(dist/2))
    out = {'chr': [chr, chr], 'start': [(a - 11), (b + 10)], 'end': [(a - 10), (b + 11)]}
    outDF = pd.DataFrame(data=out)
    ins = {'chr': [chr], 'start': [middle], 'end': [(middle + 1)]}
    insDF = pd.DataFrame(data=ins)
    bed = pd.concat([outDF, insDF])
    bedSorted = bed.sort_values('start')
    pd.DataFrame.to_csv(bedSorted, path_or_buf = outpath, sep = "\t", header = False, index = False)
elif nsites == 2 and nchr == 2:
    print("translocation_integration")
elif nsites == 1:
    print("unmatched_integration")
else:
    print("multi-breakpoint_integration")

