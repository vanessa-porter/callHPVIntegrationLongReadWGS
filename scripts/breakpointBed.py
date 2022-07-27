import pandas as pd
import os
import itertools   

bedpath = os.environ.get("BED")
outpath = os.environ.get("OUT")
sumpath = os.environ.get("SUMMARY")
reads = pd.read_csv(bedpath, sep='\t', lineterminator='\n', names = ["chr","start","end","read", "mapq", "strand"])
summary = pd.read_csv(sumpath, sep='\t', lineterminator='\n', names = ["chr","pos","HPVchr","HPVpos","HPV.site","event","VAF","read.depth"])

# select for mapq
select = reads['mapq'] == 60
reads = reads[select]

# get only the human alignments on the chromosomes with HPV integration
chr = list(set(summary.chr.tolist()))
subreads = reads[reads['chr'].isin(chr)]

# function for counting occurances
def countOccurrence(a):
  k = {}
  for j in a:
    if j in k:
      k[j] +=1
    else:
      k[j] =1
  return k

# add the chromosome name as the dict keys 
l = []
for c in chr:
  sub = reads[reads['chr'] == c]
  all = sub.start.tolist() + sub.end.tolist()
  bp = countOccurrence(all)
  bpsub = {key: value for key, value in bp.items() if value > 5}
  bpOrder = dict(sorted(bpsub.items()))
  pos = list(bpOrder.keys())
  
  # add outer values to the list to get "background" sections
  if len(pos) == 1:
    posAll = [pos[0] - 10000] + pos + [pos[0] + 10000]
  else:
    posAll = [pos[0] - 10000] + pos + [pos[len(pos)-1] + 10000]
  
  # create a bedfile between all listed breakpoints
  chrCol = list(itertools.repeat(c, len(posAll)-1))
  startCol = posAll[0:(len(posAll)-1)]
  endCol = posAll[1:(len(posAll))]
  bed = {"chr": chrCol, "start": startCol, "end": endCol}
  bedDF = pd.DataFrame(bed)
  l.append(bedDF)

# concatonate the df together
outDF = pd.concat(l)

# remove regions < 5 bp
outDF["length"] = outDF["end"] - outDF["start"]
select = outDF['length'] > 5
outDFsub = outDF[select]
del outDFsub["length"]

pd.DataFrame.to_csv(outDFsub, path_or_buf = outpath, sep = "\t", header = False, index = False)
