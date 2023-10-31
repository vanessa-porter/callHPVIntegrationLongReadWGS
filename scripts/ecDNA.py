#!/usr/bin/env python

import pandas as pd
import os

# import paths
subpath = os.environ.get("SUB")
infopath = os.environ.get("INFO")
depthpath = os.environ.get("DEPTH")
asmpath = os.environ.get("ASM")

# read in files
sub = pd.read_csv(subpath, sep='\t', lineterminator='\n', names = ["chr", "start", "end", "read", "MAPQ", "strand"])
depth = pd.read_csv(depthpath, sep='\t', lineterminator='\n', names = ["chr", "position", "depth"])
asm = pd.read_csv(asmpath, sep='\t', lineterminator='\n', names = ["chr", "start", "end", "read", "MAPQ", "strand"])
info = pd.read_csv(infopath, sep='\t', lineterminator='\n', header=0)

## -----------------------------------------------------------------------------------------
## Check for ecDNA
## -----------------------------------------------------------------------------------------

# test positions for depth
startpos = asm['start'].to_list() 
startpos = startpos + [number + 1 for number in startpos]

endpos = asm['end'].to_list() 
endpos = endpos + [number - 1 for number in endpos]
pos = startpos + endpos
dtest = depth[depth["position"].isin(pos)]

# test the depth of the edges
def condition1(x):
  return x < 2

test1 = sum(condition1(x) for x in dtest["depth"])

# test for coverage outside the asm
sub["length"] = sub["end"] - sub["start"]

def condition2(x):
  return x > 200

test2 = sum(condition2(x) for x in sub["length"])

# test for the number of contigs
test3 = info.shape[0]

# test for the circular asm
test4 = info['circ.'].to_list() 

# test for the graph path
test5 = info['graph_path'].values[0]

# print out result
if ((test1 == 0) and (test2 == 0) and (test3 == 1) and (test4 == "Y")):
    print("ecDNA_detected")
elif ((test1 == 0) and (test2 == 0)):
    print("ecDNA_conditions_partially_met")
else:
    print("ecDNA_conditions_not_met")
