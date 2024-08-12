import pandas as pd
import os

# import paths
subpath = os.environ.get("SUB")
infopath = os.environ.get("INFO")
depthpath = os.environ.get("DEPTH")
asmpath = os.environ.get("ASM")

# read in files

#subpath = "output/HTMCP-03-06-02182/intType/event1/reads_sub_asm.bed"
#infopath = "output/HTMCP-03-06-02182/asm/event1/assembly_info.txt"
#depthpath = "output/HTMCP-03-06-02182/intType/event1/reads_cov_asm.bed"
#asmpath = "output/HTMCP-03-06-02182/intType/event1/assembly.hybrid.filt.bed"


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


# Filter for rows where length > 200 and MAPQ > 10
filtered_df = sub[(sub['length'] > 200) & (sub['MAPQ'] > 10)]

# Sort the DataFrame by 'chr' and 'start'
filtered_df = filtered_df.sort_values(by=['chr', 'start']).reset_index(drop=True)

# Function to find overlapping intervals
def find_overlaps(df):
    overlaps = []
    for i in range(len(df) - 1):
        for j in range(i + 1, len(df)):
            if df.loc[i, 'chr'] != df.loc[j, 'chr']:
                continue
            if df.loc[j, 'start'] <= df.loc[i, 'end']:
                overlap_start = max(df.loc[i, 'start'], df.loc[j, 'start'])
                overlap_end = min(df.loc[i, 'end'], df.loc[j, 'end'])
                if overlap_start < overlap_end:
                    overlaps.append({
                        'chr': df.loc[i, 'chr'],
                        'interval_1': (df.loc[i, 'start'], df.loc[i, 'end']),
                        'interval_2': (df.loc[j, 'start'], df.loc[j, 'end']),
                        'overlap': (overlap_start, overlap_end)
                    })
            else:
                break
    return overlaps

# Find overlaps in the filtered DataFrame
overlaps = find_overlaps(filtered_df)

# number of overlaps
num_overlaps = len(overlaps)

# test for the number of contigs
test3 = info.shape[0]

# test for the circular asm
test4 = info['circ.'].to_list() 
test4 = test4[0]

# test for the graph path
test5 = info['graph_path'].values[0]

# print out result
if ((test1 == 0) and (num_overlaps == 0) and (test4 == 'Y')):
    print("ecDNA_detected")
else:
    print("ecDNA_conditions_not_met")
