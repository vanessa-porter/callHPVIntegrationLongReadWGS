import os
import pandas as pd
SAMPLE = os.environ.get("SAMPLE")

# get the events for the sample
intPath = "output/" + SAMPLE + "/intType/"
events = [f for f in os.listdir(intPath) if not f.startswith('.')]

# check if the events have two breaks
eventsPath = [intPath + s for s in events]

# genes
gene_gff = "/gsc/resources/annotation/snpeff/snpeff_binaries/5.0/data/GRCh38.100/genes.gtf"

# two-break events
tbe = []
for d in eventsPath:
    files = os.listdir(d)
    for f in files:
        if(f == 'twoBreak.bed'):
            e = d.replace(intPath, "")
            tbe.append(e)

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/{sample}/intType/{event}/intTest2.txt", sample=SAMPLE, event=tbe),
        expand("output/{sample}/intType/{event}/genic_test/event_location.txt" , sample=SAMPLE, event=tbe)

### -------------------------------------------------------------------
### get depth of regions before and after integration sites
### -------------------------------------------------------------------

rule int_depth:
    input:
        bam = "output/{sample}/events/{event}.sorted.bam",
        bed = "output/{sample}/intType/{event}/twoBreak.bed"
    output:
        "output/{sample}/intType/{event}/twoBreakDepth.bed"
    shell:
        "samtools depth -a -b {input.bed} {input.bam} > {output}"

### -------------------------------------------------------------------
### Intersect the events with genes
### -------------------------------------------------------------------

rule event_bed:
    input:
        "output/{sample}/intType/{event}/event_summary.txt"
    output: 
        "output/{sample}/intType/{event}/{event}_region.bed"
    run:
        df = pd.read_csv(input[0], sep='\t', lineterminator='\n', names = ["chr","pos","HPVchr","HPVpos","HPV.site","event","VAF","read.depth"])
        out = {'chr': [df["chr"][0]], 'start': [min(df["pos"])], 'end': [max(df["pos"])]}
        outDF = pd.DataFrame(data=out)
        pd.DataFrame.to_csv(outDF, path_or_buf = output[0], sep = "\t", header = False, index = False)

rule genic:
    input:
        region = "output/{sample}/intType/{event}/{event}_region.bed",
        genes = gene_gff
    output: 
        "output/{sample}/intType/{event}/genic_test/genic.bed"
    shell:
        "bedtools intersect -wa -f 0.9 -a {input.region} -b {input.genes} | sort | uniq > {output}"

rule partial_genic:
    input:
        region = "output/{sample}/intType/{event}/{event}_region.bed",
        genes = gene_gff
    output: 
        "output/{sample}/intType/{event}/genic_test/partial_genic.bed"
    shell:
        "bedtools intersect -wa -f 0.1 -a {input.region} -b {input.genes} | sort | uniq > {output}"

rule genic_test:
    input:
        genic = "output/{sample}/intType/{event}/genic_test/genic.bed",
        part_genic = "output/{sample}/intType/{event}/genic_test/partial_genic.bed"
    output: 
        "output/{sample}/intType/{event}/genic_test/event_location.txt"
    run:
        t1 = os.stat(input.genic).st_size == 0
        t2 = os.stat(input.part_genic).st_size == 0

        if (t1 == False):
            with open(output[0], 'w') as f:
                f.write('genic')
        elif (t2 == False):
            with open(output[0], 'w') as f:
                f.write('partially_genic')
        else:
            with open(output[0], 'w') as f:
                f.write('intergenic')


### -------------------------------------------------------------------
### run the python test script to call two-break event categories
### -------------------------------------------------------------------

rule intTest2:
    input:
        "output/{sample}/intType/{event}/twoBreakDepth.bed"
    output:
        "output/{sample}/intType/{event}/intTest2.txt"
    shell:
        "DEPTH={input} python scripts/intTypeTestTwoBreak.py > {output}"

