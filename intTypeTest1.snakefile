import os
import pandas as pd
SAMPLE = os.environ.get("SAMPLE")

# get the events for the sample
eventsPath = "output/" + SAMPLE + "/events/summary.txt"
summary = pd.read_csv(eventsPath, sep='\t', lineterminator='\n', header = 0)
EVENTS = list(set(summary.event))

none = EVENTS.count('none')
if (none > 0):
    EVENTS = EVENTS.remove('none')

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
            expand("output/{sample}/intType/{event}/intTest1.txt", sample=SAMPLE, event=EVENTS),
            expand("output/{sample}/cn/{event}/regionsMeanDepth.bed", sample=SAMPLE, event=EVENTS)

### -------------------------------------------------------------------
### RepeatMaster
### -------------------------------------------------------------------

rule event_fasta:
    input:
        bam="output/{sample}/events/{event}.sorted.bam"
    output:
        "output/{sample}/intTypeTest/{event}/reads.fasta"
    shell:
        "samtools fasta {input.bam} > {output}"

rule repeat_master:
    input:
        fasta="output/{sample}/intTypeTest/{event}/reads.fasta"
    output:
        "output/{sample}/intTypeTest/{event}/reads.fasta.out.gff"
    threads: 20
    shell:
        "singularity exec -B /projects,/home /projects/vporter_prj/tools/repeatmasker_4.1.2.p1--pl5321hdfd78af_1.sif RepeatMasker {input.fasta} -pa {threads} -gff -species human"

### -------------------------------------------------------------------
### get regions before and after integration sites
### -------------------------------------------------------------------

rule split_summary:
    input:
        summary = "output/{sample}/events/summary.txt"
    output:
        "output/{sample}/intType/{event}/event_summary.txt"
    shell:
        "grep {wildcards.event} {input.summary} > {output}"

rule intTest1:
    input:
        summary = "output/{sample}/intType/{event}/event_summary.txt",
        rep = "output/{sample}/intTypeTest/{event}/reads.fasta.out.gff"
    output:
        "output/{sample}/intType/{event}/intTest1.txt"
    shell:
        "SUMMARY={input.summary} REP={input.rep} OUT=output/{wildcards.sample}/intType/{wildcards.event}/twoBreak.bed python scripts/twoBreakTest.py > {output}"


### -------------------------------------------------------------------
### Make bed file of regions within the event
### -------------------------------------------------------------------

rule bamtobed:
    input:
        bam = "output/{sample}/events/{event}.sorted.bam"
    output: 
        "output/{sample}/events/{event}.bed"
    shell:
        "bedtools bamtobed -i {input.bam} > {output}"

rule breakpointBed:
    input:
        bed = "output/{sample}/events/{event}.bed",
        summary = "output/{sample}/intType/{event}/event_summary.txt"
    output:
        "output/{sample}/cn/{event}/breakpoint_regions.bed"
    shell:
        "OUT={output} BED={input.bed} SUMMARY={input.summary} python scripts/breakpointBed.py"

rule positionDepth:
    input:
        bed = "output/{sample}/cn/{event}/breakpoint_regions.bed",
        bam = "output/{sample}/bam/illumina_wgs.sorted.bam"
    output:
        "output/{sample}/cn/{event}/position_depth.bed"
    shell:
        "samtools depth {input.bam} -b {input.bed} | awk '{{print $1\"\\t\"$2\"\\t\"$2+1\"\\t\"$3}}' > {output}"

rule regionMeanDepth:
    input:
        bed = "output/{sample}/cn/{event}/breakpoint_regions.bed",
        depth = "output/{sample}/cn/{event}/position_depth.bed"
    output:
        "output/{sample}/cn/{event}/regionsMeanDepth.bed"
    shell:
        "bedtools map -a {input.bed} -b {input.depth} -c 4 -o mean > {output}"

#rule regionCN:
#    input:
#        ploidy = "../ploidetect/illumina/Ploidetect-pipeline/ploidetect_out/{sample}/*/models.txt",
#        cna = "../ploidetect/illumina/Ploidetect-pipeline/ploidetect_out/{sample}/*/cna_condensed.txt",
#        depth = "output/{sample}/cn/{event}/regionsMeanDepth.bed"
#    output: 
#        "output/{sample}/cn/{event}/region_cna.txt"
#    shell:
#        "scripts/calculateCN.R -p {input.ploidy} -c {input.cna} -d {input.depth} -o {output}"

