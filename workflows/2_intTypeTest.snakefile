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
            #expand("output/{sample}/cn/{event}/regionsMeanDepth.bed", sample=SAMPLE, event=EVENTS),
            #expand("output/{sample}/cn/{event}/breakpoint_merge.bed", sample=SAMPLE, event=EVENTS),
            expand("output/{sample}/methylation/event_methyl_freq.tsv", sample=SAMPLE),
            expand("output/{sample}/hpv_size/hpvSizeCategories.txt", sample=SAMPLE),
            expand("output/{sample}/events/{event}_SVsubset.bedpe", sample=SAMPLE, event=EVENTS),
            expand("output/{sample}/events/{event}_SVsubset.bed", sample=SAMPLE, event=EVENTS)

### -------------------------------------------------------------------
### RepeatMaster
### -------------------------------------------------------------------

rule event_fasta:
    input:
        bam="output/{sample}/events/{event}.sorted.bam"
    output:
        "output/{sample}/intType/{event}/reads.fasta"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/event_fasta.log"
    shell:
        "samtools fasta {input.bam} > {output}"

rule repeat_master:
    input:
        fasta="output/{sample}/intType/{event}/reads.fasta"
    output:
        "output/{sample}/intType/{event}/reads.fasta.out.gff"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/repeat_master.log"
    threads: 20
    shell:
        "singularity exec -B /projects,/home RepeatMasker {input.fasta} -pa {threads} -gff -species human"

### -------------------------------------------------------------------
### Subset the methylation by event 
### -------------------------------------------------------------------

rule hpv_methyl_freq:
    input:
        tsv="output/{sample}/methylation/hpv_reads_methylation.tsv"
    output:
        "output/{sample}/methylation/event_methyl_freq.tsv"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/hpv_methyl_freq.log"
    shell: 
        "scripts/eventMethylationNewChemistry.R -m {input.tsv} -d output/{wildcards.sample}/events -o {output}"

### -------------------------------------------------------------------
### Calculate the HPV size
### -------------------------------------------------------------------

rule hpv_integrant_size:
    input:
        paf="output/{sample}/bam/hpv_reads.paf",
        summary="output/{sample}/events/summary.txt"
    output:
        "output/{sample}/hpv_size/hpvSizeCategories.txt"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/hpv_integrant_size.log"
    shell: 
        "scripts/hpv_integrant_size.R -p {input.paf} -s {input.summary} -o output/{wildcards.sample}/hpv_size"

### -------------------------------------------------------------------
### SVs associateed with the event
### -------------------------------------------------------------------

rule subset_vcf:
    input:
        vcf = "output/{sample}/vcf/sniffles_output.vcf",
        reads = "output/{sample}/events/{event}.txt",
    output:
        vcf = "output/{sample}/events/{event}_SVsubset.vcf",
        bed = "output/{sample}/events/{event}_SVsubset.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/subset_vcf.log"
    shell:
        "scripts/lookForSVs.py {input.vcf} {input.reads} {output.vcf} {output.bed}"

rule bedpe:
    input:
        vcf = "output/{sample}/events/{event}_SVsubset.vcf",
    output: 
        "output/{sample}/events/{event}_SVsubset.bedpe"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/bedpe.log"
    shell:
        "SURVIVOR vcftobed {input.vcf} 5 -1 {output}"

### -------------------------------------------------------------------
### get regions before and after integration sites
### -------------------------------------------------------------------

rule split_summary:
    input:
        summary = "output/{sample}/events/summary.txt"
    output:
        "output/{sample}/intType/{event}/event_summary.txt"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/split_summary.log"
    shell:
        "grep {wildcards.event} {input.summary} > {output}"

rule intTest1:
    input:
        summary = "output/{sample}/intType/{event}/event_summary.txt",
        rep = "output/{sample}/intType/{event}/reads.fasta.out.gff"
    output:
        "output/{sample}/intType/{event}/intTest1.txt"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/intTest1.log"
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
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/bamtobed.log"
    shell:
        "bedtools bamtobed -i {input.bam} > {output}"

rule breakpointBed:
    input:
        bed = "output/{sample}/events/{event}.bed",
        summary = "output/{sample}/intType/{event}/event_summary.txt"
    output:
        "output/{sample}/cn/{event}/breakpoint_regions.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/breakpointBed.log"
    shell:
        "OUT={output} BED={input.bed} SUMMARY={input.summary} python scripts/breakpointBed.py"

rule positionDepth:
    input:
        bed = "output/{sample}/cn/{event}/breakpoint_regions.bed",
        bam = "output/{sample}/bam/illumina_wgs.sorted.bam"
    output:
        "output/{sample}/cn/{event}/position_depth.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/positionDepth.log"
    shell:
        "samtools depth {input.bam} -b {input.bed} | awk '{{print $1\"\\t\"$2\"\\t\"$2+1\"\\t\"$3}}' > {output}"

rule regionMeanDepth:
    input:
        bed = "output/{sample}/cn/{event}/breakpoint_regions.bed",
        depth = "output/{sample}/cn/{event}/position_depth.bed"
    output:
        "output/{sample}/cn/{event}/regionsMeanDepth.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/regionMeanDepth.log"
    shell:
        "bedtools map -a {input.bed} -b {input.depth} -c 4 -o mean > {output}"

rule merge:
    input:
        bed = "output/{sample}/cn/{event}/breakpoint_regions.bed"
    output: 
        "output/{sample}/cn/{event}/breakpoint_merge.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/merge.log"
    shell:
        "bedtools merge -d 500000 -i {input} > {output}"

