import os
SAMPLE = os.environ.get("SAMPLE")
asmPath = "output/" + SAMPLE + "/asm/"
EVENTS_ASM = [f for f in os.listdir(asmPath) if not f.startswith('.')]

# check if the events have an assembly
eventsPath = [asmPath + s for s in EVENTS_ASM]

# two-break events
eventsASM = []
for d in eventsPath:
    files = os.listdir(d)
    for f in files:
        if(f == 'assembly.fasta'):
            e = d.replace(asmPath, "")
            eventsASM.append(e)

## Load config values
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

REF = config["GENOME_MMI"]

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/{sample}/intType/{event}/intTest3.txt", sample=SAMPLE, event=eventsASM)

### -------------------------------------------------------------------
### Make bed files
### -------------------------------------------------------------------

rule asm_bam:
    input:
        fasta = "output/{sample}/asm/{event}/assembly.fasta",
        ref = REF
    output:
        "output/{sample}/intType/{event}/assembly.hybrid.bam"
    threads: 50
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/asm_bam.log"
    shell:
        "minimap2 -x map-ont -a -t {threads} {input.ref} {input.fasta} | samtools view -b - > {output}"

rule asm_bed:
    input:
        bam = "output/{sample}/intType/{event}/assembly.hybrid.bam"
    output:
        "output/{sample}/intType/{event}/assembly.hybrid.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/asm_bed.log"
    shell:
        "bedtools bamtobed -i {input.bam} > {output}"

rule read_bed:
    input:
        bam = "output/{sample}/events/{event}.sorted.bam"
    output:
        "output/{sample}/intType/{event}/read.hybrid.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/read_bed.log"
    shell:
        "bedtools bamtobed -i {input.bam} > {output}"

### -------------------------------------------------------------------
### filter the read bed file to only the human chromosomes with > 3 reads
### -------------------------------------------------------------------

rule filter_asm_bed:
    input:
        bed = "output/{sample}/intType/{event}/assembly.hybrid.bed"
    output:
        "output/{sample}/intType/{event}/assembly.hybrid.filt.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/filter_asm_bed.log"
    shell:
        "grep 'chr' {input.bed} > {output}"
    
rule uniq_chr:
    input:
        "output/{sample}/intType/{event}/assembly.hybrid.filt.bed"
    output:
        "output/{sample}/intType/{event}/chr.txt"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/uniq_chr.log"
    shell:
        "awk '{{print $1}}' {input} | sort | uniq > {output}"

rule filter_read_bed:
    input:
        reads = "output/{sample}/intType/{event}/read.hybrid.bed",
        chrs = "output/{sample}/intType/{event}/chr.txt"
    output:
        "output/{sample}/intType/{event}/read.hybrid.filt.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/filter_read_bed.log"
    shell:
        "grep -F -f {input.chrs} {input.reads} > {output}"

### -------------------------------------------------------------------
### Subtract
### -------------------------------------------------------------------

rule bed_subtract:
    input:
        reads = "output/{sample}/intType/{event}/read.hybrid.filt.bed",
        asm = "output/{sample}/intType/{event}/assembly.hybrid.filt.bed"
    output:
        "output/{sample}/intType/{event}/reads_sub_asm.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/bed_subtract.log"
    shell:
        "bedtools subtract -a {input.reads} -b {input.asm} > {output}"

### -------------------------------------------------------------------
### depth
### -------------------------------------------------------------------

rule asm_depth:
    input:
        bam = "output/{sample}/events/{event}.sorted.bam",
        bed = "output/{sample}/intType/{event}/assembly.hybrid.filt.bed"
    output:
        "output/{sample}/intType/{event}/reads_cov_asm.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/asm_depth.log"
    shell:
        "samtools depth -a -b {input.bed} {input.bam} > {output}"

### -------------------------------------------------------------------
### Run the ecDNA test conditions
### -------------------------------------------------------------------

rule intTypeTest:
    input:
        sub = "output/{sample}/intType/{event}/reads_sub_asm.bed",
        info = "output/{sample}/asm/{event}/assembly_info.txt",
        depth = "output/{sample}/intType/{event}/reads_cov_asm.bed",
        asm = "output/{sample}/intType/{event}/assembly.hybrid.filt.bed"
    output:
        "output/{sample}/intType/{event}/intTest3.txt"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/{event}/intTypeTest.log"
    shell:
        "SUB={input.sub} INFO={input.info} DEPTH={input.depth} ASM={input.asm} python scripts/ecDNA.py > {output}"

        
