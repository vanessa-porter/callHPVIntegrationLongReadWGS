import os
SAMPLE = os.environ.get("SAMPLE")

## Load config values
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/{sample}/bam/hpv_reads.bam.bai", sample=SAMPLE),
		expand("output/{sample}/bam/methyl_hpv_reads.bam.bai", sample=SAMPLES),
        expand("output/{sample}/methylation/hpv_reads_extract.tsv", sample=SAMPLE),
        expand("output/{sample}/methylation/hpv_reads_methylation.tsv", sample=SAMPLE),
        expand("output/{sample}/methylation/hpv_methylation_frequency.bed", sample=SAMPLE),
        expand("output/{sample}/.combined/events_combined1.txt", sample=SAMPLE),
		expand("output/{sample}/.combined/events_combined2.txt", sample=SAMPLE),
		expand("output/{sample}/.combined/events_combined3.txt", sample=SAMPLE),
        expand("output/{sample}/.combined/events_combined4.txt", sample=SAMPLE),
        expand("output/{sample}/.combined/events_combined5.txt", sample=SAMPLE)

### -------------------------------------------------------------------
### Creating the HPV-only bam files
### -------------------------------------------------------------------

rule select_HPV_reads:
    input:
        bam = lambda w: config["samples"][w.sample]["bam"]
    output:
        "output/{sample}/bam/hpv_read_names.txt"
    shell:
        "samtools view {input.bam} HPV16 HPV18 HPV26 HPV30 HPV31 HPV33 HPV35 HPV45 HPV51 HPV52 HPV58 HPV59 HPV68 HPV69 HPV73 HPV82 | cut -f1 | sort | uniq - > {output}"

rule filter_HPV_reads:
    input:
        names="output/{sample}/bam/hpv_read_names.txt",
        bam=lambda w: config["samples"][w.sample]["bam"]
    output:
        "output/{sample}/bam/hpv_reads.bam"
    shell:
        "picard FilterSamReads -I {input.bam} -O {output} -READ_LIST_FILE {input.names} -FILTER includeReadList -VALIDATION_STRINGENCY SILENT"

rule index_HPV_reads:
    input:
        bam="output/{sample}/bam/hpv_reads.bam"
    output:
        "output/{sample}/bam/hpv_reads.bam.bai"
    shell:
        "samtools index {input.bam}"


rule filter_HPV_methyl_reads:
    input:
        names="output/{sample}/bam/hpv_read_names.txt",
        bam="output/{sample}/bam/methyl_reads.sorted.bam"
    output:
        "output/{sample}/bam/methyl_hpv_reads.bam"
    shell:
        "picard FilterSamReads -I {input.bam} -O {output} -READ_LIST_FILE {input.names} -FILTER includeReadList -VALIDATION_STRINGENCY SILENT"

rule index_HPV_methyl_reads:
    input:
        bam="output/{sample}/bam/methyl_hpv_reads.bam"
    output:
        "output/{sample}/bam/methyl_hpv_reads.bam.bai"
    shell:
        "samtools index {input.bam}"

### -------------------------------------------------------------------
### Make paf of the HPV reads
### -------------------------------------------------------------------

rule HPV_fasta:
    input:
        bai="output/{sample}/bam/hpv_reads.bam.bai",
        bam="output/{sample}/bam/hpv_reads.bam"
    output:
        "output/{sample}/bam/hpv_reads.fasta"
    shell:
        "samtools fasta {input.bam} > {output}"

rule HPV_paf_reads:
    input:
        fasta="output/{sample}/bam/hpv_reads.fasta"
    output:
        "output/{sample}/bam/hpv_reads.paf"
    shell:
        "minimap2 -cx map-ont /projects/alignment_references/9606/hg38_no_alt_TCGA_HTMCP_HPVs/genome/minimap2-2.15-map-ont/hg38_no_alt_TCGA_HTMCP_HPVs_map-ont.mmi {input.fasta} > {output}"


### -------------------------------------------------------------------
### Subset the methylation information to the HPV reads (old chemistry)
### -------------------------------------------------------------------

rule methyl_hpv_reads:
    input:
        names="output/{sample}/bam/hpv_read_names.txt",
        tsv="output/{sample}/methylation/methylation.tsv"
    output:
        "output/{sample}/methylation/hpv_reads_methylation.tsv"
    shell:
        "grep -f {input.names} -F {input.tsv} > {output}"

rule methyl_freq_hpv:
    input:
        tsv="output/{sample}/methylation/methylation_frequency.tsv"
    output:
        "output/{sample}/methylation/hpv_methylation_frequency.tsv"
    shell:
        "grep HPV {input.tsv} > {output}"

### -------------------------------------------------------------------
### Run the integration event caller on the Sniffles VCF
### -------------------------------------------------------------------

rule sniffles:
    input:
        bam=lambda w: config["samples"][w.sample]["bam"]
    output:
        "output/{sample}/vcf/sniffles_output.vcf"
    threads: 30
    shell:
        "sniffles --threads {threads} --max_distance 50 --max_num_splits -1 --report_BND --num_reads_report -1 --min_support 5 --min_seq_size 500 -m {input.bam} -v {output}"

checkpoint call_integration:
    input:
        vcf="output/{sample}/vcf/sniffles_output.vcf",
	    paf="output/{sample}/bam/hpv_reads.paf"
    output:
        directory("output/{sample}/events")
    shell:
        "scripts/callHPVintegration.R -v {input.vcf} -o {output}"


### -------------------------------------------------------------------
### Make bams and fastas for individual events
### -------------------------------------------------------------------

rule event_bams:
    input:
        names="output/{sample}/events/event{i}.txt",
        bai="output/{sample}/bam/hpv_reads.bam.bai",
        bam="output/{sample}/bam/hpv_reads.bam"
    output:
        "output/{sample}/events/event{i}.bam"
    shell:
        "picard FilterSamReads -I {input.bam} -O {output} -READ_LIST_FILE {input.names} -FILTER includeReadList -VALIDATION_STRINGENCY SILENT"

rule event_bams_sort:
    input:
        bam="output/{sample}/events/event{i}.bam"
    output:
        "output/{sample}/events/event{i}.sorted.bam"
    shell:
        "sambamba sort {input.bam} > {output}"

rule event_fastq:
    input:
        bam="output/{sample}/events/event{i}.sorted.bam"
    output:
        "output/{sample}/events/event{i}.fastq"
    shell:
        "samtools fastq {input.bam} > {output}"

rule event_fasta:
    input:
        bam="output/{sample}/events/event{i}.sorted.bam"
    output:
        "output/{sample}/events/event{i}.fasta"
    shell:
        "samtools fasta {input.bam} > {output}"

rule event_depth:
    input:
        bam="output/{sample}/events/event{i}.sorted.bam"
    output:
        "output/{sample}/depth/event{i}_depth.txt"
    shell:
        "samtools depth {input.bam} > {output}"

rule site_bams:
    input:
        names="output/{sample}/events/hpv.site{i}.txt",
        bam="output/{sample}/bam/hpv_reads.bam"
    output:
        "output/{sample}/events/hpv.site{i}.bam"
    shell:
        "picard FilterSamReads -I {input.bam} -O {output} -READ_LIST_FILE {input.names} -FILTER includeReadList -VALIDATION_STRINGENCY SILENT"

rule site_bams_sort:
    input:
        bam="output/{sample}/events/hpv.site{i}.bam"
    output:
        "output/{sample}/events/hpv.site{i}.sorted.bam"
    shell:
        "sambamba sort {input.bam} > {output}"

rule site_depth:
    input:
        bam="output/{sample}/events/hpv.site{i}.sorted.bam"
    output:
        "output/{sample}/depth/hpv.site{i}_depth.txt"
    shell:
        "samtools depth {input.bam} > {output}"

### -------------------------------------------------------------------
### Create the event assemblies
### -------------------------------------------------------------------

rule flye:
    input:
        fastq="output/{sample}/events/event{i}.fastq",
        fasta="output/{sample}/events/event{i}.fasta"
    output:
        "output/{sample}/asm/event{i}/assembly.fasta"
    threads: 10
    shell:
        "flye --nano-raw {input.fastq} -i 3 -t {threads} -o output/{wildcards.sample}/asm/event{wildcards.i}"

### -------------------------------------------------------------------
### Map the assemblies
### -------------------------------------------------------------------

rule map_asm_ref_paf:
    input:
        "output/{sample}/asm/event{i}/assembly.fasta"
    output:
        "output/{sample}/asm/event{i}/assembly.hybrid.paf"
    shell:
        "minimap2 -x asm5 /projects/alignment_references/9606/hg38_no_alt_TCGA_HTMCP_HPVs/genome/minimap2-2.15-map-ont/hg38_no_alt_TCGA_HTMCP_HPVs_map-ont.mmi {input} > {output}"

rule map_reads_asm_paf:
    input:
        reads="output/{sample}/events/event{i}.fastq",
        asm="output/{sample}/asm/event{i}/assembly.fasta"
    output:
        "output/{sample}/asm/event{i}/reads.asm.paf"
    shell:
        "minimap2 -L --MD -Y -x map-ont {input.asm} {input.reads} > {output}"

### -------------------------------------------------------------------
### Call SVs within the assemblies
### -------------------------------------------------------------------

rule map_reads_asm_sam:
    input:
        reads="output/{sample}/events/event{i}.fastq",
        asm="output/{sample}/asm/event{i}/assembly.fasta"
    output:
        "output/{sample}/asm/event{i}/reads.asm.sam"
    shell:
        "minimap2 -L --MD -Y -ax map-ont {input.asm} {input.reads} > {output}"

rule view_reads_asm:
    input:
        sam="output/{sample}/asm/event{i}/reads.asm.sam"
    output:
        "output/{sample}/asm/event{i}/reads.asm.bam"
    shell:
        "samtools view -S -b {input.sam} > {output}"

rule sort_reads_asm:
    input:
        bam="output/{sample}/asm/event{i}/reads.asm.bam"
    output:
        "output/{sample}/asm/event{i}/reads.asm.sorted.bam"
    shell:
        "sambamba sort {input.bam}"

rule asm_sniffles:
    input:
        bam="output/{sample}/asm/event{i}/reads.asm.sorted.bam"
    output:
        "output/{sample}/asm/event{i}/sniffles_asm.vcf"
    threads: 5
    shell:
        "sniffles --threads {threads} --max_distance 50 --max_num_splits -1 --report_BND --num_reads_report -1 --min_support 3 --min_seq_size 500 -m {input.bam} -v {output}"


### -------------------------------------------------------------------
### Aggregate the final inputs to define the wildcards
### -------------------------------------------------------------------

def aggregate_input1(wildcards):
    checkpoint_output = checkpoints.call_integration.get(**wildcards).output[0]
    return expand('output/{sample}/asm/event{i}/assembly.hybrid.paf',
		sample=SAMPLE,
		i=glob_wildcards(os.path.join(checkpoint_output, 'event{i}.txt')).i)

def aggregate_input2(wildcards):
    checkpoint_output = checkpoints.call_integration.get(**wildcards).output[0]
    return expand('output/{sample}/asm/event{i}/reads.asm.paf',
		sample=SAMPLE,
		i=glob_wildcards(os.path.join(checkpoint_output, 'event{i}.txt')).i)

def aggregate_input3(wildcards):
    checkpoint_output = checkpoints.call_integration.get(**wildcards).output[0]
    return expand('output/{sample}/asm/event{i}/sniffles_asm.vcf',
		sample=SAMPLE,
		i=glob_wildcards(os.path.join(checkpoint_output, 'event{i}.txt')).i)

def aggregate_input4(wildcards):
    checkpoint_output = checkpoints.call_integration.get(**wildcards).output[0]
    return expand('output/{sample}/depth/event{i}_depth.txt',
		sample=SAMPLE,
		i=glob_wildcards(os.path.join(checkpoint_output, 'event{i}.txt')).i)

def aggregate_input5(wildcards):
    checkpoint_output = checkpoints.call_integration.get(**wildcards).output[0]
    return expand('output/{sample}/depth/hpv.site{i}_depth.txt',
		sample=SAMPLE,
		i=glob_wildcards(os.path.join(checkpoint_output, 'hpv.site{i}.txt')).i)

rule final1:
    input:
        aggregate_input1
    output:
        "output/{sample}/.combined/events_combined1.txt"
    shell:
        "touch {output}"

rule final2:
    input:
        aggregate_input2
    output:
        "output/{sample}/.combined/events_combined2.txt"
    shell:
        "touch {output}"

rule final3:
    input:
        aggregate_input3
    output:
        "output/{sample}/.combined/events_combined3.txt"
    shell:
        "touch {output}"

rule final4:
    input:
        aggregate_input4
    output:
        "output/{sample}/.combined/events_combined4.txt"
    shell:
        "touch {output}"

rule final5:
    input:
        aggregate_input5
    output:
        "output/{sample}/.combined/events_combined5.txt"
    shell:
        "touch {output}"