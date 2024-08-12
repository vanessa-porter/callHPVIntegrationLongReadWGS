import os
SAMPLE = os.environ.get("SAMPLE")

## Load config values
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

# get the events for the sample
#eventsPath = "output/" + SAMPLE + "/events/summary.txt"
#summary = pd.read_csv(eventsPath, sep='\t', lineterminator='\n', header = 0)
#EVENTS = list(set(summary.event))

#none = EVENTS.count('none')
#if (none > 0):
#    EVENTS = EVENTS.remove('none')

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
            expand("output/{sample}/methylation/cpg_island_methyl_average.txt",sample=SAMPLE)
            #expand("output/{sample}/methylation/dmrHotspotIntersectHPV.bed", sample=SAMPLE),
            #expand("output/{sample}/methylation/diff_meth.csv", sample=SAMPLE),
            #expand("output/{sample}/methylation/densityPlotDMRs.png", sample=SAMPLE),
            #expand("output/{sample}/methylation/dmrHotspotIntersectHPV-wa.bed", sample=SAMPLE),
            #expand("output/{sample}/events/hpv_integration_events_dist.txt", sample=SAMPLE),
            #expand("output/{sample}/methylation/dmr.sorted.bed", sample=SAMPLE)

### -------------------------------------------------------------------
### DMR Hotspots
### -------------------------------------------------------------------

#rule call_dmrs:
#    input:
#        hp1="output/{sample}/methylation/HP1_MethylFrequency.tsv",
#        hp2="output/{sample}/methylation/HP2_MethylFrequency.tsv"
#    output:
#        "output/{sample}/methylation/diff_meth.csv"
#    threads: 30
#    shell:
#        "scripts/call_dmrs_with_dss.R --methfile1 {input.hp1} --methfile2 {input.hp2} -o {output} -t {threads} -m 3"

rule hpv_int_dist_merge:
    input:
        "output/{sample}/events/hpv_integration_sites.bed"
    output:
        "output/{sample}/events/hpv_integration_events_distance.bed"
    shell:
        "bedtools sort -i {input} | bedtools merge -d 500000 -i - -c 4 -o distinct > {output}"

rule hpv_int_dist_txt:
    input:
        "output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/events/hpv_integration_events_dist.txt"
    shell:
        "scripts/distIdentifierHPVEvents.py SAMPLE={input}"

rule cp_dmr:
    input:
        dmr = lambda w: config["samples"][w.sample]["dmr"]
    output:
        "output/{sample}/methylation/diff_meth.csv"
    shell:
        """
        cp {input.dmr} output/{wildcards.sample}/methylation/diff_meth.csv.gz
        gunzip output/{wildcards.sample}/methylation/diff_meth.csv.gz
        """

rule dmr_hotspots:
    input:
        dmr = "output/{sample}/methylation/diff_meth.csv",
        hpv="output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/methylation/densityDMRHotspotsRegions.bed",
        "output/{sample}/methylation/densityPlotDMRs.png",
    shell:
        "scripts/callDMRHotspots.R -d {input.dmr} -v {input.hpv} -o output/{wildcards.sample}/methylation"

rule sort:
    input:
        dmr="output/{sample}/methylation/densityDMRHotspotsRegions.bed"
    output:
        "output/{sample}/methylation/densityDMRHotspotsRegions.sorted.bed"
    shell:
        "bedtools sort -i {input.dmr} > {output}"

rule intersect:
    input:
        dmr="output/{sample}/methylation/densityDMRHotspotsRegions.sorted.bed",
        hpv="output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/methylation/dmrHotspotIntersectHPV.bed"
    shell:
        "bedtools intersect -wao -a {input.hpv} -b {input.dmr} > {output}"

rule intersect_wa:
    input:
        dmr="output/{sample}/methylation/densityDMRHotspotsRegions.sorted.bed",
        hpv="output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/methylation/dmrHotspotIntersectHPV-wa.bed"
    shell:
        "bedtools intersect -wa -a {input.hpv} -b {input.dmr} > {output}"

rule make_bed:
    input:
        dmr = "output/{sample}/methylation/diff_meth.csv"
    output:
        "output/{sample}/methylation/dmr.sorted.bed"
    shell:
        "cat {input.dmr} | tail -n +2 | cut -f 1-3 | grep -v ""e"" | bedtools sort -i - > {output}"

rule intersect_cpg:
    input:
        bam=lambda w: config["samples"][w.sample]["bam"],
        cpg="tables/hg38_cpg_islands_bed3.bed"
    output:
        "output/{sample}/methylation/cpg_island_methyl_mod.bed"
    threads: 30
    log: "output/{sample}/log/intersect_cpg.log"
    shell:
        "modkit pileup {input.bam} {output} --include-bed {input.cpg} --log-filepath {log} -t {threads}"

rule intersect_cpg2:
    input:
        freq="output/{sample}/methylation/cpg_island_methyl_mod.bed",
        cpg="tables/hg38_cpg_islands.bed"
    output:
        "output/{sample}/methylation/cpg_island_methyl_freq.bed"
    shell:
        "cat {input.freq} | grep ""m"" | bedtools intersect -wo -a stdin -b {input.cpg} | cut -f 1-4,11,22 > {output}"

rule average_Cpg:
    input:
        freq="output/{sample}/methylation/cpg_island_methyl_freq.bed",

    output:
        "output/{sample}/methylation/cpg_island_methyl_average.txt"
    shell:
        "python scripts/cpgIslandAverageR10.py {input.freq} {output}"

#checkpoint make_individual_bed:
#    input:
#        "output/{sample}/events/hpv_integration_events_distance.bed"
#    output:
#        directory("output/{sample}/events/dist_events")
#    shell:
#        """
#        mkdir {output}
#        awk '{{filename = sprintf("hpv_dist_event%d.bed", NR); print >filename; close(filename)}}' {input}
#        mv hpv_dist_event*.bed {output}
#        """

#rule int_dmr_hotspot:
#    input:
#        dmr="output/{sample}/methylation/dmr.sorted.bed",
#        hpv="output/{sample}/events/dist_events/hpv_dist_event{i}.bed"
#    output:
#        "output/{sample}/methylation/dist_events/dmrIntersectHPV{i}.bed"
#    shell:
#        "bedtools closest -D ref -a {input.dmr} -b {input.hpv} > {output}" 

#rule int_dmr_hotspot_combine:
#    input:
#        bed="output/{sample}/methylation/dist_events/dmrIntersectHPV{i}.bed"
#    output:
#        "output/{sample}/methylation/dmrIntersectHPV.bed"
#    shell:
#        "cat {input} > {output}" 

#def aggregate_input1(wildcards):
#    checkpoint_output = checkpoints.make_individual_bed.get(**wildcards).output[0]
#    return expand('output/{sample}/methylation/dist_events/dmrIntersectHPV{i}.bed',
#		sample=SAMPLE,
#		i=glob_wildcards(os.path.join(checkpoint_output, 'hpv_dist_event{i}.bed')).i)

#def aggregate_input2(wildcards):
#    checkpoint_output = checkpoints.make_individual_bed.get(**wildcards).output[0]
#    return expand('output/{sample}/methylation/dmrIntersectHPV.bed',
#		sample=SAMPLE,
#		i=glob_wildcards(os.path.join(checkpoint_output, 'hpv_dist_event{i}.bed')).i)

#rule final1:
#    input:
#        aggregate_input1
#    output:
#        "output/{sample}/.combined/dist_combined1.txt"
#    shell:
#        "touch {output}"

#rule final2:
#    input:
#        aggregate_input2
#    output:
#        "output/{sample}/.combined/dist_combined2.txt"
#    shell:
#        "touch {output}"

#rule HP1_readnames:
#    input:
#        bed="output/{sample}/events/hpv_integration_events_distance.bed",
#        bam="output/{sample}/bam/HP1_Converted2Bisulfite.sorted.bam"
#    output:
#        "output/{sample}/event_phase/{event}/HP1_readnames.txt"
#    shell:
#        "samtools view -L {input.bed} {input.bam} | cut -f1 | sort | uniq > {output}"

#rule HP2_readnames:
#    input:
#        bed="output/{sample}/events/hpv_integration_events_distance.bed",
#        bam="output/{sample}/bam/HP2_Converted2Bisulfite.sorted.bam"
#    output:
#        "output/{sample}/event_phase/{event}/HP2_readnames.txt"
#    shell:
#        "samtools view -L {input.bed} {input.bam} | cut -f1 | sort | uniq > {output}"
