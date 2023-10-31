import os
SAMPLE = os.environ.get("SAMPLE")
CHROM_SIZES = config["CHROM_SIZE_PATH"]

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
            expand("output/{sample}/methylation/dmrHotspotIntersectHPV.bed", sample=SAMPLE),
            expand("output/{sample}/methylation/dmrIntersectHPV.bed", sample=SAMPLE),
            expand("output/{sample}/methylation/densityPlotDMRs.png", sample=SAMPLE),
            expand("output/{sample}/methylation/dmrHotspotIntersectHPV-wa.bed", sample=SAMPLE)

### -------------------------------------------------------------------
### DMR Hotspots
### -------------------------------------------------------------------

rule call_dmrs:
    input:
        hp1="output/{sample}/methylation/HP1_MethylFrequency.tsv",
        hp2="output/{sample}/methylation/HP2_MethylFrequency.tsv"
    output:
        "output/{sample}/methylation/diff_meth.csv"
    threads: 30
    shell:
        "scripts/call_dmrs_with_dss.R --methfile1 {input.hp1} --methfile2 {input.hp2} -o {output} -t {threads} -m 3"

rule hpv_int_dist_merge:
    input:
        "output/{sample}/events/hpv_integration_sites.bed"
    output:
        "output/{sample}/events/hpv_integration_events_distance.bed"
    shell:
        "bedtools sort -i {input} | bedtools merge -d 500000 -i - -c 4 -o distinct > {output}"

rule dmr_hotspots:
    input:
        dmr="output/{sample}/methylation/diff_meth.csv",
        hpv="output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/methylation/densityDMRHotspotsRegions.bed",
        "output/{sample}/methylation/densityPlotDMRs.png",
    shell:
        "scripts/callDMRHotspots.R -d {input.dmr} -v {input.hpv} -c CHROM_SIZES -o output/{wildcards.sample}/methylation"

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
        dmr="output/{sample}/methylation/diff_meth.csv"
    output:
        "output/{sample}/methylation/dmr.sorted.bed"
    shell:
        "cat {input.dmr} | tail -n +2 | cut -f 1-3 | grep -v ""e"" | bedtools sort -i - > {output}"

rule int_dmr_hotspot:
    input:
        dmr="output/{sample}/methylation/dmr.sorted.bed",
        hpv="output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/methylation/dmrIntersectHPV.bed"
    shell:
        "bedtools closest -D ref -a {input.dmr} -b {input.hpv} > {output}" 

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
