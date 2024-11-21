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
            expand("output/{sample}/methylation/dmrHotspotIntersectHPV.bed", sample=SAMPLE),
            expand("output/{sample}/methylation/diff_meth.csv", sample=SAMPLE),
            expand("output/{sample}/methylation/densityPlotDMRs.png", sample=SAMPLE),
            expand("output/{sample}/methylation/dmrHotspotIntersectHPV-wa.bed", sample=SAMPLE),
            expand("output/{sample}/events/hpv_integration_events_dist.txt", sample=SAMPLE),
            expand("output/{sample}/methylation/dmr.sorted.bed", sample=SAMPLE),
            expand("output/{sample}/methylation/cpg_island_methyl_average.txt",sample=SAMPLE)
            expand("output/{sample}/.combined/dist_combined2.txt",sample=SAMPLE),
            expand("output/{sample}/methylation/dmrIntersectHPV.bed", sample=SAMPLE)

### -------------------------------------------------------------------
### DMR Hotspots
### -------------------------------------------------------------------

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
        freq="output/{sample}/methylation/methylation_frequency.tsv.gz",
        cpg="tables/hg38_cpg_islands.bed"
    output:
        "output/{sample}/methylation/cpg_island_methyl_freq.bed"
    shell:
        "zcat {input.freq} | tail -n +2 | bedtools intersect -wo -a stdin -b {input.cpg} | cut -f 1-8,12 > {output}"

rule average_Cpg:
    input:
        freq="output/{sample}/methylation/cpg_island_methyl_freq.bed",

    output:
        "output/{sample}/methylation/cpg_island_methyl_average.txt"
    shell:
        "python scripts/cpgIslandAverageR9.py {input.freq} {output}"

