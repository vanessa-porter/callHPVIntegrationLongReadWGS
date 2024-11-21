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
            expand("output/{sample}/methylation/cpg_island_methyl_average.txt",sample=SAMPLE)
            expand("output/{sample}/methylation/dmrHotspotIntersectHPV.bed", sample=SAMPLE),
            expand("output/{sample}/methylation/diff_meth.csv", sample=SAMPLE),
            expand("output/{sample}/methylation/densityPlotDMRs.png", sample=SAMPLE),
            expand("output/{sample}/methylation/dmrHotspotIntersectHPV-wa.bed", sample=SAMPLE),
            expand("output/{sample}/events/hpv_integration_events_dist.txt", sample=SAMPLE),
            expand("output/{sample}/methylation/dmr.sorted.bed", sample=SAMPLE)

### -------------------------------------------------------------------
### DMR Hotspots
### -------------------------------------------------------------------

rule hpv_int_dist_merge:
    input:
        "output/{sample}/events/hpv_integration_sites.bed"
    output:
        "output/{sample}/events/hpv_integration_events_distance.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/hpv_int_dist_merge.log"
    shell:
        "bedtools sort -i {input} | bedtools merge -d 500000 -i - -c 4 -o distinct > {output}"

rule hpv_int_dist_txt:
    input:
        "output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/events/hpv_integration_events_dist.txt"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/hpv_int_dist_txt.log"
    shell:
        "scripts/distIdentifierHPVEvents.py SAMPLE={input}"

rule cp_dmr:
    input:
        dmr = lambda w: config["samples"][w.sample]["dmr"]
    output:
        "output/{sample}/methylation/diff_meth.csv"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/cp_dmr.log"
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
        "output/{sample}/methylation/densityPlotDMRs.png"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/dmr_hotspots.log"
    shell:
        "scripts/callDMRHotspots.R -d {input.dmr} -v {input.hpv} -o output/{wildcards.sample}/methylation"

rule sort:
    input:
        dmr="output/{sample}/methylation/densityDMRHotspotsRegions.bed"
    output:
        "output/{sample}/methylation/densityDMRHotspotsRegions.sorted.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/sort.log"
    shell:
        "bedtools sort -i {input.dmr} > {output}"

rule intersect:
    input:
        dmr="output/{sample}/methylation/densityDMRHotspotsRegions.sorted.bed",
        hpv="output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/methylation/dmrHotspotIntersectHPV.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/intersect.log"
    shell:
        "bedtools intersect -wao -a {input.hpv} -b {input.dmr} > {output}"

rule intersect_wa:
    input:
        dmr="output/{sample}/methylation/densityDMRHotspotsRegions.sorted.bed",
        hpv="output/{sample}/events/hpv_integration_events_distance.bed"
    output:
        "output/{sample}/methylation/dmrHotspotIntersectHPV-wa.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/intersect_wa.log"
    shell:
        "bedtools intersect -wa -a {input.hpv} -b {input.dmr} > {output}"

rule make_bed:
    input:
        dmr = "output/{sample}/methylation/diff_meth.csv"
    output:
        "output/{sample}/methylation/dmr.sorted.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/make_bed.log"
    shell:
        "cat {input.dmr} | tail -n +2 | cut -f 1-3 | grep -v ""e"" | bedtools sort -i - > {output}"

rule intersect_cpg:
    input:
        bam=lambda w: config["samples"][w.sample]["bam"],
        cpg="tables/hg38_cpg_islands_bed3.bed"
    output:
        "output/{sample}/methylation/cpg_island_methyl_mod.bed"
    threads: 30
    conda: "config/conda.yaml"
    log: "output/{sample}/log/intersect_cpg.log"
    shell:
        "modkit pileup {input.bam} {output} --include-bed {input.cpg} --log-filepath {log} -t {threads}"

rule intersect_cpg2:
    input:
        freq="output/{sample}/methylation/cpg_island_methyl_mod.bed",
        cpg="tables/hg38_cpg_islands.bed"
    output:
        "output/{sample}/methylation/cpg_island_methyl_freq.bed"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/intersect_cpg2.log"
    shell:
        "cat {input.freq} | grep ""m"" | bedtools intersect -wo -a stdin -b {input.cpg} | cut -f 1-4,11,22 > {output}"

rule average_Cpg:
    input:
        freq="output/{sample}/methylation/cpg_island_methyl_freq.bed"
    output:
        "output/{sample}/methylation/cpg_island_methyl_average.txt"
    conda: "config/conda.yaml"
    log: "output/{sample}/log/average_Cpg.log"
    shell:
        "python scripts/cpgIslandAverageR10.py {input.freq} {output}"

