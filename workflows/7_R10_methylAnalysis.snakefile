import os
import pandas as pd
SAMPLE = os.environ.get("SAMPLE")

## Load config values
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

# get the events for the sample
eventsPath = "output/" + SAMPLE + "/events/hpv_integration_events_dist.txt"
summary = pd.read_csv(eventsPath, sep='\t', lineterminator='\n', header = 0)
EVENTS = list(summary.identifier)

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
            expand("output/{sample}/event_phase/{event}/hpv_phaseblock_hp.txt", sample=SAMPLE, event=EVENTS),
            expand("output/{sample}/event_phase/{event}/dmr_dist_hpv.bed", sample=SAMPLE, event=EVENTS),
            expand("output/{sample}/event_phase/{event}/dmr_permute_500kbup/plottingDensityValues.txt", sample=SAMPLE, event=EVENTS),
            expand("output/{sample}/event_phase/{event}/dmr_permute_500kbdown/plottingDensityValues.txt", sample=SAMPLE, event=EVENTS),
            expand("output/{sample}/event_phase/{event}/dmr_permute_1Mb/plottingDensityValues.txt", sample=SAMPLE, event=EVENTS)

### -------------------------------------------------------------------
### DMR Hotspots
### -------------------------------------------------------------------

## 1. subset the hpv site distance file for the one event (or subset of the event) and make it a bed file
## 2. get the read names within the hpv sites at those regions (from text file in events)
## 3. create phase block region file 
## 4. intersect the hpv site region bed with phase blocks bed to get the overlapping phase block 
## 5. get the read names of the reads within the overlapping phase plocks 
## 6. count the number of HPV reads within each haplotype from 5.
## 7. the haplotype with more HPV reads is called as the haplotype
## 8. output a file that has the haplotype (HP1/HP2), the phase block (chr/start/end), and the HPV site names

rule subset_bed:
    input:
        txt="output/{sample}/events/hpv_integration_events_dist.txt"
    output:
        "output/{sample}/event_phase/{event}/event_region.bed"
    shell:
        """
        grep {wildcards.event} {input.txt} > {output}
        """

rule int_dmr_hotspot:
    input:
        dmr="output/{sample}/methylation/dmr.sorted.bed",
        hpv="output/{sample}/event_phase/{event}/event_region.bed"
    output:
        "output/{sample}/event_phase/{event}/dmr_dist_hpv.bed"
    shell:
        "bedtools closest -D ref -a {input.dmr} -b {input.hpv} > {output}" 

rule get_phaseblocks:
    input:
        vcf=lambda w: config["samples"][w.sample]["phase"]
    output:
        "output/{sample}/event_phase/phase_blocks.txt"
    shell:
        """
        whatshap stats {input.vcf} --block-list {output}
        """

rule phaseblock_bed:
    input:
        txt="output/{sample}/event_phase/phase_blocks.txt"
    output:
        "output/{sample}/event_phase/phase_blocks.bed"
    shell:
        """
        cat {input.txt} | tail -n +2 | cut -f 2,3,5 | grep 'chr' > {output}
        """

rule get_readnames:
    input:
        txt="output/{sample}/event_phase/{event}/event_region.bed"
    output:
        "output/{sample}/event_phase/{event}/hpv_read_names.txt"
    run:
        df = pd.read_csv(input[0], sep='\t', lineterminator='\n', names = ['chr', 'start', 'end','hpv_sites','identifier'])
        s = df.hpv_sites[0]
        sample = wildcards.sample
        sites = s.split(",")
        reads = []
        for s in sites:
            readPath = "output/" + sample + "/events/" + s + ".txt"
            with open(readPath, 'r') as fh:
                r = []
                for line in fh:
                    line = line.rstrip("\n")
                    r.append(line)
            reads.append(r)
        readsOut = [item for sublist in reads for item in sublist]
        df = pd.DataFrame({'col':readsOut})
        pd.DataFrame.to_csv(df, path_or_buf = output[0], sep = "\t", header = False, index = False)

rule dist_event_bam:
    input:
        bam="output/{sample}/bam/hpv_reads.bam",
        txt="output/{sample}/event_phase/{event}/hpv_read_names.txt"
    output:
        "output/{sample}/event_phase/{event}/hpv_reads_dist.bam"
    shell:
        """
        picard FilterSamReads --INPUT {input.bam} --OUTPUT {output} --QUIET true --READ_LIST_FILE {input.txt} --FILTER includeReadList --VALIDATION_STRINGENCY SILENT
        samtools index {output} 
        """

rule compare_HP:
    input:
        bam="output/{sample}/event_phase/{event}/hpv_reads_dist.bam"
    output:
        "output/{sample}/event_phase/{event}/event_haplotype.txt"
    shell:
        """
        /projects/vporter_prj/tools/miniconda3/envs/hpv_integration_events/bin/python scripts/countHPfromBAM.py {input.bam} > {output}
        """

#rule HP1_readnames:
#    input:
#        bam="output/{sample}/bam/HP1_Converted2Bisulfite.sorted.bam",
#        txt="output/{sample}/event_phase/{event}/hpv_read_names.txt"
#    output:
#        "output/{sample}/event_phase/{event}/HP1_hpv_readnames.txt"
#    shell:
#        "picard FilterSamReads --INPUT {input.bam} --OUTPUT /dev/stdout --QUIET true --READ_LIST_FILE {input.txt} --FILTER includeReadList --VALIDATION_STRINGENCY SILENT | samtools view | cut -f1 | sort | uniq > {output}"

#rule HP2_readnames:
#    input:
#        bam="output/{sample}/bam/HP2_Converted2Bisulfite.sorted.bam",
#        txt="output/{sample}/event_phase/{event}/hpv_read_names.txt"
#    output:
#        "output/{sample}/event_phase/{event}/HP2_hpv_readnames.txt"
#    shell:
#        "picard FilterSamReads --INPUT {input.bam} --OUTPUT /dev/stdout --QUIET true --READ_LIST_FILE {input.txt} --FILTER includeReadList --VALIDATION_STRINGENCY SILENT | samtools view | cut -f1 | sort | uniq > {output}"

#rule compare_HP:
#    input:
#        hp1 = "output/{sample}/event_phase/{event}/HP1_hpv_readnames.txt",
#        hp2 = "output/{sample}/event_phase/{event}/HP2_hpv_readnames.txt"
#    output:
#        "output/{sample}/event_phase/{event}/event_haplotype.txt"
 #   run:
#        with open(input.hp1) as f:
#            n1 = f.readlines()
#        with open(input.hp2) as f:
#            n2 = f.readlines()
#        if (len(n1) > len(n2)):
#            with open(output[0], 'w') as f:
#                f.write('HP1')
#        else:
#            with open(output[0], 'w') as f:
#                f.write('HP2')

rule intersect_phaseblock:
    input:
        event="output/{sample}/event_phase/{event}/event_region.bed",
        pb="output/{sample}/event_phase/phase_blocks.bed"
    output:
        "output/{sample}/event_phase/{event}/event_phase_block.bed"
    shell:
        """
        bedtools intersect -wo -a {input.event} -b {input.pb} > {output}
        """

rule make_DMR_test_regions_1Mb:
    input:
        "output/{sample}/event_phase/{event}/event_region.bed"
    output:
        df="output/{sample}/event_phase/{event}/event_region_1Mb.bed"
    run:
        df = pd.read_csv(input[0], sep='\t', lineterminator='\n', names = ['chr', 'start', 'end','hpv_sites','identifier'])
        s = df.start[0]
        e = df.end[0]
        mid = round((e + s)/2)
        # full region
        df = df
        df.start[0] = mid - 500000
        df.end[0] = mid + 500000
        pd.DataFrame.to_csv(df, path_or_buf = output[0], sep = "\t", header = False, index = False)

rule make_DMR_test_regions_500kbup:
    input:
        "output/{sample}/event_phase/{event}/event_region.bed"
    output:
        df="output/{sample}/event_phase/{event}/event_region_500kbup.bed",
    run:
        df = pd.read_csv(input[0], sep='\t', lineterminator='\n', names = ['chr', 'start', 'end','hpv_sites','identifier'])
        s = df.start[0]
        e = df.end[0]
        mid = round((e + s)/2)
        # upstream data frame
        df.start[0] = mid - 500000
        df.end[0] = mid
        pd.DataFrame.to_csv(df, path_or_buf = output[0], sep = "\t", header = False, index = False)

rule make_DMR_test_regions_500kbdown:
    input:
        "output/{sample}/event_phase/{event}/event_region.bed"
    output:
        df="output/{sample}/event_phase/{event}/event_region_500kbdown.bed",
    run:
        df = pd.read_csv(input[0], sep='\t', lineterminator='\n', names = ['chr', 'start', 'end','hpv_sites','identifier'])
        s = df.start[0]
        e = df.end[0]
        mid = round((e + s)/2)
        # downstream data frame
        df = df
        df.start[0] = mid 
        df.end[0] = mid + 500000
        pd.DataFrame.to_csv(df, path_or_buf = output[0], sep = "\t", header = False, index = False)

## 1. Intersect the HPV site midpoint with the DMR hotspots
## 2. Bedtools shuffle the DMR hotspot 1000 times (make 1000 bed files)
## 3. Subset the DMRs within the 1000 random regions

rule permute_test_upstream:
    input:
        df="output/{sample}/event_phase/{event}/event_region_500kbup.bed"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute_500kbup/densityPermuteTest.txt"
    threads: 15
    shell:
        """
        /projects/vporter_prj/tools/miniconda3/envs/hpv_integration_events/bin/python3.6 scripts/permutateDMRHotspots.py -d sample_txtfiles/dmr_file_locations.txt -b {input.df} -s {wildcards.sample} -t {threads} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute_500kbup
        """

rule permute_test_downstream:
    input:
        df="output/{sample}/event_phase/{event}/event_region_500kbdown.bed"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute_500kbdown/densityPermuteTest.txt"
    threads: 15
    shell:
        """
        /projects/vporter_prj/tools/miniconda3/envs/hpv_integration_events/bin/python3.6 scripts/permutateDMRHotspots.py -d sample_txtfiles/dmr_file_locations.txt -b {input.df} -s {wildcards.sample} -t {threads} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute_500kbdown
        """

rule permute_test_1Mb:
    input:
        df="output/{sample}/event_phase/{event}/event_region_1Mb.bed"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute_1Mb/densityPermuteTest.txt"
    threads: 15
    shell:
        """
        /projects/vporter_prj/tools/miniconda3/envs/hpv_integration_events/bin/python3.6 scripts/permutateDMRHotspots.py -d sample_txtfiles/dmr_file_locations.txt -b {input.df} -s {wildcards.sample} -t {threads} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute_1Mb
        """

rule pvalue_upstream:
    input:
        den="output/{sample}/event_phase/{event}/dmr_permute_500kbup/densityPermuteTest.txt"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute_500kbup/plottingDensityValues.txt"
    shell:
        """
        scripts/dmrHotspotPvalue.R -d {input.den} -s {wildcards.sample} -e {wildcards.event} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute_500kbup
        """

rule pvalue_downstream:
    input:
        den="output/{sample}/event_phase/{event}/dmr_permute_500kbdown/densityPermuteTest.txt"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute_500kbdown/plottingDensityValues.txt"
    shell:
        """
        scripts/dmrHotspotPvalue.R -d {input.den} -s {wildcards.sample} -e {wildcards.event} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute_500kbdown
        """

rule pvalue_1Mb:
    input:
        den="output/{sample}/event_phase/{event}/dmr_permute_1Mb/densityPermuteTest.txt"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute_1Mb/plottingDensityValues.txt"
    shell:
        """
        scripts/dmrHotspotPvalue.R -d {input.den} -s {wildcards.sample} -e {wildcards.event} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute_1Mb
        """

rule final:
    input:
        hp="output/{sample}/event_phase/{event}/event_haplotype.txt",
        pb="output/{sample}/event_phase/{event}/event_phase_block.bed"
    output:
        "output/{sample}/event_phase/{event}/hpv_phaseblock_hp.txt"
    shell:
        """
        paste {input.pb} {input.hp} > {output}
        """
