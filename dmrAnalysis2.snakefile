import os
import pandas as pd
SAMPLE = os.environ.get("SAMPLE")

# get the events for the sample
eventsPath = "output/" + SAMPLE + "/events/hpv_integration_events_dist.txt"
summary = pd.read_csv(eventsPath, sep='\t', lineterminator='\n', header = 0)
EVENTS = list(summary.identifier)

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
            expand("output/{sample}/event_phase/{event}/hpv_phaseblock_hp.txt", sample=SAMPLE, event=EVENTS)

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

rule get_phaseblocks:
    input:
        vcf="output/{sample}/vcf/phased_merge_output.vcf.gz"
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
        readsOut = reads[0]
        df = pd.DataFrame({'col':readsOut})
        pd.DataFrame.to_csv(df, path_or_buf = output[0], sep = "\t", header = False, index = False)

rule HP1_readnames:
    input:
        bam="output/{sample}/bam/HP1_Converted2Bisulfite.sorted.bam",
        txt="output/{sample}/event_phase/{event}/hpv_read_names.txt"
    output:
        "output/{sample}/event_phase/{event}/HP1_hpv_readnames.txt"
    shell:
        "picard FilterSamReads --INPUT {input.bam} --OUTPUT /dev/stdout --QUIET true --READ_LIST_FILE {input.txt} --FILTER includeReadList --VALIDATION_STRINGENCY SILENT | samtools view | cut -f1 | sort | uniq > {output}"

rule HP2_readnames:
    input:
        bam="output/{sample}/bam/HP2_Converted2Bisulfite.sorted.bam",
        txt="output/{sample}/event_phase/{event}/hpv_read_names.txt"
    output:
        "output/{sample}/event_phase/{event}/HP2_hpv_readnames.txt"
    shell:
        "picard FilterSamReads --INPUT {input.bam} --OUTPUT /dev/stdout --QUIET true --READ_LIST_FILE {input.txt} --FILTER includeReadList --VALIDATION_STRINGENCY SILENT | samtools view | cut -f1 | sort | uniq > {output}"

rule compare_HP:
    input:
        hp1 = "output/{sample}/event_phase/{event}/HP1_hpv_readnames.txt",
        hp2 = "output/{sample}/event_phase/{event}/HP2_hpv_readnames.txt"
    output:
        "output/{sample}/event_phase/{event}/event_haplotype.txt"
    run:
        with open(input.hp1) as f:
            n1 = f.readlines()
        with open(input.hp2) as f:
            n2 = f.readlines()
        if (len(n1) > len(n2)):
            with open(output[0], 'w') as f:
                f.write('HP1')
        else:
            with open(output[0], 'w') as f:
                f.write('HP2')

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
