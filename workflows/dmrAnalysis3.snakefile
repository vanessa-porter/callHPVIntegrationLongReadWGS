import os
import pandas as pd
samples_dict = config["samples"]
sample_ids = samples_dict.keys()
EVENT = os.environ.get("EVENT")

dmrPaths = '/path/to/dmr_file_locations.txt'

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
            expand("output/{sample}/event_phase/{event}/dmr_permute/plottingDensityValues.txt", sample=sample_ids, event=EVENT)

### -------------------------------------------------------------------
### DMR Hotspots
### -------------------------------------------------------------------

## 1. Intersect the HPV site midpoint with the DMR hotspots
## 2. Bedtools shuffle the DMR hotspot 1000 times (make 1000 bed files)
## 3. Subset the DMRs within the 1000 random regions

rule subset_dmr_hotspot:
    input:
        hpv="output/{sample}/event_phase/{event}/event_region.bed",
        hs="output/{sample}/methylation/densityDMRHotspotsRegions.sorted.bed"
    output:
        "output/{sample}/event_phase/{event}/hpv_DMR_hotspot.bed"
    shell:
        """
        bedtools intersect -wb -a {input.hpv} -b {input.hs} | cut -f 6,7,8 > {output}
        """

rule permute_test:
    input:
        hs="output/{sample}/event_phase/{event}/hpv_DMR_hotspot.bed"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute/densityPermuteTest.txt"
    shell:
        """
        python scripts/permutateDMRHotspots.py -d sample_txtfiles/dmr_file_locations.txt -b {input.hs} -s {wildcards.sample} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute
        """

rule pvalue:
    input:
        den="output/{sample}/event_phase/{event}/dmr_permute/densityPermuteTest.txt"
    output:
        "output/{sample}/event_phase/{event}/dmr_permute/plottingDensityValues.txt"
    shell:
        """
        scripts/dmrHotspotPvalue.R -d {input.den} -s {wildcards.sample} -e {wildcards.event} -o output/{wildcards.sample}/event_phase/{wildcards.event}/dmr_permute
        """