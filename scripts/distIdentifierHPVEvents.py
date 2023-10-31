#!/usr/bin/env python

import os
import pandas as pd

SAMPLE = os.environ.get('SAMPLE')

# get the events for the sample
eventsPath = "output/" + SAMPLE + "/events/hpv_integration_events_distance.bed"
eventsPath2 = "output/" + SAMPLE + "/events/summary.txt"
summary = pd.read_csv(eventsPath, sep='\t', lineterminator='\n', names=['chr', 'start', 'end','hpv_sites'])
summary2 = pd.read_csv(eventsPath2, sep='\t', lineterminator='\n', header = 0)

# get the events for the sample and make an identifier
id = []
for i in [*range(0,len(summary),1)]:
    site = summary.hpv_sites[i]
    sites = site.split(",")
    event = summary2.event[summary2["HPV.site"] == sites[0]]
    event = event.to_list()[0]
    chrom = summary.chr[i]
    start = summary.start[i]
    identifier = event + "_" + chrom + "_" + str(start)
    id.append(identifier)

summary['identifier'] = id

# save as a new table
eventsOutPath = "output/" + SAMPLE + "/events/hpv_integration_events_dist.bed"
pd.DataFrame.to_csv(summary, path_or_buf = eventsOutPath, sep = "\t", header = True, index = False)