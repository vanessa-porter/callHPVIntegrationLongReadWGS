import os
import pandas as pd
SAMPLE = os.environ.get("SAMPLE")

# get the events for the sample
eventsPath = "output/" + SAMPLE + "/events/summary.txt"
summary = pd.read_csv(eventsPath, sep='\t', lineterminator='\n', header = 0)
EVENTS = list(set(summary.event))

none = EVENTS.count('none')
if (none > 0):
    EVENTS = EVENTS.remove('none')

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/{sample}/intType/{event}/integration_type.txt", sample=SAMPLE, event=EVENTS)

### -------------------------------------------------------------------
### Rule
### -------------------------------------------------------------------

rule final:
    input:
        directory("output/{sample}/intType/{event}")
    output:
        "output/{sample}/intType/{event}/integration_type.txt"
    shell:
        "DIR={input} python scripts/final.py > {output}"
