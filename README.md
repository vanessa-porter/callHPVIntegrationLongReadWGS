# Calling HPV Integration Using Oxford Nanopore Long-Read Sequencing
Workflow for calling and analyzing HPV integration events in ONT long-read sequencing data. 

# Installation
This will clone the repository. You can run workflow within this directory.
```
git clone https://github.com/vanessa-porter/callONTIntegration.git
```

### Dependencies
> To run this workflow, you must have snakemake (v6.12.3) and conda. You can install snakemake using [this guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). The remaining dependencies will be downloaded automatically within the snakemake workflow.

# Input Files

### **Method 1**: R9 Flow Cells <br />
- bam: ONT alignment (bam file)
- methyl: Methyl-tagged ONT alignment from NanoPolish (bam file)
- HP1: HP1-tagged ONT alignment from NanoMethPhase (bam file)
- HP2: HP2-tagged ONT alignment from NanoMethPhase (bam file)
- dmr: Allelic differentially methylated regions (DSS callDMR.txt)
- phase: Phased SNVs from WhatsHap (vcf file)

### **Method 2**: R10 Flow Cells <br />
- bam: ONT alignment with added HP tags (bam file)
- dmr: Allelic differentially methylated regions (DSS callDMR.txt)
- phase: Phased SNVs from WhatsHap (vcf file)

# Set Up Configuration Files

### **Edit the config files**

#### **Example parameters.yaml:** <br />
Config files to specify parameters and paths needed for the workflow. The main parameter to include is the genome path.

```
genome_path: /path/to/genome/fasta
```

#### **Example samples.yaml:** <br />
Main config file to specify input files.

```
samples:
    sample_name
        bam: /path/to/bam/file
        dmr: /path/to/dmr/file
        phase: /path/to/vcf/file
```

#### Converting sample paths to yaml file
A text file can be converted to the samples.yaml file using the scripts/samplestoyaml.py script. The tsv file should have the sample name in one column and the path in another and be tab delimited (no header). 

```
scripts/samplestoyaml.py -t samples.txt -o config/samples.yaml
```

# **Run Workflow**
Each section of the workflow has to be run in succession, according to the listed number on the workflow. The correct version should be used for steps that have different workflows for R9 and R10 flow cells. Samples must be ran one at a time, with each named in the command to pull the paths from the config file.  

```
SAMPLE=sample_name snakemake -s 1_R10_callIntegration.snakefile -c 30 --use-conda
SAMPLE=sample_name snakemake -s 2_intTypeTest.snakefile -c 30 --use-conda
SAMPLE=sample_name snakemake -s 3_intTypeTest.snakefile -c 30 --use-conda
SAMPLE=sample_name snakemake -s 4_intTypeTest.snakefile -c 30 --use-conda
SAMPLE=sample_name snakemake -s 5_intTypeTest.snakefile -c 30 --use-conda
SAMPLE=sample_name snakemake -s 6_R10_methylAnalysis.snakefile -c 30 --use-conda
SAMPLE=sample_name snakemake -s 7_R10_methylAnalysis.snakefile -c 30 --use-conda
```

# Contributors
The contributors of this project are Vanessa Porter

<a href="https://github.com/vanessa-porter/illuminaCallHPVInt/graphs/contributors">
</a>

