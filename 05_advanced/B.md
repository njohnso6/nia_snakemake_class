# Input Functions

This section provides additional context from the Snakemake tutorial which I feel is unclear.

First, we must understand Snakemake operates in three phases
1. In the initialization phase, the files defining the workflow are parsed and all rules are instantiated.
2. In the DAG phase, the directed acyclic dependency graph of all jobs is built by filling wildcards and matching input files to output files.
3. In the scheduling phase, the DAG of jobs is executed, with jobs started according to the available resources.

Let's say we made a config file with samples in it:
```yaml
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
```

```snakemake
rule bwa_map:
    input:
        "data/genome.fa",
        single_fastq_to_pass_to_bwa_mem
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```
We want to pass a single fastq file to BWA mem. If we wish to retreive it from the config file, we need to know which sample letter it corresponds to. That will be determined by the wildcards.

But the wildcards are not determined until phase two, even though Snakemake wants to resolve what that input file we need is at the beginning.

If we create the following function, we can tell Snakemake to wait until other wildcards have been resolved (only possible because there is a wildcard in the output section). Then we can index into the config file.

```snakemake
def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```
Here, our input function returns the file in the config selected by whichever wildcard is passed to that job, either A or B, but not both.
This would not be necessary if we were not using a config file.

Another way to imagine this is if we had some wildcard that determined the file type. We might
```
def input_function(wildcards):
    if wildcards.condition == "csv":
        return f"data/{wildcards.sample}_A.csv"
    else:
        return f"data/{wildcards.sample}_A.tsv"
