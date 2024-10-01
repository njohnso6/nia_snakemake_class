# Being a good Biowulf Citizen

Snakemake makes launching a slew of 1000s of jobs as easy as hitting a few keys.
Unfortunately, what feels best for you is not always the best for other Biowulf users (or yourself).

<img width="999" alt="Screenshot 2024-09-26 at 1 21 01 PM" src="https://github.com/user-attachments/assets/fab0b47c-1a94-45c0-88d8-d7b4fc811d58">

<img width="1216" alt="Screenshot 2024-09-26 at 1 27 54 PM" src="https://github.com/user-attachments/assets/4e3cad63-3d9f-48a3-839a-409cd961eec7">

<img width="1224" alt="Screenshot 2024-09-26 at 1 30 00 PM" src="https://github.com/user-attachments/assets/8b733d37-a433-4f9e-af69-9be8b827ef4e">


<img width="1224" alt="Screenshot 2024-09-26 at 1 30 07 PM" src="https://github.com/user-attachments/assets/f5cc7b1d-6f50-4957-9d68-e771ff8a8719">


Please refer to this [presentation](https://hpc.nih.gov/training/handouts/effective_storage_class.pdf) for more information about this topic.

There are a number of ways to accomplish this:

## 1. Combine several jobs into one:
A similarly enchanced workflow will look something like this:

In this case we have three jobs that depend on the sorted alignments (bams) and sorted indicies (bais):
1. infer_experiment
2. geneBody_coverage check
3. Wig file generation
```snakemake
rule rseqc:
    """
    post-process the sorted bam file.
    Copy bam file to /tmp  and operate there. This is done b/c /tmp is
    a bind mount of /lscratch/$SLURM_JOB_ID.
    """
    input: bam = "02aln/{sample}.bam",
           bai = "02aln/{sample}.bam.bai",
           bed = "00ref/R64-1-1.genes.bed12",
           gs  = "00ref/chromosomes"
    output: ie  = "01qc/{sample}.infer_experiment",
            gb  = "01qc/{sample}.geneBodyCoverage.txt",
            gb2 = "01qc/{sample}.geneBodyCoverage.r",
            gb3 = "01qc/{sample}.geneBodyCoverage.curves.pdf",
            wig = "03track/{sample}.wig"
    threads: 1
    resources: mem_mb = 4096
    params: tmp_bam = lambda wc: "/tmp/{s}.bam".format(s=wc.sample)
    singularity:
        "library://wresch/classes/rnaseq:0.8"
    shell:
        """
        cp {input.bam} {input.bai} /tmp
        infer_experiment.py -i {params.tmp_bam} -r {input.bed} > {output.ie}
        geneBody_coverage.py -i {params.tmp_bam} -r {input.bed} \
                -o 01qc/{wildcards.sample}
        bam2wig.py -i {params.tmp_bam} -s {input.gs} \
                -o $(echo {output.wig} | sed 's/.wig//') -t 1000000 -u
        """
```
In this example, notice we avoided reading the bam (sorted alignment) and sorted index (bai) files multiple times.
Instead, we copied them to a temp directory.

## Grouping jobs

Imagine a workflow that looks like this. Each column represents the path of a sample.
<img width="442" alt="Screenshot 2024-10-01 at 9 29 03 AM" src="https://github.com/user-attachments/assets/5f2257bd-8b69-48a2-878e-858ceabcbc7d">

By assigning each of the three stages to the same group, we can force each sample run to run all three steps on the same node.

<img width="442" alt="Screenshot 2024-10-01 at 9 35 15 AM" src="https://github.com/user-attachments/assets/a2af3344-f4fb-40da-bbd6-60f6f2dd42bf">


You can see I have done something similar in the following Snakefile. Notice how each step is assigned the same group.

```snakemake
from snakemake.utils import read_job_properties
import os

# Define the input FASTQ files
PREFIXES = "GL142"
SAMPLES = ["1", "2"]
DEMULTIPLEXES = ["810", "811"]

rule all:
    input:
        expand("fakepkls/{prefix}_{sample}_{demultiplex}.pkl.gzip", 
                prefix=PREFIXES, sample=SAMPLES, demultiplex=DEMULTIPLEXES)


rule demultiplex:
    input:
        "data/{prefix}-{sample}_S{sample}_L001_R1_001.fastq.gz"
    output:
        temp("{prefix}_{sample}_{demultiplex}_trimmed.fastq")
    envmodules:
        "cutadapt"
    threads: 10
    group: "group1"
    resources: runtime=720, mem_mb=64000, disk_mb=64000
    log: "logs/demultiplex/{prefix}_{sample}_{demultiplex}.log"
    shell:
        """
        echo tmpdir: {resources.tmpdir}
        cutadapt --discard-untrimmed --no-indels -e 1 \
        -a 810=ATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -a 811=AGCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -a 812=CGTAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -a 813=CTAGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -m32 -M41 \
        --cores {threads} \
        -o $TMPDIR/{output} {input} &> {log} && 
        touch {output}
        """

rule remove_non_coding:
    input:
        "{prefix}_{sample}_{demultiplex}_trimmed.fastq"
    output:
        fastq=temp("{prefix}_{sample}_{demultiplex}_nomatch.fastq"),
        sam="SAM/{prefix}_{sample}_{demultiplex}.SAM"
    threads: 20
    group: "group1"
    resources: mem_mb=4000, disk_mb=64000
    envmodules: "bowtie/1.3.1"
    log: "logs/remove_non_coding/{prefix}_{sample}/{demultiplex}.log"
    shell:
        """
        bowtie -v 2 -y -S -p {threads} --trim5 2 --trim3 5 --un \
                $TMPDIR/{output.fastq} \
                bowtie/Hs_ncRNA_T_CCA \
                $TMPDIR/{input} \
                {output.sam} &> {log} &&
        touch {output.fastq}
        """

rule dedulication:
    input:
        "{prefix}_{sample}_{demultiplex}_nomatch.fastq"
    output:
        temp("{prefix}_{sample}_{demultiplex}_nomatch_dd.fastq")
    conda: 
        "environment.yml"
    group: "group1"
    resources: mem_mb=64000, disk_mb=64000
    log: "logs/dedeup/{prefix}_{sample}/{demultiplex}.log"
    shell:
        """
        python scripts/dedup_function.py $TMPDIR/{input} $TMPDIR/{output} &> {log} &&
        touch {output}
        """

rule cutadapt_trim:
    input: "{prefix}_{sample}_{demultiplex}_nomatch_dd.fastq"
    output: temp("{prefix}_{sample}_{demultiplex}_nomatch_dd_t2.fastq")
    threads: 20
    group: "group1"
    resources: mem_mb=3000, runtime=720, disk_mb=64000
    envmodules: "cutadapt"
    log: "logs/cutadapt_trim/{prefix}_{sample}/{demultiplex}.log"
    shell:
        """
        cutadapt --cores={threads} -u 2 -u -5 -o $TMPDIR/{output} $TMPDIR/{input} &> {log}
        touch {output}
        """

# Note: bowtie 1 used because it was designed for shorter reads like riboseq has. 
# Bowtie2 is for longer reads like in modern RNA-seq
rule transcriptome_alignment:
    input:
        "{prefix}_{sample}_{demultiplex}_nomatch_dd_t2.fastq"
    output:
        sam="SAM/{prefix}_{sample}_{demultiplex}_2.SAM",
        fastq="{prefix}_{sample}_{demultiplex}_bowtie2match.fastq"
    log: "logs/transcription_alignment/{prefix}_{sample}/{demultiplex}.log"
    group: "group1"
    envmodules: "bowtie/1.3.1"
    threads: 20
    resources: mem_mb=4000, disk_mb=64000
    shell:
        # -al means Write all reads for which at least one alignment was reported to a 
        # file with name <filename>
        # -S means Print alignments in SAM format
        # -v 1 means no more than 1 alignment mismatch allowed
        # -y = Try as hard as possible to find valid alignments when they exist, 
        # including paired-end alignments. This is equivalent to specifying very high 
        # values for the --maxbts and --pairtries options
        # -p number of threads. Speedup is close to linear
        """
        bowtie -S -v 1 -y -p {threads} --al {output.fastq} \
                bowtie/refseq_reduced \
                $TMPDIR/{input} \
                {output.sam} &> {log} 
        """


# NOTE: 1e6 prevents normalization
rule pkldensity:
    input: 
        sam="SAM/{prefix}_{sample}_{demultiplex}_2.SAM",
        fasta="fastas/refseqandmane.fasta_longnames.fasta"
    output:
        "fakepkls/{prefix}_{sample}_{demultiplex}.pkl.gzip"
    log: 
        "logs/pkldensity/{prefix}_{sample}/{demultiplex}.log"
    conda: 
        "environment.yml"
    group: "group1"
    shell:
        """
        python scripts/density_function.py {input.fasta} \
                {input.sam} \
                fakepkls/{wildcards.prefix}_{wildcards.sample}_{wildcards.demultiplex} \
                1000000 25 34 \
                &> {log}
        """
```
Also notice how in the shell section I refer to the TMPDIR, which is lscratch. If done correctly, because all three of these run on the same name node, they can share the lscratch space and I can pass data between them. In order to to that, I must be using the Biowulf profile and I must set aside some disk_mem

## Local rules

Lastly, you may have rules that run very fast, like a rule that makes symlinks for every file. These can be run on the host node and not need to be assigned their own node.

To do that simply add `localrules: rule_you_want_to_run_on_the_host` under your master rule.
