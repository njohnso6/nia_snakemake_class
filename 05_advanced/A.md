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

- 
