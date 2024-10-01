# Setting up Resources

**GOAL** In this section we will learn how to set up a rule to run as an independent batch job.
In order to do that, each job must be provisioned with enough resources.

We will use the material in exercise02 in our exercise folder. In all the exercises in that folder, there is a Snakefile and Snakefile.finished
Snakefile.finished is what it might look like when you are done.

The Biowulf Snakemake profile allows you to assign resource requirements to your jobs.

The following is a list of keywords and what they translate to in normal Biowulf submissions.

- The threads keyword is translated to --cpus-per-task
- mem_mb=# required memory in megabytes. Translates to --mem
- disk_mb=# translates to --gres=lscratch
- gpu=# Number of gpus needed. If no gpu_model is given translates to --gres=gpu:#
- gpu_model="MODEL" Which gpu to use. Translates to --gres=gpu:MODEL:#. If the string contains a "|" it is interpreted as a constraint and the gpu allocation is translated to --gres=gpu:# -- constraint=MODEL. Note that for the latter use the feature names for the gpus have to be used. On biowulf these start with 'gpu'. See the example below
- runtime=# Runtime of the rule. Translates to --time=#
- ntasks and nodes
- slurm_partition translates to partition

**NOTE**: If you are using a student account, all rules you want to submit to SLURM must use slurm_partitition="student"

Multiple resouces should be separated by a comma

## In rules

Let's go back to our hisat2 rule. We're using the hisat2 aligner on some fastqs.
We can add more threads (really cpus in this case) by
```snakemake
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    threads: 4
```

For threads specifically, we must also tell the action how to use the threads

```snakemake

  shell:
        """
        hisat2 {params.hisat} -x {input.idx} -U {input.fq} --threads {threads} \
          | samtools sort -T tmp/{wildcards.sample} -O BAM \
          > {output.bam}
        samtools index {output.bam}
        """
```

Go ahead and also set aside 6144 MB of RAM

**Don't forget**, if you're using a student account you must add as one of your resources slurm_partition="student"

Now, go ahead and submit:

```bash
snakemake -s Snakefile --profile ~/snakemake_profile --cores 8 --use-singularity --singularity-prefix=../00container
```
If you make the above into a script you can also submit that using sbatch if you'd like.

Congratulations! You've submitted your first slurm jobs through snakemake!

But wait!

What if something is wrong and you want to end it prematurely! Your brooms (with the Sorcer's Apprentice reference) are doing the wrong thing!
If you 
- didn't use `sbatch`
- are on an interactive node
- used the Biowulf profile (which includes the `--cluster-cancel scancel` automatically)
  just type cntrl-c to end the head job and all others!

Try it!
Then use `sacct` to see your cancelled jobs. (Doesn't always work, seems to be a little buggy)

If you need to cancel everything and you used `sbatch` job your best bet is to do scancel -u your-user-name after preparing to lose your current interactive session.

Then start over.

After canceling a snakemake run partway through, you must run the following in order to continue

```bash
snakemake --unlock
```

## Benchmarking

Often, we are not sure what resources will be required when we first set up a pipeline. 

One way to deal with that is with the benchmarking capability and give a surplus of resources for the first run.

Simply add 
```snakemake
benchmark:
    "benchmarks/{sample}.hisat2.benchmark.txt"
```

## Dynamic resource allocation

It is possible to dynamically assign resources based on the size of your files.

In your resources section, you can put something like
```snakemake
resources:
        mem_mb=lambda wc, input: max(2.5 * input.size_mb, 2000)
```
This gives either gives 2000 MB (the minimum required for a job), or 2.5 * the size of the input files. Whichever is greater.


