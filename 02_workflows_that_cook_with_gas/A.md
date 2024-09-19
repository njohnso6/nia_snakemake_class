# Wildcards,Pattern Rules, and the Master Rule

Here is where we will begin to talk about scalability. One of the most powerful features of Snakemake
is the wildcard feature. It allows us to avoid repeating ourselves and to write D.R.Y. (Do Not Repeat Yourself) code.

**Task**: Process several files in the data folder D.R.Y.ly and then combine them into a single combine_results.txt file.

## Introducing the master rule

Snakemake always runs the first rule in the Snakefile, unless you specify otherwise.
If the required files are not present, it searches for rules that tell it how to make them.
Until now, we have only had one rule. When there are more rules present, we need a "master rule"
to tell Snakemake which files we ultimately want. It will then look for rules that allow the creation
of those files. This rule is conventionally called "all". It only has inputs. This will be important in 
the next steps.

First we'll specify our final file. 
```snakemake
# Rule all defines the final output of the pipeline.
rule all:
    input:
        "combined_results.txt"  # Final output, combining all processed samples
```

The input of that rule may require the output of another rule. Thus, Snakemake is commonly
understood to "work backwards". 

**Important**: in the chain of input to output files towards 
the target rule, every output file path must be unique. In other words, you 
cannot modify a file and have the output have the same name. You will get
some very mysterious errors and skipped rules if you try. This is because Snakemake works
by creating a DAG -- a directed acyclic graph. This is our DAG so far.

Having an input and output with the same name introduces a cycle into the graph and
will confuse Snakemake.

If you insist on having input and output files with the same names, you must 
have the input and output go into different folders. Here is an example from the online
tutorial:

```snakemake
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}
```
We'll talk more about dealing with temporary and intermediate files tomorrow.
But it's good to point out now that it is better to avoid having intermediate files written to disk
for Biowulf's sake as well as for efficiency.

## Wildcards

We said we will "work backwards" from the master rule. But it is conventional to 
place the beginning of the chain immediately after the master rule.

But how can we avoid repeating ourselves if we want to process three files?

```snakemake
# First, define the samples you want to process
SAMPLES = ["A", "B", "C"]

# Rule to process each sample
rule process_sample:
    input:
        "data/{sample}.txt"  # Input is a text file corresponding to each sample
    output:
        "results/{sample}_processed.txt"  # Output is the processed file for each sample
    shell:
        """
        # Simulate processing (for example, copying the input to output)
        echo "Processing {wildcards.sample}" > {output}
        cat {input} >> {output}
        """

# Rule to combine all the processed results into a single file
rule combine_results:
    input:
        expand("results/{sample}_processed.txt", sample=SAMPLES)  # Use expand to get all processed sample files
    output:
        "combined_results.txt"  # The final combined file
    shell:
        """
        # Concatenate all processed sample files into one combined file
        cat {input} > {output}
        """
```

## IDEAS
- General idea, we'll go by themes?
- Documentation, automation
- Have two versions. One with full text and one with bullets for people who prefer to listen
- Have sections and try best to have all the basic stuff out of the way on Monday so veterans can peace out and learn cool stuff Tuesday.

