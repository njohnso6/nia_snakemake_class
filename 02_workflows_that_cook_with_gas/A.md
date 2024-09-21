# Wildcards, Pattern Rules, and the Master Rule

Here is where we will begin to talk about scalability. One of the most powerful features of Snakemake
is the wildcard feature. It allows us to avoid repeating ourselves and to write D.R.Y. (Do Not Repeat Yourself) code.

**Task**: Process several files in the data folder D.R.Y.ly and then combine them into a single combine_results.txt file.

## Wildcards

How can we avoid repeating ourselves if we want to process three files?

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

Don't worry. The above won't work just yet. We still need a target rule definition. We'll cover that next.

**Notice**:
- when referred to within a rule action, wildcards are treated the same as either "input" or "output" and must be referenced through an attribute.
- the expand function actually permutes the path string and all strings in SAMPLES.
- Instead of listing every letter of SAMPLES, we can do `SAMPLES = glob_wildcards("data/{sample}.txt").sample`



## The rule to rule them *all*

Snakemake always runs the first rule in the Snakefile, unless you specify otherwise.
If the required files are not present, it searches for rules that tell it how to make them.
Until now, we have only had one rule. When there are more rules present, or rules have wildcars, we need a "master rule"
to tell Snakemake which files we ultimately want. It will then look for rules that allow the creation
of those files. This rule is conventionally called "all". It only has inputs. 

The master rule is conventionally the first rule of the Snakefile, but likely not the first line.

```snakemake
# Rule all defines the final output of the pipeline.
rule all:
    input:
        "combined_results.txt"  # Final output, combining all processed samples
```

The input of that rule may require the output of another rule. Thus, Snakemake is commonly
understood to "work backwards". Working backwards allows it to disambiguate wildcards.

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

## The DAG

Create your own with `snakemake --forceall --rulegraph | dot -Tpdf > dag_rulegraph.pdf`
May need to install graphviz

![dag_rulegraph](https://github.com/user-attachments/files/17081234/dag_rulegraph.pdf)




## IDEAS
- General idea, we'll go by themes?
- Documentation, automation
- Have two versions. One with full text and one with bullets for people who prefer to listen
- Have sections and try best to have all the basic stuff out of the way on Monday so veterans can peace out and learn cool stuff Tuesday.
- For the second half, just to save some time on making things just have links to Ryan's stuff to show it in action, and then have them build examples

