# Wildcards and Pattern Rules

Here is where we will begin to talk about scalability. One of the most powerful features of Snakemake
is the wildcard feature. It allows us to avoid repeating ourselves and to write D.R.Y. (Do Not Repeat Yourself) code.


## Let's introduce the master rule

Snakemake always runs the first rule in the Snakefile, unless you specify otherwise.
If the required files are not present, it searches for rules that tell it how to make them.
Until now, we have only had one rule. When there are more rules present, we need a "master rule"
to tell Snakemake which files we ultimately want. It will then look for rules that allow the creation
of those files. This rule is conventionally called "all". It only has inputs. This will be important in 
the next steps.

```snakemake
rule modify_and_reverse:
    input:
        hello="input.txt",
        name="my_name.txt"
    output:
        "output.txt"
    log: "logs/modify.log"
    shell:
        """
        # reverse the text, and save to output
        rev {input.hello} > {output} &&
        cat {input.name} >> {output} &&
        echo " I made my first rule!" >> {output}
        echo !! &> {log} &&
        time >> {log}
        """

```


The input of that rule may require the output of another rule. Thus, Snakemake is commonly
understood to "work backwards". 
Very important to know: in the chain of input to output files towards 
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
But it's good to point out now that is better to avoid having intermediate files written to disk
for Biowulf's sake as well as for efficiency.

  General idea, we'll go by themes?
Documentation, automation 
