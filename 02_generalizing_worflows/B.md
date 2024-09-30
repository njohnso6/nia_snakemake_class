# The DAG

Create your own with `snakemake --forceall --rulegraph | dot -Tpdf > dag_rulegraph.pdf`
May need to install graphviz

![dag_rulegraph](https://github.com/user-attachments/assets/bee97cbd-049a-4889-9ed4-3da6b5d790c3)

## Disambiguation

Sometimes we run into issues likes this where there are two filetypes that need to be processed differently, but will result in a similarly useable file.

It would be unwieldly to append the type of original file to the new file. But it must be differentiated somehow. 

Try creating a new Snakefile in this directory with this code and see what happens.

```snakemake
SAMPLES = ['sample1', 'sample2']

rule all:
    input:
        expand("results/{sample}_processed.txt", sample=SAMPLES)

rule process_csv:
    input:
        "data/pre_samples/{sample}.csv"
    output:
        "results/{sample}_processed.txt"
    shell:
        "echo 'Processing CSV file: {input}' > {output}"

rule process_tsv:
    input:
        "data/pre_samples/{sample}.tsv"
    output:
        "results/{sample}_processed.txt"
    shell:
        "echo 'Processing TSV file: {input}' > {output}"
```

There are various methods for constraining how wildcards are parsed. But that is not sufficient in this instance.

There are two options:
1. Add a `ruleorder` directive at the top to prefer one original filetype over the other.
```snakemake
ruleorder: process_csv > process_tsv

SAMPLES = ['sample1', 'sample2']

rule all:
    input:
        expand("results/{sample}_processed.txt", sample=SAMPLES)

rule process_csv:
    input:
        "data/pre_samples/{sample}.csv"
    output:
        "results/{sample}_processed.txt"
    shell:
        "echo 'Processing CSV file: {input}' > {output}"

rule process_tsv:
    input:
        "data/pre_samples/{sample}.tsv"
    output:
        "results/{sample}_processed.txt"
    shell:
        "echo 'Processing TSV file: {input}' > {output}"
```
  
2. Send the results of each process to different output folders.

```snakemake
rule all:
    input:
        "results/csv/sample1_processed.txt",
        "results/tsv/sample2_processed.txt"

rule process_csv:
    input:
        "data/pre_samples/{sample}.csv"
    output:
        "results/csv/{sample}_processed.txt"
    shell:
        "echo 'Processing CSV file: {input}' > {output}"

rule process_tsv:
    input:
        "data/pre_samples/{sample}.tsv"
    output:
        "results/tsv/{sample}_processed.txt"
    shell:
        "echo 'Processing TSV file: {input}' > {output}"
```
Generally, the second option is better.

In other circumstances it may be wise to constrain wildcards. The following example is taken from the Snakemake tutorial:
Consider the output file {sample}.{group}.txt and assume that the target file is A.1.normal.txt. It is not clear whether dataset="A.1" and group="normal" or dataset="A" and group="1.normal" is the right assignment. Here, constraining the dataset wildcard by {sample,[A-Z]+}.{group} solves the problem.


## Running workflow subsets

Sometimes running the entire workflow isn't desired. For example, while debugging. 

If no wildcards are present, one can request the output of a particular rule simply by
```snakemake
snakemake --cores 2 my_rule
```

If there are wildcards, there are two options:

1. ```snakemake
   snakemake --cores 2 some_output_file_i_want
   ```
2. ```snakemake
   snakemake --cores 2 --until my_rule
   ```
