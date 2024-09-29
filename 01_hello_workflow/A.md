# Creating a rule

Snakemake is based around rules. The first rule in a snakefile is always the first to run. That will be important later.

Let's start by creating our very first rule!
As mentioned earlier, all rules have an input, and most have an output and an action.
There are multiple types of actions that accomodate script files, python code, or shell commands.
Here, we'll start with the shell command to run some bash commands.

```snakemake
rule modify_and_reverse:
    input:
        "input.txt"
    output:
        "output.txt"
    shell:
        """
        # reverse the text, and save to output
        rev "input.txt" > {output}
        """
```

The input.txt file contains "Hello World!" but backwards! This rule will quickly right it.

You may notice "input.txt" is in quotes, whereas output is in brackets. This is simply to show that 
variables can be hardcoded into a script or they can reused. A variable can be referred to anywhere in
the shell section by enclosing it in curly braces. Please fix that now.

## Chaining commands

Workflows in Snakemake can get complicated and debugging them can become difficult
if the proper precautions are not taken. Therefore, always remember to chain commands with `&&` 
to prevent subsequent commands from running should any command in the chain fail -- a series of failed commands
can complicate debugging.
And, as usual in bash, each line in a multiline command should end in `\`

```snakemake
    shell:
        """
        # reverse the text, and save to output
        rev {input} > {output} &&
        echo " I made my first rule!" >> {output}
        """
```

## Named inputs or outputs

It is possible to have multiple input or output files. They can even be named! (Pretty much any parameter can be 
named). Let's do that here. We'll augment our output with a file that uses your name, just for fun.
Go ahead and make a new input file called "my_name.txt" that contains the text "My name is (your name) and ".

```snakemake
rule modify_and_reverse:
    input:
        hello="input.txt",
        name="my_name.txt"
    output:
        "output.txt"
    shell:
        """
        # reverse the text, and save to output
        rev {input.hello} > {output} &&
        cat {input.name} >> {output} &&
        echo " I made my first rule!" >> {output}
        """
```
Notice how each input "object" now has an attribute that we can reference by name.

Borrowing from the Snakemake tutorial is the following rule:
```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```
Notice how all the input files can be provided at once to bwa mem simply with {input}

NOTE:
It is best practice to have subsequent steps of a workflow in separate, unique, output folders. This keeps the working directory structured. Further, such unique prefixes allow Snakemake to quickly discard most rules in its search for rules that can provide the requested input. This accelerates the resolution of the rule dependencies in a workflow.
