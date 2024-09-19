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
