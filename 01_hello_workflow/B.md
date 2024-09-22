Do multiple inputs
Then logs

# Log Files

These early lessons are particularly about establishing good Snakemake habits
and learning how to get back up when something goes wrong. Snakemake keeps
its own detailed logs about each run, which can be accessed in `.snakemake/`. 
However, they are often incomplete. Therefore, it is always good practice to
create our own log files. You can either redirect standard output and error to the 
log file `&>`, or if you want you can use `tee` to make sure some output blazes by while
you watch.

Let's modify our modify and reverse rule. We'll make some basic output
since there isn't much that can be logged at this point.

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
Notice there must be a comma between items under the same parameter.

Also notice how we aren't just telling Snakemake where to make the log and it makes it there.
That's because every parameter we give to Snakemake merely tells Snakemake about the existence of something.
Some other input from us will be required to
make that thing happen. In this case, we need to tell Snakemake how to make the logs.

NOTE: These must use the same wildcards as the output files. We'll cover wildcards momentarily.

## Let's make our first error

Try commenting out all of the actions in the shell so no output is created. 
Then run snakemake and see what happens.

  
