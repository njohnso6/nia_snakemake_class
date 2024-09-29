# Common Errors

As mentioned, because Snakemake builds a DAG, it's important for every file to have a unique name and location. 
Let's demonstrate that.

Commonly, people want to make a modification to a file and leave it essentially the same.
Therefore, they don't change the name.
```
rule modify_and_reverse:
    input:
        "input.txt"
    output:
        "input.txt"
    shell:
        """
        # reverse the text, and save to output
        rev "input.txt" > {output}
        """
```
What is the result of this?

Another common error is to forget the comma between the input or output items. Since Python concatenates subsequent strings, this can lead to unexpected behavior.
