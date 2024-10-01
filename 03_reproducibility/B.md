# Config Files

## Using config files for sample metadata and other configuration

Config files are helpful when you want to pass your Snakefile to someone else but they have different data.
That way, they can avoid editing your Snakefile and focus on what is important.

For this section, we will work entirely in exercise 4 of our exercise folder. Let's move there now.

## Using config files to reduce duplication

If you look at most of our complete Snakefiles, you will notice that most file paths are duplicated at least twice: in the output of one rule and in the input of the following rule.

But the same path could be used dozens of times! That's bad! What if I need to edit the path? Then I have to edit each one correctly. That could be a huge source of error.

Most people recommend abstracting these to a global variable, but that clutters the Snakefile.

It is considered best practice not to use the config file for this and to find other ways, but I like to use the config file.

Let's consider exercise06 in the exercises folder.

Before continuing, we will need to fix one error in the Snakefile. 

You will need to change `expand("05salmon/{s}", s=samples)` to `expand("05salmon/{s}/quant.sf", s=samples)`
<img width="429" alt="Screenshot 2024-10-01 at 9 00 08 AM" src="https://github.com/user-attachments/assets/5f591f99-7c34-4483-ac50-ad987b0d5192">

Take a look at the number of times "02aln/{sample}.bam" is repeated.

I'll wait a minute.

We can replace that by making an entry in the config file (config.yaml)

```yaml
paths:
  aln_bam: "02aln/{sample}.bam"
```
And we can address this anywhere in the Snakefile using `config['paths']['aln_bam']`

Go ahead and replace a few of those bam entries with the above.

To check if everything is kosher, simply run
```bash
snakemake -n --cores 2
```

If everything is good, you should see
<img width="898" alt="Screenshot 2024-10-01 at 9 16 29 AM" src="https://github.com/user-attachments/assets/694a6835-0b73-4b3a-9f91-63bf30776062">


