# Helpful Commands

When running and debugging Snakemake, there are a number of commands that will make your life much easier.

1. `-n` Never run a workflow for the first time without running `snakemake -n`
2. `-p`: Use this when debugging to print out all the jobscripts Snakemake will submit on your behalf. Especially useful for cluster work, which we'll talk about later.
3. `--debug-dag` This will help you see what is going wrong with wildcards, which we'll talk about later.
