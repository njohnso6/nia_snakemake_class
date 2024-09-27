Before submitting anything to Biowulf, it's important to be aware that all environment variables, modules, etc., are passed 
to the Snakemake master job. Therefore, whatever is needed to parse and run the Snakefile must be loaded.
```
module load snakemake
module load singularity
module load conda
```

