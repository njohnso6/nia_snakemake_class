Before submitting anything to Biowulf, it's important to be aware that all environment variables, modules, etc., are passed 
to the Snakemake master job. Therefore, whatever is needed to parse and run the Snakefile must be loaded.
```
module load snakemake
module load singularity
module load conda
```

Need to be able to get the below correct for the singularity container.
 snakemake --cores=4 --use-singularity \
     --singularity-args '-B $PWD:/data --pwd /data' \
     --singularity-prefix=00container
