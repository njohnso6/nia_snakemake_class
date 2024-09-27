# Scaling on Biowulf

It's finally time to really kick Snakemake into high gear! One of its great strengths is the ease with which it deploys large-scale
workflows to HPC clusters. Let's take a look at the important pieces!

## Profiles
Using profiles is critical to making Snakemake run large-scale jobs on Biowulf.
In your home directory, go ahead and clone Biowulf's Snakemake profile
```
cd ~
`git clone https://github.com/NIH-HPC/snakemake_profile.git
cd -
```

This is greatly reduces the work needed to scale jobs on Biowulf. The profile contains:
- Resource tags you can use in your rules to specify resource requirements
- Automatic assignment of the tmpdir (useful if you request /lscratch space)
- Settings to aid in job canceling
- Restrictions on job submission numbers so you don't go too crazy and overwhelm the system
- Sensible default settings for Singularity
- Automatic mandates to continue where Snakemake left off

Before submitting anything to Biowulf, it's important to be aware that all environment variables, modules, etc., are passed 
to the Snakemake master job. Therefore, whatever is needed to parse and run the Snakefile must be loaded.
```
module load snakemake
module load singularity
module load conda
```

