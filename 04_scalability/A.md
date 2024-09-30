# Scaling on Biowulf

It's finally time to really kick Snakemake into high gear! One of its great strengths is the ease with which it deploys large-scale
workflows to HPC clusters. Let's take a look at the important pieces!

## Profiles

It can get pretty pretty cumbersome to re-type all the parameters for snakemake. Snakemake profiles can set default values for command line parameters. A profile is a directory containing at a minimum a config.yaml with keys corresponding to command line flags. A profile specified by name on the command line is search for in $HOME/.config/snakemake. Alternatively profiles can be specified by path: 
`snakemake --cores 8 --profile ~/snakemake_profile`

There is a Biowulf-specific profile that contains all the biolerplate code to make deploying large-scale workflows to Biowulf simple.

In your home directory, go ahead and clone Biowulf's Snakemake profile
```
cd ~
`git clone https://github.com/NIH-HPC/snakemake_profile.git
cd -
```

This is greatly reduces the work needed to scale jobs on Biowulf. The profile contains:
- Resource tags you can use in your rules to specify resource requirements (mem, lscratch, threads, partition, etc.)
- Automatic assignment of the tmpdir (useful if you request /lscratch space)
- Settings to aid in job canceling
- Restrictions on job submission numbers so you don't go too crazy and overwhelm the system
- Sensible default settings for Singularity
- Automatic mandates to continue where Snakemake left off

Before submitting anything to Biowulf, it's important to be aware that all environment variables, modules, etc., are passed 
to the Snakemake master job. Therefore, whatever is needed to parse and run the Snakefile must be loaded.
```snakemake
module load snakemake
module load singularity
module load conda
```

## Submission scripts

I like to creat submission scripts, as there are often multiple arguments in addition to the profile I want to pass. These can be put in the profile as well, but sometimes I don't want to have something apply to everything. 

```bash
snakemake --use-singularity --use-envmodules --use-conda -S testing_snakefile
```

