

On the NIH HPC systems, start an interactive session and load the modules we will need for this class.
Then, clone this repository:
```console
sinteractive --cpus-per-task=12 --mem=24g --gres=lscratch:20
user@cn1234> module load git
user@cn1234> module load python
user@cn1234> module load snakemake
user@cn1234> module load singularity 
user@cn1234> ## change to a suitable directory somewhere in /data:
user@cn1234> cd /data/$USER
user@cn1234> module load git
user@cn1234> git clone https://github.com/njohnso6/nia_snakemake_class.git`
user@cn1234> cd nia_snakemake_class
```

Normally, I would say it's not cool to hog all those resources for this but it will make our class easier when things come up and so you can experiment during class.

Then, create a separate folder called exercises. Enter it and
clone the exercises there. These will supplement the exercises 
we run in class.

```console
user@cn1234> mkdir exercises
user@cn1234> cd exercises/
user@cn1234> git clone https://github.com/NIH-HPC/snakemake-class.git
user@cn1234> cd snakemake-class
```
We will use some of these items for exercises during class.

Setup will be different on other systems.

When you are ready run the setup script to fetch data and create all files
necessary for the exercies.

```console
user@cn1234> nohup ./setup.sh
```


