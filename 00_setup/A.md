On the NIH HPC systems, start an interactive session 

```console
sinteractive --cpus-per-task=12 --mem=24g --gres=lscratch:20
```

and load the modules we will need for this class.

```console
module load git
module load python
module load snakemake
module load singularity 
```

Then, clone this repository after changing to a suitable directory somewhere in /data:

```console
cd /data/$USER
module load git
git clone https://github.com/njohnso6/nia_snakemake_class.git
cd nia_snakemake_class
```

Normally, I would say it's not cool to hog all those resources for this but it will make our class easier when things come up and so you can experiment during class.

Then, create a separate folder called exercises. Enter it and
clone the exercises there. These will supplement the exercises 
we run in class.

```console
mkdir exercises
cd exercises/
git clone https://github.com/NIH-HPC/snakemake-class.git
cd snakemake-class
```

We will use some of these items for exercises during class.

Setup will be different on other systems.

When you are ready run the setup script to fetch data and create all files necessary for the exercise.

```console
nohup ./setup.sh
```
