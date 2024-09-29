On the NIH HPC systems, start an interactive session, load the snakemake and
singularity modules, and clone this repository:

```console
user@headnode> sinteractive --cpus-per-task=12 --mem=24g --gres=lscratch:20
...
user@cn1234> ## change to a suitable directory somewhere in /data
user@cn1234> cd /data/$USER
user@cn1234> module load git
user@cn1234> git clone https://github.com/NIH-HPC/snakemake-class.git
user@cn1234> cd snakemake-class
```



We will use some of these items for exercises during class.

Setup will be different on other systems.

When you are ready run the setup script to fetch data and create all files
necessary for the exercies.

```console
user@cn1234> nohup ./setup.sh
...
