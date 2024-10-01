# Environments Support Reproducibility

Suppose you’re conducting an RNA-Seq project using the ROSMAP dataset. Reproducibility means that:

The data set you used is publicly available or properly documented.
You share the exact scripts that you used for processing the files (e.g., removing adapters, running DESeq2).
You document your software environment (e.g., R version 4.4, specific libraries such DESeq2 Release (3.19)).
When someone else runs your code with the same data and environment, they get the same results.
By doing this, your analysis is not just a "black box" but something others can verify, extend, or adapt.

## Why are environments critical to reproducibility?

Bioinformatics projects depend on multiple interconnected libraries. Sometimes these libraries have conflicting requirements (e.g., one library needs a specific version of Python or another library). Without a controlled environment, managing these dependencies manually can be a nightmare, leading to version conflicts or broken installations. The trick is to isolate each project's dependencies. This isolation prevents issues like "dependency hell," where upgrading one package for one project breaks another project that relies on an older version of the same package.

It's not just the code or the libraries that matter but the entire software, and even hardware, stack, including compilers, interpreters (e.g., specific versions of Python or R), system libraries, and even hardware dependencies (like CUDA for GPU-based computation). 

Especially given the hardware and firmware dependency, some environment types are better than others. Software environments allow users to specify and capture the entire stack required for the project. This ensures that the code runs in the same environment, even if your local machine changes or the software ecosystem evolves.

Over time, software ecosystems change. Libraries get deprecated, new versions come out, or the original environment may no longer be supported. Without environments, if someone tries to reproduce an analysis or run code years later, they may run into compatibility issues. With a defined environment, the code can still run, as it encapsulates all the necessary software in a "time capsule" of sorts.

For instance, by storing the exact environment configuration file (e.g., environment.yml for Conda or requirements.txt for Python), someone can recreate the same setup years later assuming the machine is capable of running those older versions.

## Reproducibility using Envmodules, Conda and Singularity

There are several ways to integrate reproducible environments into Snakemake:
1. Use Singularity containers
2. Use Biowulf modules
3. The Conda package manager
4. Other package managers (less explicit support for this)

## Singularity Containers

![image](https://github.com/user-attachments/assets/bfa52abe-aa4c-4baa-8dd7-08967aeb5165)
Image taken from: https://stephen-odaibo.medium.com/docker-containers-python-virtual-environments-virtual-machines-d00aa9b8475

In most cases, Singularity environments are the ideal environment on Biowulf. Whereas a virtual environment like Conda
must be used on the same operating system on which it's built, a container can 
be used on different operating systems and is therefore more flexible. 
The exact version of R might change slightly when installed directly on a machine to match the operating system, but a container will contain the same version regardless of the machine it's on.
Advantages:
1. Can be passed to collaborators on other computer systems with the exact same version
2. More or less immutable--prevents unexpected changes

Using publicly made singularity containers, or importing Docker containers and running with singularity
is highly recommended for the sake of reprodcibility, especially if you plan to share a pipeline. If one
is available, it should generally be the first choice.

Before using Singularity with Snakemake, you must have it in your path either in your environment
or using `module load singularity`

Let's take a look at our exercises in exercise01.

To use a singularity container simply add it to your rule:

You should see the following rule in 
```
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    singularity:
        "library://wresch/classes/rnaseq:0.8"
    shell: "Lorum ipsum code do not run"
```
If it is a Docker container you found somewhere online, simply look for the pull command and take the user and package name listed. Prepend them with `'docker://'` and add the result to the singularity item in your rule.


<img width="572" alt="Screenshot 2024-09-30 at 8 36 50 PM" src="https://github.com/user-attachments/assets/553bb8c2-4a81-453d-aec5-6df5379787d4">

This is an example. It is not the correct version for our samples so do not use this docker.
```
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    singularity: "docker://dceoy/hisat2"
    shell: "Lorum ipsum code do not run"
```

Before running run:
```bash
source /usr/local/current/singularity/app_conf/sing_binds
```

This should prevent problems by giving the container the interface points it needs with your operating system and your file system.

However, depending on how the container is organized, you may have to fiddle with how exactly the directories inside
the container are bound in relation to directories outside. That can be addressed with the `--singularity-args` when running snakemake, as in the following code.
`--singularity-prefix` simply determines where the container that is downloaded is stored or where to search for a container if you already have one downloaded.
Here is an example.
 ```
 snakemake --cores=4 --use-singularity \
     --singularity-args '-B $PWD:/data --pwd /data' \
     --singularity-prefix=00container
```

However, generally `--singularity-args` shouldn't be necessary.

Let's do it together!
Go to our `exercises` folder and exercise01.

Take a look at Snakefile.finished

Notice how the container is listed with `library://`. That means it came from the [sylabs library](https://cloud.sylabs.io/library/wresch/classes/rnaseq)

**To run**

```bash
 snakemake --cores=4 --use-singularity  --singularity-prefix=00container -s Snakefile.finished
```

Now go ahead and run 
```bash
 snakemake --cores=4 -s Snakefile.finished clean
```
so we can use the Snakefile again.


## Using Biowulf modules
Using system modules helps maintain a low system footprint and makes things easy.
These are, in fact, singularity modules under the hood.

To use the Biowulf modules, simply add an envmodules section to your rule:
```
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    envmodules:
            "hisat2/2.2.1"
    shell: "Loreum ispum code do not run"
```
Notice the exact version is specified -- Biowulf may change the default version at any time, so it's smart to specify the version so it will break if it's not there.

Unfortunately, workflows using that use Biowulf modules cannot be published to the wider community. Further, Biowulf may not have the tool
you need.

Let's modify our Snakefile.final in exercise02 again.
Replace the singularity requirement with envmodules.

```bash
 snakemake --cores=4 --use-singularity  --singularity-prefix=00container -s Snakefile.finished
```

## Conda

Conda/Mamba is a package manager that allows reproducible implementations of a workflow

It is not as ideal as Singularity because:
- an environment specified in conda becomes less portable as it becomes more specific (e.g. Python=3.8 instead of Python = 3)
- the workflow may work differently because of less specific package versions

Let's take a look at the file in this folder called fully_pinned_example.yml

Such a spcific environment would be unlikely to build on a machine other than Biowulf

Snakemake allows environments for each rule to be specified, if desired.
There are other ways to use Conda in Snakemake but this is the most recommended way becuase:
- a written record exists for each rule
- this allows the environment specifications to be version-controlled
- including only the necessary packages for each rule increases the liklihoodmore specific versions can be used

Go back to exercise01 in the exercise folder

Replace `envmodules` with `conda: "hisat-environment.yml"`

Then create a file called "hisat-environment.yml"

Place the following in it:
```yml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3
  - hisat2
  - samtools
```

Run with 
```bash
snakemake -s Snakefile.finished --use-conda --cores 2 --conda-prefix hisat2-env
```

Notice the conda prefix in this case. This tells snakemake not to use the system conda modules but rather to download its own and place them in our current folder.


For even greater reproducibility it is always recommended to pair whatever you do with git.
Any environment specification you use, or an output file detailing the exact version of your current 
environment is recommended to be committed and labeled alongside each publication-worthy analysis you
do. So essentially any time you make a change, you should commit any config files or version outputs to git.








