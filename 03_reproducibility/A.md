# Environments Support Reproducibility

Suppose you’re conducting an RNA-Seq project using the ROSMAP dataset. Reproducibility means that:

The data set you used is publicly available or properly documented.
You share the exact scripts that you used for processing the files (e.g., removing adapters, running DESeq2).
You document your software environment (e.g., R version 4.4, specific libraries such DESeq2 Release (3.19)).
When someone else runs your code with the same data and environment, they get the same results.
By doing this, your analysis is not just a "black box" but something others can verify, extend, or adapt.

## Why are environments critical to reproducibility?

Bioinformatics projects depend on multiple interconnected libraries. Sometimes these libraries have conflicting requirements (e.g., one library needs a specific version of Python or another library). Without a controlled environment, managing these dependencies manually can be a nightmare, leading to version conflicts or broken installations. Environments like Conda The trick is to isolate each project's dependencies. This isolation prevents issues like "dependency hell," where upgrading one package for one project breaks another project that relies on an older version of the same package.

It's not just the code or the libraries that matter but the entire software, and even hardware, stack, including compilers, interpreters (e.g., specific versions of Python or R), system libraries, and even hardware dependencies (like CUDA for GPU-based computation). 

Especially given the hardware and firmware dependency, some environment types are better than others. Software environments allow users to specify and capture the entire stack required for the project. This ensures that the code runs in the same environment, even if your local machine changes or the software ecosystem evolves.

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

To use a singularity container simply add it to your rule:
```
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    singularity:
        "library://wresch/classes/rnaseq:0.6"
    shell: "code"
```
If it is a Docker container you found somewhere online, simply take the user and package name listed in the pull command and prepend
the path with 'docker://'
<img width="572" alt="Screenshot 2024-09-30 at 8 36 50 PM" src="https://github.com/user-attachments/assets/553bb8c2-4a81-453d-aec5-6df5379787d4">
```
rule hisat2:
    input: fq = "00fastq/{sample}.fastq.gz",
           idx = "00ref/hisat_index/R64-1-1"
    output: bam = "02aln/{sample}.bam",
            bai = "02aln/{sample}.bam.bai"
    singularity: "docker://dceoy/hisat2"
    shell: "code"
```

Depending on how the container is organized, you may have to fiddle with how exactly the directories inside
the container are bound in relation to directories outside. That can be addressed with the `--singularity-args` when running snakemake, as in the following code.
`--singularity-prefix` simply determines where the container that is downloaded is stored or where to search for a container if you already have one downloaded.
 ```
 snakemake --cores=4 --use-singularity \
     --singularity-args '-B $PWD:/data --pwd /data' \
     --singularity-prefix=00container
```

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
    shell: "code"
```
Notice the exact version is specified -- Biowulf may change the default version at any time.

Unfortunately, workflows using that use Biowulf modules cannot be published to the wider community. Further, Biowulf may not have the tool
you need.

## Conda

Conda/Mamba is a package manager that allows reproducible implementations of a workflow
  - on the same machine
  - if the user has been diligent enough to make a fully pinned env.yml and not to 
     publish environment states with each run. 

Conda does some of the work for you. It tracks changes to the environment, so you can see when it has changed if you think to look.

However, there are a couple big problems with this system:
    1. it is easy to make changes to the environment without noticing or remembering, especially years later.
    2. exact versions can only be guaranteed on the same system

Using Snakemake helps with both of these problems. 

Snakemake provides an easy way to deal with the first problem. 
The second problem can only be dealt with using Singularity containers.


For even greater reproducibility it is always recommended to pair whatever you do with git.
Any environment specification you use, or an output file detailing the exact version of your current 
environment is recommended to be committed and labeled alongside each publication-worthy analysis you
do. So essentially any time you make a change, you should commit any config files or version outputs to git.



IDEA: show an example of a fully pinned env specification.

To do that, consider making a rule at the end that triggers a git commit any time one of your environments changes,
or simply make sure those files are tracked by git.







