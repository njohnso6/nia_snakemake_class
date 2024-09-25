**Reproducibility using Envmodules, Conda and Singularity**

There are several ways to ensure reproducibility of workflows. Some are better than others.
1. Use the system modules. We'll go over how to do that.
2. Use a package manager such as Conda (please use the mamba solver, for your own sanity)
3. Use Singularity

Using system modules helps maintain a low system footprint and makes things easy.
However, you never know when Biowulf will change versions. If using the system versions,
be sure to specify the version you want to maintain reproducibility. 

Any workflows using these cannot be published for the wider community, and Biowulf may not have the tool
you need.

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

Singularity containers are fully contained systems including the operating system and allow 
virtualization of basically any computation. Containers can be made with all sorts of software,
and they can be guaranteed not to change without explicit updating and therefore allow better versioning.

Using publicly made singularity containers, or importing Docker containers and running with singularity
is highly recommended for the sake of reprodcibility, especially if you plan to share a pipeline. If one
is available, it should generally be the first choice.

For even more reproducibility it is always recommended to pair whatever you do with git.
Any environment specification you use, or an output file detailing the exact version of your current 
environment is recommended to be committed and labeled alongside each publication-worthy analysis you
do. So essentially any time you make a change, you should commit any config files or version outputs to git.

IDEA: show an example of a fully pinned env specification.

To do that, consider making a rule at the end that triggers a git commit any time one of your environments changes,
or simply make sure those files are tracked by git.

Exercise:
  1. Write the following workflow with three steps. I have given a template.
  2. 
  It should be a matter of adding the correct container, envmodule, and conda yaml requirements.
  I've thrown in a couple wrenches that you will need to fix.
  IDEA: Make the conda one require python 2 instead of 3, so they will have to specify which version.
  Perhaps do a quick cutadapt module one. Just have it check whether the version output is correct. 





