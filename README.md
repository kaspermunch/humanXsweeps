# Code used the for the Cell Genomics paper

**Extraordinary selection on the human X chromosome associated with archaic admixture**

Analysis of X chromosomes from the Simons Diversity Panel. Once cloned the repository folder needs to be renamed to `kmt` to fit file paths used in the workflows. Also the following subfolders must be created in the repository folder:

    mkdir steps results figures 

## Conda environments

### simons

This is the main py3 environment for the workflows. You can create it like this:

    conda env create -f simons.yml
    
### simons_jupyter

This is the main py3 environment for jupyter analysis. You can create it like this:

    conda env create -f simons_jupyter.yml    

### argweaver

Separate py2 environment that is only used from within the GWF workflow (see below). You can create it like this:

You can create it like this:

    conda env create -f argweaver.yml
    
This has argweaver executables and argweaver python lib (must be downloaded and installed using pip). It requires compbio which is also installed the same way.

For some reason I hade to compile these executables manually to keep them from coredumping:

    g++ argweaver/src/argweaver/*.cpp argweaver/src/smc2bed.cpp -o ~/anaconda3/envs/argweaver/bin/smc2bed
    g++ argweaver/src/argweaver/*.cpp argweaver/src/arg-summarize.cpp -o ~/anaconda3/envs/argweaver/bin/arg-summarize
    cd argweaver ; pip install .
    
Note that you need to chagne the option in the workflow from --resume to --overwrite if you want to redo the analysis

If resampling is aborted the stats file in each sampling dir may be corrupt. Delete the last sample (here 2000):

    find steps/argweaver/samples -name '*2000.smc.gz' -exec rm {} \;

and roll back the stats file to the correct iteration (here 1900):

    find steps/argweaver/samples -name '*[01].stats' | python scripts roll_back_argweaver_stats_files.py 1900    
    
## GWF workflows

Main analysis of SGDP:

    gwf -f workflow_simons.py run

Additional analysis of 1000 genomes

    gwf -f workflow_1000genomes.py run

## Data analysis and visualization using jupyter notebooks

All notebooks are listed in the folder notebooks. They should be run in the order they are numbered. Some of them require a *lot* of memory. For easy use of jupyter analysis on a slurm computing cluster use [slurm_jupyter](https://github.com/kaspermunch/slurm-jupyter).
