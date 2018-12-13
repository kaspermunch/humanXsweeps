# humanXsweeps

Code for analysis of X chromosomes from the Simons Diversity Panel


## Environments

### simons

This is the main py3 environment for scripts and jupyter analysis. You can create it like this:

    conda env create -f simons.yml

### argweaver

Separate py2 environment that is only used from within the GWF workflow (see below). You can create it like this:

    You can create it like this:
    
This has argweaver executables and argweaver python lib (must be downloaded and installed using pip). It requires compbio which is also installed the same way.

For some reason I hade to compile these executables manually to keep them from coredumping:

    g++ argweaver/src/argweaver/*.cpp argweaver/src/smc2bed.cpp -o ~/anaconda2/envs/argweaver/bin/smc2bed
    g++ argweaver/src/argweaver/*.cpp argweaver/src/arg-summarize.cpp -o ~/anaconda2/envs/argweaver/bin/arg-summarize



## GWF workflow

    gwf run


## Data analysis using jupyter notebooks

