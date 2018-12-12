# humanXsweeps
Code for analysis of X chromosomes from the Simons Diversity Panel



# Environments

## simons

This is the main py3 environment for scripts and jupyter analysis

## argweaver

Separate py2 environment with argweaver executables and argweaver python lib (downloaded and installed using pip). This requires compbio which is also installed the same way.

For some reason I hade to compile these executables manually to keep them from coredumping:

    g++ argweaver/src/argweaver/*.cpp argweaver/src/smc2bed.cpp -o ~/anaconda2/envs/argweaver/bin/smc2bed
    g++ argweaver/src/argweaver/*.cpp argweaver/src/arg-summarize.cpp -o ~/anaconda2/envs/argweaver/bin/arg-summarize


