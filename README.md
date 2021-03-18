# humanXsweeps

Code for analysis of X chromosomes from the Simons Diversity Panel

Create directory structure for analysis:

    mkdir steps results figures 

## Environments

### simons

This is the main py3 environment for scripts and jupyter analysis. You can create it like this:

    conda create --name simons -c gwforg -c micknudsen -c etetoolkit -c anaconda -c conda-forge -c bioconda python=3 gwf gwf-utilization ete3 biopython ipyparallel jupyter jupyterlab matplotlib mpld3 nbconvert numpy pandas pytables scikit-learn scipy seaborn scikit-learn statsmodels pyfaidx scikit-bio mygene msprime openblas descartes basemap-data-hires basemap cartopy networkx mygene psutil bitarray

    conda activate simons

    pip install pyslim

    git clone ssh://git@github.com/kaspermunch/ChromosomeWindows.git
    cd ChromosomeWindows ; pip install .

    git clone ssh://git@github.com/kaspermunch/GenomicWindows.git
    cd GenomicWindows ; pip install .

The libraries genomicintervals and genomicwindows are my own and have their own git repositories.

### argweaver

Separate py2 environment that is only used from within the GWF workflow (see below). You can create it like this:

You can create it like this:

    conda create --name argweaver -c anaconda python=2.7 libgfortran numpy scipy

This has argweaver executables and argweaver python lib (must be downloaded and installed using pip). It requires compbio which is also installed the same way.

For some reason I hade to compile these executables manually to keep them from coredumping:

    g++ argweaver/src/argweaver/*.cpp argweaver/src/smc2bed.cpp -o ~/anaconda3/envs/argweaver/bin/smc2bed
    g++ argweaver/src/argweaver/*.cpp argweaver/src/arg-summarize.cpp -o ~/anaconda3/envs/argweaver/bin/arg-summarize
    cd argweaver ; pip install .

## GWF workflows

Main analysis of SGDP

    gwf -f workflow_simons.py run

Additional analysis of 1000 genomes

    gwf -f workflow_1000genomes.py run

Additional analysis of WestEurasian individuals using Clues

    gwf -f workflow_clues.py run

## Data analysis using jupyter notebooks

All notebooks are listed in the folder notebooks. They should be run in the order they are numbered. 

For easy use of notebook on a slurm computing cluster use [slurm_jupyter]()
