
from pathlib import Path


## This file contains paths and values that apply globally across notebooks ###


## Directories ################################################################

# Analysis root dir
root_dir = Path('/home/kmt/simons/faststorage/people/kmt')

# scripts dir
scripts_dir = root_dir / 'scripts'

# meta data dir
meta_data_dir = Path('/home/kmt/simons/faststorage/data/metadata')

# gwf steps dir
steps_dir = root_dir / 'steps'

# argweaver raw results
argweaver_dir = steps_dir / 'argweaver/output'

# results dir (for hdf files from nobooks and such)
results_dir = root_dir / 'results'

# figures from notebooks
figures_dir = root_dir / 'figures'

# analysis data (not produced in this analysis)
data_dir = root_dir / 'data'

# results from great ape ils analysis
great_ape_ils_analysis_dir = Path('/home/kmt/projects/great_ape_ils_maps/analyses')

# pairwise distances between individuals from the same populations (pi), all chromosomes
pi_dir = steps_dir / 'pi_stores'

# pairwise distances between 
dist_dir = steps_dir / 'dist_stores'

# pairwise distances between male chrX haplotypes
male_x_haploid_dir = steps_dir / 'male_x_haploids'

## Constants ##################################################################

# # exclude sub-Saharan African populations:
# excluded_pops = ['Masai', 'Somali', # show some non-African component in Structure plot in main paper (Laurits does that too)
#                  'Mozabite', 'Saharawi', # known non-african ancestry and Neanderthal admxiture
#                  'Luhya', 'Luo', 'BantuKenya', # 1000 genomes study sugggest that Kenyans (Luhya) have some backflow. (Laurits info)
#                  'Gambian'] # Not sure what the argument is here, but Laurits excludes it too...

# maximum uncalled bases in a 100kb window
max_uncalled_bases = 50000
ust_ishim_max_uncalled_bases = 70000

# pairwise distance between two male windows to call it low distance
#pwdist_cutoff = 6e-05
pwdist_cutoff = 5e-05

# maximum proportion of individual haplotypes missing in a window
max_prop_indiv_missing = 0.1

# minimum size a clade that a haplotype is part of.
#min_sweep_clade_size = 35 # corresponds to 20% of all males in data set
min_sweep_clade_size = 42 # corresponds to 25% of all males in data set

# minimum number of consequtive windows 100kb scale
min_run_length = 5

# minimum number of consequtive windows 1Mb scale
min_run_length_Mb_scale = 1

# minimum proportion of individuals swept for a region to be conservatively called as swept
min_prop_swept = 0.5

# minimum proportion of individuals swept for a region to be included in peaks
peak_min_prop_swept = 0.25

