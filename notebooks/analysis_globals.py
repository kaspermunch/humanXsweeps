
from pathlib import Path
import pandas as pd

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


# mutation rate per year
mut_per_year = 4.3e-10 

# generation time
gen_time = 29

par1_start, par1_end = 60001, 2699520
par2_start, par2_end = 154931044, 155260560

#############################

# # populations (read from pop_names.tsv)
# g1000_populations = [
#     'CHS', #	Southern Han Chinese, China
#     'MXL', #	Mexican Ancestry in Los Angeles, California
#     'ACB', #	African Caribbean in Barbados
#     'CHB', #	Han Chinese in Bejing, China
#     'IBS', #	Iberian populations in Spain
#     'KHV', #	Kinh in Ho Chi Minh City, Vietnam
#     'ASW', #	African Ancestry in Southwest US
#     'ITU', #	Indian Telugu in the UK
#     'MSL', #	Mende in Sierra Leone
#     'TSI', #	Toscani in Italy
#     'GBR', #	British in England and Scotland
#     'BEB', #	Bengali in Bangladesh
#     'STU', #	Sri Lankan Tamil in the UK
#     'CEU', #	Utah residents with Northern and Western European ancestry
#     'CDX', #	Chinese Dai in Xishuangbanna, China
#     'JPT', #	Japanese in Tokyo, Japan
#     'ESN', #	Esan in Nigeria
#     'PJL', #	Punjabi in Lahore,Pakistan
#     'CLM', #	Colombian in Medellin, Colombia
#     'FIN', #	Finnish in Finland
#     'LWK', #	Luhya in Webuye, Kenya
#     'PEL', #	Peruvian in Lima, Peru
#     'YRI', #	Yoruba in Ibadan, Nigeria
#     'PUR', #	Puerto Rican in Puerto Rico
#     'GIH', #	Gujarati Indian in Houston,TX
#     'GWD', #	Gambian in Western Division, The Gambia
# ]

g1000_pop_info = pd.DataFrame.from_records([
    ('CHB',	'Han Chinese in Beijing, China', 'EAS'),
    ('JPT',	'Japanese in Tokyo, Japan', 'EAS'),
    ('CHS',	'Southern Han Chinese', 'EAS'),
    ('CDX',	'Chinese Dai in Xishuangbanna, China', 'EAS'),
    ('KHV',	'Kinh in Ho Chi Minh City, Vietnam', 'EAS'),
    ('CEU',	'Utah Residents (CEPH) with Northern and Western European Ancestry', 'EUR'),
    ('TSI',	'Toscani in Italia', 'EUR'),
    ('FIN',	'Finnish in Finland', 'EUR'),
    ('GBR',	'British in England and Scotland', 'EUR'),
    ('IBS',	'Iberian Population in Spain', 'EUR'),
    ('YRI',	'Yoruba in Ibadan, Nigeria', 'AFR'),
    ('LWK',	'Luhya in Webuye, Kenya', 'AFR'),
    ('GWD',	'Gambian in Western Divisions in the Gambia', 'AFR'),
    ('MSL',	'Mende in Sierra Leone', 'AFR'),
    ('ESN',	'Esan in Nigeria', 'AFR'),
    ('ASW',	'Americans of African Ancestry in SW USA', 'AFR'),
    ('ACB',	'African Caribbeans in Barbados', 'AFR'),
    ('MXL',	'Mexican Ancestry from Los Angeles USA', 'AMR'),
    ('PUR',	'Puerto Ricans from Puerto Rico', 'AMR'),
    ('CLM',	'Colombians from Medellin, Colombia', 'AMR'),
    ('PEL',	'Peruvians from Lima, Peru', 'AMR'),
    ('GIH',	'Gujarati Indian from Houston, Texas', 'SAS'),
    ('PJL',	'Punjabi from Lahore, Pakistan', 'SAS'),
    ('BEB',	'Bengali from Bangladesh', 'SAS'),
    ('STU',	'Sri Lankan Tamil from the UK', 'SAS'),
    ('ITU',	'Indian Telugu from the UK', 'SAS'),
], columns=['population', 'description', 'superpop'])
g1000_pop_info['superpop'] = pd.Categorical(g1000_pop_info.superpop, 
                                            categories=['AFR', 'EUR', 'SAS', 'EAS', 'AMR'], 
                                            ordered=True)
g1000_pop_info.sort_values('superpop', inplace=True)
g1000_pop_info.reset_index(inplace=True)

g1000_min_sweep_clade_proportion = 0.3
# could also use 0.3 which would correspond to 42 / 140 (min_clade_size of 42 in simons out of 140 nr non-Africans in simons)