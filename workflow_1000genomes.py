
from gwf import Workflow
import sys, os, glob, itertools, re
from collections import defaultdict
#from pathlib import Path

from random import seed
seed(42)

sys.path.append('./scripts')
sys.path.append('./notebooks')

import simons_meta_data
import hg19_chrom_sizes

import analysis_globals

from templates import *

gwf = Workflow(defaults={'account': 'simons'})


#################################################################################
# Dirs
#################################################################################

faststorage = '/home/kmt/simons/faststorage'
mydir = os.path.join(faststorage, 'people', 'kmt')



#################################################################################
# Dict of 1000 genome VCF files
#################################################################################

g1000_vcf_files = dict()

for chrom in map(str, range(1, 23)):
    g1000_vcf_files[chrom] = os.path.join(faststorage, 'data', '1000Genomes',
                                'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom))

g1000_vcf_files['X'] = os.path.join(faststorage, 'data', '1000Genomes',
                                'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz')
for f in g1000_vcf_files.values():
    assert os.path.exists(f), f

#################################################################################

# Extracted sample names from VCF (sample_names.txt)
#     zcat data/1000Genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz | head -n 10000 | grep CHROM | perl -pe 's/\s+/\n/g' > sample_names.txt

# Converted metainfo excel file to csv (sample_info.csv)

# write files with sample names devided by population and sex (also writes pop_names.tsv):
# python scripts/write_1000gen_meta_info.py ../../data/1000Genomes/metainfo/sample_names.txt  ../../data/1000Genomes/metainfo/sample_info.csv  ../../data/1000Genomes/metainfo

# populations (read from pop_names.tsv)
g1000_populations = [
    'CHS', #	Southern Han Chinese, China
    'MXL', #	Mexican Ancestry in Los Angeles, California
    'ACB', #	African Caribbean in Barbados
    'CHB', #	Han Chinese in Bejing, China
    'IBS', #	Iberian populations in Spain
    'KHV', #	Kinh in Ho Chi Minh City, Vietnam
    'ASW', #	African Ancestry in Southwest US
    'ITU', #	Indian Telugu in the UK
    'MSL', #	Mende in Sierra Leone
    'TSI', #	Toscani in Italy
    'GBR', #	British in England and Scotland
    'BEB', #	Bengali in Bangladesh
    'STU', #	Sri Lankan Tamil in the UK
    'CEU', #	Utah residents with Northern and Western European ancestry
    'CDX', #	Chinese Dai in Xishuangbanna, China
    'JPT', #	Japanese in Tokyo, Japan
    'ESN', #	Esan in Nigeria
    'PJL', #	Punjabi in Lahore,Pakistan
    'CLM', #	Colombian in Medellin, Colombia
    'FIN', #	Finnish in Finland
    'LWK', #	Luhya in Webuye, Kenya
    'PEL', #	Peruvian in Lima, Peru
    'YRI', #	Yoruba in Ibadan, Nigeria
    'PUR', #	Puerto Rican in Puerto Rico
    'GIH', #	Gujarati Indian in Houston,TX
    'GWD', #	Gambian in Western Division, The Gambia
]


#################################################################################
# Write gzipped phased haplotypes for all indiviudals 
#################################################################################

# steps dir for all 1000genomes:
g1000_dir = os.path.join(mydir, 'steps', '1000genomes')
if not os.path.exists(g1000_dir):
    os.makedirs(g1000_dir)

autosomes = list(map(str, range(1, 23)))

#################################################################################
# Write a maked href
#################################################################################

g1000_callability_mask_dir = os.path.join(faststorage, 'data/1000Genomes/callabilitymask')
g1000_callability_mask_files = [os.path.join(g1000_callability_mask_dir,
                '20141020.chr{}.strict_mask.fasta.gz'.format(chrom)) for chrom in autosomes + ['X']]
human_reference = os.path.join(faststorage, 'data/cteam_lite_public3/FullyPublic/Href.fa')

# dir for male haplotypes
g1000_masked_ref_dir = os.path.join(mydir, 'steps', '1000genomes', 'masked_ref')
if not os.path.exists(g1000_masked_ref_dir):
    os.makedirs(g1000_masked_ref_dir)

for i, mask_file in enumerate(g1000_callability_mask_files):

    masked_ref = modpath(mask_file, parent=g1000_masked_ref_dir, suffix='.fa')

    #g1000_masked_reference = os.path.join(faststorage, 'steps', '1000genomes', 'masked_ref', 'masked_reference.fa')

    g = gwf.target("mask_reference_g1000_{}".format(i), inputs=[human_reference, mask_file], 
        outputs=[masked_ref], 
        memory='16g', walltime='01:00:00') << """

        source activate simons
        python scripts/1000gen_masked_href.py {href} {mask} {output}

    """.format(href=human_reference, mask=mask_file, output=masked_ref)

#################################################################################
# Write phased haplotypes for all male X 
#################################################################################

# read male sample names:
g1000_males_file = os.path.join(faststorage, 'data', '1000Genomes', 'metainfo', 'all_males.txt')
with open(g1000_males_file) as f:
    g1000_males = f.read().strip().split()


# dir for male haplotypes
g1000_male_haplo_dir = os.path.join(mydir, 'steps', '1000genomes', 'haplotypes', 'males')
if not os.path.exists(g1000_male_haplo_dir):
    os.makedirs(g1000_male_haplo_dir)

g1000_male_x_haplotype_files = list()
for chrom in ['X']:
    chrom_dir = os.path.join(g1000_male_haplo_dir, chrom)
    if not os.path.exists(chrom_dir): 
        os.makedirs(chrom_dir)
    for sample_id in g1000_males:
        haplo_file1 = modpath("{}_{}-A.fa".format(sample_id, chrom), parent=chrom_dir)
        g1000_male_x_haplotype_files.append(haplo_file1)

        gwf.target_from_template("vcf2haplo_male_{}_{}".format(chrom, sample_id),
            vcf2haplo(vcf_file=g1000_vcf_files[chrom], masked_ref=g1000_masked_reference,
            sample_id=sample_id, 
            out_file1=haplo_file1))


#################################################################################
# Write phased haplotypes for all female X and auto
#################################################################################

# read male sample names:
g1000_females_file = os.path.join(faststorage, 'data', '1000Genomes', 'metainfo', 'all_females.txt')
with open(g1000_females_file) as f:
    g1000_females = f.read().strip().split()


# dir for male haplotypes
g1000_female_haplo_dir = os.path.join(mydir, 'steps', '1000genomes', 'haplotypes', 'females')
if not os.path.exists(g1000_female_haplo_dir):
    os.makedirs(g1000_female_haplo_dir)

g1000_female_haplotype_files = list()
for chrom in autosomes + ['X']:
    chrom_dir = os.path.join(g1000_female_haplo_dir, chrom)
    if not os.path.exists(chrom_dir): 
        os.makedirs(chrom_dir)
    for sample_id in g1000_females:
        haplo_file1 = modpath("{}_{}-A.fa".format(sample_id, chrom), parent=chrom_dir)
        g1000_female_haplotype_files.append(haplo_file1)
        haplo_file2 = modpath("{}_{}-B.fa".format(sample_id, chrom), parent=chrom_dir)
        g1000_female_haplotype_files.append(haplo_file2)

        # assert haplo_file1 != '/' and haplo_file2 != '/'
        # print(haplo_file1, haplo_file2)

        gwf.target_from_template("vcf2haplo_female_{}_{}".format(chrom, sample_id),
            vcf2haplo(vcf_file=g1000_vcf_files[chrom], masked_ref=g1000_masked_reference,
            sample_id=sample_id, 
            out_file1=haplo_file1, out_file2=haplo_file2))








#################################################################################
# 
#################################################################################







#################################################################################
# 
#################################################################################
