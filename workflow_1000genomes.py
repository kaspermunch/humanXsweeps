
from gwf import Workflow
import sys, os, glob, itertools, re
from collections import defaultdict
import numpy as np

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
# Meta information
#################################################################################

"""
Extracted sample names from VCF (sample_names.txt)

zcat data/1000Genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz | head -n 10000 | grep CHROM | perl -pe 's/\s+/\n/g' > sample_names.txt

Converted metainfo excel file to csv (sample_info.csv)

write files with sample names devided by population and sex (also writes pop_names.tsv):
python scripts/write_1000gen_meta_info.py ../../data/1000Genomes/metainfo/sample_names.txt  ../../data/1000Genomes/metainfo/sample_info.csv  ../../data/1000Genomes/metainfo
"""

# read male sample names:
g1000_males_file = os.path.join(faststorage, 'data', '1000Genomes', 'metainfo', 'all_males.txt')
with open(g1000_males_file) as f:
    g1000_males = f.read().strip().split()

# read male sample names:
g1000_females_file = os.path.join(faststorage, 'data', '1000Genomes', 'metainfo', 'all_females.txt')
with open(g1000_females_file) as f:
    g1000_females = f.read().strip().split()

# mappings between male and pop
g1000_males_by_pop = defaultdict(list)
g1000_male_pop_files = dict()
for pop in analysis_globals.g1000_pop_info.population:
    pop_file = os.path.join(faststorage, 'data', '1000Genomes', 'metainfo', '{}_male.txt'.format(pop))
    g1000_male_pop_files[pop] = pop_file
    with open(pop_file) as f:
        g1000_males_by_pop[pop] = f.read().split()
g1000_pops_by_male = dict()
for pop, indivs in g1000_males_by_pop.items():
    for i in indivs:
        g1000_pops_by_male[i] = pop


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
g1000_callability_mask_files = dict([(chrom, os.path.join(g1000_callability_mask_dir,
                '20140520.chr{}.strict_mask.fasta.gz'.format(chrom))) for chrom in autosomes + ['X']])
g1000_callability_mask_files['X'] = os.path.join(g1000_callability_mask_dir, '20141020.chrX.strict_mask.fasta.gz')
assert len(g1000_callability_mask_files) == 23
                
human_reference = os.path.join(faststorage, 'data/cteam_lite_public3/FullyPublic/Href.fa')

# dir for male haplotypes
g1000_masked_ref_dir = os.path.join(mydir, 'steps', '1000genomes', 'masked_ref')
if not os.path.exists(g1000_masked_ref_dir):
    os.makedirs(g1000_masked_ref_dir)

g1000_masked_ref_files = dict()
for chrom, mask_file in g1000_callability_mask_files.items():

    masked_ref = modpath(mask_file, parent=g1000_masked_ref_dir, suffix='.fa')
    g1000_masked_ref_files[chrom] = masked_ref
    #g1000_masked_reference = os.path.join(faststorage, 'steps', '1000genomes', 'masked_ref', 'masked_reference.fa')

    g = gwf.target("mask_reference_g1000_{}".format(chrom), inputs=[human_reference, mask_file], 
        outputs=[masked_ref], 
        memory='16g', walltime='01:00:00') << """

        source activate simons
        python scripts/1000gen_masked_href.py {href} {mask} {output}

    """.format(href=human_reference, mask=mask_file, output=masked_ref)

#################################################################################
# Write phased haplotypes for all males
#################################################################################

# dir for male haplotypes
g1000_male_haplo_dir = os.path.join(mydir, 'steps', '1000genomes', 'haplotypes', 'males')
if not os.path.exists(g1000_male_haplo_dir):
    os.makedirs(g1000_male_haplo_dir)

g1000_male_haplotype_files = defaultdict(list)
for chrom in ['2'] + ['X']:
    chrom_dir = os.path.join(g1000_male_haplo_dir, chrom)
    if not os.path.exists(chrom_dir): 
        os.makedirs(chrom_dir)
    for sample_id in g1000_males:
        # only haplotype A for males:
        haplo_file1 = modpath("{}_{}-A.fa".format(sample_id, chrom), parent=chrom_dir)
        g1000_male_haplotype_files[chrom].append(haplo_file1)

        # get single haplotypes from only the haploid part of the male X
        gwf.target_from_template("vcf2haplo_male_{}_{}".format(chrom, sample_id),
            vcf2haplo(vcf_file=g1000_vcf_files[chrom], masked_ref=g1000_masked_ref_files[chrom],
            sample_id=sample_id, 
            haploid=chrom=='X', # only haploid if chrom is X
            out_file1=haplo_file1))

#################################################################################
# Write phased haplotypes for all female X and auto
#################################################################################

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

        # get two haplotypes from the entire female X
        gwf.target_from_template("vcf2haplo_female_{}_{}".format(chrom, sample_id),
            vcf2haplo(vcf_file=g1000_vcf_files[chrom], masked_ref=g1000_masked_ref_files[chrom],
            sample_id=sample_id, 
            out_file1=haplo_file1, out_file2=haplo_file2))


# #################################################################################
# # mask admxiture segments in male x chromosomes
# #################################################################################

# g1000_admix_masked_male_haplotype_files = dict()
# g1000_admix_masked_male_haplotype_files['X'] = g1000_male_haplotype_files['X']
g1000_admix_masked_male_haplotype_files = g1000_male_haplotype_files

# g1000_admix_masked_male_haploids_dir = 

# # dir for files
# admix_masked_male_x_haploids_dir = os.path.join(mydir, 'steps', 'male_x_haploids_admix_masked')
# if not os.path.exists(admix_masked_male_x_haploids_dir):
#     os.makedirs(admix_masked_male_x_haploids_dir)

# admix_masked_male_x_haploids = [modpath(x, parent=admix_masked_male_x_haploids_dir, suffix='.fa') for x in male_x_haploids]

# min_admix_post_prob = 0.8

# laurits_admix_pred_file = os.path.join(mydir, 'data/laurits_data/RestofworldHMMHaploid_samePAR.txt')

# for i, (unmasked, masked) in enumerate(zip(male_x_haploids, admix_masked_male_x_haploids)):
#     gwf.target_from_template("admixmask1_x_{}".format(i), 
#         admix_mask(unmasked_file=str(unmasked), masked_file=str(masked), 
#         admix_pred_file=laurits_admix_pred_file, min_post_prob=min_admix_post_prob))


#################################################################################
# pairwise differences for male admix-masked haplotypes
#################################################################################

dist_binsize = 100000

g1000_male_admix_masked_dist_file_names = dict()

# root dir for distance data for each population
g1000_male_admix_masked_dist_dir = os.path.join(mydir, 'steps', '1000genomes', 'male_haploid_dist_admix_masked')

for chrom in ['2'] + ['X']:

    chrom_dir = os.path.join(g1000_male_haplo_dir, chrom)
    if not os.path.exists(chrom_dir): 
        os.makedirs(chrom_dir)

    chrom_dir = os.path.join(g1000_male_admix_masked_dist_dir, chrom)
    if not os.path.exists(chrom_dir):
        os.makedirs(chrom_dir)

    # dict with male x admix masked haplotype files by population
    pop_file_dict = defaultdict(list)
    for file_name in g1000_admix_masked_male_haplotype_files[chrom]:
        indiv, chrom, hap = re.search(r'/([^/]+)_(\S+)-([AB]).fa', str(file_name)).groups() 
        pop = g1000_pops_by_male[indiv]
        pop_file_dict[pop].append(file_name)

    # dict with distance data files
    g1000_male_admix_masked_dist_file_names[chrom] = defaultdict(list)

    for pop, pop_files in pop_file_dict.items():

        # dir for files
        # pop_dist_dir = os.path.join(mydir, 'steps', '1000genomes', 'male_x_haploid_dist_admix_masked', pop)
        pop_dist_dir = os.path.join(g1000_male_admix_masked_dist_dir, chrom, pop)
        if not os.path.exists(pop_dist_dir):
            os.makedirs(pop_dist_dir)

        i = 0

        # iter male haploid pairs
        for file1, file2 in itertools.combinations(sorted(pop_files), 2):

            indiv1, chrom1, hap1 = re.search(r'/([^/]+)_(\S+)-([AB]).fa', str(file1)).groups() 
            indiv2, chrom2, hap2 = re.search(r'/([^/]+)_(\S+)-([AB]).fa', str(file2)).groups() 

            assert chrom1 == chrom and chrom2 == chrom

            # we do not compare chromosome from the same 
            # individul to avoid inbreeding arterfcts
            if indiv1 == indiv2:
                continue

            output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, hap1, indiv2, hap2, bp2str(dist_binsize))
            out_file_name = modpath(output_base_name, parent=pop_dist_dir)

            g1000_male_admix_masked_dist_file_names[chrom][pop].append(out_file_name)

            gwf.target_from_template('g1000_male_dist_admix_masked_windows1_{}_{}_{}'.format(pop, chrom, i), 
                admix_masked_dist_for_x_pair_template(str(file1), str(file2), 
                dist_binsize, 'NA', indiv1, hap1, indiv2, hap2, str(out_file_name)))

            i += 1


#################################################################################
# Build male distance data set for each chromsome for each population
#################################################################################

g1000_male_dist_admix_masked_store_files = defaultdict(dict)

# dir for files
g1000_male_dist_admix_masked_store_dir = os.path.join(mydir, 'steps', '1000genomes', 'male_dist_admix_masked_stores')
if not os.path.exists(g1000_male_dist_admix_masked_store_dir):
    os.makedirs(g1000_male_dist_admix_masked_store_dir)

for chrom in ['2'] + ['X']:

    chrom_dir = os.path.join(g1000_male_dist_admix_masked_store_dir, chrom)
    if not os.path.exists(chrom_dir): 
        os.makedirs(chrom_dir)

    for pop, pop_files in g1000_male_admix_masked_dist_file_names[chrom].items():

        pop_dist_dir = os.path.join(g1000_male_admix_masked_dist_dir, pop)

        pop_store_dir = os.path.join(chrom_dir, pop)
        if not os.path.exists(pop_store_dir):
            os.makedirs(pop_store_dir)

        pop_store_base_name = "male_dist_data_chr{}_{}_{}".format(chrom, bp2str(dist_binsize), pop)
        pop_store_file = modpath(pop_store_base_name, parent=pop_store_dir, suffix='.hdf')

        g1000_male_dist_admix_masked_store_files[chrom][pop] = pop_store_file

        g = gwf.target("g1000_build_dist_datasets_{}_{}".format(chrom, pop), inputs=pop_files, outputs=[pop_store_file], 
            memory='8g', walltime='1:00:00') << """

            source activate simons
            python scripts/g1000_build_male_dist_admix_masked_datasets.py \
                --dist-dir {dist_dir} \
                --out-file {out_file}

        """.format(dist_dir=pop_dist_dir, out_file=pop_store_file)


#################################################################################
# Call sweeps on the distance data with given pwdist_cutoff and min_sweep_clade_size
#################################################################################

min_sweep_clade_percent = int(analysis_globals.g1000_min_sweep_clade_proportion * 100)

g1000_male_dist_admix_masked_sweep_data_files = defaultdict(dict)

for chrom in ['2'] + ['X']:

    for pop, pop_store_file in g1000_male_dist_admix_masked_store_files[chrom].items():

        pop_sweep_data_file = modpath("sweep_data_{}_{}%.hdf".format(analysis_globals.pwdist_cutoff, min_sweep_clade_percent),
                                        parent=os.path.dirname(pop_store_file))

        pop_store_base_name = modpath(pop_store_file, parent='', suffix='')
        pop_dist_twice_file = modpath(pop_store_base_name + '_twice', parent=os.path.dirname(pop_store_file), suffix='.hdf')

        g1000_male_dist_admix_masked_sweep_data_files[chrom][pop] = pop_sweep_data_file

        gwf.target_from_template('g1000_male_sweep_data_{}_{}'.format(chrom, pop), 
            g1000_sweep_data(pop_store_file, pop_sweep_data_file, dump_dist_twice=pop_dist_twice_file))




#################################################################################
# Fst in 100kb windows for all females for all chromosomes
#################################################################################

# dir for files
g1000_fst_dir = os.path.join(mydir, 'steps', '1000genomes', 'fst')
if not os.path.exists(g1000_fst_dir):
    os.makedirs(g1000_fst_dir)

fst_pop_sets = [
    ['YRI', 'MSL', 'ESN'],
    ['FIN', 'GBR', 'IBS'],
    ['BEB', 'PJL', 'STU'],
    ['KHV', 'CDX', 'CHS'],
    ['PUR', 'CLM', 'PEL'],
]
fst_pop_sets += np.array(fst_pop_sets).T.tolist()

for fst_pops in fst_pop_sets:

    fst_pop_files = [g1000_male_pop_files[p] for p in fst_pops]

    for chrom in autosomes + ['X']:

        chrom_dir = os.path.join(g1000_fst_dir, chrom)
        if not os.path.exists(chrom_dir): 
            os.makedirs(chrom_dir)

        out_file = modpath('weir_fst_{}_{}.txt'.format(chrom, '_'.join(fst_pops)), parent=chrom_dir)

        gwf.target_from_template('fst_{}_{}'.format(chrom, '_'.join(fst_pops)),
            g1000_fst(g1000_vcf_files[chrom], fst_pop_files, out_file))
