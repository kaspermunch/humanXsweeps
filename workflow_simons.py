
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
# Load meta data
#################################################################################

individuals, populations, regions = simons_meta_data.get_meta_data()

#################################################################################
# Project root dir
#################################################################################

#faststorage = '../../'
faststorage = '/home/kmt/simons/faststorage'
#faststorage = '/project/simons/faststorage'
mydir = os.path.join(faststorage, 'people', 'kmt')
#mydir = '.'


#################################################################################
# simons input files
#################################################################################

# reference sequence file
reference_file_name = os.path.join(faststorage, 'data', 'cteam_lite_public3', 'FullyPublic', 'Href.fa')

ust_ishim_sample_id = 'Ust_Ishim'
altai_sample_id = 'Altai'
denisova_sample_id = 'Denisova'

orig_sample_files = list()
orig_mask_files = list()
for sample_id in sorted(individuals):
    orig_sample_files.append(os.path.join(faststorage, 
        'data', 'cteam_lite_public3', 
        'FullyPublic', '{}.ccomp.fa.rz'.format(sample_id)))
    orig_mask_files.append(os.path.join(faststorage, 
        'data', 'cteam_lite_public3', 
        'FullyPublic', '{}.ccompmask.fa.rz'.format(sample_id)))

# ust ishim:
orig_ust_ishim_sample_file = os.path.join(faststorage, 
                                        'data', 'cteam_lite_public3', 
                                        'FullyPublic', '{}.ccomp.fa.rz'.format(ust_ishim_sample_id))
orig_ust_ishim_mask_file = os.path.join(faststorage, 
                                        'data', 'cteam_lite_public3', 
                                        'FullyPublic', '{}.ccompmask.fa.rz'.format(ust_ishim_sample_id))

# altai:
orig_altai_sample_file = os.path.join(faststorage, 
                                        'data', 'cteam_lite_public3', 
                                        'FullyPublic', '{}.ccomp.fa.rz'.format(altai_sample_id))
orig_altai_mask_file = os.path.join(faststorage, 
                                        'data', 'cteam_lite_public3', 
                                        'FullyPublic', '{}.ccompmask.fa.rz'.format(altai_sample_id))

# denisova:
orig_denisova_sample_file = os.path.join(faststorage, 
                                        'data', 'cteam_lite_public3', 
                                        'FullyPublic', '{}.ccomp.fa.rz'.format(denisova_sample_id))
orig_denisova_mask_file = os.path.join(faststorage, 
                                        'data', 'cteam_lite_public3', 
                                        'FullyPublic', '{}.ccompmask.fa.rz'.format(denisova_sample_id))
#################################################################################
# turn rz files into gzip files
#################################################################################

# I do this because there is some trailing grabage in 
# the rz files that python gzip cannot handle

# dir for files
sample_dir = os.path.join(mydir, 'steps', 'gziped_samples')
if not os.path.exists(sample_dir):
    os.makedirs(sample_dir)

sample_files = [modpath(x, parent=sample_dir, suffix='.gz') for x in orig_sample_files]

mask_files = [modpath(x, parent=sample_dir, suffix='.gz') for x in orig_mask_files]

for i, (orig_sample_file, sample_file) in enumerate(zip(orig_sample_files, sample_files)):
    gwf.target_from_template('smpl_rz2gz_{}'.format(i), rz2gz(rz_file=str(orig_sample_file), gz_file=str(sample_file)))

for i, (orig_mask_file, mask_file) in enumerate(zip(orig_mask_files, mask_files)):
    gwf.target_from_template('mask_rz2gz_{}'.format(i), rz2gz(rz_file=str(orig_mask_file), gz_file=str(mask_file)))

# ust_ishim
ust_ishim_sample_file = modpath(orig_ust_ishim_sample_file, parent=sample_dir, suffix='.gz')
ust_ishim_mask_file = modpath(orig_ust_ishim_mask_file, parent=sample_dir, suffix='.gz')
gwf.target_from_template('smpl_rz2gz_{}'.format('ust_ishim'), rz2gz(rz_file=str(orig_ust_ishim_sample_file), gz_file=str(ust_ishim_sample_file)))
gwf.target_from_template('mask_rz2gz_{}'.format('ust_ishim'), rz2gz(rz_file=str(orig_ust_ishim_mask_file), gz_file=str(ust_ishim_mask_file)))

# altai
altai_sample_file = modpath(orig_altai_sample_file, parent=sample_dir, suffix='.gz')
altai_mask_file = modpath(orig_altai_mask_file, parent=sample_dir, suffix='.gz')
gwf.target_from_template('smpl_rz2gz_{}'.format('altai'), rz2gz(rz_file=str(orig_altai_sample_file), gz_file=str(altai_sample_file)))
gwf.target_from_template('mask_rz2gz_{}'.format('altai'), rz2gz(rz_file=str(orig_altai_mask_file), gz_file=str(altai_mask_file)))

# denisova
denisova_sample_file = modpath(orig_denisova_sample_file, parent=sample_dir, suffix='.gz')
denisova_mask_file = modpath(orig_denisova_mask_file, parent=sample_dir, suffix='.gz')
gwf.target_from_template('smpl_rz2gz_{}'.format('denisova'), rz2gz(rz_file=str(orig_denisova_sample_file), gz_file=str(denisova_sample_file)))
gwf.target_from_template('mask_rz2gz_{}'.format('denisova'), rz2gz(rz_file=str(orig_denisova_mask_file), gz_file=str(denisova_mask_file)))


# Hack to make the seq and mask same length for Ust Ishim (trimmed to shortest one of the two):
trimmed_ust_ishim_sample_file = modpath(ust_ishim_sample_file, suffix=('.fa.gz', '.trimmed.fa.gz'))
trimmed_ust_ishim_mask_file = modpath(ust_ishim_mask_file, suffix=('.fa.gz', '.trimmed.fa.gz'))

g = gwf.target("trim_ust_ishim", inputs=[ust_ishim_sample_file, ust_ishim_mask_file],
                      outputs=[trimmed_ust_ishim_sample_file, trimmed_ust_ishim_mask_file], 
                      memory='15g', walltime='11:00:00') << """

    source activate simons
    python scripts/trim_ust_ishim.py {} {} {} {}

""".format(ust_ishim_sample_file, ust_ishim_mask_file, trimmed_ust_ishim_sample_file, trimmed_ust_ishim_mask_file)


# Hack to make the archaic sequences same length as all other:
# random sample file as template:
template_sample_file = sample_files[0]

padded_altai_sample_file = modpath(altai_sample_file, suffix=('.fa.gz', '.padded.fa.gz'))
padded_altai_mask_file = modpath(altai_mask_file, suffix=('.fa.gz', '.padded.fa.gz'))
gwf.target_from_template("pad_altai", pad_archaic_files(template_file=template_sample_file,
           input_file=altai_sample_file, pad_char='N', output_file=padded_altai_sample_file))
gwf.target_from_template("pad_altai_mask", pad_archaic_files(template_file=template_sample_file,
           input_file=altai_mask_file, pad_char='0', output_file=padded_altai_mask_file))

padded_denisova_sample_file = modpath(denisova_sample_file, suffix=('.fa.gz', '.padded.fa.gz'))
padded_denisova_mask_file = modpath(denisova_mask_file, suffix=('.fa.gz', '.padded.fa.gz'))
gwf.target_from_template("pad_denisova", pad_archaic_files(template_file=template_sample_file,
           input_file=denisova_sample_file, pad_char='N', output_file=padded_denisova_sample_file))
gwf.target_from_template("pad_denisova_mask", pad_archaic_files(template_file=template_sample_file,
           input_file=denisova_mask_file, pad_char='0', output_file=padded_denisova_mask_file))

# g = gwf.target("pad_archaic", inputs=[altai_sample_file, altai_mask_file],
#                       outputs=[padded_altai_sample_file, padded_altai_mask_file], 
#                       memory='15g', walltime='11:00:00') << """

#     source activate simons
#     python scripts/pad_archaic_genome.py {template} {input_seq} N {output_seq}
#     python scripts/pad_archaic_genome.py {template} {input_mask} 0 {output_mask}

# """.format(template=template_sample_file,
#            input_seq=altai_sample_file, output_seq=padded_altai_sample_file,
#            input_mask=altai_mask_file, output_mask=padded_altai_mask_file)


#################################################################################
# mask samples
#################################################################################

mask_level = 1

# dir for files
masked_sample_dir = os.path.join(mydir, 'steps', 'masked_samples')
if not os.path.exists(masked_sample_dir):
    os.makedirs(masked_sample_dir)

masked_sample_files = [modpath(x, parent=masked_sample_dir) for x in sample_files]

for i, (unmasked, mask, masked) in enumerate(zip(sample_files, mask_files, masked_sample_files)):
    gwf.target_from_template("masking_{}".format(i), mask_sample(unmasked_file=str(unmasked), 
        mask_file=str(mask), masked_file=str(masked), mask_level=mask_level))

# ust_ishim
ust_ishim_masked_sample_file = modpath(trimmed_ust_ishim_sample_file, parent=masked_sample_dir)
gwf.target_from_template("masking_{}".format('ust_ishim'), mask_sample(unmasked_file=str(trimmed_ust_ishim_sample_file), 
        mask_file=str(trimmed_ust_ishim_mask_file), masked_file=str(ust_ishim_masked_sample_file), 
        mask_level=mask_level,
        skip=['Y'])) 

# altai
altai_masked_sample_file = modpath(padded_altai_sample_file, parent=masked_sample_dir)
gwf.target_from_template("masking_{}".format('altai'), mask_sample(unmasked_file=str(padded_altai_sample_file), 
        mask_file=str(padded_altai_mask_file), masked_file=str(altai_masked_sample_file), 
        mask_level=mask_level,
        skip=['Y']))

# denisova
denisova_masked_sample_file = modpath(padded_denisova_sample_file, parent=masked_sample_dir)
gwf.target_from_template("masking_{}".format('denisova'), mask_sample(unmasked_file=str(padded_denisova_sample_file), 
        mask_file=str(padded_denisova_mask_file), masked_file=str(denisova_masked_sample_file), 
        mask_level=mask_level,
        skip=['Y']))


#################################################################################
# generate pseudohaploids
#################################################################################

# dir for files
pseudohaploid_dir = os.path.join(mydir, 'steps', 'pseudohaploid_genomes')
if not os.path.exists(pseudohaploid_dir):
    os.makedirs(pseudohaploid_dir)

# Build targets for generating pseudhaploids
pseudohaploid_file_names = defaultdict(list)
for i, sample_file_name in enumerate(masked_sample_files):
    basename = os.path.basename(sample_file_name).split('.')[0]
    output_file_name1 = modpath('{}-A.fa.gz'.format(basename), parent=pseudohaploid_dir)
    output_file_name2 = modpath('{}-B.fa.gz'.format(basename), parent=pseudohaploid_dir)

    pop = individuals[basename]['Population ID']

    pseudohaploid_file_names[pop].extend([os.path.join(output_file_name1), os.path.join(output_file_name2)])

    # # NB: if male only the pseudohaploid to list of pesudohaploid files for downstread analysis
    # is_female = individuals[basename]['Genetic sex assignment'] == 'XX'
    # pseudohaploid_file_names[pop].append(os.path.join(output_file_name1))
    # if is_female:
    #     pseudohaploid_file_names[pop].append(os.path.join(output_file_name2))

    gwf.target_from_template('psudohaploids_{}'.format(i), pseudohaploids(input_file=sample_file_name, 
        ref_file_name=str(reference_file_name),
        output_file1=str(output_file_name1),
        output_file2=str(output_file_name2)))


# ust_ishim:
basename = os.path.basename(ust_ishim_masked_sample_file).split('.')[0]
ust_ishim_output_file_name1 = modpath('{}-A.fa.gz'.format(basename), parent=pseudohaploid_dir)
ust_ishim_output_file_name2 = modpath('{}-B.fa.gz'.format(basename), parent=pseudohaploid_dir)

ust_ishim_pseudohaploid_file_names = [os.path.join(ust_ishim_output_file_name1), 
                                        os.path.join(ust_ishim_output_file_name2)]

gwf.target_from_template('psudohaploids_{}'.format('ust_ishim'), pseudohaploids(input_file=ust_ishim_masked_sample_file, 
    ref_file_name=str(reference_file_name),
    output_file1=str(ust_ishim_output_file_name1),
    output_file2=str(ust_ishim_output_file_name2)))


# altai:
basename = os.path.basename(altai_masked_sample_file).split('.')[0]
altai_output_file_name1 = modpath('{}-A.fa.gz'.format(basename), parent=pseudohaploid_dir)
altai_output_file_name2 = modpath('{}-B.fa.gz'.format(basename), parent=pseudohaploid_dir)

altai_pseudohaploid_file_names = [os.path.join(altai_output_file_name1), 
                                  os.path.join(altai_output_file_name2)]

gwf.target_from_template('psudohaploids_{}'.format('altai'), pseudohaploids(input_file=altai_masked_sample_file, 
    ref_file_name=str(reference_file_name),
    output_file1=str(altai_output_file_name1),
    output_file2=str(altai_output_file_name2)))


# denisova:
basename = os.path.basename(denisova_masked_sample_file).split('.')[0]
denisova_output_file_name1 = modpath('{}-A.fa.gz'.format(basename), parent=pseudohaploid_dir)
denisova_output_file_name2 = modpath('{}-B.fa.gz'.format(basename), parent=pseudohaploid_dir)

denisova_pseudohaploid_file_names = [os.path.join(denisova_output_file_name1), 
                                  os.path.join(denisova_output_file_name2)]

gwf.target_from_template('psudohaploids_{}'.format('denisova'), pseudohaploids(input_file=denisova_masked_sample_file, 
    ref_file_name=str(reference_file_name),
    output_file1=str(denisova_output_file_name1),
    output_file2=str(denisova_output_file_name2)))


archaic_pseudohaploid_file_names = altai_pseudohaploid_file_names + denisova_pseudohaploid_file_names


#################################################################################
# compute pi in windows over masked genomes
#################################################################################

# size of windows for computing pi
pi_binsize = 100000
    
pi_dir = os.path.join(mydir, 'steps', 'population_pi')
if not os.path.exists(pi_dir):
    os.makedirs(pi_dir)

pi_file_names = list()

#arg_list = list()
# iter populations
i = 0
for pop, pop_samples in sorted(pseudohaploid_file_names.items()):
    
    # iter pseudohaploid pairs
    for file1, file2 in itertools.combinations(sorted(pop_samples), 2):
        
        indiv1, pseud1 = re.search(r'/([^/]+)-([AB]).fa.gz', str(file1)).groups() 
        indiv2, pseud2 = re.search(r'/([^/]+)-([AB]).fa.gz', str(file2)).groups() 

        # we do not compare chromosome from the same 
        # individul to avoid inbreeding arterfcts
        if indiv1 == indiv2:
            continue
            
        # open files for the pair of pseudohaploids to compare
        f1 = modpath(file1, parent=pseudohaploid_dir)
        f2 = modpath(file2, parent=pseudohaploid_dir)

        output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, pseud1, indiv2, pseud2, bp2str(pi_binsize))
        out_file_name = modpath(output_base_name, parent=pi_dir)

        pi_file_names.append(out_file_name)

        gwf.target_from_template('pi_windows_{}'.format(i), pi_for_pair_template(str(f1), str(f2), 
            pi_binsize, pop, indiv1, pseud1, indiv2, pseud2, str(out_file_name)))

        i += 1



# ust_ishim (in this case we compute the heterozygosity):
file1, file2 = ust_ishim_pseudohaploid_file_names
indiv1, pseud1 = re.search(r'/([^/]+)-([AB]).fa.gz', str(file1)).groups() 
indiv2, pseud2 = re.search(r'/([^/]+)-([AB]).fa.gz', str(file2)).groups() 

f1 = modpath(file1, parent=pseudohaploid_dir)
f2 = modpath(file2, parent=pseudohaploid_dir)

output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, pseud1, indiv2, pseud2, bp2str(pi_binsize))

ust_ishim_pi_file_name = modpath(output_base_name, parent=pi_dir)

gwf.target_from_template('pi_windows_{}'.format('ust_ishim'), pi_for_pair_template(str(f1), str(f2), 
    pi_binsize, pop, indiv1, pseud1, indiv2, pseud2, str(ust_ishim_pi_file_name)))



#################################################################################
# compute diffs in windows over masked chromosomes
# NB: only compares X pseudohaploids between pairs of which at least one is African
#################################################################################

dist_binsize = 100000

# dir for files
dist_dir = os.path.join(mydir, 'steps', 'afr_nonafr_x_pseudohap_dist')
if not os.path.exists(dist_dir):
    os.makedirs(dist_dir)

dist_file_names = list()

# iter populations
i = 0

x_pseudohaploids = list()
for pop in sorted(pseudohaploid_file_names):
    x_pseudohaploids.extend(pseudohaploid_file_names[pop])

# iter pseudohaploid pairs
for file1, file2 in itertools.combinations(sorted(x_pseudohaploids), 2):
    
    indiv1, pseud1 = re.search(r'/([^/]+)-([AB]).fa.gz', str(file1)).groups() 
    indiv2, pseud2 = re.search(r'/([^/]+)-([AB]).fa.gz', str(file2)).groups() 

    # we do not compare chromosome from the same 
    # individul to avoid inbreeding arterfcts
    if indiv1 == indiv2:
        continue

    # only compare two individuals if one is an African:
    if not (indiv1 in regions['Africa'] or indiv2 in regions['Africa']):
        continue
        
    # open files for the pair of pseudohaploids to compare
    f1 = modpath(file1, parent=pseudohaploid_dir)
    f2 = modpath(file2, parent=pseudohaploid_dir)

    output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, pseud1, indiv2, pseud2, bp2str(dist_binsize))
    out_file_name = modpath(output_base_name, parent=dist_dir)

    dist_file_names.append(out_file_name)

    gwf.target_from_template('dist_windows_{}'.format(i), dist_for_x_pair_template(str(f1), str(f2), 
        dist_binsize, 'NA', indiv1, pseud1, indiv2, pseud2, str(out_file_name)))

    i += 1


#################################################################################
# Build pi data sets for each chromosome with added meta info
#################################################################################

# dir for files
pi_store_dir = os.path.join(mydir, 'steps', 'pi_stores')
if not os.path.exists(pi_store_dir):
    os.makedirs(pi_store_dir)

metadata_dir = '/home/kmt/simons/faststorage/data/metadata'

pi_store_base_names = ["pi_data_{}_{}".format(x, bp2str(pi_binsize)) for x in hg19_chrom_sizes.hg19_chrom_sizes.keys()]
pi_store_files = [modpath(x, parent=pi_store_dir, suffix='.store') for x in pi_store_base_names]

g = gwf.target("build_pi_datasets", inputs=pi_file_names, outputs=pi_store_files, 
    memory='60g', walltime='11:00:00') << """

    source activate simons
    python scripts/build_pi_datasets.py \
        --pi-dir {pi_dir} \
        --result-dir {pi_store_dir} \
        --meta-data-dir {metadata_dir}

""".format(pi_dir=pi_dir, pi_store_dir=pi_store_dir, metadata_dir=metadata_dir)


#################################################################################
# Build distance data sets for each chromosome with added meta info
# NB: only X pseudohaploids between pairs of which at least one is African
#################################################################################

# dir for files
dist_store_dir = os.path.join(mydir, 'steps', 'dist_stores')
if not os.path.exists(dist_store_dir):
    os.makedirs(dist_store_dir)

#dist_store_base_names = ["dist_data_{}_{}".format(x, bp2str(dist_binsize)) for x in hg19_chrom_sizes.hg19_chrom_sizes.keys()]
dist_store_base_names = ["dist_data_chrX_{}".format(bp2str(dist_binsize))]
dist_store_files = [modpath(x, parent=dist_store_dir, suffix='.store') for x in dist_store_base_names]

g = gwf.target("build_dist_datasets", inputs=dist_file_names, outputs=dist_store_files, 
    memory='150g', walltime='11:00:00') << """

    source activate simons
    python scripts/build_dist_datasets.py \
        --dist-dir {dist_dir} \
        --result-dir {dist_store_dir} \
        --meta-data-dir {metadata_dir}

""".format(dist_dir=dist_dir, dist_store_dir=dist_store_dir, metadata_dir=metadata_dir)


#################################################################################
# extract male x chromosomes
#################################################################################

male_subset = list()
for pop in sorted(pseudohaploid_file_names):
    for file_name in pseudohaploid_file_names[pop]:
        basename = os.path.basename(file_name).split('.')[0]
        if basename.endswith('-A'):
            if individuals[basename.replace('-A', '')]['Genetic sex assignment'] == 'XY':
                male_subset.append(file_name)

# spike in ust ishim
file_name = ust_ishim_pseudohaploid_file_names[0]
basename = os.path.basename(file_name).split('.')[0]
assert basename.endswith('-A')
male_subset.append(file_name)

# dir for files
male_x_haploids_dir = os.path.join(mydir, 'steps', 'male_x_haploids')
if not os.path.exists(male_x_haploids_dir):
    os.makedirs(male_x_haploids_dir)

male_x_haploids = [modpath(x, parent=male_x_haploids_dir, suffix='') for x in male_subset]

for i, (full_genome, only_x) in enumerate(zip(male_subset, male_x_haploids)):
    gwf.target_from_template("extract_x_{}".format(i), 
        extract_x(full_genome=str(full_genome), only_x=str(only_x)))


#################################################################################
# extract x pseudohaploids for altai and denisova
#################################################################################

# dir for files
archaic_x_pseudohaploids_dir = os.path.join(mydir, 'steps', 'archaic_x_pseudohaploids')
if not os.path.exists(archaic_x_pseudohaploids_dir):
    os.makedirs(archaic_x_pseudohaploids_dir)

archaic_x_pseudohaploids = [modpath(x, parent=archaic_x_pseudohaploids_dir, suffix='') for x in archaic_pseudohaploid_file_names]

for i, (full_genome, only_x) in enumerate(zip(archaic_pseudohaploid_file_names, archaic_x_pseudohaploids)):
    gwf.target_from_template("extract_archaic_x_{}".format(i), 
        extract_x(full_genome=str(full_genome), only_x=str(only_x)))


##################################################################################
## mask out ampliconic regions (replace with N) from extracted male X chromosomes
##################################################################################
#
#male_x_haploids_ampl_masked_dir = os.path.join(mydir, 'steps', 'male_x_haploids_ampl_masked')
#if not os.path.exists(male_x_haploids_ampl_masked_dir):
#    os.makedirs(male_x_haploids_ampl_masked_dir)
#
#male_x_haploids_ampl_masked = [modpath(x, parent=male_x_haploids_ampl_masked_dir, suffix='.fa') for x in male_x_haploids]
#
#ampl_regions_file = '/home/kmt/simons/faststorage/people/kmt/data/coordinates_hg18_hg19_hg38_Amplicons_Gap.txt'
#
#for i, (unmasked, masked) in enumerate(zip(male_x_haploids, male_x_haploids_ampl_masked)):
#    gwf.target_from_template("mask_ampliconic_regions_{}".format(i), 
#        mask_ampliconic_regions(unmasked_file=str(unmasked), masked_file=str(masked), ampl_regions_file=ampl_regions_file))


#################################################################################
# mask out ampliconic regions (replace with N) from extracted male X chromosomes
#################################################################################

# dir for files
male_x_haploids_ampl_masked_dir = os.path.join(mydir, 'steps', 'male_x_haploids_ampl_masked')
if not os.path.exists(male_x_haploids_ampl_masked_dir):
    os.makedirs(male_x_haploids_ampl_masked_dir)

male_x_haploids_ampl_masked = [modpath(x, parent=male_x_haploids_ampl_masked_dir, suffix='.fa') for x in male_x_haploids]

ampl_regions_file = '/home/kmt/simons/faststorage/people/kmt/data/coordinates_hg18_hg19_hg38_Amplicons_Gap.txt'

for i, (unmasked, masked) in enumerate(zip(male_x_haploids, male_x_haploids_ampl_masked)):
    gwf.target_from_template("mask_ampliconic_regions_{}".format(i), 
        mask_ampliconic_regions(unmasked_file=str(unmasked), masked_file=str(masked), ampl_regions_file=ampl_regions_file))


#################################################################################
# mask admxiture segments in male x chromosomes
#################################################################################

# dir for files
admix_masked_male_x_haploids_dir = os.path.join(mydir, 'steps', 'male_x_haploids_admix_masked')
if not os.path.exists(admix_masked_male_x_haploids_dir):
    os.makedirs(admix_masked_male_x_haploids_dir)

admix_masked_male_x_haploids = [modpath(x, parent=admix_masked_male_x_haploids_dir, suffix='.fa') for x in male_x_haploids]

min_admix_post_prob = 0.8

laurits_admix_pred_file = os.path.join(mydir, 'data/laurits_data/RestofworldHMMHaploid_samePAR.txt')

for i, (unmasked, masked) in enumerate(zip(male_x_haploids, admix_masked_male_x_haploids)):
    gwf.target_from_template("admixmask1_x_{}".format(i), 
        admix_mask(unmasked_file=str(unmasked), masked_file=str(masked), 
        admix_pred_file=laurits_admix_pred_file, min_post_prob=min_admix_post_prob))


# #################################################################################
# # same but for haplotypes with ampliconic regions masked out
# #################################################################################

# # dir for files
# ampl_and_admix_masked_male_x_haploids_dir = os.path.join(mydir, 'steps', 'male_x_haploids_ampl_and_admix_masked')
# if not os.path.exists(ampl_and_admix_masked_male_x_haploids_dir):
#     os.makedirs(ampl_and_admix_masked_male_x_haploids_dir)

# ampl_and_admix_masked_male_x_haploids = [modpath(x, parent=ampl_and_admix_masked_male_x_haploids_dir, suffix='.fa') for x in male_x_haploids_ampl_masked]

# for i, (unmasked, masked) in enumerate(zip(male_x_haploids_ampl_masked, ampl_and_admix_masked_male_x_haploids)):
#     gwf.target_from_template("admixmask2_x_{}".format(i), 
#         admix_mask(unmasked_file=str(unmasked), masked_file=str(masked), 
#         admix_pred_file=laurits_admix_pred_file, min_post_prob=min_admix_post_prob))


# #################################################################################
# # compute diffs in windows over all male x haplotypes
# #################################################################################

# # dir for files
# male_dist_dir = os.path.join(mydir, 'steps', 'male_x_haploid_dist')
# if not os.path.exists(male_dist_dir):
#     os.makedirs(male_dist_dir)

# male_dist_file_names = list()

# i = 0

# # iter male haploid pairs
# for file1, file2 in itertools.combinations(sorted(male_x_haploids), 2):
    
#     indiv1, hap1 = re.search(r'/([^/]+)-([AB]).fa', str(file1)).groups() 
#     indiv2, hap2 = re.search(r'/([^/]+)-([AB]).fa', str(file2)).groups() 

#     # we do not compare chromosome from the same 
#     # individul to avoid inbreeding arterfcts
#     if indiv1 == indiv2:
#         continue

#     # open files for the pair of pseudohaploids to compare
#     f1 = modpath(file1, parent=male_x_haploids_dir)
#     f2 = modpath(file2, parent=male_x_haploids_dir)

#     output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, hap1, indiv2, hap2, bp2str(dist_binsize))
#     out_file_name = modpath(output_base_name, parent=male_dist_dir)

#     male_dist_file_names.append(out_file_name)

#     gwf.target_from_template('male_dist_windows1_{}'.format(i), dist_for_x_pair_template(str(f1), str(f2), 
#         dist_binsize, 'NA', indiv1, hap1, indiv2, hap2, str(out_file_name)))

#     i += 1

# #################################################################################
# # same but for male haplotypes with ampliconic regions masked out
# #################################################################################

# # dir for files
# male_dist_dir_ampl_masked = os.path.join(mydir, 'steps', 'male_x_haploid_dist_ampl_masked')
# if not os.path.exists(male_dist_dir_ampl_masked):
#     os.makedirs(male_dist_dir_ampl_masked)

# male_ampl_masked_dist_file_names = list()

# i = 0

# # iter male haploid pairs
# for file1, file2 in itertools.combinations(sorted(male_x_haploids_ampl_masked), 2):
    
#     indiv1, hap1 = re.search(r'/([^/]+)-([AB]).fa', str(file1)).groups() 
#     indiv2, hap2 = re.search(r'/([^/]+)-([AB]).fa', str(file2)).groups() 

#     # we do not compare chromosome from the same 
#     # individul to avoid inbreeding arterfcts
#     if indiv1 == indiv2:
#         continue

#     # open files for the pair of pseudohaploids to compare
#     f1 = modpath(file1, parent=male_x_haploids_ampl_masked_dir)
#     f2 = modpath(file2, parent=male_x_haploids_ampl_masked_dir)

#     output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, hap1, indiv2, hap2, bp2str(dist_binsize))
#     out_file_name = modpath(output_base_name, parent=male_dist_dir_ampl_masked)

#     male_ampl_masked_dist_file_names.append(out_file_name)

#     gwf.target_from_template('male_dist_windows2_{}'.format(i), dist_for_x_pair_template(str(f1), str(f2), 
#         dist_binsize, 'NA', indiv1, hap1, indiv2, hap2, str(out_file_name)))

#     i += 1

#################################################################################
# same but for the admix-masked haplotypes (produces separate stats for admix masked)
#################################################################################

# dir for files
male_admix_masked_dist_dir = os.path.join(mydir, 'steps', 'male_x_haploid_dist_admix_masked')
if not os.path.exists(male_admix_masked_dist_dir):
    os.makedirs(male_admix_masked_dist_dir)

male_admix_masked_dist_file_names = list()

i = 0

# iter male haploid pairs
for file1, file2 in itertools.combinations(sorted(admix_masked_male_x_haploids), 2):
    
    indiv1, hap1 = re.search(r'/([^/]+)-([AB]).fa', str(file1)).groups() 
    indiv2, hap2 = re.search(r'/([^/]+)-([AB]).fa', str(file2)).groups() 

    # we do not compare chromosome from the same 
    # individul to avoid inbreeding arterfcts
    if indiv1 == indiv2:
        continue

    # open files for the pair of pseudohaploids to compare
    f1 = modpath(file1, parent=admix_masked_male_x_haploids_dir)
    f2 = modpath(file2, parent=admix_masked_male_x_haploids_dir)

    output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, hap1, indiv2, hap2, bp2str(dist_binsize))
    out_file_name = modpath(output_base_name, parent=male_admix_masked_dist_dir)

    male_admix_masked_dist_file_names.append(out_file_name)

    gwf.target_from_template('male_dist_admix_masked_windows1_{}'.format(i), admix_masked_dist_for_x_pair_template(str(f1), str(f2), 
        dist_binsize, 'NA', indiv1, hap1, indiv2, hap2, str(out_file_name)))

    i += 1

# #################################################################################
# # same but for the admix-masked haplotypes WITH ampliconic regions masked out (produces separate stats for admix masked)
# #################################################################################

# # dir for files
# male_ampl_and_admix_masked_dist_dir = os.path.join(mydir, 'steps', 'male_x_haploid_dist_ampl_and_admix_masked')
# if not os.path.exists(male_ampl_and_admix_masked_dist_dir):
#     os.makedirs(male_ampl_and_admix_masked_dist_dir)

# male_ampl_and_admix_masked_dist_file_names = list()

# i = 0

# # iter male haploid pairs
# for file1, file2 in itertools.combinations(sorted(ampl_and_admix_masked_male_x_haploids), 2):
    
#     indiv1, hap1 = re.search(r'/([^/]+)-([AB]).fa', str(file1)).groups() 
#     indiv2, hap2 = re.search(r'/([^/]+)-([AB]).fa', str(file2)).groups() 

#     # we do not compare chromosome from the same 
#     # individul to avoid inbreeding arterfcts
#     if indiv1 == indiv2:
#         continue

#     # open files for the pair of pseudohaploids to compare
#     f1 = modpath(file1, parent=ampl_and_admix_masked_male_x_haploids_dir)
#     f2 = modpath(file2, parent=ampl_and_admix_masked_male_x_haploids_dir)

#     output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, hap1, indiv2, hap2, bp2str(dist_binsize))
#     out_file_name = modpath(output_base_name, parent=male_ampl_and_admix_masked_dist_dir)

#     male_ampl_and_admix_masked_dist_file_names.append(out_file_name)

#     gwf.target_from_template('male_dist_admix_masked_windows2_{}'.format(i), admix_masked_dist_for_x_pair_template(str(f1), str(f2), 
#         dist_binsize, 'NA', indiv1, hap1, indiv2, hap2, str(out_file_name)))

#     i += 1

# #################################################################################
# # Build distance data sets and call sweeps for each male chromosome with added meta info
# #################################################################################

# # dir for files
# male_dist_store_dir = os.path.join(mydir, 'steps', 'male_dist_stores')
# if not os.path.exists(male_dist_store_dir):
#     os.makedirs(male_dist_store_dir)

# #male_dist_store_base_names = ["male_dist_data_{}_{}".format(x, bp2str(dist_binsize)) for x in hg19_chrom_sizes.hg19_chrom_sizes.keys()]
# male_dist_store_base_name = "male_dist_data_chrX_{}".format(bp2str(dist_binsize))
# male_dist_store_file = modpath(male_dist_store_base_name, parent=male_dist_store_dir, suffix='.hdf')

# g = gwf.target("build_male_dist_datasets1", inputs=male_dist_file_names, outputs=[male_dist_store_file], 
#     memory='80g', walltime='11:00:00') << """

#     source activate simons
#     python scripts/build_male_dist_datasets.py \
#         --dist-dir {dist_dir} \
#         --meta-data-dir {metadata_dir} \
#         --out-file {out_file}

# """.format(dist_dir=male_dist_dir, out_file=male_dist_store_file, metadata_dir=metadata_dir)

# #################################################################################
# # Call sweeps on the distance data with given pwdist_cutoff and min_sweep_clade_size
# #################################################################################

# male_dist_sweep_data_file = os.path.join(male_dist_store_dir, "sweep_data_{}_{}.hdf".format(analysis_globals.pwdist_cutoff, 
#                                                                                      analysis_globals.min_sweep_clade_size))
# gwf.target_from_template('male_dist_sweep_data', sweep_data(male_dist_store_file, male_dist_sweep_data_file))


# #################################################################################
# # same but with ampliconic regions masked out
# #################################################################################

# # dir for files
# male_dist_ampl_masked_store_dir = os.path.join(mydir, 'steps', 'male_dist_ampl_masked_stores')
# if not os.path.exists(male_dist_ampl_masked_store_dir):
#     os.makedirs(male_dist_ampl_masked_store_dir)

# male_dist_ampl_masked_store_base_name = "male_dist_data_chrX_{}".format(bp2str(dist_binsize))
# male_dist_ampl_masked_store_file = modpath(male_dist_ampl_masked_store_base_name, parent=male_dist_ampl_masked_store_dir, suffix='.hdf')

# g = gwf.target("build_male_dist_datasets2", inputs=male_ampl_masked_dist_file_names, outputs=[male_dist_ampl_masked_store_file], 
#     memory='80g', walltime='11:00:00') << """

#     source activate simons
#     python scripts/build_male_dist_datasets.py \
#         --dist-dir {dist_dir} \
#         --meta-data-dir {metadata_dir} \
#         --out-file {out_file}

# """.format(dist_dir=male_dist_dir_ampl_masked, out_file=male_dist_ampl_masked_store_file, metadata_dir=metadata_dir)

# #################################################################################
# # Call sweeps on the distance data with given pwdist_cutoff and min_sweep_clade_size
# #################################################################################

# male_dist_ampl_masked_sweep_data_file = os.path.join(male_dist_ampl_masked_store_dir, "sweep_data_{}_{}.hdf".format(analysis_globals.pwdist_cutoff, 
#                                                                                      analysis_globals.min_sweep_clade_size))
# gwf.target_from_template('male_dist_ampl_masked_sweep_data', sweep_data(male_dist_ampl_masked_store_file, male_dist_ampl_masked_sweep_data_file))


#################################################################################
# same but for the admix-masked haplotypes
#################################################################################


# dir for files
male_dist_admix_masked_store_dir = os.path.join(mydir, 'steps', 'male_dist_admix_masked_stores')
if not os.path.exists(male_dist_admix_masked_store_dir):
    os.makedirs(male_dist_admix_masked_store_dir)

male_dist_admix_masked_store_base_name = "male_dist_data_chrX_{}".format(bp2str(dist_binsize))
male_dist_admix_masked_store_file = modpath(male_dist_admix_masked_store_base_name, parent=male_dist_admix_masked_store_dir, suffix='.hdf')

g = gwf.target("build_male_dist_admix_masked_datasets1", inputs=male_admix_masked_dist_file_names, outputs=[male_dist_admix_masked_store_file], 
    memory='80g', walltime='11:00:00') << """

    source activate simons
    python scripts/build_male_dist_admix_masked_datasets.py \
        --dist-dir {dist_dir} \
        --meta-data-dir {metadata_dir} \
        --out-file {out_file}

""".format(dist_dir=male_admix_masked_dist_dir, out_file=male_dist_admix_masked_store_file, metadata_dir=metadata_dir)

#################################################################################
# Call sweeps on the distance data with given pwdist_cutoff and min_sweep_clade_size
#################################################################################

male_dist_admix_masked_sweep_data_file = \
    os.path.join(male_dist_admix_masked_store_dir, "sweep_data_{}_{}.hdf".format(analysis_globals.pwdist_cutoff, 
                                                                                 analysis_globals.min_sweep_clade_size))

male_dist_admix_masked_dist_twice_file = modpath(male_dist_admix_masked_store_base_name + '_twice', parent=male_dist_admix_masked_store_dir, suffix='.hdf')

gwf.target_from_template('male_dist_admix_masked_sweep_data', sweep_data(male_dist_admix_masked_store_file, 
                                                                         male_dist_admix_masked_sweep_data_file, 
                                                                         dump_dist_twice=male_dist_admix_masked_dist_twice_file))


# #################################################################################
# # same but for the ampliconic region masked AND admix-masked haplotypes
# #################################################################################


# # dir for files
# male_dist_ampl_and_admix_masked_store_dir = os.path.join(mydir, 'steps', 'male_dist_ampl_and_admix_masked_stores')
# if not os.path.exists(male_dist_ampl_and_admix_masked_store_dir):
#     os.makedirs(male_dist_ampl_and_admix_masked_store_dir)

# male_dist_ampl_and_admix_masked_store_base_name = "male_dist_data_chrX_{}".format(bp2str(dist_binsize))
# male_dist_ampl_and_admix_masked_store_file = modpath(male_dist_ampl_and_admix_masked_store_base_name, parent=male_dist_ampl_and_admix_masked_store_dir, suffix='.store')

# g = gwf.target("build_male_dist_admix_masked_datasets2", inputs=male_ampl_and_admix_masked_dist_file_names, outputs=[male_dist_ampl_and_admix_masked_store_file], 
#     memory='80g', walltime='11:00:00') << """

#     source activate simons
#     python scripts/build_male_dist_admix_masked_datasets.py \
#         --dist-dir {dist_dir} \
#         --meta-data-dir {metadata_dir} \
#         --out-file {out_file}

# """.format(dist_dir=male_ampl_and_admix_masked_dist_dir, out_file=male_dist_ampl_and_admix_masked_store_file, metadata_dir=metadata_dir)


# #################################################################################
# # Call sweeps on the distance data with given pwdist_cutoff and min_sweep_clade_size
# #################################################################################

# male_dist_ampl_and_admix_masked_sweep_data_file = os.path.join(male_dist_ampl_and_admix_masked_store_dir, "sweep_data_{}_{}.hdf".format(analysis_globals.pwdist_cutoff, 
#                                                                                      analysis_globals.min_sweep_clade_size))
# gwf.target_from_template('male_dist_ampl_and_admix_masked_sweep_data', sweep_data(male_dist_ampl_and_admix_masked_store_file, male_dist_ampl_and_admix_masked_sweep_data_file))


#################################################################################
# compute additional diffs between archaic female pseudohaplotids and all male x haplotypes
#################################################################################

# dir for files
archaic_dist_dir = os.path.join(mydir, 'steps', 'archaic_x_pseudohaploid_dist')
if not os.path.exists(archaic_dist_dir):
    os.makedirs(archaic_dist_dir)

archaic_dist_file_names = list()

i = 0

for file1, file2 in itertools.product(sorted(male_x_haploids), archaic_x_pseudohaploids):

    indiv1, hap1 = re.search(r'/([^/]+)-([AB]).fa', str(file1)).groups() 
    indiv2, hap2 = re.search(r'/([^/]+)-([AB]).fa', str(file2)).groups() 

    output_base_name = '{}_{}_{}_{}_{}.pickle' .format(indiv1, hap1, indiv2, hap2, bp2str(dist_binsize))
    out_file_name = modpath(output_base_name, parent=archaic_dist_dir)

    archaic_dist_file_names.append(out_file_name)

    gwf.target_from_template('archaic_dist_windows_{}'.format(i), dist_for_x_pair_template(str(file1), str(file2), 
        dist_binsize, 'NA', indiv1, hap1, indiv2, hap2, str(out_file_name)))

    i += 1


#################################################################################
# Build distance data sets for archaic pseudohaploids and male x chromosomes
#################################################################################

# dir for files
archaic_dist_store_dir = os.path.join(mydir, 'steps', 'archaic_dist_stores')
if not os.path.exists(archaic_dist_store_dir):
    os.makedirs(archaic_dist_store_dir)

#male_dist_store_base_names = ["male_dist_data_{}_{}".format(x, bp2str(dist_binsize)) for x in hg19_chrom_sizes.hg19_chrom_sizes.keys()]
archaic_dist_store_base_name = "archaic_dist_data_chrX_{}".format(bp2str(dist_binsize))
archaic_dist_store_file = modpath(archaic_dist_store_base_name, parent=archaic_dist_store_dir, suffix='.hdf')

g = gwf.target("build_archaic_dist_datasets", inputs=archaic_dist_file_names, outputs=[archaic_dist_store_file], 
    memory='80g', walltime='11:00:00') << """

    source activate simons
    python scripts/build_archaic_dist_datasets.py \
        --dist-dir {dist_dir} \
        --meta-data-dir {metadata_dir} \
        --out-file {out_file}

""".format(dist_dir=archaic_dist_dir, out_file=archaic_dist_store_file, metadata_dir=metadata_dir)


#################################################################################
# make fasta alignments of male x chromosomes
#################################################################################

argweaver_binsize = 100000

def argweaver_input_targets(region_label, male_x_haploids):

    # dir for files
    argweaver_input_dir = os.path.join(mydir, 'steps', 'argweaver', 'input', region_label)
    if not os.path.exists(argweaver_input_dir):
        os.makedirs(argweaver_input_dir)

    # make list of expected argweaver fasta input files:
    argweaver_input_files = list()
    chrom_len = hg19_chrom_sizes.hg19_chrom_sizes['chrX']
    for i in range(0, chrom_len - argweaver_binsize, argweaver_binsize):
        file_name = 'X-{:09d}-{:09d}.fa'.format(i, min(chrom_len, i+argweaver_binsize)) # NB: next time I should format the input files like the stores...
        #file_path = os.path.join(mydir, 'steps', 'argweaver', 'input', file_name)
        #file_path = argweaver_input_dir / file_name
        file_path = modpath(file_name, parent=argweaver_input_dir)
        argweaver_input_files.append(file_path)

    # gwf.target('fasta_align_{}'.format(region_label)) << fasta_alignments(list(map(str, male_x_haploids)), 
    #     list(map(str, argweaver_input_files)),
    #     argweaver_binsize, str(argweaver_input_dir))
    gwf.target_from_template('fasta_align_{}'.format(region_label), fasta_alignments(list(map(str, male_x_haploids)), 
        list(map(str, argweaver_input_files)),
        argweaver_binsize, str(argweaver_input_dir)))

    return argweaver_input_files

argweaver_input_files = dict()
#male_x_haploids_subset = [x for x in male_x_haploids if x.name.replace('-A.fa', '') in individuals]
male_x_haploids_subset = [x for x in male_x_haploids if os.path.basename(x).replace('-A.fa', '') in individuals]
# (for world set we also filter to only get individuals not filtereded out of the meta data info)
argweaver_input_files['World'] = argweaver_input_targets('World', male_x_haploids_subset)
for region_label in list(regions.keys()):
    male_x_haploids_subset = [x for x in male_x_haploids if os.path.basename(x).replace('-A.fa', '') in regions[region_label]]
    argweaver_input_files[region_label] = argweaver_input_targets(region_label, male_x_haploids_subset)


#################################################################################
# argweaver
#################################################################################

region_labels = ['World'] + list(regions.keys())

def run_argweaver_analysis(region_label, argweaver_input_files):

    argweaver_output_dir = os.path.join(mydir, 'steps', 'argweaver', 'output', region_label)

    if not os.path.exists(argweaver_output_dir):
        os.makedirs(argweaver_output_dir)

    argweaver_output_files = list()

    for i, input_file in enumerate(argweaver_input_files):

        output_file = modpath(input_file, parent=argweaver_output_dir, suffix='.tsv.gz')
        argweaver_output_files.append(output_file)

        # gwf.target('argweaver_{}_{}'.format(region_label, i)) << argeaver_window_analysis(input_fasta_file=input_file, output_hdf_file=output_file)
        gwf.target_from_template('argweaver_{}_{}'.format(region_label, i),
         argeaver_window_analysis(input_fasta_file=input_file, output_file=output_file))

    return argweaver_output_files

argweaver_output_files = dict()
for region_label in region_labels:
    argweaver_output_files[region_label] = run_argweaver_analysis(region_label, argweaver_input_files[region_label])


#################################################################################
# For each analysis window, extract mean tmrca tmrca_half for each chain and group
#################################################################################

for region_label in region_labels:
    tmrca_dir = os.path.join(mydir, 'steps', 'argweaver', 'tmrca', region_label)
    if not os.path.exists(tmrca_dir):
        os.makedirs(tmrca_dir)
    for i, input_table_file in enumerate(argweaver_output_files[region_label]):
        output_tmrca_file = modpath(input_table_file, parent=tmrca_dir, suffix=('.tsv.gz', '.hdf'))
        gwf.target_from_template('tmrca_{}_{}'.format(region_label, i), 
            compute_tmrca_window_stats(input_table_file, output_tmrca_file))



#################################################################################
# Compute extra tmrca statistics from pruned argweaver trees for World analysis
#################################################################################

excluded_pops = ','.join(simons_meta_data.excluded_populations)
excluded_indivs = ','.join(simons_meta_data.excluded_individuals)

# dir for files
annotated_output_dir = os.path.join(mydir, 'steps', 'argweaver', 'annotated_output', 'World')
if not os.path.exists(annotated_output_dir):
    os.makedirs(annotated_output_dir)

argweaver_annotated_output_files = defaultdict(list)
for region_label in ['World']:
    for i, input_file_name in enumerate(argweaver_output_files[region_label]):
        output_extra_file = modpath(input_file_name, parent=annotated_output_dir)
        gwf.target_from_template('argweaver_extra_stats_{}_{}'.format(region_label, i), 
            argweaver_extra_stats(input_file_name, output_extra_file, excluded_pops, excluded_indivs))
        argweaver_annotated_output_files[region_label].append(output_extra_file)


#################################################################################
# For the extra stats or pruned trees, for each analysis window, extract mean tmrca tmrca_half for each chain and group
#################################################################################

for region_label in ['World']:
    tmrca_extra_dir = os.path.join(mydir, 'steps', 'argweaver', 'tmrca_extra', region_label)
    if not os.path.exists(tmrca_extra_dir):
        os.makedirs(tmrca_extra_dir)
    for i, input_table_file in enumerate(argweaver_annotated_output_files[region_label]):
        output_tmrca_extra_file = modpath(input_table_file, parent=tmrca_extra_dir, suffix=('.tsv.gz', '.hdf'))
        gwf.target_from_template('tmrca_extra_{}_{}'.format(region_label, i), 
            compute_extra_tmrca_window_stats(input_table_file, output_tmrca_extra_file))


#################################################################################
# additional summary stats for each sampled tree
#################################################################################

# for region_label in region_labels:

#     stats_dir = os.path.join(mydir, 'steps', 'argweaver', 'stats', region_label)
#     if not os.path.exists(stats_dir):
#         os.makedirs(stats_dir)

#     for i, input_hdf_file in enumerate(argweaver_output_files[region_label]):
#         output_hdf_file = modpath(input_hdf_file, parent=stats_dir, suffix=('.tsv.gz', '.hdf'))
#         component_hdf_file = modpath(input_hdf_file, parent=stats_dir, suffix=('.tsv.gz', '.comp.hdf'))
#         component_stats_hdf_file = modpath(input_hdf_file, parent=stats_dir, suffix=('.tsv.gz', '.compstats.hdf'))
#         sweep_sister_clade_hdf_file = modpath(input_hdf_file, parent=stats_dir, suffix=('.tsv.gz', '.sweepsister.hdf'))
#         nonsweep_sister_clade_hdf_file = modpath(input_hdf_file, parent=stats_dir, suffix=('.tsv.gz', '.nonsweepsister.hdf'))

#         gwf.target_from_template('treestats_{}_{}'.format(region_label, i), 
#             compute_tree_stats(input_hdf_file, output_hdf_file, 
#                 component_hdf_file, component_stats_hdf_file,
#                 sweep_sister_clade_hdf_file, nonsweep_sister_clade_hdf_file))

#################################################################################
# liftovers
#################################################################################


def reciprocal_liftover(intervals_files, forwards_chain_file, backwards_chain_file, 
                        slurm_tag, steps_dir, target_chromosomes):
    """
    Does reciprocal lift over of a set of intervals.
    """

    if not steps_dir.exists():
        os.makedirs(str(steps_dir))

    # output files
    mapped_files= [steps_dir / x.with_suffix('.mapped').name for x in intervals_files]
    unmapped_files = [x.with_suffix('.unmapped') for x in mapped_files]
    backmapped_files = [x.with_suffix('.backmapped') for x in mapped_files]
    unbackmapped_files = [x.with_suffix('.nobackmapped') for x in mapped_files]
    filtered_files = [x.with_suffix('.filtered') for x in mapped_files]

    lifted_files = [steps_dir / 'sorted' / "{}.bed".format(x) for x in target_chromosomes]

    for i in range(len(intervals_files)):

        # lift over intervals
        gwf.target_from_template('{}_lift_{}'.format(slurm_tag, i),
            liftover(bed_file=intervals_files[i], chain_file=forwards_chain_file,
            mapped_file=mapped_files[i], unmapped_file=unmapped_files[i]))

        # lift back to orginal coordinates to ensure one to one correspondence
        gwf.target_from_template('{}_liftback_{}'.format(slurm_tag, i),
            liftover(bed_file=mapped_files[i], chain_file=backwards_chain_file,
            mapped_file=backmapped_files[i], unmapped_file=unbackmapped_files[i]))

        # filter out intervals that does not map both ways
        gwf.target_from_template('{}_filter_{}'.format(slurm_tag, i),
            bed_difference(bed_file1=mapped_files[i], bed_file2=unbackmapped_files[i],
            output_file=filtered_files[i]))

    # filter out intervals that does not map both ways
    gwf.target_from_template('{}_merge_and_split'.format(slurm_tag),
        bed_merge_and_split(input_files=filtered_files, output_files=lifted_files))

    return lifted_files

# split map file per chromosome......



chains_dir = Path('data/chain_files')

# decode hg38 map files
hg38_map_files = []

# chromosomes we are interested in (not other random contigs)
target_chromosomes = ['chr{}'.format(x) for x in list(range(1,23)) + ['X']]

# lift decode map from hg38 to hg19
hg19_map_files = reciprocal_liftover(hg38_map_files,
    forwards_chain_file=chains_dir/'hg38ToHg19.over.chain', 
    backwards_chain_file=chains_dir/'hg19ToHg38.over.chain',
    slurm_tag='liftover',
    steps_dir=Path(os.getcwd(), 'steps', 'decode_liftover'),
    target_chromosomes=target_chromosomes)


#################################################################################
# slim simulations
#################################################################################

# slim_file_list = ['scripts/bottle.slim', 
#                   'scripts/complete_before_bottle.slim',                  
#                   'scripts/complete_in_bottle.slim',                  
#                   'scripts/complete_after_bottle.slim',
#                   'scripts/partial_before_bottle.slim',                  
#                   'scripts/partial_in_bottle.slim',                  
#                   'scripts/partial_after_bottle.slim',
#                   'scripts/complete_after_no_bottle.slim']

slim_tree_files = list()
slim_dist_files = list()

# how many to simulate:
nr_non_africans = sum(x['Region'] != 'Africa' and x['Genetic sex assignment'] == 'XY' for x in individuals.values())

simulations_dir = os.path.join(mydir, 'steps', 'slim', 'simulations')

slim_output_dir = simulations_dir

total_sim_generations = 100000

# pasted fro nb_25_slim_simulations notebook:
standard_demography = \
[(1, 20000),
 (191379, 12000),
 (194827, 6000),
 (196551, 4000),
 (197586, 3000),
 (198275, 4000),
 (198793, 6000),
 (199137, 12000),
 (199482, 20000),
 (199793, 100000)]

# I SHOULD PRODUCE THE DEMOGRAPHIES HERE
# IT WOULD BE NICER WITH YEARS IN THE NAMES OF FILES...

# named autosomal population size demographies:
demographies = [('standard', standard_demography)]

# make SLiM model X chromosome
chrom = 'X'

# x/autosome ratios:
#x_auto_ratios = [0.55 * x for x in [1, 0.73]] # <- recuction below one reflects region pi vs. global pi in Africans.
x_auto_ratios = [0.66 * x for x in [1, 0.73]]

# mean per generation recombination rate in regions (new decode map):
#rec_rates_per_gen = [0.465e-8 * x for x in [1, 2/4, 2/6]] # <- reductions below one reflect male/famale retios of 2 and 4
rec_rates_per_gen = [0.465e-8, 1.16e-8]

for x_auto_ratio in x_auto_ratios:
    for rec_rate_per_gen in rec_rates_per_gen:

        assert chrom == 'X'
        meiosis_rec_rate =  rec_rate_per_gen * 3 / 2
        xreduction = x_auto_ratio / 0.75
        
        for demog_name, demog in demographies:
            for sweep_type in ['nosweep']: #['complete', 'partial', 'nosweep']:
                for sweep_start in [198275]: #[198965, 198275, 197586, 196551]: # pasted fro nb_25_slim_simulations notebook
                    for selcoef in [0.1]: #[0.01, 0.05, 0.1, 0.2]: 

                        id_str = '{}_{}_{}_{}_{}_{}_{}'.format(demog_name,
                        round(x_auto_ratio*100), round(rec_rate_per_gen * 1e12),
                        chrom, sweep_type, sweep_start, int(selcoef*100))

                        slim_output_dir = os.path.join(simulations_dir, id_str.replace('_', '/'))
                        if not os.path.exists(slim_output_dir): os.makedirs(slim_output_dir)

                        for i in range(15):
                            sim_output_prefix = os.path.join(slim_output_dir, "{}_{}".format(id_str, i))
                            slim_tree_file = sim_output_prefix + '.trees'
                            slim_tree_files.append(slim_tree_file)
                            slim_dist_file = sim_output_prefix + '.hdf'
                            slim_dist_files.append(slim_dist_file)

                            gwf.target_from_template("{}_{}".format(id_str, i),
                                slim_sim(selcoef, analysis_globals.gen_time, 
                                '{:.12f}'.format(analysis_globals.mut_per_year), 
                                meiosis_rec_rate,
                                nr_non_africans,
                                sweep_type, sweep_start, demog, 
                                chrom, xreduction, 
                                total_sim_generations,
                                slim_tree_file, slim_dist_file))

##################################################################################
# calling sweeps on slim simulations
#################################################################################

slim_sweep_data_dir = os.path.join(mydir, 'steps', 'slim', 'sweep_data')
if not os.path.exists(slim_sweep_data_dir):
    os.makedirs(slim_sweep_data_dir)

slim_sweep_data_files = list()
for i, slim_dist_file in enumerate(slim_dist_files):
    slim_sweep_data_file = modpath(slim_dist_file, parent=slim_sweep_data_dir)
    slim_sweep_data_files.append(slim_sweep_data_file)
    gwf.target_from_template('slim_sweep_data_{}'.format(i), sweep_data(slim_dist_file, slim_sweep_data_file, 
                                                                        cores=1, memory='8g', walltime='00:30:00'))


#################################################################################
# compute prop swept for all slim sweep data in a file that includes simulation info
#################################################################################

slim_summary_file = os.path.join(mydir, 'steps', 'slim', 'slim_summary.hdf')

g = gwf.target("slim_summary", 
    inputs=slim_sweep_data_files, outputs=[slim_summary_file], 
    memory='10g', walltime='00:30:00') << """

    source activate simons
    python scripts/slim_summary.py {slim_sweep_data_dir} {out_file}

""".format(slim_sweep_data_dir=slim_sweep_data_dir, out_file=slim_summary_file)


#################################################################################

if __name__ == "__main__":
    print(len(gwf.targets))