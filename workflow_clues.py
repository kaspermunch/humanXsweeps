from gwf import Workflow, AnonymousTarget
import os, re, sys
import numpy as np
os.environ['NUMEXPR_MAX_THREADS'] = '8' # to suppress a warning
import pandas as pd
from subprocess import PIPE, Popen


sys.path.append('./scripts')
sys.path.append('./notebooks')

import simons_meta_data
import hg19_chrom_sizes


gwf = Workflow(defaults={'account': 'simons'})


#################################################################################
# Load meta data
#################################################################################

individuals, populations, regions = simons_meta_data.get_meta_data()


gwf = Workflow()

def modpath(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


# def arg_sites_file(start, end, sites_file, fasta_files):
def arg_sites_file(start, end, sites_file, fasta_files):

    inputs = fasta_files

    outputs = {'sites_file': sites_file}
    options = {
        'memory': '4g',
        'walltime': '01:00:00'
    }

    spec = f'''
    mkdir -p {os.path.dirname(sites_file)}
    python scripts/argsample_sites_file.py X {start} {end} {sites_file} {" ".join(fasta_files)}
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def arg_recomb_file(recomb_map, start, end, recomb_file):
    inputs = [recomb_map]
    outputs = {'recomb_file': recomb_file}
    options = {
        'memory': '4g',
        'walltime': '01:00:00'
    }

    spec = f'''
    mkdir -p {os.path.dirname(recomb_file)}
    python scripts/argsample_rec_window.py {recomb_map} chrX {start} {end} {recomb_file}
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def argsample(sites_file, times_file, popsize_file, recomb_file, bed_file):

    output_dir = os.path.dirname(bed_file)
    arg_sample_base_name = modpath(bed_file, suffix='')
    # TODO should be: arg_sample_base_name = modpath(bed_file, suffix=('.bed.gz', '')
    log_file = modpath(arg_sample_base_name, suffix='.log')
    tabix_file = modpath(arg_sample_base_name, suffix='.bed.gz.tbi')

    inputs = {'sites_file': sites_file, 'recomb_file': recomb_file}
    outputs = {'bed_file': bed_file, 'log_file': log_file, 'tabix_file': tabix_file}
    options = {
        'memory': '40g',
        'walltime': '14-00:00:00'
    }

    spec = f'''
    mkdir -p {output_dir}
    arg-sample -s {sites_file} \
            --times-file {times_file} \
            --popsize-file {popsize_file} \
            --recombmap {recomb_file} \
            -m 1.247e-08 \
            -c 25 \
            -n 30000 \
            --overwrite \
            -o {arg_sample_base_name} \
    && \
    ./argweaver/bin/smc2bed-all {arg_sample_base_name}
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_transition_matrices_from_argweaver(selection_coef, arg_weaver_log_file):

    output_file = f'trans.s_{selection_coef}.h5'

    path = 'steps/trans/matrices'

    inputs=[arg_weaver_log_file]
    outputs=[os.path.join(path, output_file)]
    options = {'memory': '16g', 'walltime': '5-00:00:00'}

    spec = f'''
    mkdir -p {path}
    ORIGDIR=`pwd`
    cd {path}

    python $ORIGDIR/clues-v0/make_transition_matrices_from_argweaver.py 10000 {selection_coef} \
        $ORIGDIR/{arg_weaver_log_file} {output_file} --breaks 0.95 0.025 --debug    
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def clues(bed_file, sites_file, clues_output_file, cond_trans_matrix_file, snp_pos, chrom, win_start, win_end, derived_allele, derived_freq):

    log_file = modpath(bed_file, suffix=('.bed.gz', '.log'))
    trees_file = modpath(clues_output_file, suffix=('.h5', '.trees'))
    clues_output_base_name=modpath(clues_output_file, suffix=('.h5', ''))
    
    inputs = [bed_file]
    outputs = [clues_output_file]
    options = {
        'memory': '8g',
        'walltime': '1-00:00:00'
    }

    spec = f'''        
    mkdir -p {os.path.dirname(clues_output_file)}
    arg-summarize -a {bed_file} -r {chrom}:{snp_pos}-{snp_pos} -l {log_file} -E > {trees_file} \
    && \
    python ./clues-v0/clues.py {trees_file} {cond_trans_matrix_file} {sites_file} {derived_freq} --posn {snp_pos} \
        --derivedAllele {derived_allele} --noAncientHap --approx 1000 --thin 10 --burnin 100 --output {clues_output_base_name} --debug
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def execute(cmd, stdin=None):
    process = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate(stdin)
    if not process.returncode == 0:
        print(cmd)
        print(stderr.decode())
    return stdout, stderr


def read_snp_info(snp_file):
    snp_list = list()
    with open('snps.txt', 'r') as snp_file:
        for line in snp_file:
            chrom, snp_pos, derived_allele, derived_freq, prop_missing = line.split()
            snp_pos = int(snp_pos)
            derived_freq = float(derived_freq)
            snp_list.append((chrom, snp_pos, derived_allele, derived_freq))
    return snp_list


def get_single_snp(freq_data_file, chrom, pop, snp_pos):
    snp_file_name = 'snps.txt'
    if os.path.exists(snp_file_name):
        os.remove(snp_file_name)  
    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {snp_file_name} --snppos {snp_pos}")
    snp_list = read_snp_info(snp_file_name)
    return snp_list


def get_snps(freq_data_file, chrom, window_start, window_end, min_freq, nr_snps):
    snp_file_name = 'snps.txt'
    if os.path.exists(snp_file_name):
        os.remove(snp_file_name)  
    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {snp_file_name} --start {window_start} --end {window_end} --minfreq {min_freq} --maxsnps {nr_snps}")
    snp_list = read_snp_info(snp_file_name)
    return snp_list


arg_sample_times_file = 'data/tennessen_times_fine.txt'
arg_sample_popsize_file = 'data/tennessen_popsize_fine.txt'
decode_recomb_map_file = 'data/decode_hg38_sexavg_per_gen_lifted_tohg19_chrX.tsv'

freq_data_file = 'steps/clues-v0/derived_freq_info_WestEurasia.h5'

#individuals, populations, regions = simons_meta_data.get_meta_data()

sample_fasta_files = ['steps/male_x_haploids/{}-A.fa'.format(x) for x in regions['WestEurasia'] if individuals[x]['Genetic sex assignment'] == 'XY']





# # set of europeans to include (must be same as that in compile_clues_snp_info.py)
# subset_of_europeans = ['B_Crete-2', 'B_French-3', 'B_Sardinian-3', 'S_Basque-1', 'S_Bulgarian-1', 'S_Bulgarian-2', 
#                        'S_Czech-2', 'S_English-1', 'S_Estonian-1', 'S_Estonian-2', 'S_Finnish-3', 'S_French-1', 
#                        'S_Greek-1', 'S_Greek-2', 'S_Hungarian-2', 'S_Polish-1', 'S_Saami-2', 'S_Sardinian-1', 
#                        'S_Spanish-1', 'S_Tuscan-2']
# sample_fasta_files = [f for f in sample_fasta_files if f in ]






# small run to get 
window_start, window_end = 29500000, 30500000

# make sites file
dummy_sites_task = gwf.target_from_template(
    name=f'dummy_sites_{window_start}_{window_end}',
    template=arg_sites_file(
        start=window_start,
        end=window_end,
        sites_file=f'steps/clues/dummy/dummy_{window_start}_{window_end}.sites',
        fasta_files=sample_fasta_files
    )
)

# get recombination file
dummy_recomb_task = gwf.target_from_template(
    name=f'dummy_rec_{window_start}_{window_end}',
    template=arg_recomb_file(
        recomb_map=decode_recomb_map_file,
        start=window_start,
        end=window_end,
        recomb_file=f'steps/clues/dummy/dummy_{window_start}_{window_end}.rec'
    )
)

# run argsample && smc2bed-all
dummy_argsample_target = gwf.target_from_template(
    name=f'dummy_args_{window_start}_{window_end}',
    template=argsample(
        sites_file=dummy_sites_task.outputs['sites_file'], 
        times_file=arg_sample_times_file, 
        popsize_file=arg_sample_popsize_file, 
        recomb_file=dummy_recomb_task.outputs['recomb_file'], 
        bed_file=f'steps/clues/dummy/dummy_{window_start}_{window_end}.bed.bed.gz'
        # TODO: should be bed_file=f'steps/clues/dummy/dummy_{window_start}_{window_end}.bed.gz'
    )
)


# make transition matrices for clues for a range of selection coeficients:
selection_coeficients = np.logspace(-5, np.log10(0.5), 50)
trans_task_list = list()
for i, sel_coef in enumerate(selection_coeficients):
    task = gwf.target_from_template(f'trans_mat_{i}', 
        make_transition_matrices_from_argweaver(sel_coef, dummy_argsample_target.outputs['log_file']))
    trans_task_list.append(task)

# all transition matrix files:
trans_mat_files = [output for task in trans_task_list for output in task.outputs]

# make transition matrices conditional on snp sampling frequencies:
#freqs = ' '.join([str(round(x, 2)) for x in np.linspace(0.02, 0.98, 49)])
freqs = '0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.42 0.44 0.46 0.48 0.5 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98'
cond_trans_matrix_file = 'steps/clues/trans/trans_tennesen_fine.hdf5'
cond_trans_matrix_file_no_suffix = modpath(cond_trans_matrix_file, suffix='')

gwf.target('cond_matrices', inputs=trans_mat_files, outputs=[cond_trans_matrix_file], walltime='4-00:00:00', memory='36g') << f"""
mkdir -p {os.path.dirname(cond_trans_matrix_file)}
python ./clues/conditional_transition_matrices.py {dummy_argsample_target.outputs['log_file']} steps/trans/matrices/ \
    --listFreqs {freqs} -o {cond_trans_matrix_file_no_suffix} --debug

"""

extended_peak_regions_90 = pd.read_hdf('results/extended_peak_regions_5e-05_25%_90%.hdf')

arg_win_size = 1500000
center_analyzed = 500000
flank = int((arg_win_size - center_analyzed) / 2)

clues_windows = [(int(pos - arg_win_size/2), int(pos + arg_win_size/2)) for pos in extended_peak_regions_90.pos]

for window_start, window_end in clues_windows:

    # make sites file
    sites_task = gwf.target_from_template(
        name=f'sites_{window_start}_{window_end}',
        template=arg_sites_file(
            start=window_start,
            end=window_end,
            sites_file=f'steps/clues/argsample/{window_start}_{window_end}/{window_start}_{window_end}.sites',
            fasta_files=sample_fasta_files
        )
    )

    # get recombination file
    recomb_task = gwf.target_from_template(
        name=f'rec_{window_start}_{window_end}',
        template=arg_recomb_file(
            recomb_map=decode_recomb_map_file,
            start=window_start,
            end=window_end,
            recomb_file=f'steps/clues/argsample/{window_start}_{window_end}/{window_start}_{window_end}.rec'
        )
    )

    min_freq = 0.25
    nr_snps = 20 # when made larger this needs to be multipla of the prev value to reuse prev snps

    snps_start, snps_end = window_start + flank, window_end - flank
    snp_list = get_snps(freq_data_file, 'X', snps_start, snps_end, min_freq, nr_snps)
    
    for chain in range(1, 3):

        argsample_target = gwf.target_from_template(
            name=f'args_{window_start}_{window_end}_{chain}',
            template=argsample(
                sites_file=sites_task.outputs['sites_file'], 
                times_file=arg_sample_times_file, 
                popsize_file=arg_sample_popsize_file, 
                recomb_file=recomb_task.outputs['recomb_file'], 
                bed_file=f'steps/clues/argsample/{window_start}_{window_end}/{window_start}_{window_end}_{chain}.bed.bed.gz'
                # TODO: should be bed_file=f'steps/clues/argsample/{window_start}_{window_end}/{window_start}_{window_end}_{chain}.bed.gz'
            )
        )

        clues_task_list = list()
        for chrom, snp_pos, derived_allele, derived_freq in snp_list:
            clues_task = gwf.target_from_template(
                name=f'clues_{window_start}_{window_end}_{chain}_{snp_pos}',
                template=clues(
                    bed_file=argsample_target.outputs['bed_file'], 
                    sites_file=sites_task.outputs['sites_file'], 
                    clues_output_file = modpath(argsample_target.outputs['bed_file'], 
                        suffix=('.bed.gz', f'_{snp_pos}.h5'), parent=f'steps/clues/clues'),
                    cond_trans_matrix_file=cond_trans_matrix_file, 
                    snp_pos=snp_pos, chrom=chrom, win_start=window_start, win_end=window_end, 
                    derived_allele=derived_allele, derived_freq=derived_freq                )
            )
            clues_task_list.append(clues_task)
        clues_files = [output for task in clues_task_list for output in task.outputs]


        # Extract info from all clues files for each window
        clues_csv_file_name = f'steps/clues/csv/clues_{window_start}_{window_end}_{chain}.csv'

        clues_file_base_names = ' '.join([modpath(f, parent='') for f in clues_files])

        gwf.target(f'clues_{window_start}_{window_end}_{chain}_csv', 
            inputs=clues_files, outputs=[clues_csv_file_name], walltime='1:00:00', memory='1g') << f"""

        mkdir -p {os.path.dirname(clues_csv_file_name)}
        python scripts/extract_clues_info.py {clues_csv_file_name} steps/clues/clues {clues_file_base_names}
        """

# These clues runs kept dying. So I had to just touch them and handle empty h5 files when producing csv files:

# clues_20500000_22000000_1_21036206       shouldrun      75.00% [1/0/0/3]
# clues_20500000_22000000_2_21036206       shouldrun      75.00% [1/0/0/3]
# clues_20500000_22000000_2_21479192       shouldrun      75.00% [1/0/0/3]
# clues_35550000_37050000_1_36519431       shouldrun      75.00% [1/0/0/3]
# clues_35550000_37050000_2_36519431       shouldrun      75.00% [1/0/0/3]
# clues_72400000_73900000_1_73298899       shouldrun      75.00% [1/0/0/3]
# clues_72400000_73900000_2_73298899       shouldrun      75.00% [1/0/0/3]
# clues_76300000_77800000_1_77124388       shouldrun      75.00% [1/0/0/3]
# clues_76300000_77800000_2_77124388       shouldrun      75.00% [1/0/0/3]
# clues_132000000_133500000_1_132633596    shouldrun      75.00% [1/0/0/3]
# clues_132000000_133500000_2_132633596    shouldrun      75.00% [1/0/0/3]

# touch steps/clues/clues/20500000_22000000_1.bed_21036206.h5
# touch steps/clues/clues/20500000_22000000_2.bed_21036206.h5
# touch steps/clues/clues/20500000_22000000_2.bed_21479192.h5
# touch steps/clues/clues/35550000_37050000_1.bed_36519431.h5
# touch steps/clues/clues/35550000_37050000_2.bed_36519431.h5
# touch steps/clues/clues/72400000_73900000_1.bed_73298899.h5
# touch steps/clues/clues/72400000_73900000_2.bed_73298899.h5
# touch steps/clues/clues/76300000_77800000_1.bed_77124388.h5
# touch steps/clues/clues/76300000_77800000_2.bed_77124388.h5
# touch steps/clues/clues/132000000_133500000_1.bed_132633596.h5
# touch steps/clues/clues/132000000_133500000_2.bed_132633596.h5        
        
        
        
        


# srun -A simons --mem-per-cpu 40g --time 4-00:00:00 arg-sample -s steps/sites_files/3_54250000_55750000_CEU.sites --times-file data/tennessen_times_fine.txt --popsize-file data/tennessen_popsize_fine.txt --recombmap steps/recomb_file/3_54250000_55750000.rec -m 1.247e-08 -c 25 -n 30000 -o steps/argsample/3_54250000_55750000_CEU_1/3_54250000_55750000_CEU  --resume


# ################################################################################################
# # Example run of ARGweaver to check that it works:
# # We also need the resulting log file to build transition matrices.
# ################################################################################################

# # wondow to analyse around lactase gene:
# chrom = '2'
# window_start = 136000000
# window_end = 137000000
# pop = 'CEU'
# sites_file = '../kmt/argweaver-clues/steps/sites/chr2.sites' # I just use the one Janne made

# # single lactase snp position:
# snp_pos = 136608646
# argweaver_bed_file = 'steps/arg_sample/arg_sample.bed.gz'
# argweaver_log_file = 'steps/arg_sample/arg_sample.log'
# argweaver_trees_file = 'steps/arg_sample/arg_sample.trees'
# argweaver_base_name = modpath(argweaver_log_file, suffix='')

# # run arg-sample (sampling 30,000 args), smc2bed-all, and arg-summarize:
# gwf.target('argweaver', 
#     inputs=[], 
#     outputs=[argweaver_bed_file, argweaver_log_file, argweaver_trees_file],
#     walltime='5-00:00:00', memory='10g') << f"""

# mkdir -p steps/arg_sample

# arg-sample -s {sites_file} --times-file data/tennessen_times_fine.txt \
#     --popsize-file data/tennessen_popsize_fine.txt -r 1e-8 -m 1.2e-8 -c 25 -n 10000 --overwrite -o {argweaver_base_name}

# ../../software/argweaver/bin/smc2bed-all {argweaver_base_name}

# arg-summarize -a {argweaver_bed_file} -r {chrom}:{snp_pos}-{snp_pos} \
#     -l {argweaver_log_file} -E > {argweaver_trees_file}

# """

# #############################################################################################
# # Build transition matrices using tennensen_fine demography:
# ################################################################################################

# def make_transition_matrices_from_argweaver(selection_coef, arg_weaver_log_file):

#     output_file = f'trans.s_{selection_coef}.h5'

#     path = 'steps/trans/matrices'

#     inputs=[arg_weaver_log_file]
#     outputs=[os.path.join(path, output_file)]
#     options = {'memory': '16g', 'walltime': '5-00:00:00'}

#     spec = f'''

#     mkdir -p {path}
#     ORIGDIR=`pwd`
#     cd {path}

#     python $ORIGDIR/../../software/clues/make_transition_matrices_from_argweaver.py 10000 {selection_coef} \
#         $ORIGDIR/{arg_weaver_log_file} {output_file} --breaks 0.95 0.025 --debug

#     # python $ORIGDIR/../../software/clues/make_transition_matrices_from_argweaver.py 10000 {selection_coef} \
#     #     $ORIGDIR/{arg_weaver_log_file} {output_file} --debug        

#     '''
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# # make transition matrices for clues for a range of selection coeficients:
# selection_coeficients = np.logspace(-5, np.log10(0.5), 50)
# trans_task_list = list()
# for i, sel_coef in enumerate(selection_coeficients):
#     task = gwf.target_from_template(f'trans_mat_{i}', make_transition_matrices_from_argweaver(sel_coef, argweaver_log_file))
#     trans_task_list.append(task)

# # all transition matrix files:
# trans_mat_files = [output for task in trans_task_list for output in task.outputs]

# # make transition matrices conditional on snp sampling frequencies:
# # freqs = ' '.join(map(str, 
# # [0.,     0.001,  0.002,  0.003,  0.004,  0.005,  0.006,  0.007,  0.0085, 0.0105,
# #  0.0125, 0.0145, 0.017,  0.02,   0.023,  0.0265, 0.0305, 0.035,  0.04,   0.0455,
# #  0.0515, 0.058,  0.065,  0.072,  0.0795, 0.0875, 0.096,  0.105,  0.1145, 0.1245,
# #  0.135,  0.146,  0.157,  0.1685, 0.1805, 0.1925, 0.205,  0.218,  0.231,  0.244,
# #  0.2575, 0.2715, 0.2855, 0.2995, 0.314,  0.329,  0.344,  0.359,  0.374,  0.389,
# #  0.4045, 0.4205, 0.4365, 0.4525, 0.4685, 0.4845, 0.496,  0.5,    0.504,  0.5155,
# #  0.5315, 0.5475, 0.5635, 0.5795, 0.5955, 0.611,  0.626,  0.641,  0.656,  0.671,
# #  0.686,  0.7005, 0.7145, 0.7285, 0.7425, 0.756,  0.769,  0.782,  0.795,  0.8075,
# #  0.8195, 0.8315, 0.843,  0.854,  0.865,  0.8755, 0.8855, 0.895,  0.904,  0.9125,
# #  0.9205, 0.928,  0.935,  0.942,  0.9485, 0.9545, 0.96,   0.965,  0.9695, 0.9735,
# #  0.977,  0.98,   0.983,  0.9855, 0.9875, 0.9895, 0.9915, 0.993,  0.994,  0.995,
# #  0.996,  0.997,  0.998,  0.999,  1.]))

# #freqs = ' '.join([str(round(x, 2)) for x in np.linspace(0.02, 0.98, 49)])
# freqs = '0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24 0.26 0.28 0.3 0.32 0.34 0.36 0.38 0.4 0.42 0.44 0.46 0.48 0.5 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98'

# #freqs = ' '.join([str(round(x, 2)) for x in np.linspace(0.05, 1, 20)])
# #freqs = '0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0'
# #freqs = '0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9'
# cond_trans_matrix_file = 'steps/trans/trans_tennesen_fine.hdf5'
# cond_trans_matrix_file_no_suffix = modpath(cond_trans_matrix_file, suffix='')
# gwf.target('cond_matrices', inputs=trans_mat_files, outputs=[cond_trans_matrix_file], walltime='4-00:00:00', memory='36g') << f"""

# python ../../software/clues/conditional_transition_matrices.py {argweaver_log_file} steps/trans/matrices/ \
#     --listFreqs {freqs} -o {cond_trans_matrix_file_no_suffix} --debug

# """

# ################################################################################################
# # Run clues on one lactase snp:
# ################################################################################################

# clues_output_file = 'steps/clues/chr2_136608646.h5'
# clues_output_base_name = modpath(clues_output_file, suffix='')
# derived_lactase_allele = 'G'
# derived_lactase_freq = 75e-2
# lactase_snp_pos = 136608646

# gwf.target('clues', inputs=[cond_trans_matrix_file], outputs=[clues_output_file], walltime='10:00:00', memory='36g') << f"""

# mkdir -p steps/clues

# python ../../software/clues/clues.py {argweaver_trees_file} {cond_trans_matrix_file} {sites_file} {derived_lactase_freq} \
#     --posn {lactase_snp_pos} --derivedAllele {derived_lactase_allele} --noAncientHap --approx 10000 \
#     --thin 10 --burnin 100 --output {clues_output_base_name}

# """

# ################################################################################################
# # Run clues on many SNPs in a window:
# ################################################################################################

# # file with data for looking up information about snps
# freq_data_file = 'data/derived_pop_freqs.h5'

# # we only want to analyze snps with some minimum derived frequency:
# min_freq = 0.25

# # how many snps we want to analyze:
# nr_snps = 5

# def execute(cmd, stdin=None):
#     process = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
#     stdout, stderr = process.communicate(stdin)
#     return stdout, stderr

# def read_snp_info(snp_file):
#     snp_list = list()
#     with open('snps.txt', 'r') as snp_file:
#         for line in snp_file:
#             chrom, snp_pos, derived_allele, derived_freq = line.split()
#             snp_pos = int(snp_pos)
#             derived_freq = float(derived_freq)
#             snp_list.append((chrom, snp_pos, derived_allele, derived_freq))
#     return snp_list
    
# def get_single_snp(freq_data_file, chrom, pop, snp_pos):
#     snp_file_name = 'snps.txt'
#     execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} {snp_file_name} --snppos {snp_pos}")
#     snp_list = read_snp_info(snp_file_name)
#     return snp_list

# def get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps):
#     snp_file_name = 'snps.txt'
#     execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} {snp_file_name} --start {window_start} --end {window_end} --minfreq {min_freq} --nrsnps {nr_snps}")
#     snp_list = read_snp_info(snp_file_name)
#     return snp_list

# clues_task_list = list()

# # # selected snps in window:
# # overlapping_of_window = 2e5 # in case you do 1.2 Mb windows with 200kb overlap
# # snp_list = get_snps(freq_data_file, chrom, pop, window_start+overlapping_of_window, window_end, min_freq, nr_snps)

# # single snp (the lactase snp):
# snp_list = get_single_snp(freq_data_file, chrom, pop, 136608646)

# for chrom, snp_pos, derived_allele, derived_freq in snp_list:

#     clues_output_file = f'steps/clues/clues_{chrom}_{snp_pos}_{pop}.h5'
#     clues_output_base_name = modpath(clues_output_file, suffix='')

#     clues_trees_file = modpath(clues_output_file, suffix='.trees')
# #    clues_trees_file = modpath(clues_output_file, suffix='.trees', parent='/scratch/$GWF_JOBID')

#     clues_task = gwf.target(f'clues_{chrom}_{snp_pos}_{pop}', inputs=[cond_trans_matrix_file], outputs=[clues_output_file], walltime='10:00:00', memory='36g') << f"""

#     source ./scripts/conda_init.sh
#     conda activate clues

#     mkdir -p steps/clues

#     arg-summarize -a {argweaver_bed_file} -r {chrom}:{snp_pos}-{snp_pos} \
#         -l {argweaver_log_file} -E > {clues_trees_file}

#     python ../../../software/clues/clues.py {clues_trees_file} {cond_trans_matrix_file} {sites_file} {derived_freq} \
#         --posn {snp_pos} --derivedAllele {derived_allele} --noAncientHap --approx 10000 \
#         --thin 10 --burnin 100 --output {clues_output_base_name}
#     """
#     clues_task_list.append(clues_task)

# clues_files = [output for task in clues_task_list for output in task.outputs]

# ################################################################################################
# # Extract info from all clues files
# ################################################################################################

# clues_csv_file_name = f'steps/clues/clues_{chrom}_{window_start}_{window_end}_{pop}.csv'

# clues_file_base_names = ' '.join([modpath(f, parent='', suffix='') for f in clues_files])

# gwf.target(f'clues_{chrom}_{window_start}_{window_end}_{pop}_csv', inputs=clues_files, outputs=[clues_csv_file_name], walltime='1:00:00', memory='1g') << f"""

# python scripts/extract_clues_info.py {clues_csv_file_name} steps/clues {clues_file_base_names}
# """












# ################################################################################################
# # Old version:
# ################################################################################################



# # clues_output_file = 'steps/clues/clues'
# # clues_output_base_name = modpath(clues_output_file, suffix='')

# # gwf.target('clues', inputs=[cond_trans_matrix_file], outputs=[clues_output_file], walltime='10:00:00', memory='36g') << f"""

# # mkdir -p steps/clues

# # JOBDIR=/scratch/$GWF_JOBID

# # JOBDIR=.

# # python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} $JOBDIR/snps.txt \
# #     --start {window_start} --end {window_end} --minfreq {min_freq} --nrsnps {nr_snps} 


# # python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} $JOBDIR/snps.txt \
# #     --snppos {snp_pos} 

# # while IFS= read -r line
# # do
# #     POS=$(cut -f1 -d " " <<<$line)
# #     BASE=$(cut -f2 -d " " <<<$line)
# #     FREQ=$(cut -f3 -d " " <<<$line)

# #     ~/anaconda3/envs/clues/bin/arg-summarize -a {argweaver_bed_file} -r {chrom}:$POS-$POS \
# #         -l {argweaver_log_file} -E > $JOBDIR/$POS.trees

# #     python ../../software/clues/clues.py $JOBDIR/$POS.trees {cond_trans_matrix_file} {sites_file} $FREQ \
# #         --posn $POS --derivedAllele $BASE --noAncientHap \
# #         --thin 10 --burnin 100 --output {clues_output_base_name}_{}_$POS

# # done < "$JOBDIR/snps.txt"
# # """ 
