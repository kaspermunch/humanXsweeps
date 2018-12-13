

import os, sys
#from gwf import template

##############################################################################
# Templates
##############################################################################


# rz2gz = template(inputs=['{rz_file}'], outputs=['{gz_file}']) << """

#     zcat {rz_file} | gzip > {gz_file}

#     """
def rz2gz(rz_file, gz_file):

    spec = """

    zcat {rz_file} | gzip > {gz_file}

    """.format(rz_file=rz_file, gz_file=gz_file)

#    return (options, shell_spec)
    return [rz_file], [gz_file], {}, spec


# mask_sample = template(inputs=['{unmasked_file}', '{mask_file}'],
#                           outputs=['{masked_file}'],
#                           memory='4g',
#                           walltime='11:00:00') << """

#     source activate simons
#     python scripts/mask_sample.py {mask_level} {unmasked_file} {mask_file} {masked_file}

#     """
def mask_sample(unmasked_file, mask_file, mask_level, masked_file, skip=[]):

    skipped = ''
    for s in skip:
        skipped += ' --skip {}'.format(s)

    options = {'memory': '4g',
               'walltime': '11:00:00'
              }

    spec = """

    source activate simons
    python scripts/mask_sample.py {skipped} --mask-level {mask_level} {unmasked_file} {mask_file} {masked_file}

    """.format(mask_file=mask_file, unmasked_file=unmasked_file, 
        mask_level=mask_level, masked_file=masked_file, skipped=skipped)

    return [unmasked_file, mask_file], [masked_file], options, spec


# # generate two pesudohaploids from one pseudodiploid
# pseudohaploids = template(inputs=['{input_file}'],
#                           outputs=['{output_file1}', '{output_file2}'],
#                           memory='15g',
#                           walltime='11:00:00') << """

#     source activate simons
#     python scripts/pseudohaploids.py {input_file} {ref_file_name} {output_file1} {output_file2}

#     """
# generate two pesudohaploids from one pseudodiploid
def pseudohaploids(input_file, ref_file_name, output_file1, output_file2):

    options = {'memory': '15g',
               'walltime': '11:00:00'
              } 

    spec = """

    source activate simons
    python scripts/pseudohaploids.py {input_file} {ref_file_name} {output_file1} {output_file2}

    """.format(input_file=input_file, ref_file_name=ref_file_name, 
        output_file1=output_file1, output_file2=output_file2)

    return [input_file], [output_file1, output_file2], options, spec


# # compute pi for a pair of pseudohaploids
# def pi_for_pair_template(*args):

#     f1, f2, binsize, pop, indiv1, pseud1, indiv2, pseud2, output_file_name = args

#     options = {'inputs': [f1, f2],
#             'outputs': [output_file_name],
#             'memory': '4g',
#             'walltime': '11:00:00'} 

#     shell_spec = '''

#         source activate simons
#         python scripts/pi_windows.py {args}

#     '''.format(args=' '.join(map(str, args)))

#     return (options, shell_spec)   
# compute pi for a pair of pseudohaploids
def pi_for_pair_template(*args):

    f1, f2, binsize, pop, indiv1, pseud1, indiv2, pseud2, output_file_name = args

    options = {'memory': '4g',
               'walltime': '11:00:00'
               } 

    spec = '''

        source activate simons
        python scripts/pi_windows.py {args}

    '''.format(args=' '.join(map(str, args)))

    return [f1, f2], [output_file_name], options, spec

# compute pi for a pair of pseudohaploids
def dist_for_x_pair_template(*args):

    # This is identical to pi_for_pair_template and the script called is the identical
    # pi_windows.py except that it skips all but the X chromosome

    f1, f2, binsize, pop, indiv1, pseud1, indiv2, pseud2, output_file_name = args

    options = {'memory': '4g',
               'walltime': '11:00:00'
               } 

    spec = '''

        source activate simons
        python scripts/x_diff_windows.py {args}

    '''.format(args=' '.join(map(str, args)))

    return [f1, f2], [output_file_name], options, spec


def admix_masked_dist_for_x_pair_template(*args):

    # This is identical to pi_for_pair_template and the script called is the identical
    # pi_windows.py except that it skips all but the X chromosome

    f1, f2, binsize, pop, indiv1, pseud1, indiv2, pseud2, output_file_name = args

    options = {'memory': '4g',
               'walltime': '11:00:00'
               } 

    spec = '''

        source activate simons
        python scripts/x_admix_masked_diff_windows.py {args}

    '''.format(args=' '.join(map(str, args)))

    return [f1, f2], [output_file_name], options, spec


# def make_fasta_index(fasta_file, copied_fasta_file, index_file):

#     options = {'memory': '4g',
#                'walltime': '01:00:00'
#                } 

#     spec = '''

#         source activate simons
#         cat {fasta_file} | gzip -d > {copied_fasta_file}
#         python scripts/index_fasta.py {copied_fasta_file}

#     '''.format(fasta_file=fasta_file, copied_fasta_file=copied_fasta_file)

#     return [fasta_file], [copied_fasta_file, index_file], options, spec


# def dist_matrix(fasta_files, fasta_index_files, output_file):

#     options = {'memory': '4g',
#                'walltime': '11:00:00'
#                } 

#     spec = '''

#         source activate simons
#         python scripts/dist_matrix.py --output {output_file} {fasta_files}

#     '''.format(output_file=output_file, ' '.join(fasta_files))

#     return fasta_files + fasta_index_files, [output_file], options, spec



# # argweaver: simulate args
# simulate_sites = template(inputs=[], outputs=['{output_prefix}.sites']) << """

#     source activate argweaver
#     arg-sim --model=coal_recomb -k 100 -L 1000000 -N 10000 -r 1.6e-8 -m 1.8e-8 -o {output_prefix}

#     """

# # extract the pseudohaploid x chromosomes from males
# extract_male_x = template(inputs=['{full_genome}'], outputs=['{only_x}'], 
#     memory='10g') << """

#     source activate simons
#     python scripts/extract_male_x.py {full_genome} {only_x}

#     """
# extract the pseudohaploid x chromosomes from males
def extract_x(full_genome, only_x): 

    options = {'memory': '10g'}

    spec = """

    source activate simons
    python scripts/extract_x.py {full_genome} {only_x}

    """.format(full_genome=full_genome, only_x=only_x)

    return [full_genome], [only_x], options, spec


def mask_ampliconic_regions(unmasked_file, masked_file, ampl_regions_file):

    options = {'memory': '1g'}

    spec = """

    source activate simons
    python scripts/mask_ampl_regions.py {ampl_regions_file} {unmasked_file} {masked_file}

    """.format(unmasked_file=str(unmasked_file), masked_file=str(masked_file),
               ampl_regions_file=str(ampl_regions_file))

    return [unmasked_file, ampl_regions_file], [masked_file], options, spec


def admix_mask(unmasked_file, masked_file, admix_pred_file, min_post_prob): 

    options = {'memory': '10g'}

    spec = """

    source activate simons
    python scripts/admix_mask.py {admix_pred_file} {min_post_prob} {unmasked_file} {masked_file}

    """.format(admix_pred_file=admix_pred_file, min_post_prob=min_post_prob,
     unmasked_file=unmasked_file, masked_file=masked_file)

    return [unmasked_file, admix_pred_file], [masked_file], options, spec

# # generate alignments for argweaver
# def fasta_alignments(male_x_file_names, argweaver_input_files,
#     argweaver_binsize, argweaver_input_dir):

#     options = {'inputs': male_x_file_names,
#                'outputs': argweaver_input_files,
#                'memory': '4g',
#                'walltime': '11:00:00' }

#     shell_spec = """

#     source activate simons
#     python scripts/generate_argweaver_input.py \
#         -w {argweaver_binsize} -o {argweaver_input_dir} {male_x_files}

#     """.format(argweaver_binsize=argweaver_binsize,
#                argweaver_input_dir=argweaver_input_dir,
#                male_x_files=" ".join(male_x_file_names))

#     return options, shell_spec)
# generate alignments for argweaver
def fasta_alignments(male_x_file_names, argweaver_input_files,
    argweaver_binsize, argweaver_input_dir):

    options = {'inputs': male_x_file_names,
               'outputs': argweaver_input_files,
               'memory': '4g',
               'walltime': '11:00:00' }

    spec = """

    source activate simons
    python scripts/generate_argweaver_input.py \
        -w {argweaver_binsize} -o {argweaver_input_dir} {male_x_files}

    """.format(argweaver_binsize=argweaver_binsize,
               argweaver_input_dir=argweaver_input_dir,
               male_x_files=" ".join(male_x_file_names))

    return male_x_file_names, argweaver_input_files, options, spec


def argeaver_window_analysis(input_fasta_file, output_file):

    options = {'memory': '4g',
               'walltime': '48:00:00'}

    spec = """

    ITERATIONS=5000

    JOBDIR=/scratch/$GWF_JOBID

    source activate argweaver

    for CHAIN in 1 2; do

        INPUTFILE={input_fasta_file}
        INPUTBASENAME=${{INPUTFILE##*/}}
        INPUTBASENAMENOSUFFIX=${{INPUTBASENAME%.*}}
        OUTPUTPREFIX=$JOBDIR/${{INPUTBASENAMENOSUFFIX}}_samples_$CHAIN

        arg-sample --fasta $INPUTFILE \
            --popsize 10000 --recombrate 1.6e-8 \
            --mutrate 1.8e-8 --ntimes 20 --maxtime 200e3 \
            --randseed $CHAIN \
            --iters $ITERATIONS --sample-step 10 \
            --output $OUTPUTPREFIX \
            --overwrite

        if [ ! -f $OUTPUTPREFIX.$ITERATIONS.smc.gz ]; then
            >&2 echo "ARGWEAVER SAMPLING DID NOT COMPLETE FOR CHAIN $CHAIN, MISSING $OUTPUTPREFIX.$ITERATIONS.smc.gz"
        else

            # extract times from log file and put in seperate file:
            grep 'times = \[' $OUTPUTPREFIX.log | tail -n 1 | cut -c12- | awk '{{ gsub(/[,\]]/,"\\n") }}; 1' > $OUTPUTPREFIX.times

            for SMCFILE in $OUTPUTPREFIX.*.smc.gz ; do
                PARTS=($(echo `basename $SMCFILE` | tr "." "\n"))
                SAMPLENR=${{PARTS[1]}} 
                BEDFILE=${{SMCFILE%.*}}.bed.gz
                STATSFILE=${{SMCFILE%.*}}.stats
                smc2bed --times $OUTPUTPREFIX.times --sample $SAMPLENR $SMCFILE | sort-bed - | gzip > $BEDFILE
                arg-summarize -a $BEDFILE --time-file $OUTPUTPREFIX.times -E -T -B -R -K -H -F -P -C > $STATSFILE
            done

        fi
    done

    ALLSTATSFILES=($JOBDIR/*.$ITERATIONS.smc.stats)
    NRSTATSFILES=${{#ALLSTATSFILES[@]}}
    if [ $NRSTATSFILES -eq 2 ]; then
        source deactivate
        source activate simons
        #python scripts/write_argweaver_to_hdf.py --store {output_file} --sample-dir $JOBDIR --glob-pattern '*.smc.stats'        
        python scripts/write_argweaver_to_table.py --sample-dir $JOBDIR --glob-pattern '*.smc.stats' | gzip > {output_file}
    fi

    """.format(input_fasta_file=input_fasta_file, output_file=output_file)

    return [input_fasta_file], [output_file], options, spec



# def compute_tree_stats(input_table_file, output_hdf_file):

#     options = {'memory': '4g',
#                'walltime': '24:00:00'}

#     spec = """

#     python ./scripts/tree_stats_hdf.py {input_table_file} {output_hdf_file}

#     """.format(input_table_file=input_table_file, output_hdf_file=output_hdf_file)

#     return [input_table_file], [output_hdf_file], options, spec



# # argweaver: sample args
# def sample_args(input_file, output_prefix, output_files, 
#     sample_iterations, sample_steps, input_format, random_seed):

#     options = {'inputs': [input_file], 
#         'outputs': output_files + ['{}.times'.format(output_prefix)],
#         'memory': '4g',
#         'walltime': '24:00:00'}

#     shell_spec = '''

#     source activate argweaver
#     arg-sample --{input_format} {input_file} \
#         --popsize 10000 --recombrate 1.6e-8 \
#         --mutrate 1.8e-8 --ntimes 20 --maxtime 200e3 \
#         --randseed {random_seed} \
#         --iters {sample_iterations} --sample-step {sample_steps} \
#         --output {output_prefix} \
#         --overwrite

#     # extract times from log file and put in seperate file:
#     grep 'times = \[' {output_prefix}.log | tail -n 1 | cut -c12- | \
#     awk '{{ gsub(/[,\]]/,"\\n") }}; 1' > {output_prefix}.times

#     '''.format(input_file=input_file, 
#         sample_iterations=sample_iterations, sample_steps=sample_steps,
#         output_prefix=output_prefix, input_format=input_format, random_seed=random_seed)

#     return (options, shell_spec)


# argweaver: extract TMRCA
# extract_tmrca = template(inputs=['{smc_file}'], outputs=['{tmrca_file}'], 
#     walltime='00:59:00', memory='4g') << """

#     source activate argweaver
#     arg-extract-tmrca {smc_file} > {tmrca_file}

#     """

# # argweaver: extract TMRCA
# def extract_tmrca(smc_files, tmrca_files):

#     options = {'inputs': smc_files,
#         'outputs': tmrca_files,
#         'memory': '4g',
#         'walltime': '48:00:00'}

#     input_files_prefix = os.path.commonprefix(smc_files)
#     output_files_prefix = os.path.commonprefix(tmrca_files)

#     input_files_suffixes = [x[len(input_files_prefix):] for x in smc_files]
#     output_files_suffixes = [x[len(output_files_prefix):] for x in tmrca_files]

#     shell_spec = """

#     source activate argweaver
#     INPUTFILESUFFIXES=({input_suffixes})
#     OUTPUTFILESUFFIXES=({output_suffixes})
#     MAXINDEX=$((${{#INPUTFILESUFFIXES[@]}}-1))
#     for NR in `seq 0 $MAXINDEX` ; do
#         arg-extract-tmrca {input_prefix}${{INPUTFILESUFFIXES[$NR]}} > {output_prefix}${{OUTPUTFILESUFFIXES[$NR]}}
#     done

#     """.format(input_prefix=input_files_prefix, output_prefix=output_files_prefix,
#         input_suffixes=' '.join(input_files_suffixes), output_suffixes=' '.join(output_files_suffixes))

#     return (options, shell_spec)


# # argweaver: turn arweaver output into bed files
# smc2bed = template(inputs=['{smc_file}'], outputs=['{bed_file}'], 
#     walltime='00:59:00', memory='4g') << """

#     source activate argweaver
#     smc2bed --times {times_file} --sample {sample_nr} {smc_file} | sort-bed - | gzip > {bed_file}

#     """

# # argweaver: turn arweaver output into bed files
# def smc2bed(smc_files, bed_files, sample_nr, times_files):

#     options = {'inputs': smc_files,
#         'outputs': bed_files,
#         'memory': '4g',
#         'walltime': '48:00:00'}

#     input_files_prefix = os.path.commonprefix(smc_files)
#     output_files_prefix = os.path.commonprefix(bed_files)

#     input_files_suffixes = [x[len(input_files_prefix):] for x in smc_files]
#     times_files_suffixes = [x[len(input_files_prefix):] for x in times_files]
#     output_files_suffixes = [x[len(output_files_prefix):] for x in bed_files]

#     shell_spec = """

#     source activate argweaver
#     INPUTFILESUFFIXES=({input_suffixes})
#     TIMESFILESUFFIXES=({times_suffixes})
#     OUTPUTFILESUFFIXES=({output_suffixes})
#     MAXINDEX=$((${{#INPUTFILESUFFIXES[@]}}-1))
#     for NR in `seq 0 $MAXINDEX` ; do
#         smc2bed --times {input_prefix}${{TIMESFILESUFFIXES[$NR]}} --sample {sample_nr} {input_prefix}${{INPUTFILESUFFIXES[$NR]}} | sort-bed - | gzip > {output_prefix}${{OUTPUTFILESUFFIXES[$NR]}}
#     done

#     """.format(input_prefix=input_files_prefix, 
#         output_prefix=output_files_prefix,
#         sample_nr=sample_nr,
#         input_suffixes=' '.join(input_files_suffixes), 
#         times_suffixes=' '.join(times_files_suffixes), 
#         output_suffixes=' '.join(output_files_suffixes))

#     return (options, shell_spec)


# # argweaver: compute various summaries from the sampled args
# summarize_arg = template(inputs=['{bed_file}'], outputs=['{summary_file}'], 
#     walltime='00:59:00', memory='4g') << """

#     source activate argweaver
#     arg-summarize -a {bed_file} --time-file {times_file} -E -T -B -R -K -H -F -P -C > {summary_file}

#     """

# # argweaver: compute various summaries from the sampled args
# def summarize_arg(bed_files, summary_files, times_files):

#     options = {'inputs': bed_files,
#         'outputs': summary_files,
#         'memory': '4g',
#         'walltime': '48:00:00'}

#     input_files_prefix = os.path.commonprefix(bed_files)
#     output_files_prefix = os.path.commonprefix(summary_files)

#     input_files_suffixes = [x[len(input_files_prefix):] for x in bed_files]
#     times_files_suffixes = [x[len(input_files_prefix):] for x in times_files]
#     output_files_suffixes = [x[len(output_files_prefix):] for x in summary_files]

#     shell_spec = """

#     source activate argweaver
#     INPUTFILESUFFIXES=({input_suffixes})
#     TIMESFILESUFFIXES=({times_suffixes})
#     OUTPUTFILESUFFIXES=({output_suffixes})
#     MAXINDEX=$((${{#INPUTFILESUFFIXES[@]}}-1))
#     for NR in `seq 0 $MAXINDEX` ; do
#         arg-summarize -a {input_prefix}${{INPUTFILESUFFIXES[$NR]}} --time-file {input_prefix}${{TIMESFILESUFFIXES[$NR]}} -E -T -B -R -K -H -F -P -C > {output_prefix}${{OUTPUTFILESUFFIXES[$NR]}}
#     done

#     """.format(input_prefix=input_files_prefix, 
#         output_prefix=output_files_prefix,
#         input_suffixes=' '.join(input_files_suffixes), 
#         times_suffixes=' '.join(times_files_suffixes), 
#         output_suffixes=' '.join(output_files_suffixes))

#     return (options, shell_spec)


# compute additional tree statistics and write hdf5 stores
# summary_hdf_store = template(inputs=['{sample_dir}'], outputs=['{store_file}'], 
#     walltime='00:59:00', memory='4g') << """

#     python ./scripts/tree_statistics.py {sample_dir} {store_file}

#     """


def argweaver_extra_stats(input_file_name, output_file_name, excluded_pops, excluded_indivs):

    options = {'memory': '500m',
               'walltime': '24:00:00'}

    spec = """

    source activate simons
    python ./scripts/argweaver_subset_tmrca_stats.py --exclpops {pops} --exclindivs {indivs} {input_file} {output_file}

    """.format(input_file=input_file_name, output_file=output_file_name, 
        pops=excluded_pops, indivs=excluded_indivs)

    return [input_file_name], [output_file_name], options, spec


def compute_tmrca_window_stats(input_table_file, output_tmrca_file):

    options = {'memory': '10g',
               'walltime': '00:10:00'}

    spec = """

    source activate simons
    python ./scripts/tmrca_in_windows.py {input_table_file} {output_tmrca_file}

    """.format(input_table_file=input_table_file, output_tmrca_file=output_tmrca_file)

    return [input_table_file], [output_tmrca_file], options, spec



def compute_extra_tmrca_window_stats(input_table_file, output_tmrca_file):

    options = {'memory': '10g',
               'walltime': '00:10:00'}

    spec = """

    source activate simons
    python ./scripts/tmrca_extra_in_windows.py {input_table_file} {output_tmrca_file}

    """.format(input_table_file=input_table_file, output_tmrca_file=output_tmrca_file)

    return [input_table_file], [output_tmrca_file], options, spec



def compute_tree_stats(input_table_file, output_hdf_file, component_hdf_file, component_stats_hdf_file):

    # FIXME: find a way to get discretization without hardcoding it here
    # best way would be to pass a times file to the template as argument.
    # I should make the argweaver template keep the times file and make it an output file.
    # That way the --discretization argument to tree_statistics.py could be the times file 
    # and not the long list of values.

    discretization = [0.000000,49.193481,122.586947,232.085215,395.449492,639.178343,
                      1002.805899,1545.314509,2354.701987,3562.255340,5363.846221,
                      8051.702368,12061.808515,18044.625462,26970.598323,40287.567936,
                      60155.618452,89797.454603,134021.141756,200000.000000]

    options = {'memory': '50g',
               'walltime': '24:00:00'}

    spec = """
    source activate simons
    python ./scripts/tree_statistics.py \
    --pop-size 10000 \
    --min-post-prob 0.9 \
    --min-clade-score 3 \
    --min-clade-size 5 \
    --discretization {discretization} \
    {input_table_file} {output_hdf_file} {component_hdf_file} {component_stats_hdf_file} \
    {sweep_sister_clade_hdf_file} {nonsweep_sister_clade_hdf_file}

    """.format(discretization=','.join(map(str, discretization)), 
               input_table_file=input_table_file, 
               output_hdf_file=output_hdf_file,
               component_hdf_file=component_hdf_file,
               component_stats_hdf_file=component_stats_hdf_file,
               sweep_sister_clade_hdf_file=sweep_sister_clade_hdf_file,
               nonsweep_sister_clade_hdf_file=nonsweep_sister_clade_hdf_file)

    return [input_table_file], [output_hdf_file, component_hdf_file, component_stats_hdf_file], options, spec



def slim_sim(selcoef, slim_file, output_prefix):

    options = {'memory': '18g',
               'walltime': '00:30:00'
              } 

    tree_output_file = output_prefix + '.trees'
    hdf_output_file = output_prefix + '.hdf'

    spec = """

    source activate simons
    python scripts/slim_trees.py --selcoef {selcoef} --samples 100 \
        --window 100000 {slim_file} {tree_output} {hdf_output}

    """.format(selcoef=selcoef, slim_file=slim_file, 
            tree_output=tree_output_file, hdf_output=hdf_output_file)

    return [], [tree_output_file, hdf_output_file], options, spec


# def summary_hdf_store(summary_files, store_file):

#     options = {'inputs': summary_files,
#         'outputs': [store_file],
#         'memory': '4g',
#         'walltime': '24:00:00'}

#     shell_spec = """

#     python ./scripts/tree_statistics.py --store {store_file} {summary_files}

#     """.format(store_file=store_file, summary_files=' '.join(summary_files))

#     return (options, shell_spec)