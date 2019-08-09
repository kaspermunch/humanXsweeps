

import os, sys, re
#from gwf import template

from pathlib import Path

#################################################################################
# Utility functions
#################################################################################

def modpath(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
#    if suffix is not None:
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    # return os.path.join(par, name_no_suffix + suf)
    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path

def get_basename(p):
    return os.path.splitext(os.path.basename(p))[0]


def bp2str(n):
    """
    Convert a number of bases to a human readable string
    """
    if n < 1000:
        return '{}bp'.format(n)
    elif n < 1000000:
        if n % 1000:
            return '{}kb'.format(n/1000)
        else:
            return '{}kb'.format(int(n/1000))
    else:
        if n % 1000000:
            return '{}Mb'.format(n/1000000)
        else:
            return '{}Mb'.format(int(n/1000000))

##############################################################################
# Templates
##############################################################################

def bed_difference(bed_file1, bed_file2, output_file): 

    options = {'memory': '4g',
               'walltime': '11:00:00'
               }
    shell_spec = """

    source activate simons
    python ./scripts/liftover_diff.py {bed_file1} {bed_file2} | sort -k1,1 -k2,2n -k3,3n > {output_file}

    """.format(bed_file1=bed_file1, bed_file2=bed_file2, output_file=output_file)

    return [bed_file1, bed_file2], [output_file], options, shell_spec


def bed_merge_and_split(input_files, output_files):

    output_base_names = [Path(x).name for x in output_files]

    # get output dir and make sure it is unambiguous:
    output_dirs = [Path(x).parent for x in output_files]
    assert len(set(output_dirs)) == 1
    output_dir = output_dirs[0]

    input_files = [str(x) for x in input_files]
    output_files = [str(x) for x in output_files]

    options = {'memory': '10g',
               'walltime': '11:00:00'}

    shell_spec = """

    sort -k1,1 -k2,2n -k3,3n --merge {input_files} -T /scratch/$GWF_JOBID | python ./scripts/bed_split.py {output_dir} {output_base_names}

    """.format(input_files=" ".join(input_files),
               output_dir=output_dir,
               output_base_names=" ".join(output_base_names))

    return input_files, output_files, options, shell_spec


def liftover(bed_file, chain_file, mapped_file, unmapped_file):

    options = {'memory': '4g',
               'walltime': '36:00:00'
               }

    spec = """

    source activate simons
    liftOver -bedPlus=3 {bed_file} {chain_file} {mapped_file} {unmapped_file}

    """.format(bed_file=bed_file, chain_file=chain_file, mapped_file=mapped_file, unmapped_file=unmapped_file)

    return [bed_file, chain_file], [mapped_file, unmapped_file], options, spec


def vcf2haplo(vcf_file, sample_id, masked_ref, out_file1, haploid=False, out_file2=None):

    options = {'memory': '8g',
               'walltime': '3:00:00'
              }

    outputs = [out_file1]
    out_args = '--out1 ' + out_file1
    if out_file2 is not None:
        outputs += [out_file2]
        out_args += ' --out2 ' + out_file2

    haploid_arg = haploid and '--haploid' or ''

    spec = """

    source activate simons
    source /com/extra/vcftools/0.1.14/load.sh
    vcftools --gzvcf {vcf_file} --remove-indels --remove-filtered-all --max-alleles 2 \
        --recode --recode-INFO-all --indv {sample_id} --stdout | \
    python scripts/phased_vcf_to_hapltype.py --sample {sample_id} --masked_ref {ref} {haploid_arg} {out_args}

    """.format(vcf_file=vcf_file, sample_id=sample_id, ref=masked_ref, 
                haploid_arg=haploid_arg, out_args=out_args)

    return [vcf_file, masked_ref], outputs, options, spec


def rz2gz(rz_file, gz_file):

    spec = """

    zcat {rz_file} | gzip > {gz_file}

    """.format(rz_file=rz_file, gz_file=gz_file)

#    return (options, shell_spec)
    return [rz_file], [gz_file], {}, spec

def pad_archaic_files(template_file, input_file, pad_char, output_file):

    options = {'memory': '8g',
               'walltime': '1:00:00'
              }

    spec = """
    source activate simons
    python scripts/pad_archaic_genome.py {template} {input} {pad_char} {output}

    """.format(template=template_file, input=input_file, pad_char=pad_char, output=output_file)

    return [input_file], [output_file], options, spec

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


def admix_masked_dist_for_pair_template(*args):
    """
    same as above but not skipping all but X... assumes only one chromosome...
    """

    f1, f2, binsize, pop, indiv1, pseud1, indiv2, pseud2, output_file_name = args

    options = {'memory': '4g',
               'walltime': '11:00:00'
               } 

    spec = '''

        source activate simons
        python scripts/admix_masked_diff_windows.py {args}

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



# def slim_sim(selcoef, gen_time, mut_per_year, samples, slim_file, slim_tree_file, slim_dist_file):

#     options = {'memory': '18g',
#                'walltime': '00:30:00'
#               } 

#     spec = """

#     source activate simons
#     python scripts/slim_trees.py --selcoef {selcoef} --samples {samples} \
#         --window 100000 --generationtime {gen_time} --mutationrate {mut_per_year} {slim_file} {tree_output} {hdf_output}

#     """.format(selcoef=selcoef, gen_time=gen_time, mut_per_year=mut_per_year, samples=samples,
#             slim_file=slim_file, tree_output=slim_tree_file, hdf_output=slim_dist_file)

#     return [], [slim_tree_file, slim_dist_file], options, spec



def slim_sim(selcoef, gen_time, mut_per_year, rec_rate, samples, sweep_type, sweep_start, demography,
 chrom, xreduction, total_sim_generations, slim_tree_file, slim_dist_file):

    options = {'memory': '8g',
               'walltime': '05:00:00'
              } 

    spec = """

    source activate simons
    python scripts/slim_trees.py --selcoef {selcoef} --samples {samples} \
        --window 100000 --generationtime {gen_time} \
        --mutationrate {mut_per_year} --recrate {rec_rate} \
        --sweep {sweep_type} --sweepstart {sweep_start} \
        --demography {demography} \
        --chrom {chrom} --x-size-reduction {xreduction} \
        --totalgenerations {total_sim_generations} \
        {tree_output} {hdf_output}
            
    """.format(selcoef=selcoef, gen_time=gen_time, 
            mut_per_year=mut_per_year, rec_rate=rec_rate,
            samples=samples, sweep_type=sweep_type, sweep_start=sweep_start, 
            total_sim_generations=total_sim_generations,
            demography=' '.join("{}:{}".format(*pair) for pair in demography),
            chrom=chrom, xreduction=xreduction,
            tree_output=slim_tree_file, hdf_output=slim_dist_file)

    return [], [slim_tree_file, slim_dist_file], options, spec


def slim_dist_twice(dist_file, dist_twice_file):

    options = {'memory': '8g',
               'walltime': '01:00:00'
              } 

    spec = """

    source activate simons
    python scripts/slim_dist_to_dist_twice.py {dist_file} {dist_twice_file}

    """.format(dist_file=dist_file, dist_twice_file=dist_twice_file)

    return [dist_file], [dist_twice_file], options, spec


# def sweep_data(dist_file, sweep_data_file, dump_dist_twice=None, cores=1, memory='100g', walltime='48:00:00'):
#     """
#     Make and dump data frame with all distances in both directions (indiv_1, indiv2 and indiv_2, indiv_1). 
#     Then calls sweeps.
#     """

#     options = {'memory': memory,
#                'walltime': walltime,
#                'cores': cores,
#               } 

#     out_files = [sweep_data_file]
#     extra_arg = ''
#     if dump_dist_twice:
#         extra_arg = '--dump-dist-twice {}'.format(dump_dist_twice)
#         out_files += [dump_dist_twice]
    
#     spec = """

#     source activate simons
#     python scripts/sweep_calling.py {extra_arg} {dist_file} {sweep_data_file}

#     """.format(extra_arg=extra_arg, dist_file=dist_file, sweep_data_file=sweep_data_file)

#     return [dist_file], out_files, options, spec


def sweep_data(dist_file, sweep_data_file, 
               min_sweep_clade_percent, pwdist_cutoff, 
               cores=1, memory='60g', walltime='4:00:00'):
    """
    Make and dump data frame with all distances in both directions (indiv_1, indiv2 and indiv_2, indiv_1). 
    Then calls sweeps.
    """

    options = {'memory': memory,
               'walltime': walltime,
               'cores': cores,
              } 

    spec = """

    source activate simons
    python scripts/sweep_calling.py \
            --pwdist-cutoff {pwdist_cutoff} \
            --min-sweep-clade-percent {min_sweep_clade_percent} \
            {dist_file} {sweep_data_file}

    """.format(dist_file=dist_file, sweep_data_file=sweep_data_file,
              min_sweep_clade_percent=min_sweep_clade_percent, pwdist_cutoff=pwdist_cutoff)

    return [dist_file], [sweep_data_file], options, spec

# def g1000_sweep_data(dist_file, sweep_data_file, dump_dist_twice):
#     """
#     Same as above but without adding meta data/
#     """

#     options = {'memory': '10g',
#                'walltime': '5:00:00',
#                'cores': 1,
#               } 

#     spec = """

#     source activate simons
#     python scripts/g1000_sweep_calling.py --dump-dist-twice {dump_dist_twice} {dist_file} {sweep_data_file}

#     """.format(dump_dist_twice=dump_dist_twice, dist_file=dist_file, sweep_data_file=sweep_data_file)

#     return [dist_file], [sweep_data_file, dump_dist_twice], options, spec

def g1000_sweep_data(dist_file, sweep_data_file, min_sweep_clade_percent, pwdist_cutoff):
    """
    Same as above but without adding meta data/
    """

    options = {'memory': '16g',
               'walltime': '2:00:00',
               'cores': 1,
              } 

    spec = """

    source activate simons
    python scripts/g1000_sweep_calling.py \
            --pwdist-cutoff {pwdist_cutoff} \
            --min-sweep-clade-percent {min_sweep_clade_percent} \
            {dist_file} {sweep_data_file}

    """.format(dist_file=dist_file, sweep_data_file=sweep_data_file, 
               min_sweep_clade_percent=min_sweep_clade_percent, pwdist_cutoff=pwdist_cutoff)

    return [dist_file], [sweep_data_file], options, spec


def g1000_fst(vcf_file, pop_files, out_file):
    """
    Computes weir fst between populations
    """

    options = {'memory': '10g',
               'walltime': '5:00:00',
              } 

    fst_args = " ".join("--weir-fst-pop " + pop_file for pop_file in pop_files)

    spec = """
    source activate simons
    source /com/extra/vcftools/0.1.14/load.sh
    vcftools --gzvcf {vcf_file} --remove-indels --remove-filtered-all --max-alleles 2 \
        --fst-window-size 100000 {fst_args} --out {out_file}
    """.format(vcf_file=vcf_file, fst_args=fst_args, out_file=out_file.replace('.windowed.weir.fst', ''))

    return [vcf_file] + pop_files, [out_file], options, spec



# def summary_hdf_store(summary_files, store_file):

#     options = {'inputs': summary_files,
#         'outputs': [store_file],
#         'memory': '4g',
#         'walltime': '24:00:00'}

#     shell_spec = """

#     python ./scripts/tree_statistics.py --store {store_file} {summary_files}

#     """.format(store_file=store_file, summary_files=' '.join(summary_files))

#     return (options, shell_spec)
