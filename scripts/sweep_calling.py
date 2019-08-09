import os, sys
import numpy
import pandas
from pandas import DataFrame, Series
import argparse
from multiprocessing import Pool, cpu_count

import gc
import psutil
process = psutil.Process(os.getpid())

import simons_meta_data
import hg19_chrom_sizes

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_dir + '/../notebooks')
import analysis_globals


def call_rolling_windows(df, pwdist_cutoff, min_sweep_clade_size):
    """
    Takes a df with all pwdiffs in a 500kb rolling window between 
    one indiv and all other individuals. Returns all nan if one or more 100kb 
    windows are without any data. Computes mean pwdist across 
    the five 100kb windows between each pair. Calls sweep_clade as number 
    of other indivisuals with a pwdist smaller than cutoff. Calls as
    swept if this number is above cutoff. Computes mean_clade_dist as 
    mean pwdist in sweep clade.
    """

    if len(df) == nr_wins or numpy.isnan(df.groupby('start')['dist'].mean()).any():
        # there is not five windows or one or more 100kb has no dist data
        called, clade_size, mean_clade_dist = numpy.nan, numpy.nan, numpy.nan
    else:
        # mean distance between indiv_1 and each indiv_2 for the 500kb window
        pwdiffs = df.groupby(['indiv_2']).dist.mean()
        
        # mean distance within clade
        mean_clade_dist = pwdiffs.loc[pwdiffs <= pwdist_cutoff].mean() 

        # number of indiv_2 closer to indiv_1 than cutoff across the 500kb window
        clade_size = (pwdiffs <= pwdist_cutoff).sum() 

        # call if clade size is larger then cutoff
        called = clade_size >= min_sweep_clade_size

        
    if 'dist_af' in df.columns:
        # this is not a simulation
        
        if len(df) == nr_wins or numpy.isnan(df.groupby('start')['dist_af'].mean()).any():
            # there is not five windows or one or more 100kb has no dist data
            called_af, clade_size_af, mean_clade_dist_af = numpy.nan, numpy.nan, numpy.nan
        else:
            # mean distance between indiv_1 and each indiv_2 for the 500kb window
            pwdiffs_af = df.groupby(['indiv_2']).dist_af.mean()

            # mean distance within clade
            mean_clade_dist_af = pwdiffs_af.loc[pwdiffs_af <= pwdist_cutoff].mean() 

            # number of indiv_2 closer to indiv_1 than cutoff across the 500kb window
            clade_size_af = (pwdiffs_af <= pwdist_cutoff).sum() 

            # call if clade size is larger then cutoff
            called_af = clade_size_af >= min_sweep_clade_size
                
        return df[['indiv_1', 'start', 'off', 'dist', 'dist_af']].assign(called=called, clade_size=clade_size, mean_clade_dist=mean_clade_dist,
                                called_af=called_af, clade_size_af=clade_size_af, mean_clade_dist_af=mean_clade_dist_af)
        # # Try to also downcast off
        # _df['off'] = pandas.to_numeric(_df.off, downcast='unsigned')
        # return _df
    else:
        return df[['indiv_1', 'start', 'off', 'dist']].assign(called=called, clade_size=clade_size, mean_clade_dist=mean_clade_dist) 


def call_swept(df):
    """
    Takes a df with all rolling window data for an indivisual for one 100kb window.
    Call each 100kb window as sweept if any overlapping rolling window is called as swept.
    Compute clade size and mean clade dist as from the rolling window with the largest clade size.
    """
    max_clade_size = df.clade_size.max()
    largest_clade_offsets = (df.groupby('off')
                             .filter(lambda df: (df.clade_size == max_clade_size).all() and df.called.all())
                            )
    if 'dist_af' in df.columns:
        max_clade_size_af = df.clade_size_af.max()
        largest_clade_offsets_af = (df.groupby('off')
                                .filter(lambda df: (df.clade_size_af == max_clade_size_af).all() and df.called_af.all())
                                )
        return DataFrame(dict(called=[df.called.any()], 
                            clade_size=[max_clade_size],
                            clade_mean_dist=[largest_clade_offsets['dist'].mean()],
                            called_af=[df.called_af.any()], 
                            clade_size_af=[max_clade_size_af],
                            clade_mean_dist_af=[largest_clade_offsets_af['dist_af'].mean()]))
    else:
        return DataFrame(dict(called=[df.called.any()], 
                            clade_size=[max_clade_size],
                            clade_mean_dist=[largest_clade_offsets['dist'].mean()]))


def _window_stats(df, pwdist_cutoff, min_sweep_clade_size):
    if 'dist_af' in df.columns:
        return pandas.DataFrame({
                                'mean_dist': [df.dist.mean()],
                                'mean_dist_to_africans': [df.loc[df.region_2 == 'Africa', 'dist'].mean()],
                                'mean_dist_af': [df.dist_af.mean()],
                                'mean_dist_to_africans_af': [df.loc[df.region_2 == 'Africa', 'dist_af'].mean()],
                                'win_swept': (df.dist <= pwdist_cutoff).sum() >= \
                                                min_sweep_clade_size,
                                'win_swept_af': (df.dist_af <= pwdist_cutoff).sum() >= 
                                                min_sweep_clade_size,
                                'prop_indivs_missing': [numpy.isnan(df.dist).sum() / df.dist.size]
                                })
    else:
        return pandas.DataFrame({
                        'mean_dist': [df.dist.mean()],
                        'win_swept': (df.dist <= pwdist_cutoff).sum() >= \
                                        min_sweep_clade_size,
                        'prop_indivs_missing': [numpy.isnan(df.dist).sum() / df.dist.size]
                        })


parser = argparse.ArgumentParser()
parser.add_argument('--nr_wins', type=int, default=5, help='')
parser.add_argument('--offset', type=int, default=100000, help='')
parser.add_argument('--cpus', type=int, help='')
parser.add_argument('--min-sweep-clade-percent', dest='min_sweep_clade_percent', type=int)
parser.add_argument('--pwdist-cutoff', dest='pwdist_cutoff', type=float)                
parser.add_argument('dist_file_name', type=str, help='')
parser.add_argument('sweep_data_file_name', type=str, help='')
args = parser.parse_args()

if args.cpus:
    nr_cpu = args.cpus
else:
    nr_cpu = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE'))

nr_wins = args.nr_wins
offset = args.offset
offsets = [x * offset for x in range(nr_wins)]
window_size = len(offsets) * offset

all_male_dist_twice = pandas.read_hdf(args.dist_file_name)

# remove unused columns
to_keep = ['indiv_1', 'indiv_2', 'start', 'end', 'dist', 'dist_af', 
           'region_1', 'region_2', 'pop_1', 'region_label_1', 'region_id_1']
to_drop = [x for x in all_male_dist_twice.columns if x not in to_keep]
all_male_dist_twice.drop(to_drop, axis=1, inplace=True)
gc.collect()

nr_indiv = all_male_dist_twice.indiv_1.unique().size

MIN_SWEEP_CLADE_SIZE = round(nr_indiv * args.min_sweep_clade_percent / 100)

# partial version of window_stats
import functools
window_stats = functools.partial(_window_stats, pwdist_cutoff=args.pwdist_cutoff, min_sweep_clade_size=MIN_SWEEP_CLADE_SIZE)

# merge window sweep info with distance data
if 'dist_af' in all_male_dist_twice.columns:
    # this is not a simulation
    gr_cols = ['indiv_1', 'start', 'end', 'pop_1', 'region_label_1', 'region_id_1', 'region_1']
else:
    gr_cols = ['indiv_1', 'start', 'end']
stats_data = (all_male_dist_twice
        .groupby(gr_cols)
        .apply(window_stats)
        .reset_index(level=gr_cols)
        )

lst = list()
# loop over five offsets of 500kb windows
for off in offsets:
    groups = (all_male_dist_twice
                .assign(off=off, # keep offset
                        roll_win = lambda df: (off + df.start) // window_size) # label for rolling 500kb window
                .groupby(['indiv_1', 'roll_win', 'off'], as_index=False)
                )
    # with Pool(nr_cpu) as p:
    #     df = pandas.concat(p.map(call_rolling_windows, [group for name, group in groups]))
    ##### added pwdist arg to function call

    lst.append(groups.apply(call_rolling_windows, args.pwdist_cutoff, MIN_SWEEP_CLADE_SIZE))

del all_male_dist_twice
gc.collect()

# concatenate data frames for each offset and call windows as swept
sweep_calls = (pandas.concat(lst)
                .groupby(['indiv_1', 'start'])
                .apply(call_swept)
                .reset_index(level=['indiv_1', 'start'])
                )

del lst[:]
gc.collect()

sweep_data = stats_data.merge(sweep_calls, on=['indiv_1', 'start'])

del stats_data
del sweep_calls
gc.collect()


# get run length of swept windows and call windows as part of sweeps (ECH haplotypes)
def run_id(sr):
    return (sr != sr.shift()).cumsum()

sweep_data['run_id'] = (sweep_data
                        .groupby('indiv_1')['called']
                        .apply(run_id)
                       )
sweep_data['run_length'] = (sweep_data
                            .groupby(['indiv_1', 'run_id'])['run_id']
                            .transform(numpy.size)
                           )
sweep_data['swept'] = numpy.bitwise_and(sweep_data['called'], 
                                        sweep_data['run_length'] >= analysis_globals.min_run_length)

if 'called_af' in sweep_data.columns:
    # this is not a simulation
    sweep_data['run_id_af'] = (sweep_data
                            .groupby('indiv_1')['called_af']
                            .apply(run_id)
                        )
    sweep_data['run_length_af'] = (sweep_data
                                .groupby(['indiv_1', 'run_id_af'])['run_id_af']
                                .transform(numpy.size)
                            )
    sweep_data['swept_af'] = numpy.bitwise_and(sweep_data['called_af'], 
                                            sweep_data['run_length_af'] >= analysis_globals.min_run_length)

gc.collect()

# write to hdf output file
sweep_data.to_hdf(args.sweep_data_file_name, 'df', mode='w', format='table')