import os, sys
import numpy
import pandas
from pandas import DataFrame, Series
import argparse
from multiprocessing import Pool, cpu_count

import hg19_chrom_sizes

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_dir + '/../notebooks')
import analysis_globals

def dist_twice(dist_data):

    dist_data.drop('pop_label', axis=1, inplace=True)

    # dict for swapping columns
    swap_dict = dict()
    for colname in dist_data.columns.values:
        if colname.endswith('_1'):
            swap_dict[colname] = colname[:-2] + '_2'
        if colname.endswith('_2'):
            swap_dict[colname] = colname[:-2] + '_1'

    cols = ['start', 'end', 'indiv_1', 'indiv_2', 
            'dist', 'mismatch', 'match', 
            'dist_af', 'mismatch_af', 'match_af',
            'uncalled']

    dist_data_twice = (pandas.concat([dist_data[cols],
                                      dist_data[cols].rename(columns=swap_dict)])
        .sort_values(['indiv_1', 'start'])
        .reset_index(drop=True)
        )
    
    # mask uncalled windows:
    dist_data_twice.dist.where(dist_data_twice.uncalled <= analysis_globals.max_uncalled_bases, 
                               inplace=True)
    dist_data_twice.dist_af.where(dist_data_twice.uncalled <= analysis_globals.max_uncalled_bases, 
                                inplace=True)

    return dist_data_twice


def call_rolling_windows(df):#, pwdist_cutoff, min_sweep_clade_size):
    """
    Takes a df with all pwdiffs in a 500kb rolling window between 
    one indiv and all other individuals. Returns all nan if one or more 100kb 
    windows are without any data. Computes mean pwdist across 
    the five 100kb windows between each pair. Calls sweep_clade as number 
    of other indivisuals with a pwdist smaller than cutoff. Calls as
    swept if this number is above cutoff. Computes mean_clade_dist as 
    mean pwdist in sweep clade.
    """
    
    def mean_indiv_dist(df, col):
        """
        Compute mean across 100kb windows in in 500kb window.
        """
        if len(df) != nr_wins:
            return numpy.nan
        return df[col].mean()
    
    if numpy.isnan(df.groupby('start')['dist'].mean()).any():
        # one or more 100kb has no dist data
        called, clade_size, mean_clade_dist = numpy.nan, numpy.nan, numpy.nan
    else:
        # mean distance between indiv_1 and each indiv_2 for the 500kb window
        pwdiffs = df.groupby(['indiv_2']).apply(mean_indiv_dist, 'dist')

        mean_clade_dist = pwdiffs.loc[pwdiffs <= analysis_globals.pwdist_cutoff].mean() 
        
        # number of indiv_2 closer to indiv_1 than cutoff across the 500kb window
        clade_size = (pwdiffs <= analysis_globals.pwdist_cutoff).sum() 
        
        # call if clade size is larger then cutoff
        called = clade_size >= analysis_globals.min_sweep_clade_size
    

    if numpy.isnan(df.groupby('start')['dist_af'].mean()).any():     
        called_af, clade_size_af, mean_clade_dist_af = numpy.nan, numpy.nan, numpy.nan
    else:
        pwdiffs_af = df.groupby(['indiv_2']).apply(mean_indiv_dist, 'dist_af')

        mean_clade_dist_af = pwdiffs.loc[pwdiffs_af <= analysis_globals.pwdist_cutoff].mean() 

        clade_size_af = (pwdiffs_af <= analysis_globals.pwdist_cutoff).sum()

        called_af = clade_size_af >= analysis_globals.min_sweep_clade_size

    return df.copy().assign(called=called, clade_size=clade_size, mean_clade_dist=mean_clade_dist,
                            called_af=called_af, clade_size_af=clade_size_af, mean_clade_dist_af=mean_clade_dist_af)


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


def window_stats(df):

    return pandas.DataFrame({
                            'mean_dist': [df.dist.mean()],
                            'mean_dist_af': [df.dist_af.mean()],
                            'win_swept': (df.dist <= analysis_globals.pwdist_cutoff).sum() >= \
                                            analysis_globals.min_sweep_clade_size,
                            'win_swept_af': (df.dist_af <= analysis_globals.pwdist_cutoff).sum() >= 
                                            analysis_globals.min_sweep_clade_size,
                            'prop_indivs_missing': [numpy.isnan(df.dist).sum() / df.dist.size]
                            })


parser = argparse.ArgumentParser()
parser.add_argument('--nr_wins', type=int, default=5, help='')
parser.add_argument('--offset', type=int, default=100000, help='')
parser.add_argument('--cpus', type=int, help='')
parser.add_argument('--dump-dist-twice', dest='dump_dist_twice', type=str, default=None, help='')
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

male_dist_data = pandas.read_hdf(args.dist_file_name)

all_male_dist_twice = dist_twice(male_dist_data)

lst = list()
# loop over five offsets of 500kb windows
for off in offsets:
    groups = (all_male_dist_twice
                .assign(off=off, # keep offset
                        roll_win = lambda df: (off + df.start) // window_size) # label for rolling 500kb window
                .groupby(['indiv_1', 'roll_win', 'off'])
                )
    # with Pool(nr_cpu) as p:
    #     df = pandas.concat(p.map(call_rolling_windows, [group for name, group in groups]))
    df = pandas.concat([call_rolling_windows(group) for name, group in groups])

    lst.append(df)

# concatenate data frames for each offset and call windows as swept
sweep_calls = (pandas.concat(lst)
                .groupby(['indiv_1', 'start'])
                .apply(call_swept)
                .reset_index(level=['indiv_1', 'start'])
                )

# merge window sweep info with distance data
gr_cols = ['indiv_1', 'start', 'end']
df = (all_male_dist_twice
        .groupby(gr_cols)
        .apply(window_stats)
        .reset_index(level=gr_cols)
        )
sweep_data = df.merge(sweep_calls, on=['indiv_1', 'start'])


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


if args.dump_dist_twice:
    all_male_dist_twice.to_hdf(args.dump_dist_twice, 'df', 
                               data_columns=['start', 'end', 
                                         'indiv_1', 'indiv_2'],
                               format='table', mode='w')

    
# write to hdf output file
sweep_data.to_hdf(args.sweep_data_file_name, 'df', format='table', mode='w')
