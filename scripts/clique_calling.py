import os, sys
import numpy
import pandas
from pandas import DataFrame, Series
import argparse
from multiprocessing import Pool, cpu_count

import gc

import simons_meta_data
import hg19_chrom_sizes

import networkx as nx
import itertools

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_dir + '/../notebooks')
import analysis_globals


def call_rolling_windows(df, pwdist_cutoff, min_sweep_clade_size):
    """
    Input data frame holds windows for all individuals in a 500kb window.
    """

#     """
#     dist values may be nan if there are too many uncalled bases, but may also
#     be nan if there is no unmasked bases left after admixture masking. We only
#     want to call the 500kb (and the constituent 100kb windows) as swept if all
#     five 100kb windows are not uncalled. But we want to allow nans if these are
#     induced by admixture masking.
#     """
    
#     def mean_dist(df, dist_col):
#         if (df.uncalled > analysis_globals.max_uncalled_bases).any():
#             # there are uncalled windows in this 500kb so we report 
#             # it as nan to indicdate it cannot be called as swept
#             return numpy.nan
#         else:
#             # there are no uncalled windows so we ingnore any windows 
#             # with nan dist that are induced by missing data after
#             # admixture filtering
#             if df[dist_col].isnull().sum() > 2:
#                 # but return nan if there are more than TWO such windows
#                 return numpy.nan
#             return df[dist_col].mean()


    # build data frame to return
    cols = ['indiv_1', 'start', 'end', 'off']
    result_df = pandas.DataFrame(list(df.groupby(cols).groups.keys()), columns=cols)

    ## NEW: #############################

#     # Get mean dist between pairs in the 500kb window:
#     win_dist = df.groupby(['indiv_1', 'indiv_2']).apply(mean_dist, 'dist')

    ###############################

    # Get mean dist between pairs in the 500kb window (make mean nan if any 100kb is nan):
#     win_dist = df.groupby(['indiv_1', 'indiv_2']).dist.agg(numpy.mean)
    win_dist = df.groupby(['indiv_1', 'indiv_2']).dist.agg(lambda sr: numpy.mean(sr.values))
    ###############################

    # Build graph
    graph = nx.Graph()
    for tup in win_dist[win_dist <= pwdist_cutoff].reset_index().itertuples():
        graph.add_edge(tup.indiv_1, tup.indiv_2)
        
    cliques = sorted(nx.algorithms.clique.find_cliques(graph), key=len)

    # find swept
    result_df['called'] = False
    result_df['mean_clade_dist'] = numpy.nan
    result_df['clade_size'] = numpy.nan
    
    for clique in cliques:

        if len(clique) >= min_sweep_clade_size:
            
            # get indiv pairs. (in just one orientation). that is fine because df has both orientations

            called = result_df.indiv_1.isin(clique)
            #result_df['called'] = called | result_df.called # do not make it False if it is already True   
            
            result_df.loc[called, 'called'] = True
            
            # clade_size is the size of the clique
            result_df.loc[called, 'clade_size'] = len(clique) 
            
            # mean_clade_dist is the mean dist of the clique
            indiv_pairs = list(itertools.combinations(sorted(clique), 2))
            result_df.loc[called, 'mean_clade_dist'] = win_dist.loc[indiv_pairs].mean()                


    if 'dist_af' in df.columns:

        #### NEW: #########################

#         # Get mean dist between pairs in the 500kb window (make mean nan if any 100kb is nan):
#         win_dist = df.groupby(['indiv_1', 'indiv_2']).apply(mean_dist, 'dist_af')

        #############################
        
        # Get mean dist between pairs in the 500kb window (make mean nan if any 100kb is nan):
#         win_dist = df.groupby(['indiv_1', 'indiv_2']).dist_af.agg(numpy.mean)
        win_dist = df.groupby(['indiv_1', 'indiv_2']).dist_af.agg(lambda sr: numpy.mean(sr.values))
        #############################

        # Build graph
        graph = nx.Graph()
        for tup in win_dist[win_dist <= pwdist_cutoff].reset_index().itertuples():
            graph.add_edge(tup.indiv_1, tup.indiv_2)

        cliques = sorted(nx.algorithms.clique.find_cliques(graph), key=len)

        # find swept
        result_df['called_af'] = False
        result_df['mean_clade_dist_af'] = numpy.nan
        result_df['clade_size_af'] = numpy.nan

        for clique in cliques:

            if len(clique) >= min_sweep_clade_size:
                # get indiv pairs. (in just one orientation). that is fine because df has both orientations

                called = result_df.indiv_1.isin(clique)

                result_df.loc[called, 'called_af'] = True

                # clade_size is the size of the clique
                result_df.loc[called, 'clade_size_af'] = len(clique) 

                # mean_clade_dist is the mean dist of the clique
                indiv_pairs = list(itertools.combinations(sorted(clique), 2))
                result_df.loc[called, 'mean_clade_dist_af'] = win_dist.loc[indiv_pairs].mean()

    return result_df


def call_swept(df):
    """
    Takes a df with all rolling window data for an indivisual for one 100kb window.
    Call each 100kb window as sweept if any overlapping rolling window is called as swept.
    Compute clade size and mean clade dist as from the rolling window with the largest clade size.
    """
    max_clade_size = df.clade_size.max()
    offset_with_largest_clique = (df.groupby('off')
                             .filter(lambda df: (df.clade_size == max_clade_size).all() and df.called.all())
                            )

    assert len(df) == len(df.off.unique()) # one row for each offset

    if 'called_af' in df.columns:
        max_clade_size_af = df.clade_size_af.max()
        offset_with_largest_clique_af = (df.groupby('off')
                                .filter(lambda df: (df.clade_size_af == max_clade_size_af).all() and df.called_af.all())
                                )

        return DataFrame(dict(called=[df.called.any()], 
                            clade_size=[max_clade_size],
                            clade_mean_dist=[offset_with_largest_clique.mean_clade_dist.mean()], # actually mean over identical numbers
                            called_af=[df.called_af.any()], 
                            clade_size_af=[max_clade_size_af],
                            clade_mean_dist_af=[offset_with_largest_clique_af.mean_clade_dist_af.mean()])) # actually mean over identical numbers
    else:
        return DataFrame(dict(called=[df.called.any()], 
                            clade_size=[max_clade_size],
                            clade_mean_dist=[offset_with_largest_clique.mean_clade_dist.mean()])) # actually mean over identical numbers
    


def _window_stats(df, pwdist_cutoff, min_sweep_clade_size):
    if 'dist_af' in df.columns:
        if set([ 'pop_1', 'region_label_1', 'region_id_1', 'region_1']).issubset(all_male_dist_twice.columns):
            return pandas.DataFrame({
                                    'mean_dist': [df.dist.mean()],
                                    'mean_dist_to_africans': [df.loc[df.region_2 == 'Africa', 'dist'].mean()],
                                    'mean_dist_af': [df.dist_af.mean()],
                                    'mean_dist_to_africans_af': [df.loc[df.region_2 == 'Africa', 'dist_af'].mean()],
                                    'win_swept': (df.dist <= pwdist_cutoff).sum() >= min_sweep_clade_size,
                                    'win_swept_af': (df.dist_af <= pwdist_cutoff).sum() >= 
                                                    min_sweep_clade_size,
                                    'prop_indivs_missing': [numpy.isnan(df.dist).sum() / df.dist.size]
                                    })
        else:
            # 1000 genomes
            return pandas.DataFrame({
                                    'mean_dist': [df.dist.mean()],
                                    'mean_dist_af': [df.dist_af.mean()],
                                    'win_swept': (df.dist <= pwdist_cutoff).sum() >= min_sweep_clade_size,
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

if __name__ == "__main__":

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

    # print('all_male_dist_twice', "Ust_Ishim" in all_male_dist_twice.indiv_1.unique())

    # remove unused columns
    to_keep = ['indiv_1', 'indiv_2', 'start', 'end', 'dist', 'dist_af', 
               'region_1', 'region_2', 'pop_1', 'region_label_1', 'region_id_1', 'uncalled']
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
        if set([ 'pop_1', 'region_label_1', 'region_id_1', 'region_1']).issubset(all_male_dist_twice.columns):
            # this is not 1000 genomes
            gr_cols = ['indiv_1', 'start', 'end', 'pop_1', 'region_label_1', 'region_id_1', 'region_1']
        else:
            gr_cols = ['indiv_1', 'start', 'end']
    else:
        # simulation
        gr_cols = ['indiv_1', 'start', 'end']
    stats_data = (all_male_dist_twice
            .groupby(gr_cols)
            .apply(window_stats)
            .reset_index(level=gr_cols)
            )

    # print('stats_data', "Ust_Ishim" in stats_data.indiv_1.unique())

    lst = list()
    # loop over five offsets of 500kb windows
    for off in offsets:
        print(off)
        groups = (all_male_dist_twice
                    .assign(off=off, # keep offset
                            roll_win = lambda df: (off + df.start) // window_size) # label for rolling 500kb window
                    .groupby(['roll_win', 'off'], as_index=False)

                  ### ADDED THIS TO ENSURE NO CALLS ARE MADE ON WINDOWS SMALLER THAN 500KB:
                    .filter(lambda df: df.start.unique().size == nr_wins)             
                    .groupby(['roll_win', 'off'], as_index=False)

                    )
        # with Pool(nr_cpu) as p:
        #     df = pandas.concat(p.map(call_rolling_windows, [group for name, group in groups]))
        ##### added pwdist arg to function call

        lst.append(groups.apply(call_rolling_windows, args.pwdist_cutoff, MIN_SWEEP_CLADE_SIZE))

    del all_male_dist_twice
    gc.collect()

    # concatenate data frames for each offset and call windows as swept
    sweep_calls = (pandas.concat(lst)
                    .groupby(['indiv_1', 'start', 'end'])
                    .apply(call_swept)
                    .reset_index(level=['indiv_1', 'start', 'end'])
                    )

    del lst[:]
    gc.collect()

    # print('sweep_calls', "Ust_Ishim" in sweep_calls.indiv_1.unique())

    sweep_data = stats_data.merge(sweep_calls, on=['indiv_1', 'start', 'end'])

    del stats_data
    del sweep_calls
    gc.collect()

    # print('sweep_data', "Ust_Ishim" in sweep_data.indiv_1.unique())

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

    # print('sweep_data', "Ust_Ishim" in sweep_data.indiv_1.unique())

    # write to hdf output file
    sweep_data.to_hdf(args.sweep_data_file_name, 'df', mode='w', format='table')


