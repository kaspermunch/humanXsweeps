
import sys, os
from pathlib import Path
import argparse
import pickle
import pandas
from pandas import DataFrame, Series

from hg19_chrom_sizes import hg19_chrom_sizes as chromosome_lengths

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_dir + '/../notebooks')
import analysis_globals

parser = argparse.ArgumentParser()
parser.add_argument("--dist-dir", dest="dist_dir", type=Path)
parser.add_argument("--out-file", dest="out_file", type=Path)
parser.add_argument("--dist-twice-out-file", dest="dist_twice_out_file", type=Path)
# parser.add_argument("--result-dir", dest="result_dir", type=Path)
# parser.add_argument("--result-file-prefix", dest="result_file_prefix", type=str, default='dist_data')
args = parser.parse_args()

def read_dist_table(file_name):

    col_names = ['chrom', 'start', 'end', 'pop_label',
             'indiv_1', 'pseud_1', 'indiv_2', 'pseud_2',
             'dist', 'mismatch', 'match', 
             'dist_af', 'mismatch_af', 'match_af', 
             'uncalled']
    with open(str(file_name), 'rb') as f:
        table = pickle.load(f)
    df = DataFrame(table, columns=col_names)
    df.indiv1 = [Path(x).name for x in df.indiv_1]
    df.indiv2 = [Path(x).name for x in df.indiv_2]
    return df

    
dist_data = (pandas.concat(map(read_dist_table, args.dist_dir.iterdir())) # read and concat pi tables
        #    # filter individuals
        #    .loc[indiv_filter]
           .reset_index(drop=True)
              )

# sort
dist_data.sort_values(by=['chrom', 'indiv_1', 'start'], inplace=True)

# # write stores for each chromosome
# groups = dist_data.groupby('chrom')
# for chrom, group in groups:
# #    store_path = args.result_dir / 'dist_data_chr{}_100kb.store'.format(chrom)
#     store_path = args.result_dir / '{}_chr{}_100kb.store'.format(args.result_file_prefix, chrom)
#     group.to_hdf(str(store_path), 'df',  mode='w', format="table", data_columns=['indiv_1', 'start', 'end']) # we index all data columns 

dist_data.to_hdf(str(args.out_file), 'df',  mode='w', format="table", data_columns=['indiv_1', 'start', 'end']) # we index all data columns 

##### copied this over from g1000_sweep_calling.py
def dist_twice(dist_data):

    dist_data.drop('pop_label', axis=1, inplace=True)

    # no filering of individuals, we use all them.

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


##### apply function and write hdf
all_male_dist_twice = dist_twice(dist_data)

all_male_dist_twice.to_hdf(str(args.dist_twice_out_file), 'df', 
                           data_columns=['start', 'end', 
                                     'indiv_1', 'indiv_2'],
                           format='table', mode='w')
