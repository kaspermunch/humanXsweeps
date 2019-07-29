import sys, os
from pathlib import Path
import argparse
import pandas
from pandas import DataFrame, Series

import simons_meta_data

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_dir + '/../notebooks')
import analysis_globals

parser = argparse.ArgumentParser()
parser.add_argument("dist_file", type=Path)
parser.add_argument("dist_twice_file", type=Path)
args = parser.parse_args()

dist_data = pd.read_hdf(args.dist_file)

individuals, populations, regions = simons_meta_data.get_meta_data(meta_data_dir=analysis_globals.meta_data_dir)

# dict for swapping columns
swap_dict = dict()
for colname in dist_data.columns.values:
    if colname.endswith('_1'):
        swap_dict[colname] = colname[:-2] + '_2'
    if colname.endswith('_2'):
        swap_dict[colname] = colname[:-2] + '_1'

cols =['start', 'end', 'indiv_1', 'indiv_2', 'dist']

dist_data_twice = (pandas.concat([dist_data[cols],
                                    dist_data[cols].rename(columns=swap_dict)])
    .sort_values(['indiv_1', 'start'])
    .reset_index(drop=True)
    )

dist_data_twice.to_hdf(str(args.dist_twice_file), 'df', 
                           data_columns=['start', 'end', 
                                     'indiv_1', 'indiv_2'],
                           format='table', mode='w'
)