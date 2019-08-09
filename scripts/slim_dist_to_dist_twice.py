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

dist_data = pandas.read_hdf(args.dist_file)

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

def optimize_data_frame(df, down_int='integer'):
    # down_int can be 'unsigned'
    
    converted_df = pandas.DataFrame()

    floats_optim = (df
                    .select_dtypes(include=['float'])
                    .apply(pandas.to_numeric,downcast='float')
                   )
    converted_df[floats_optim.columns] = floats_optim

    ints_optim = (df
                    .select_dtypes(include=['int'])
                    .apply(pandas.to_numeric,downcast=down_int)
                   )
    converted_df[ints_optim.columns] = ints_optim

    for col in df.select_dtypes(include=['object']).columns:
        num_unique_values = len(df[col].unique())
        num_total_values = len(df[col])
        if num_unique_values / num_total_values < 0.5:
            converted_df[col] = df[col].astype('category')
        else:
            converted_df[col] = df[col]

    unchanged_cols = df.columns[~df.columns.isin(converted_df.columns)]
    converted_df[unchanged_cols] = df[unchanged_cols]

    # keep columns order
    converted_df = converted_df[df.columns]      
            
    return converted_df

dist_data_twice = optimize_data_frame(dist_data_twice, down_int='unsigned')

dist_data_twice.to_hdf(str(args.dist_twice_file), 'df', 
                           data_columns=['start', 'end', 
                                     'indiv_1', 'indiv_2'],
                           format='table', mode='w'
)
