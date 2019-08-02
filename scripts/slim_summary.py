
import os, sys
import numpy
import pandas
from pandas import DataFrame, Series
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('sweep_data_dir', type=Path, help='')
parser.add_argument('output_file', type=Path, help='')
args = parser.parse_args()

# lst = list()
# for file_path in args.sweep_data_dir.glob('**/*.hdf'):
#     #sweep_type, pop_size, bottle_start, bottle_end, bottle_pop_size, sweep_start, int(selcoef*100), replication
#     simulation, selection_percent, replication = \
#         file_path.with_suffix('').name.rsplit('_', maxsplit=2)
#     # simulation, replication, selection_percent = \
#     #     file_path.with_suffix('').name.rsplit('_', maxsplit=2)

#     sweep_data = pandas.read_hdf(str(file_path))

#     df = (sweep_data
#            .groupby(['start', 'end'])['swept']
#            .aggregate(['sum', 'size'])
#            .rename(columns={'sum': 'nr_swept', 'size': 'total'})
#            .reset_index(level=['start', 'end'])
#           )
#     df['simulation'] = simulation
#     df['replication'] = replication
#     df['selection_coef'] = int(selection_percent) / 100

#     lst.append(df)

# pandas.concat(lst).to_hdf(args.output_file, 'df', format='table', mode='w')



lst = list()
for dir_path in args.sweep_data_dir.iterdir():

    if not dir_path.is_dir():
        continue

    #sweep_type, pop_size, bottle_start, bottle_end, bottle_pop_size, sweep_start, int(selcoef*100), replication
    simulation, selection_percent, replication =  dir_path.name.rsplit('_', maxsplit=2)

    for file_path in dir_path.glob('*.hdf'):

        pwdist_cutoff = float(file_path.name.split('_')[-1])

        sweep_data = pandas.read_hdf(str(file_path))

        df = (sweep_data
            .groupby(['start', 'end'])['swept']
            .aggregate(['sum', 'size'])
            .rename(columns={'sum': 'nr_swept', 'size': 'total'})
            .reset_index(level=['start', 'end'])
            )
        df['simulation'] = simulation
        df['replication'] = replication
        df['selection_coef'] = int(selection_percent) / 100
        df['pwdist_cutoff'] = pwdist_cutoff

        lst.append(df)

pandas.concat(lst).to_hdf(args.output_file, 'df', format='table', mode='w')



