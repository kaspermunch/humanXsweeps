
import os, sys, re
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
#    simulation, selection_percent, replication =  dir_path.name.rsplit('_', maxsplit=2)


    demography, size_reduction, rec_mode, rec_val, chrom, sweep_type, sweep_start, selcoef, replication = \
            dir_path.name.split('_')

    size_reduction = float(size_reduction) / 100
    # rec_rate_per_gen = float(rec_rate_per_gen) / 1e12
    rec_val = float(rec_val) / 1e12
    sweep_start = int(sweep_start)
    selcoef = float(selcoef) / 100
    replication = int(replication)
        
    for sub_dir_path in dir_path.iterdir():
            
        pwdist_cutoff = float(sub_dir_path.name)
            
        for file_path in sub_dir_path.glob('clique_data*.hdf'):
            
            min_clade_percent = int(re.search(r'(\d+)%.hdf', file_path.name).group(1))
         
            sweep_data = pandas.read_hdf(str(file_path))

            df = (sweep_data
                    .assign(swept=lambda df: df.swept.astype(int)) # groupby does not reliably sum booleans (return False when less than two are True)
                    .groupby(['start', 'end'])#['swept']
                    .swept
                    .aggregate(['sum', 'size'])
                    .rename(columns={'sum': 'nr_swept', 'size': 'total'})
                    .reset_index()
                    )
            df['demography'] = demography
            df['chrom'] = chrom
            df['size_reduction'] = size_reduction
            df['rec_mode'] = rec_mode
            df['rec_val'] = rec_val
            df['sweep_start'] = sweep_start
            df['sweep_type'] = sweep_type
            df['selcoef'] = selcoef
            df['replication'] = replication
            df['min_clade_percent'] = min_clade_percent
            
            lst.append(df)

pandas.concat(lst).to_hdf(args.output_file, 'df', format='table', mode='w')



