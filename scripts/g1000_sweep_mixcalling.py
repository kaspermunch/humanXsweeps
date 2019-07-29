import os, sys
import re
import argparse
from pathlib import Path
import numpy
import pandas as pd
import numpy as np
from pandas import DataFrame, Series

numpy.random.seed(7)

parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=str, help='')
parser.add_argument('--pop', type=str, help='')
parser.add_argument('input_dir_path', type=Path, help='')
parser.add_argument('output_file_path', type=Path, help='')
args = parser.parse_args()

# read in calls with different min clade sizes
df_list = list()
for path in args.input_dir_path.glob('sweep_data_*.hdf'):
    dist_cut, clade_cut = map(float, re.search(r'sweep_data_([^_]+)_([^_]+)%.hdf', path.name).groups())
    clade_cut /= 100
    df = pd.read_hdf(path)
    df['dist_cut'] = dist_cut
    df['clade_cut'] = clade_cut
    df['chrom'] = args.chrom
    df['population'] = args.pop
    df_list.append(df)    
df = pd.concat(df_list)

def get_best_sweep_call(grp):

    if grp.swept.sum():
        # extract the call with the largets min clade size
        return grp.loc[(grp.swept == True) & (grp.clade_cut == grp.clade_cut.max())]
    else:
        # if no sweeps are calle for any clade size, we use the extract calls for the smallest min clade size
        return grp.loc[grp.clade_cut == grp.clade_cut.min()]


sweep_data = df.groupby(['start']).apply(get_best_sweep_call).reset_index(drop=True)

sweep_data.to_hdf(str(args.output_file_path), 'df', format='table', mode='w')
