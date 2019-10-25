import sys, os
import gc
from pathlib import Path
import argparse
import pickle
import pandas
from pandas import DataFrame, Series
import numpy as np

import simons_meta_data
from hg19_chrom_sizes import hg19_chrom_sizes as chromosome_lengths

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_dir + '/../notebooks')
import analysis_globals

parser = argparse.ArgumentParser()
parser.add_argument("--dist-dir", dest="dist_dir", type=Path)
parser.add_argument("--meta-data-dir", dest="meta_data_dir", type=Path)
parser.add_argument("--out-file", dest="out_file", type=Path)
parser.add_argument("--dist-twice-out-file", dest="dist_twice_out_file", type=Path)
parser.add_argument("--include-ust-ishim", dest="include_ust_ishim", action='store_true', default=False)
# parser.add_argument("--result-dir", dest="result_dir", type=Path)
# parser.add_argument("--result-file-prefix", dest="result_file_prefix", type=str, default='dist_data')
args = parser.parse_args()

# easy loading of meta data in a consistent manner across code
individuals, populations, regions = simons_meta_data.get_meta_data(
    meta_data_dir=args.meta_data_dir,
    include_ust_ishim=args.include_ust_ishim)


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


def read_dist_table(file_name):

    col_names = ['chrom', 'start', 'end', 'pop_label',
             'indiv_1', 'pseud_1', 'indiv_2', 'pseud_2',
             'dist', 'mismatch', 'match', 
             'dist_af', 'mismatch_af', 'match_af', 
             'uncalled']
    with open(str(file_name), 'rb') as f:
        table = pickle.load(f)
    df = DataFrame(table, columns=col_names)
    df['indiv1'] = [Path(x).name for x in df.indiv_1]
    df['indiv2'] = [Path(x).name for x in df.indiv_2]
    return df

def indiv_filter(df):
    indiv1_in_meta = [x in individuals for x in df.indiv_1]
    indiv2_in_meta = [x in individuals for x in df.indiv_2]
    return([x and y for x, y in zip(indiv1_in_meta, indiv2_in_meta)])
    
dist_data = (pandas.concat(map(read_dist_table, args.dist_dir.iterdir())) # read and concat pi tables
           # filter individuals
           .loc[indiv_filter]
           # add meta data
           .assign(longitude_1 = lambda df: [individuals[x]['Longitude'] for x in df.indiv_1],
               latitude_1 = lambda df: [individuals[x]['Latitude'] for x in df.indiv_1],
               sex_1 = lambda df: [individuals[x]['Genetic sex assignment'] for x in df.indiv_1],
               #country_1 = lambda df: [individuals[x]['Country'] for x in df.indiv_1],
               #y_haplogroup_1 = lambda df: [individuals[x]['Y haplogroup'] for x in df.indiv_1],
               #coverage_1 = lambda df: [individuals[x]['Coverage (mean)'] for x in df.indiv_1],
               longitude_2 = lambda df: [individuals[x]['Longitude'] for x in df.indiv_2],
               latitude_2 = lambda df: [individuals[x]['Latitude'] for x in df.indiv_2],
               sex_2 = lambda df: [individuals[x]['Genetic sex assignment'] for x in df.indiv_2],
               #country_2 = lambda df: [individuals[x]['Country'] for x in df.indiv_2],
               #y_haplogroup_2 = lambda df: [individuals[x]['Y haplogroup'] for x in df.indiv_2],
               #coverage_2 = lambda df: [individuals[x]['Coverage (mean)'] for x in df.indiv_2]
          )
   .reset_index(drop=True)
  )


def swap_indiv_to_put_non_african_first(s1, s2):
    l = list()
    for indiv1, indiv2 in zip(s1, s2):
        if indiv1 in regions['Africa']:
            l.append((indiv2, indiv1))
        else:
            l.append((indiv1, indiv2))
    return l


# correct for branch shortening of Ust Ishim individual:
if args.include_ust_ishim:
    correction = analysis_globals.mut_per_year * 45000
    dist_data.loc[dist_data.indiv_1 == 'Ust_Ishim', 'dist'] += correction
    dist_data.loc[dist_data.indiv_2 == 'Ust_Ishim', 'dist'] += correction


# swap indivs to put non-African first:
dist_data['indiv_1'], dist_data['indiv_2'] = zip(*swap_indiv_to_put_non_african_first(dist_data['indiv_1'], dist_data['indiv_2']))


# Add categorical annotation of regions and populations:

# add region name
dist_data['region_label_1'] = [individuals[x]['Region'] for x in dist_data.indiv_1]
dist_data['region_label_2'] = [individuals[x]['Region'] for x in dist_data.indiv_2]

# map regions to integers
region_id = {'Africa': 0, 'WestEurasia': 1, 'SouthAsia': 2,
             'CentralAsiaSiberia': 3, 'Oceania': 4, 'EastAsia': 5, 'America':6,
             'Ust_Ishim': 7}
if args.include_ust_ishim:
    region_id['Ust_Ishim'] = 7
dist_data['region_id_1'] = [region_id[x] for x in dist_data.region_label_1]
dist_data['region_id_2'] = [region_id[x] for x in dist_data.region_label_2]

# ordered list of region categories
region_categories = ['Africa', 'WestEurasia', 'SouthAsia', 'CentralAsiaSiberia',
                    'Oceania', 'EastAsia', 'America']
if args.include_ust_ishim:
    region_categories.append('Ust_Ishim')
dist_data['region_1'] = pandas.Categorical(dist_data.region_label_1, categories=region_categories, ordered=True)
dist_data['region_2'] = pandas.Categorical(dist_data.region_label_2, categories=region_categories, ordered=True)

## NB: in generating pickle files a nonsensical pop field was genereted in workflow.py
## It is denoted pop_label in col here. We just omit any population information here:

# sort
dist_data.sort_values(by=['chrom', 'region_id_1', 'indiv_1', 'start'], inplace=True)

# # write stores for each chromosome
# groups = dist_data.groupby('chrom')
# for chrom, group in groups:
# #    store_path = args.result_dir / 'dist_data_chr{}_100kb.store'.format(chrom)
#     store_path = args.result_dir / '{}_chr{}_100kb.store'.format(args.result_file_prefix, chrom)
#     group.to_hdf(str(store_path), 'df',  mode='w', format="table", data_columns=['indiv_1', 'start', 'end']) # we index all data columns 

dist_data = optimize_data_frame(dist_data, down_int='unsigned')
gc.collect()

dist_data.to_hdf(str(args.out_file), 'df',  mode='w', format="table", data_columns=['indiv_1', 'start', 'end']) # we index all data columns 





# def dist_twice(dist_data):

#     #we use the global individuals defined above...
#     #individuals, populations, regions = simons_meta_data.get_meta_data(meta_data_dir=analysis_globals.meta_data_dir)

#     if 'dist_af' in dist_data.columns:
#         # this is not a simulation

#         dist_data.drop('pop_label', axis=1, inplace=True)

#         dist_data['pop_1'] = [individuals[x]['Population ID'] for x in dist_data.indiv_1]
#         dist_data['pop_2'] = [individuals[x]['Population ID'] for x in dist_data.indiv_2]

#     # dict for swapping columns
#     swap_dict = dict()
#     for colname in dist_data.columns.values:
#         if colname.endswith('_1'):
#             swap_dict[colname] = colname[:-2] + '_2'
#         if colname.endswith('_2'):
#             swap_dict[colname] = colname[:-2] + '_1'

#     if 'dist_af' in dist_data.columns:
#         # this is not a simulation
#         cols = ['start', 'end', 'indiv_1', 'indiv_2', 
#                 'dist', 'mismatch', 'match', 
#                 'dist_af', 'mismatch_af', 'match_af',
#                 'uncalled', 'pop_1', 'pop_2',
#                 'region_label_1', 'region_label_2', 
#                 'region_id_1', 'region_id_2', 'region_1', 'region_2']
#     else:
#         cols =['start', 'end', 'indiv_1', 'indiv_2', 'dist']

#     dist_data_twice = (pandas.concat([dist_data[cols],
#                                       dist_data[cols].rename(columns=swap_dict)])
#         .sort_values(['indiv_1', 'start'])
#         .reset_index(drop=True)
#         )
    
#     # mask uncalled windows:
#     dist_data_twice.dist.where(dist_data_twice.uncalled <= analysis_globals.max_uncalled_bases, 
#                                inplace=True)
#     if 'dist_af' in dist_data.columns:
#         dist_data_twice.dist_af.where(dist_data_twice.uncalled <= analysis_globals.max_uncalled_bases, 
#                                    inplace=True)

#     return dist_data_twice


# ##### apply function and write hdf
# all_male_dist_twice = dist_twice(dist_data)


###################

if 'dist_af' in dist_data.columns:
    # this is not a simulation

    dist_data.drop('pop_label', axis=1, inplace=True)

    dist_data['pop_1'] = [individuals[x]['Population ID'] for x in dist_data.indiv_1]
    dist_data['pop_2'] = [individuals[x]['Population ID'] for x in dist_data.indiv_2]

# dict for swapping columns
swap_dict = dict()
for colname in dist_data.columns.values:
    if colname.endswith('_1'):
        swap_dict[colname] = colname[:-2] + '_2'
    if colname.endswith('_2'):
        swap_dict[colname] = colname[:-2] + '_1'

if 'dist_af' in dist_data.columns:
    # this is not a simulation
    cols = ['start', 'end', 'indiv_1', 'indiv_2', 
            'dist', 'mismatch', 'match', 
            'dist_af', 'mismatch_af', 'match_af',
            'uncalled', 'pop_1', 'pop_2',
            'region_label_1', 'region_label_2', 
            'region_id_1', 'region_id_2', 'region_1', 'region_2']
else:
    cols =['start', 'end', 'indiv_1', 'indiv_2', 'dist']

all_male_dist_twice = (pandas.concat([dist_data[cols],
                                  dist_data[cols].rename(columns=swap_dict)])
    .sort_values(['indiv_1', 'start'])
    .reset_index(drop=True)
    )


if args.include_ust_ishim:
    # compute an adjusted max_uncalled_bases for ust ishim that gives pairs includeing ust ishim
    # the same number of non-missing windows as the average pair not including ust ishim:
    ust_ishim_uncalled = all_male_dist_twice.loc[all_male_dist_twice.indiv_1 == 'Ust_Ishim'].uncalled
    other_indiv_uncalled = all_male_dist_twice.loc[all_male_dist_twice.indiv_1 != 'Ust_Ishim'].uncalled
    quant = (other_indiv_uncalled <= analysis_globals.max_uncalled_bases).sum() / other_indiv_uncalled.size
    ust_ishim_max_uncalled_bases = int(np.quantile(ust_ishim_uncalled, quant))

    # mask uncalled windows using the two different cutoffs:
    mask = (all_male_dist_twice.uncalled <= analysis_globals.max_uncalled_bases) | \
        ((all_male_dist_twice.indiv_1 == 'Ust_Ishim') | (all_male_dist_twice.indiv_2 == 'Ust_Ishim')) & \
        (all_male_dist_twice.uncalled <= ust_ishim_max_uncalled_bases)
    
else:
    mask = all_male_dist_twice.uncalled <= analysis_globals.max_uncalled_bases

# mask uncalled windows:
all_male_dist_twice.dist.where(mask, inplace=True)
if 'dist_af' in dist_data.columns:
    all_male_dist_twice.dist_af.where(mask, inplace=True)

######################

#all_male_dist_twice = optimize_data_frame(all_male_dist_twice, down_int='unsigned')
all_male_dist_twice = optimize_data_frame(all_male_dist_twice, down_int='unsigned')
gc.collect()

all_male_dist_twice.to_hdf(str(args.dist_twice_out_file), 'df', 
                           data_columns=['start', 'end', 
                                     'indiv_1', 'indiv_2'],
                           format='table', mode='w')
