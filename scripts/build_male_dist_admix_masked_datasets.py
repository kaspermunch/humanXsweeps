
from pathlib import Path
import argparse
import pickle
import pandas
from pandas import DataFrame, Series

import simons_meta_data
from hg19_chrom_sizes import hg19_chrom_sizes as chromosome_lengths

parser = argparse.ArgumentParser()
parser.add_argument("--dist-dir", dest="dist_dir", type=Path)
parser.add_argument("--meta-data-dir", dest="meta_data_dir", type=Path)
parser.add_argument("--out-file", dest="out_file", type=Path)
# parser.add_argument("--result-dir", dest="result_dir", type=Path)
# parser.add_argument("--result-file-prefix", dest="result_file_prefix", type=str, default='dist_data')
args = parser.parse_args()

# easy loading of meta data in a consistent manner across code
individuals, populations, regions = simons_meta_data.get_meta_data(meta_data_dir=args.meta_data_dir)

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

dist_data['indiv_1'], dist_data['indiv_2'] = zip(*swap_indiv_to_put_non_african_first(dist_data['indiv_1'], dist_data['indiv_2']))


# Add categorical annotation of regions and populations:

# add region name
dist_data['region_label_1'] = [individuals[x]['Region'] for x in dist_data.indiv_1]
dist_data['region_label_2'] = [individuals[x]['Region'] for x in dist_data.indiv_2]

# map regions to integers
region_id = {'Africa': 0, 'WestEurasia': 1, 'SouthAsia': 2,
             'CentralAsiaSiberia': 3, 'Oceania': 4, 'EastAsia': 5, 'America':6}
dist_data['region_id_1'] = [region_id[x] for x in dist_data.region_label_1]
dist_data['region_id_2'] = [region_id[x] for x in dist_data.region_label_2]

# ordered list of region categories
region_categories = ['Africa', 'WestEurasia', 'SouthAsia', 'CentralAsiaSiberia',
                     'Oceania', 'EastAsia', 'America']
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

dist_data.to_hdf(str(args.out_file), 'df',  mode='w', format="table", data_columns=['indiv_1', 'start', 'end']) # we index all data columns 








