
from pathlib import Path
import argparse
import pickle
import pandas
from pandas import DataFrame, Series

import simons_meta_data
from hg19_chrom_sizes import hg19_chrom_sizes as chromosome_lengths

parser = argparse.ArgumentParser()
parser.add_argument("--pi-dir", dest="pi_dir", type=Path)
parser.add_argument("--meta-data-dir", dest="meta_data_dir", type=Path)
parser.add_argument("--result-dir", dest="result_dir", type=Path)
args = parser.parse_args()

# easy loading of meta data in a consistent manner across code
individuals, populations, regions = simons_meta_data.get_meta_data(meta_data_dir=args.meta_data_dir)

def read_pi_table(file_name):
    col_names = ['chrom', 'start', 'end', 'pop_label',
             'indiv_1', 'pseud_1', 'indiv_2', 'pseud_2',
             'pi', 'mismatch', 'match', 'uncalled']
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
    
pi_data = (pandas.concat(map(read_pi_table, args.pi_dir.iterdir())) # read and concat pi tables
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

# Add categorical annotation of regions and populations:

# add region name
pi_data['region_label'] = [individuals[x]['Region'] for x in pi_data.indiv_1]

# map regions to integers
region_id = {'Africa': 0, 'WestEurasia': 1, 'SouthAsia': 2,
             'CentralAsiaSiberia': 3, 'Oceania': 4, 'EastAsia': 5, 'America':6}
pi_data['region_id'] = [region_id[x] for x in pi_data.region_label]

# ordered list of region categories
region_categories = ['Africa', 'WestEurasia', 'SouthAsia', 'CentralAsiaSiberia',
                     'Oceania', 'EastAsia', 'America']
pi_data['region'] = pandas.Categorical(pi_data.region_label, categories=region_categories, ordered=True)

# ordered list of population categories
pop_categories = (pi_data[['region', 'pop_label', 'longitude_1']]
                  .sort_values(by=['region', 'pop_label', 'longitude_1'])
                  .loc[:, 'pop_label']
                  .drop_duplicates()
                 )                  
pi_data['population'] = pandas.Categorical(pi_data.pop_label, categories=pop_categories, ordered=True)

# sort
pi_data.sort_values(by=['chrom', 'region_id', 'pop_label', 'indiv_1', 'start'], inplace=True)


# Store population and region categories for later use:
store_path = args.result_dir / 'population_categories.store'
pop_categories.to_hdf(str(store_path), 'sr',  mode='w', format="table")

store_path = args.result_dir / 'region_categories.store'
pandas.Series(region_categories).to_hdf(str(store_path), 'sr',  mode='w', format="table")


# write stores for each chromosome
groups = pi_data.groupby('chrom')
for chrom, group in groups:
    store_path = args.result_dir / 'pi_data_chr{}_100kb.store'.format(chrom)
    group.to_hdf(str(store_path), 'df',  mode='w', format="table")


