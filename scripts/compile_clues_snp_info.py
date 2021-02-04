
import pandas as pd
import numpy as np
import re
import sys
import os
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict
from genome_window_iter import genome_window_iter

import simons_meta_data

############################################################

target_region = 'WestEurasia'

meta_data_dir = Path(os.environ['HOME'], 'simons/faststorage/data/metadata')
individuals, populations, regions = simons_meta_data.get_meta_data(meta_data_dir=meta_data_dir)
included_individuals = []
for region, indivs in regions.items():
    if region  == target_region:
        for indiv in indivs:
            if individuals[indiv]['Genetic sex assignment'] == 'XY':
                included_individuals.append(indiv)



# rename target region
target_region = 'EuropeSubset'
# only include this subset of individuals
subset_of_europeans = ['B_Crete-2', 'B_French-3', 'B_Sardinian-3', 'S_Basque-1', 'S_Bulgarian-1', 'S_Bulgarian-2', 
                       'S_Czech-2', 'S_English-1', 'S_Estonian-1', 'S_Estonian-2', 'S_Finnish-3', 'S_French-1', 
                       'S_Greek-1', 'S_Greek-2', 'S_Hungarian-2', 'S_Polish-1', 'S_Saami-2', 'S_Sardinian-1', 
                       'S_Spanish-1', 'S_Tuscan-2']
included_individuals = [i for i in included_individuals if i in subset_of_europeans]


print(target_region)
print('\n'.join(included_individuals))
print(len(included_individuals))

fasta_files = ['steps/male_x_haploids/{}-A.fa'.format(x) for x in included_individuals]
chimp_file = '../../data/cteam_lite_public3/FullyPublic/Chimp.fa'

hdf_file_name = f'steps/clues/derived_freq_info_{target_region}.h5'
############################################################

# get chimp seq
for record in SeqIO.parse(chimp_file, "fasta"):
    if record.id == 'X':
        chimp_seq = str(record.seq)
        break

names = [re.search('([^/]+)-A.fa$', f).group(1) for f in fasta_files]

four_bases = ['A', 'T', 'G', 'C']

df_list = list()
for window in genome_window_iter(*fasta_files, window_size=1000000):

    chroms, starts, ends, seqs = list(zip(*window))
    assert chroms[1:] == chroms[:-1] and starts[1:] == starts[:-1] and ends[1:] == ends[:-1]
    chrom, start, end = names[0], starts[0], ends[0]

    # create a dataframe with all the bases
    df = pd.DataFrame()
    for name, seq in zip(names, seqs):
        df[name] = list(seq)
    
    # offset positions window by start so they become genomic coordinates
    df.index += start
    
    # remove monomorphic sites (we do this first becuase it is fast)
    df = df.loc[df.ne(df.iloc[:, 0], axis=0).any(axis=1)]

    # keep only bialleleic SNPs (may seem like an odd way, but it is super fast)
    df = df.loc[sum([df.eq('G').any(axis=1), df.eq('C').any(axis=1), df.eq('T').any(axis=1), df.eq('A').any(axis=1)]) == 2]

    # get counts of each base in each row
    base_counts = df.apply(pd.value_counts, axis=1)
    base_counts.replace(np.nan, 0, inplace=True)

    records = list()
    for idx in base_counts.index.values:

        row = base_counts.loc[idx]
        bases = row[row > 0].index.values
        bases = [b for b in bases if b in four_bases]

        ancestral_base = chimp_seq[idx].upper()

        if len(bases) != 2 or ancestral_base not in bases:
            continue

        if bases[0] == ancestral_base:
            ancestral_base, derived_base = bases
        else:
            derived_base, ancestral_base = bases

        derived_count = row.loc[derived_base]
        ancestral_count = row.loc[ancestral_base]
        n_count = row.loc['N']

        assert n_count + derived_count + ancestral_count == len(names)

        derived_freq = derived_count / (derived_count + ancestral_count)

        records.append((idx+1, derived_base, derived_freq, n_count/len(names))) # add one to idx to make positions one based

    derived_info = pd.DataFrame.from_records(records, columns=['pos', 'derived_base', 'derived_freq', 'prop_missing'])

    print(derived_info, len(derived_info))

    if len(derived_info):
        df_list.append(derived_info)


pd.concat(df_list).to_hdf(hdf_file_name, key='df', mode='w', format='table', 
    data_columns=['pos', 'derived_freq', 'prop_missing'], complevel=9, complib='blosc')
    