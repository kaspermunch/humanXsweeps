

import numpy as numpy
import pandas as pd

fixed_som_table_name = '/home/kmt/simons/faststorage/data/metadata/nature18964-s2-fixed-genders.csv'
elise_table_name = '/home/kmt/simons/faststorage/data/metadata/Simons_meta_accessionnb_update2.txt'
merged_file_name = '/home/kmt/simons/faststorage/data/metadata/meta_data.tsv'

fixed_som_table = pd.read_csv(fixed_som_table_name, sep=';')
fixed_som_table = fixed_som_table.set_index('Sample ID (SGDP)')

elise_table = pd.read_csv(elise_table_name, sep='\t')
sub_table = elise_table[['ID','ENA-RUN', 'Illumina_ID', 'Sample_ID', 'Sample_ID(Aliases)', 'CTeam_ID']]
sub_table = sub_table.rename(columns = {'CTeam_ID': 'Sample ID (SGDP)'})
sub_table.drop_duplicates(subset='Sample ID (SGDP)', keep='last', inplace=True)
sub_table = sub_table.set_index('Sample ID (SGDP)')

result = pd.concat([fixed_som_table, sub_table], axis=1, verify_integrity=True)
result.reset_index().to_csv(merged_file_name, sep='\t', index=False)
