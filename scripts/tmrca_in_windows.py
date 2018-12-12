

import sys
import pandas
_, input_table_file, output_hdf_file = sys.argv

def tmrca_means(df):
    segment_lengths = df.end - df.start
    return pandas.DataFrame({'start': [df.start.min()],
    						 'end':	 [df.end.max()],
    	                     'tmrca_half': [(df.tmrca_half * segment_lengths).sum() / segment_lengths.sum()],
                             'tmrca': [(df.tmrca * segment_lengths).sum() / segment_lengths.sum()]})

table_df = pandas.read_table(input_table_file, compression='gzip')
window_df = table_df.groupby(['chain', 'MCMC_sample']).apply(tmrca_means).reset_index(level=['chain', 'MCMC_sample'])
window_df.to_hdf(output_hdf_file, 'df', mode='w')
