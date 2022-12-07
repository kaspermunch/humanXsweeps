import sys
import pandas as pd
from chromwindow import window

_, decode_map_file, output_dir = sys.argv

win_size = 10000000
chrom_start = 4000000
chrom_end = 154000000
# read map
decode_sexavg = pd.read_table(decode_map_file, comment='#')
# extract non-par chrX (15 x 10mb) between 
decode_sexavg_chrX = (decode_sexavg
                      .loc[(decode_sexavg['Chr'] == 'chrX') & \
                           (decode_sexavg['Begin'] > chrom_start) & (decode_sexavg['End'] < chrom_end)]
                     )
decode_sexavg_chrX = decode_sexavg_chrX.rename(columns={'Chr': 'chrom', 'Begin':'start', 'End':'end'})

# shift coordinates so that start of first segmnet of map is 0
min_start = decode_sexavg_chrX.start.min()
decode_sexavg_chrX['start'] -= min_start
decode_sexavg_chrX['end'] -= min_start

# add bin identifiers
decode_sexavg_chrX['bin'] = decode_sexavg_chrX['start'] // win_size

# duplicate all rows that span transitions between bins:
idx = (decode_sexavg_chrX.start // win_size) < (decode_sexavg_chrX.end // win_size)
extra_rows = decode_sexavg_chrX.loc[idx].copy(deep=True)
decode_sexavg_chrX.loc[idx, 'end'] = (decode_sexavg_chrX.bin + 1) * win_size 
extra_rows.loc[:, 'start'] = (extra_rows.bin+1) * win_size 
extra_rows['bin'] += 1
decode_sexavg_chrX_with_split_rows = pd.concat([decode_sexavg_chrX, extra_rows]).sort_values(['start', 'end'])

# take case of start of first bin and end of last bin
decode_sexavg_chrX_with_split_rows.loc[lambda df: df.start == df.start.min(), 'start'] = 0
decode_sexavg_chrX_with_split_rows.loc[lambda df: df.end == df.end.max(), 'end'] = chrom_end - chrom_start

# make relative coordinates in each bin:
decode_sexavg_chrX_with_split_rows['rel_start'] = decode_sexavg_chrX_with_split_rows.start - decode_sexavg_chrX_with_split_rows.bin * win_size
decode_sexavg_chrX_with_split_rows['rel_end'] = decode_sexavg_chrX_with_split_rows.end - decode_sexavg_chrX_with_split_rows.bin * win_size

# SLiM needs a rate for when recombination can physically occur (i.e. in the female
# between the Xs). To get that from the sex averated recombination rate, we need to
# account for hte fact that only 2/3 of X chromosomes have the oportunity to combine 
# in each generation (assuming even sex ratios).
decode_sexavg_chrX_with_split_rows['cMperMb'] *= 3 / 2

# write to tsv files for reading from slim:
gr = decode_sexavg_chrX_with_split_rows.groupby('bin')
for name in gr.groups:
    df = gr.get_group(name).sort_values(['start', 'end'])
    df['rel_start'] += 1
    df[['rel_start', 'cMperMb']].to_csv(f'{output_dir}/{name}.tsv', sep="\t", header=False, index=False)
    # print(name)
    # print(df[['rel_start', 'rel_end', 'cMperMb']].head())
    # print(df[['rel_start', 'rel_end', 'cMperMb']].tail())