
import sys
import pandas as pd
import numpy as np
from chromwindow import window

_, rec_map_file, chrom, start, end, output_file = sys.argv
start, end = int(start), int(end)

df = pd.read_csv(rec_map_file, sep='\t', names=['chrom', 'start', 'end', 'rate', 'maplen'])

df = df.loc[(df.chrom == chrom) & (df.end > start) & (df.start < end)]
df.loc[df.index[0], 'start'] = start
df.loc[df.index[-1], 'end'] = end

# offset to zero before dowing windows
df['start'] -= start
df['end'] -= start

window_size = 10000

@window(size=window_size)
def mean_rate(df):
    "Compute mean rate in window but return nan if coverage is below half"
    if (df.end - df.start).sum() < window_size / 2:
        return np.nan
    return np.average(df.rate, weights=df.end - df.start)

means_df = mean_rate(df)

print(means_df.loc[means_df.mean_rate.isnull()])

nr_nans = means_df.mean_rate.isnull().sum()
assert nr_nans / len(means_df) <= 0.2, nr_nans

# 
means_df.loc[means_df.mean_rate.isnull()] = means_df.mean_rate.mean()

# offset back
means_df['start'] += start
means_df['end'] += start

means_df['start'] = means_df['start'].astype('int32')
means_df['end'] = means_df['end'].astype('int32')

means_df.to_csv(output_file, sep='\t', index=False, header=False, na_rep='nan')

# WHAT DO DO WITH REGIONS THAT DOES MAP TO X IN HG38