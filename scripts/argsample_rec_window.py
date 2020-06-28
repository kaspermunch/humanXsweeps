
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

nr_nans = means_df.mean_rate.isnull().sum()
assert nr_nans == 0, nr_nans

means_df.to_csv(output_file, sep='\t', index=False, header=False, na_rep='nan')

