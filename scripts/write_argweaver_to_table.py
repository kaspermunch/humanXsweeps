
import re
from pathlib import Path
import pandas as pd
import argparse

def write_sample_stats_to_table(stats_files, store_path):
    """
    Reads the sample stats for one sample produced by arg-summarize,
    adds a few of my own tree statistics and writes the whole thing
    to a HDF5 file for quick reading.
    """
    # coord_rgx = re.compile(r'([^-]+)-(\d+)-(\d+)_samples')
    # chain_smpl_rgx = re.compile(r'samples_(\d+)\.(\d+)')
    
    # chrom, start, end = coord_rgx.search(analysis_dir.name).groups()
    # start, end = int(start), int(end)

    # df_list = list()
    # for stats_file in analysis_dir.glob('*.bed.stats'):

    regex = re.compile(r'([^-]+)-(\d+)-(\d+)_samples_(\d+)\.(\d+)')

    printed_header = False

    df_list = list()
    for stats_file in stats_files:
        chrom, start, end, chain, sample = regex.search(stats_file.name).groups()
        start, end, chain, sample = int(start), int(end), int(chain), int(sample)

        df = pd.read_table(stats_file, skiprows=[0, 1])
        df = df.rename(columns={'#chrom': 'chrom', 'chromStart': 'start', 'chromEnd': 'end'})
        df['chain'] = chain
        df['chrom'] = chrom
        df['sample'] = sample
        df['start'] += start
        df['end'] += start
       
        if not printed_header:
            header = list(df)
            print(*header, sep='\t')
            printed_header = True

        for row in df.itertuples():
            print(*row, sep='\t')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path)
    parser.add_argument("--sample-dir", dest='sample_dir', type=Path)
    parser.add_argument("--glob-pattern", dest='glob_pattern', type=str)
    args = parser.parse_args()

    stats_files = args.sample_dir.glob(args.glob_pattern)

    write_sample_stats_to_table(stats_files, args.table)





