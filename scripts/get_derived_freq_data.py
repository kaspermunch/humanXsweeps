

import pandas as pd
import sys
import argparse

description = """
Extract derived frequency information of one or many SNPs
"""

parser = argparse.ArgumentParser(description=description)


parser.add_argument("--start",
                  type=int,
                  help="Start of a genomic window")
parser.add_argument("--end",
                  type=int,
                  help="End of a genomic window")
parser.add_argument("--nrsnps",
                  type=int,
                  help="Number of SNPs to pick in window")
parser.add_argument("--maxsnps",
                  type=int,
                  help="Pick at most this number of SNPs to pick in window")
parser.add_argument("--minfreq",
                  type=float,
                  default=0,
                  help="Minimum frequency of selected SNPs")
                  
parser.add_argument("--snppos",
                  type=int,
                  help="Position of a single SNP")

parser.add_argument("derived_info_file",
                  type=str,
                  help="HDF file with derived variant info.")   
parser.add_argument("chrom",
                  type=str,
                  help="chromosome name (without 'chr' in front)")
# parser.add_argument("pop",
#                   type=str,
#                   help="Abreviation for 1000 genomes population.")
parser.add_argument("output_file_name",
                  type=str,
                  help="Output file")

args = parser.parse_args()

if args.snppos and args.minfreq:
    print("Do not use --snppos with --minfreq")
    sys.exit()
if args.nrsnps and args.maxsnps:
    print("Do not use --nrsnps with other --maxsnps")
    sys.exit()
if args.snppos and any([args.start, args.end, args.nrsnps, args.maxsnps]):
    print("Do not use --snppos with other arguments")
    sys.exit()
if bool(args.start) != bool(args.end):
    print("Specify either both start end end or none of them.")
    sys.exit()    

if args.snppos:
    query = f'pos=={args.snppos}'
else:
    query = f'(pos>={args.start})&(pos<{args.end})&(derived_freq>{args.minfreq})'

df = pd.read_hdf(args.derived_info_file, key='df', where=[query])
df.insert(0, 'chrom', args.chrom)

# no missing data allowed
df = df.loc[df.prop_missing == 0]

if args.maxsnps:
    all_snps = len(df)
    step = all_snps // args.maxsnps
    if step > 0:
        df = df.iloc[0:step*args.maxsnps:step]
    assert len(df) <= args.maxsnps
elif args.nrsnps:
    all_snps = len(df)
    step = all_snps // args.nrsnps
    assert step > 0
    df = df.iloc[0:step*args.nrsnps:step]
    assert len(df) == args.nrsnps

df.to_csv(args.output_file_name, index=False, header=False, sep='\t')


