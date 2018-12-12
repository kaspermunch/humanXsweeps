
import sys, os, gzip, pickle, argparse
from pathlib import Path
from random import seed
seed(42)


from genome_window_iter import genome_window_iter

parser = argparse.ArgumentParser()
parser.add_argument("--skip", dest="skip", type=str, action='append', default=[])
parser.add_argument("--mask-level", dest='mask_level', type=int)
parser.add_argument("unmasked_file_name", type=Path)
parser.add_argument("mask_file_name", type=Path)
parser.add_argument("output_file_name", type=Path)
args = parser.parse_args()

mask_level = args.mask_level

outfile = gzip.open(str(args.output_file_name), 'wb')

prev_chrom = None

for window in genome_window_iter(str(args.unmasked_file_name), str(args.mask_file_name), window_size=1000000, skip=args.skip):

    names, starts, ends, seqs = list(zip(*window))

    assert names[1:] == names[:-1], names
    assert starts[1:] == starts[:-1], starts
    assert ends[1:] == ends[:-1], ends

    chrom, start, end = names[0], starts[0], ends[0]

    seq, mask = seqs

    if not prev_chrom or chrom != prev_chrom:
        header_template = prev_chrom and '\n>{}\n' or '>{}\n'
        outfile.write(header_template.format(chrom).encode()) 

    # N and '-' chars are not masked
    mask_bool = (x.isdigit() and int(x) >= mask_level for x in mask)

    masked_seq = ''.join(m and x or 'N' for x, m in zip(seq, mask_bool))

    outfile.write(''.join(masked_seq).encode())

    prev_chrom = chrom
            
outfile.write(b'\n')



# from genome_window_iter import genome_window_iter

# mask_level, unmasked_file_name, mask_file_name, output_file_name = sys.argv[1:]
# mask_level = int(mask_level)

# outfile = gzip.open(output_file_name, 'wb')

# prev_chrom = None

# for window in genome_window_iter(unmasked_file_name, mask_file_name, window_size=1000000):

#     names, starts, ends, seqs = list(zip(*window))

#     assert names[1:] == names[:-1]
#     assert starts[1:] == starts[:-1]
#     assert ends[1:] == ends[:-1]

#     chrom, start, end = names[0], starts[0], ends[0]
#     seq, mask = seqs

#     if not prev_chrom or chrom != prev_chrom:
#         header_template = prev_chrom and '\n>{}\n' or '>{}\n'
#         outfile.write(header_template.format(chrom).encode()) 

#     # N and '-' chars are not masked
#     mask_bool = (x.isdigit() and int(x) >= mask_level for x in mask)

#     masked_seq = ''.join(m and x or 'N' for x, m in zip(seq, mask_bool))

#     outfile.write(''.join(masked_seq).encode())

#     prev_chrom = chrom
            
# outfile.write(b'\n')

