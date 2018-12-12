
import sys, os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

ampl_regions_file, input_file_name, output_file_name = sys.argv[1:]

# read single fasta record
record = SeqIO.read(input_file_name, "fasta")

# get sequence as list
seq = list(str(record.seq))

orig_len = len(seq)

# loop ampl regions and mask as X
with open(ampl_regions_file) as f:

    skip_first = next(f)

    for line in f:
        lst = line.split('\t')
        start, end = int(lst[6]), int(lst[7])
        seq[start:end] = 'N' * (end - start)

assert len(seq) == orig_len

# add seq back to fasta record
record = SeqRecord(Seq(''.join(seq)), id='X', description='')
#record.seq.seq = ''.join(seq)

with open(output_file_name, "w") as f:
    SeqIO.write([record], f, "fasta")
