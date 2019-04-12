
import sys
import os
from pathlib import Path
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--sample", type=str)
parser.add_argument("--out1", type=Path)
parser.add_argument("--masked_ref", type=Path)
group = parser.add_mutually_exclusive_group()
group.add_argument("--haploid", action='store_true')
group.add_argument("--out2", type=Path)
args = parser.parse_args()

out1 = open(str(args.out1), 'w')
if args.out2:
    out2 = open(str(args.out2), 'w')
else:
    out2 = None

with open(str(args.masked_ref)) as f:
    a_haplo = list(next(SeqIO.parse(f, "fasta")).seq)
    b_haplo = a_haplo[:]


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19010

haploid = args.haploid

for line in sys.stdin:
    if line.startswith('#'):
       continue
    chrom, pos, snpid, ref, alt, qual, filt, info, form, genotype = line.split()

    if haploid:
        if genotype == '1':
            a_haplo[int(pos)-1] = alt
        elif genotype != '0':   # not haploid
            a_haplo[int(pos)-1] = 'N'
    else:
        assert haploid or '|' in genotype

        a_call, b_call = genotype.split('|')
        if a_call == '1':
            a_haplo[int(pos)-1] = alt
        if out2:
            if b_call == '1':
                b_haplo[int(pos)-1] = alt

out1.write(">{}\n{}\n".format(args.sample, ''.join(a_haplo)))
if out2:
    out2.write(">{}\n{}\n".format(args.sample, ''.join(b_haplo)))


