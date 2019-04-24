
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


with open(str(args.masked_ref)) as f:
    a_haplo = list(next(SeqIO.parse(f, "fasta")).seq)
    b_haplo = a_haplo[:]


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19010

haploid = args.haploid

prev_chrom = None
for line in sys.stdin:
    if line.startswith('#'):
       continue
    chrom, pos, snpid, ref, alt, qual, filt, info, form, genotype = line.split()
    assert chrom == prev_chrom or prev_chrom is None

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
        if args.out2:
            if b_call == '1':
                b_haplo[int(pos)-1] = alt

out1 = open(str(args.out1), 'w')
out1.write(">{}\n{}\n".format(chrom, ''.join(a_haplo)))
if args.out2:
    out2 = open(str(args.out2), 'w')
    out2.write(">{}\n{}\n".format(chrom, ''.join(b_haplo)))


