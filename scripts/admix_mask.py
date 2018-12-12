
import sys, os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

admix_pred_file, min_post_prob, input_file_name, output_file_name = sys.argv[1:]

min_post_prob = float(min_post_prob)

# get name of individual from file fasta name
indiv = os.path.splitext(os.path.basename(input_file_name))[0].replace('-A', '')

# read single fasta record
record = SeqIO.read(input_file_name, "fasta")

# get sequence as list
seq = list(str(record.seq))

tot_masked = 0
orig_len = len(seq)

# loop admix segments and mask as X
with open(admix_pred_file) as f:

	skip_first = next(f)

	for line in f:
		lst = line.split('\t')
		name = lst[0]
		start = int(lst[2]) - 1 # laurits snps are one-based 
		end = int(lst[3]) - 1
		postprob = float(lst[10])

		# filter for individual and to only use 0.8 postprob
		if indiv == name and postprob >= min_post_prob:
			seq[start:end] = [x.lower() for x in seq[start:end]]

			tot_masked += end-start

assert len(seq) == orig_len
assert tot_masked == sum([int(x.islower()) for x in seq]), (tot_masked, sum([int(x.islower()) for x in seq]))

# add seq back to fasta record
record = SeqRecord(Seq(''.join(seq)), id='X', description='')
#record.seq.seq = ''.join(seq)

with open(output_file_name, "w") as f:
    SeqIO.write([record], f, "fasta")
