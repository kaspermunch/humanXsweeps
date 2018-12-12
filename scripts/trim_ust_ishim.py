
import sys, gzip

_, seq_input_file, mask_input_file, seq_output_file, mask_output_file = sys.argv

from Bio import SeqIO
from Bio.Seq import Seq

from hg19_chrom_sizes import hg19_chrom_sizes

seqin = gzip.open(seq_input_file, "rt")
maskin = gzip.open(mask_input_file, "rt")
seqout = gzip.open(seq_output_file, "wt")
maskout = gzip.open(mask_output_file, "wt")

seq_iter = SeqIO.parse(seqin, "fasta")
mask_iter = SeqIO.parse(maskin, "fasta")


seq_records = list()
mask_records = list()

for seq_record in seq_iter:
	mask_record = next(mask_iter)

	min_len = min(len(seq_record.seq), len(mask_record.seq))

	missing_len = hg19_chrom_sizes['chr' + mask_record.id] - min_len

	s = str(seq_record.seq)
	seq_record.seq = Seq(s[:min_len] + 'N' * missing_len)

	s = str(mask_record.seq)
	mask_record.seq = Seq(s[:min_len] + 'N' * missing_len)

	print(len(seq_record.seq), hg19_chrom_sizes['chr' + mask_record.id])

	seq_records.append(seq_record)
	mask_records.append(mask_record)
#	break
SeqIO.write(seq_records, seqout, "fasta")
SeqIO.write(mask_records, maskout, "fasta")
