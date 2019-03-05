import sys, gzip
from Bio import SeqIO
from Bio.Seq import Seq

_, template_file_name, input_file_name, pad_char, output_file_name = sys.argv

template_file = gzip.open(template_file_name, 'rt')
input_file = gzip.open(input_file_name, 'rt')
output_file = gzip.open(output_file_name, 'wt')

template_records = SeqIO.parse(template_file, "fasta")
for fasta_record in SeqIO.parse(input_file, "fasta"):
    template_record = next(template_records)
    if fasta_record.id == 'MT':
        continue
    assert fasta_record.id == template_record.id
    diff = len(template_record) - len(fasta_record)
    assert diff >= 0
    fasta_record.seq += Seq(pad_char * diff)
    SeqIO.write(fasta_record, output_file, 'fasta')