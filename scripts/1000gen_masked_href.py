import sys
from Bio import SeqIO
import gzip

_, href_file, mask_file, out_file = sys.argv

with gzip.open(mask_file, 'rt') as f:
    record = next(SeqIO.parse(f, "fasta"))
    chrom = record.id
    mask = list(record.seq)

with open(href_file) as f:
    for record in SeqIO.parse(f, "fasta"):
        if record.id == chrom:
            href = list(record.seq)
            break

assert len(href) == len(mask), (len(href), len(mask))

for i in range(len(href)):
    if mask[i] != 'P':
        href[i] = 'N'

with open(out_file, 'w') as f:
    f.write(">masked_href\n{}\n".format(''.join(href)))



