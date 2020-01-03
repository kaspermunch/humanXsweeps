import sys, re
from Bio import SeqIO
import gzip
from pathlib import Path
from pandas import DataFrame
from collections import defaultdict


_, input_dir, out_file = sys.argv

bases = set(['A', 'T', 'C', 'G'])
records = list()
for path in Path(input_dir).glob('**/*'):

    if not path.suffix == '.fa':
        continue
    
    name, chrom, haplo = re.search(r'(\w+)_(\w+)-(\w+).fa', path.name).groups()
        
    with open(path) as f:
        record = next(SeqIO.parse(f, "fasta"))

    counts = defaultdict(int)
    for i, char in enumerate(record.seq):
        if char in bases:
            counts[i // 100000 * 100000] += 1

    for start in range(0, max(counts.keys())+100000, 100000):
        records.append((name, chrom, haplo, start, counts[start]))

df = DataFrame().from_records(records, columns=['name', 'chrom', 'haplo', 'start', 'nr_missing'])
df.to_hdf(out_file, 'df', format='table', mode='w')
