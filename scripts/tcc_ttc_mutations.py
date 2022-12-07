
import sys, os, gzip, pickle
import numpy as np
from random import seed
seed(42)


from genome_window_iter import genome_window_iter

f1, f2, binsize, chromosome, indiv1, pseud1, indiv2, pseud2, output_file_name = sys.argv[1:]
binsize = int(binsize)
result_list = list()

for window in genome_window_iter(f1, f2, window_size=binsize):

    names, starts, ends, seqs = list(zip(*window))

    assert names[1:] == names[:-1]
    assert starts[1:] == starts[:-1]
    assert ends[1:] == ends[:-1]

    chrom, start, end = names[0], starts[0], ends[0]

    if chrom != chromosome:
        continue

    # int repr of chars
    uncalled, match, mismatch = 0, 0, 0

    # count site classes
    
    for i, (a, b) in enumerate(zip(*seqs)):
        if a == 'N' or a == '-' or b == 'N' or b == '-':
            uncalled += 1
        else:
            mismatch += 1

#             check for cpg

            if a == 'C' and b == 'T':
                c_t_mutations += 1

            triplet, chimp_triplet = seqs[0][i:i+3], seqs[1][i:i+3]
        
            if chimp_triplet == 'TCC' and triplet == 'TTC':
                tcc_ttc_mutations += 1
            
            
            
    # add to result table
    result_list.append([chrom, start, start+binsize, 'global', 
                     indiv1, pseud1, indiv2, pseud2,
                     pi, mismatch, match, uncalled])


turn into data frame instead

# with open(output_file_name, 'wb') as f:
#     pickle.dump(result_list, f)

  
# python scripts/pi_windows.py steps/pseudohaploid_genomes/S_Even-1-A.fa.gz steps/pseudohaploid_genomes/S_Even-2-A.fa.gz 100000 Even S_Even-1 A S_Even-2 A tmp.out
