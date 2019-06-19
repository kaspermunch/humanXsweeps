
import sys, os, gzip, pickle
import numpy as np
from random import seed
seed(42)


from genome_window_iter import genome_window_iter

f1, f2, binsize, pop, indiv1, pseud1, indiv2, pseud2, output_file_name = sys.argv[1:]
binsize = int(binsize)
result_list = list()

for window in genome_window_iter(f1, f2, window_size=binsize):

    names, starts, ends, seqs = list(zip(*window))

    assert names[1:] == names[:-1], names
    assert starts[1:] == starts[:-1]
    assert ends[1:] == ends[:-1]

    chrom, start, end = names[0], starts[0], ends[0]

    # #########################################
    # # only difference from pi_windows.py:
    # if chrom != 'X':
    #     continue
    # #########################################

    # int repr of chars
    uncalled, match, mismatch = 0, 0, 0
    match_admix_masked, mismatch_admix_masked = 0, 0

    # count site classes
    for a, b in zip(*seqs):
        if a in 'Nn-' or b in 'Nn-': #a == 'N' or a == '-' or b == 'N' or b == '-':
            uncalled += 1
        elif a.upper() == b.upper():
            match += 1
            if not a.islower() and not b.islower():
                match_admix_masked += 1
        else:
            mismatch += 1
            if not a.islower() and not b.islower():
                mismatch_admix_masked += 1

    if mismatch + match:
        pi = mismatch / (mismatch + match)
    else:
        pi = np.nan

    if mismatch_admix_masked + match_admix_masked:
        pi_admix_masked = mismatch_admix_masked / (mismatch_admix_masked + match_admix_masked)
    else:
        pi_admix_masked = np.nan
            
    # add to result table
    result_list.append([chrom, start, start+binsize, pop, 
                     indiv1, pseud1, indiv2, pseud2,
                     pi, mismatch, match, pi_admix_masked, mismatch_admix_masked, match_admix_masked, uncalled])


with open(output_file_name, 'wb') as f:
    pickle.dump(result_list, f)

  
# python scripts/pi_windows.py steps/pseudohaploid_genomes/S_Even-1-A.fa.gz steps/pseudohaploid_genomes/S_Even-2-A.fa.gz 100000 Even S_Even-1 A S_Even-2 A tmp.out
