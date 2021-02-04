
import sys
import re
import os
import pandas as pd
import h5py

_, output_file_name, steps_dir, *clues_file_names = sys.argv # pylint: disable=unbalanced-tuple-unpacking

# open output file:
output_file = open(output_file_name, 'w')

# loop over base names of 
for clues_file_name in clues_file_names:
#    98000000_99500000_1.bed_98614821.h5
    start, end, chain, pos = re.search(r'(\d+)_(\d+)_(\d+).bed_(\d+).h5', clues_file_name).groups()
    h5_path = os.path.join(steps_dir, clues_file_name)
    if os.path.getsize(h5_path) == 0:
        log_likelihood_ratio = 'NA'
        selection_coef = 'NA'
    else:
        h5 = h5py.File(h5_path, 'r')
        log_likelihood_ratio = h5['logLikelihoodRatios'][h5.attrs['iHat'], h5.attrs['jHat']]
        selection_coef = h5.attrs['sHat']

    print(start, end, pos, chain, log_likelihood_ratio, selection_coef, sep='\t', file=output_file)
