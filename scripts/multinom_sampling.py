


import os
import sys
import numpy as np
import pandas as pd
import itertools
from collections import defaultdict
from math import exp
from multiprocessing import Pool

script, output_file_name = sys.argv

np.random.seed(7)

# african X/A ratio is 0.65 but is further reduced inside regions (lower african X pi in regions):
nonafr_x_auto_ratios = [0.65 * x for x in [1, 0.71]] # outside and inside regions

# mean per generation recombination rate in regions (new decode map):
sexavg_rec_rates_per_gen = [0.46e-8, # mean in regions
                            1.16e-8] # global for chrX
autosomal_bottle_N = 3000

# parameters for computation
haploid_N = [int(2 * autosomal_bottle_N * x) for x in nonafr_x_auto_ratios]
target_freqs = np.linspace(0.1, 1, 10)
g = [int(10000 / 29), int(15000 / 29)]
r = sexavg_rec_rates_per_gen
L = [5e5, 1e6]
n_samples = 10000
parameters = list(itertools.product(haploid_N, target_freqs, g, r, L))

def get_gens_to_freq(N, target_freqs):
    trajectory = []
    counts = np.array(np.ones(N))
    trajectory.append(counts)
    while np.count_nonzero(counts) > 1:
        counts = np.random.multinomial(counts.size, counts/counts.size)
        trajectory.append(counts)
    trajectory = np.array(trajectory)

    results = []
    for f in target_freqs:
        max_freqs = np.amax(trajectory, axis=1) / N
        nr_gens = np.argmax(max_freqs > f)
        if nr_gens != 0: # freq reached
            results.append((N, f, nr_gens))
    return results

# compute nr of gens until one allele reaches frequency f in a population of size N
gens_freq_reached = defaultdict(list)
for N in haploid_N:
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        with Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as p:
            results = list(p.starmap(get_gens_to_freq, [(N, target_freqs)]*n_samples))
    else:
        results = list(itertools.starmap(get_gens_to_freq, [(N, target_freqs)]*n_samples))

    for batch in results:
        for N, f, nr_gens in batch:
            gens_freq_reached[(N, f)].append(nr_gens) 

records = []
for N, f, g, r, L in parameters:
    total_prob = 0
    total_freqs_reached = 0
    for nr_gens in gens_freq_reached[(N, f)]:
        if nr_gens < g:
            total_freqs_reached += 1
            total_prob += exp(-r*L*nr_gens) + (1-exp(-r*L*nr_gens))/2 

    p = total_prob / n_samples
    records.append((N, f, g, r, L, p, total_freqs_reached, n_samples))

df = pd.DataFrame.from_records(records, columns=['haploid_N', 'freq', 'max_gens', 'rec', 'length', 'prob', 'tot_reached', 'n_samples'])

df.to_hdf(output_file_name, 'df', format='table', mode='w')
