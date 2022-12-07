


import os
import sys
import numpy as np
import pandas as pd
import itertools
from collections import defaultdict
from math import exp
from multiprocessing import Pool, cpu_count
import argparse

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
# evenly spaced target freqs:

def compute(parameters, target_freqs, haploid_N, output_file_name):

    
    gens_freq_reached = defaultdict(list)
    for N in haploid_N:
        if 'SLURM_CPUS_PER_TASK' in os.environ:
            with Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as p:
                results = list(p.starmap(get_gens_to_freq, [(N, target_freqs)]*n_samples))
        # else:
        #     with Pool(cpu_count()) as p:
        #         results = list(p.starmap(get_gens_to_freq, [(N, target_freqs)]*n_samples))
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--years", dest="years", type=int, action='append')
    parser.add_argument('output_file_name_spaced_freqs', type=str, help='')
    parser.add_argument('output_file_name_ech_freqs', type=str, help='')
    args = parser.parse_args()

#    script, output_file_name_spaced_freqs, output_file_name_ech_freqs = sys.argv

    np.random.seed(7)

    # african X/A ratio is 0.65 but is further reduced inside regions (lower african X pi in regions):
    # nonafr_x_auto_ratios = [0.65 * x for x in [1, 0.71]] # outside and inside regions
    nonafr_x_auto_ratios = [0.65, 0.51] # African and non-African

    # mean per generation recombination rate in regions (new decode map):
    # sexavg_rec_rates_per_gen = [0.46e-8, # mean in regions
    #                             1.16e-8] # global for chrX
    # sexavg_rec_rates_per_gen = [1.16e-8] # global for chrX
    sexavg_rec_rates_per_gen = [0.77e-8] # sex-avaraged global for chrX (2/3*1.16e-8)
    
    autosomal_bottle_N = [3000, 1500]

    # parameters for computation
    haploid_N = [int(2*x*y) for (x, y) in itertools.product(autosomal_bottle_N, nonafr_x_auto_ratios)]
    #haploid_N = [int(2 * autosomal_bottle_N * x) for x in nonafr_x_auto_ratios]
    spaced_target_freqs = np.linspace(0.05, 1, 20)
    g = [int(x / 29) for x in args.years]
    # g = [int(10000 / 29), int(15000 / 29)]
    r = sexavg_rec_rates_per_gen
    # L = [5e5, 1e6]
    # n_samples = 1000000
    L = [5e5]
    n_samples = 500000

    ech_target_freqs =  pd.read_hdf('results/extended_peak_regions_5e-05_25%_90%.hdf').prop_swept


    spaced_parameters = list(itertools.product(haploid_N, spaced_target_freqs, g, r, L))
    ech_parameters = list(itertools.product(haploid_N, ech_target_freqs, g, r, L))

    compute(ech_parameters, ech_target_freqs, haploid_N, args.output_file_name_ech_freqs)
    compute(spaced_parameters, spaced_target_freqs, haploid_N, args.output_file_name_spaced_freqs)

