from multiprocessing import Pool
from math import exp
import os, sys
import itertools
import numpy as np
import pandas as pd

script, output_file_name = sys.argv
        
def freq_trajectory(N, n=1, s=0):
    yield n
    while N > n > 0:
        n = np.random.binomial(N, n*(1+s)/((N-n)+(n)*(1+s)), 1)[0]
        yield n   

def get_trajectory(N, f):
    allele_count = f * N
    trajectory = list(freq_trajectory(N, n=1, s=0))
    return trajectory

def prob_norec_freq_reached(N, f, g, r, L, n_samples=None):
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        # in parallel:
        with Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as p:
            trajectories = list(p.starmap(get_trajectory, [(N, f)]*n_samples))
    else:
        # single core:
        trajectories = list(itertools.starmap(get_trajectory, [(N, f)]*n_samples))
    tot_prob = 0
    tot_freq_reached = 0
    for trajectory in trajectories:
        max_freq = max(trajectory[:g])
        if max_freq >= allele_count:
            nr_gens = trajectory[:g].index(max_freq) + 1
            tot_freq_reached += 1
            tot_prob += exp(-r*L*nr_gens) + (1-exp(-r*L*nr_gens))/2  # prob of no recombination + prob of recombination onto the same kind of haplotype 
    return tot_freq_reached, tot_prob / n_samples


#########################

def freq_trajectory(N, n=1, s=0):
    yield n
    while N > n > 0:
        n = np.random.binomial(N, n*(1+s)/((N-n)+(n)*(1+s)), 1)[0]
        yield n            

def work(N, f, g, r, L):
    allele_count = f * N
    trajectory = list(freq_trajectory(N, n=1, s=0))
    max_freq = max(trajectory[:g])
    if max_freq >= allele_count:
        nr_gens = trajectory[:g].index(max_freq) + 1
        return exp(-r*L*nr_gens) + (1-exp(-r*L*nr_gens))/2  # prob of no recombination + prob of recombination onto the same kind of haplotype 
                                                            # (assuming a symemtric logistic frequency trajectory)
    return 0

def prob_norec_freq_reached(N, f, g, r, L, n_samples=None):
    if 'SLURM_CPUS_PER_TASK' in os.environ:
        # in parallel:
        with Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as p:
            results = list(p.starmap(work, [(N, f, g, r, L)]*n_samples))
    else:
        # single core:
        results = list(itertools.starmap(work, [(N, f, g, r, L)]*n_samples))
    return n_samples - results.count(0), sum(results) / n_samples

###################


np.random.seed(7)

# african X/A ratio is 0.65 but is further reduced inside regions (lower african X pi in regions):
nonafr_x_auto_ratios = [0.65 * x for x in [1, 0.71]] # outside and inside regions

# mean per generation recombination rate in regions (new decode map):
sexavg_rec_rates_per_gen = [0.46e-8, # mean in regions
                            1.16e-8] # global for chrX
autosomal_bottle_N = 3000
generations = int(10000 / 29)

# parameters for computation
N = [2 * autosomal_bottle_N * x for x in nonafr_x_auto_ratios]
target_freqs = np.linspace(0.1, 1, 10)
g = [generations, generations*2]
r = sexavg_rec_rates_per_gen
L = [5e5, 1e6]
n_samples = 1000000
parameters = list(itertools.product(N, target_freqs, g, r, L))

# comopute
records = list()
for args in parameters:
    nr_freq_reached, p = prob_norec_freq_reached(*args, n_samples=n_samples)
    records.append((*args, p, nr_freq_reached))
df = pd.DataFrame.from_records(records, columns=['N', 'f', 'g', 'r', 'L', 'p', 'n'])


df['n_samples'] = n_samples

df.to_hdf(output_file_name, 'df', format='table', mode='w')



# def norec_freq_reached_by_any(N, allele_count, g, r, L, n_samples=None):
#     min_allele_count_reached = 0
#     if 'SLURM_CPUS_PER_TASK' in os.environ:
#         # in parallel:
#         with Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as p:
#             min_allele_count_reached = sum(p.starmap(work, [(N, allele_count, g, r, L)]*n_samples))
#     else:
#         # single core:
#         min_allele_count_reached = sum(itertools.starmap(work, [(N, allele_count, g, r, L)]*n_samples))            
#     prob_min_allele_count_reached = min_allele_count_reached / n_samples 
#     # prob_min_allele_count_reached_by_any = prob_min_allele_count_reached * (1 - (1-prob_min_allele_count_reached)**N)
#     prob_min_allele_count_reached_by_any = (1 - (1-prob_min_allele_count_reached)**N)
#     return prob_min_allele_count_reached_by_any


# # probability that it happened to any of the L-sized windows on the chromosome
# # df['prob_none_on_chr'] = df.p * 153e6 / df.L
# df['prob_none_on_chr'] = 1 - (1 - df.p)**(153e6 / df.L)

# # upper bound on probability above probability if it is reported as zero do to insuficcient sampling
# #df['limit_prob_none_on_chr'] = (1/n_samples) * ( np.exp(-df.r*df.L*df.g) + (1-np.exp(-df.r*df.L*df.g))/2 )  * (1 - (1-(1/n_samples))**df.N) * 153e6 / df.L
# limit_prob = (1/n_samples) * ( np.exp(-df.r*df.L*df.g) + (1-np.exp(-df.r*df.L*df.g))/2 )
# limit_prob_none = 1 - (1 - limit_prob**df.N)
# limit_prob_none_on_chr = 1 - (1 - limit_prob_none)**(153e6 / df.L)
# df['limit_prob_none_on_chr'] = limit_prob_none_on_chr
