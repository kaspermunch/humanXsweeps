from itertools import combinations
import argparse
import random
import subprocess
import re, os, sys
import tempfile
import numpy as np
import pandas as pd
from pandas import DataFrame
import msprime, pyslim

import psutil
process = psutil.Process(os.getpid())


initialization = r'''
initialize() {
	initializeTreeSeq();
	initializeMutationRate(0);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeMutationType("m2", 1.0, "f", SEL_COEF);        // introduced
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, 10e6-1);
	initializeRecombinationRate(REC_RATE);
	initializeSex("CHROM");
}
'''

first = r'''
GENERATION {
	defineConstant("simID", getSeed());
	sim.addSubpop("p1", SIZE);
}
'''

size_change = r'''
GENERATION {
	p1.setSubpopulationSize(SIZE);
}
'''

finish = r'''
TOTAL_GEN {
	sim.treeSeqOutput("OUTPUT_FILE");
	sim.simulationFinished();
}
'''

complete_sweep = r'''
SWEEP_START late() {
	target = sample(p1.genomes, 1);
	target.addNewDrawnMutation(m2, 5e6);
	sim.treeSeqOutput("TMPDIR/slim_" + simID + ".trees");
}
SWEEP_START: late() {
	if (sim.countOfMutationsOfType(m2) == 0) {
		if (sum(sim.substitutions.mutationType == m2) == 1) {
			cat(simID + ": FIXED\n");
			sim.deregisterScriptBlock(self);
		} else {
			cat(simID + ": LOST - RESTARTING\n");
			sim.readFromPopulationFile("TMPDIR/slim_" + simID + ".trees");
			setSeed(rdunif(1, 0, asInteger(2^32) - 1));
		}
	}
}
'''

partial_sweep = r'''
SWEEP_START late() {
	target = sample(p1.genomes, 1);
	target.addNewDrawnMutation(m2, 5e6);
	sim.treeSeqOutput("TMPDIR/slim_" + simID + ".trees");
}
SWEEP_START: late() {
	mut = sim.mutationsOfType(m2);
	if (size(mut) == 1)
	{
		if (sim.mutationFrequencies(NULL, mut) > 0.5)
		{
			cat(simID + ": ESTABLISHED\n");
			sim.deregisterScriptBlock(self);
		}
	}
	else
	{
		cat(simID + ": LOST – RESTARTING\n");
		sim.readFromPopulationFile("TMPDIR/slim_" + simID + ".trees");
		setSeed(rdunif(1, 0, asInteger(2^32) - 1));
		
	}
}
'''

drive = '''
// only positive selection in males
fitness(m3) {
	if (individual.sex == 'M') {
		return 1.0 + mut.selectionCoeff;
	} else {
		return 1.0;
	}
}
'''




random.seed(7)

parser = argparse.ArgumentParser()
parser.add_argument("--selcoef", type=float)
parser.add_argument("--window", type=int)
parser.add_argument("--samples", type=int)
parser.add_argument("--mutationrate", type=float)
parser.add_argument("--recrate", type=float)
parser.add_argument("--generationtime", type=int)
parser.add_argument("--sweep", type=str, choices=['partial', 'complete', 'nosweep'])
parser.add_argument("--sweepstart", type=int)
parser.add_argument("--xdrive")
parser.add_argument('--demography', nargs='+', type=str)
parser.add_argument("--popsize", type=int)
parser.add_argument("--chrom", type=str, choices=['X', 'A'])
parser.add_argument("--size-reduction", dest='size_reduction', type=float)
parser.add_argument("--totalgenerations", type=int)
parser.add_argument("--dumpscript", action='store_true')
#parser.add_argument("slurm_script", type=str)
parser.add_argument("trees_file", type=str)
parser.add_argument("dist_file", type=str)
parser.add_argument("sites_file", type=str)
args = parser.parse_args()

window_size = args.window

# slim needs output file to be absolute
if not os.path.isabs(args.trees_file):
    args.trees_file = os.path.abspath(args.trees_file)

# initialization
slurm_script = (initialization.replace('SEL_COEF', str(args.selcoef))
							  .replace('CHROM', args.chrom)
							  .replace('REC_RATE', str(args.recrate))
				)

# end of simulation
slurm_script += finish.replace('TOTAL_GEN', str(args.totalgenerations)).replace('OUTPUT_FILE', str(args.trees_file))

# add size changes
for pair in args.demography:
    gen, size = pair.split(':')

    assert int(gen) <= args.totalgenerations
    
    size = str(int(int(size) * args.size_reduction))

    if gen == '1':
        slurm_script += first.replace('GENERATION', gen).replace('SIZE', size)
    else:
        slurm_script += size_change.replace('GENERATION', gen).replace('SIZE', size)

# add complete or partial sweep
if args.sweep == 'partial':
    slurm_script += partial_sweep.replace('SWEEP_START', str(args.sweepstart))
elif args.sweep == 'complete':
    slurm_script += complete_sweep.replace('SWEEP_START', str(args.sweepstart))


# make positive selection act on X only in males
if args.xdrive:
    assert args.chrom == 'X' and args.sweep != 'nosweep'
    slurm_script += drive

if 'GWF_JOBID' in os.environ:
	slurm_script = slurm_script.replace('TMPDIR', '/scratch/' + os.environ['GWF_JOBID'])
else:
	slurm_script = slurm_script.replace('TMPDIR', '/tmp')

if args.dumpscript:
	print(slurm_script)
	sys.exit()

# write slim script file with the right output name
slurm_script_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
slurm_script_file.write(slurm_script)
slurm_script_file.close()

print('Before running slim:', process.memory_info().rss/1024**3)

# run slim
cmd = './slim {}'.format(slurm_script_file.name)
p = subprocess.Popen(cmd.split(), 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
print(stdout)
print(stderr)
# args.trees_file = '/home/kmt/simons/faststorage/people/kmt/steps/slim/simulations/standard/66/4500/X/nosweep/0/0/standard_66_4500_X_nosweep_0_0_6.trees'

# ts = pyslim.load("./sweep.trees").simplify()
# mutated = msprime.mutate(ts, rate=1e-7, random_seed=1, keep=True)
# mutated.dump("./sweep_overlaid.trees") 

# load trees from slim (simplify removes founder invididual in each sub population and thus support for fixed mutations)
ts = pyslim.load(args.trees_file).simplify()

print('Loaded trees:', process.memory_info().rss/1024**3)
sys.stdout.flush()

# get nodes for female individuals:
female_nodes = list()
for ind in ts.individuals():
    if pyslim.decode_individual(ind.metadata).sex == 0:
        female_nodes.extend(ind.nodes)

# prune tree seqeuence to only retain information about the sampled famale nodes:
subset_ts = ts.simplify(samples=random.sample(female_nodes, k=args.samples))

# overlay mutations
mutated_ts = msprime.mutate(sample_ts, rate=args.mutationrate*args.generationtime, random_seed=7)

# # haplotypes are ordered the same way as samples
# samples = list(ts.samples())
# samples_set = set(samples)
# # so we get the sample indexes of the females to know which haplotypes to consider:
# famale_haplo_idx = list()
# for fem in females:
#     if fem in samples_set:
#         famale_haplo_idx.append(samples.index(fem))

# # random sample among females:
# sample_idx = set(random.sample(famale_haplo_idx, k=args.samples))

# print('Computed indexes:', process.memory_info().rss/1024**3)
# sys.stdout.flush()

# # random indexes for samples
# #sample_idx = set(random.sample(range(ts.num_individuals), args.samples))

# # overlay mutations
# mutated_ts = msprime.mutate(ts, rate=args.mutationrate*args.generationtime, random_seed=7)

print('Mutated trees:', process.memory_info().rss/1024**3)
sys.stdout.flush()

# get the corresponding sample haplotypes
sample = list()
for i, hap in enumerate(mutated_ts.haplotypes()):
	if i in sample_idx:
		print('loading haplotype', i, process.memory_info().rss/1024**3)
		sys.stdout.flush()
		sample.append(hap)

assert len(sample) == args.samples

print('Loaded haplotypes:', process.memory_info().rss/1024**3)
sys.stdout.flush()

# get the positions of each segregating site
positions = [site.position for site in mutated_ts.sites()]  

# make table with sampled haplotypes
table = np.array([list(map(np.int8, hap)) for hap in sample]).transpose()

print('Made table:', process.memory_info().rss/1024**3)
sys.stdout.flush()

# turn table into dataframe with positions
df = DataFrame(table, dtype='int8')
df['pos'] = positions
# write sites to hdf
df.to_hdf(args.sites_file, 'df', format='table', mode='w')

print('Dumped positions to hdf:', process.memory_info().rss/1024**3)
sys.stdout.flush()

# add a row with zeros for the start of each window so there is at least
# one row in each window
zeros = dict((x, 0) for x in range(args.samples))
extra_df = pd.DataFrame({'pos': range(0, int(mutated_ts.sequence_length), window_size), **zeros})
df = df.append(extra_df)

# make a start column grouping all rows in same window
df['start'] = ((df.pos // window_size) * window_size).astype('uint32')
df.drop('pos', axis=1, inplace=True)
df.set_index('start', inplace=True)

def pw_dist(df):
    "computes differences bewteen all pairs in a window"
    pairs = list(combinations(df.columns, 2))
    site_diffs = [np.bitwise_xor(df[p[0]], df[p[1]]) for p in pairs]
    return pd.concat(site_diffs, axis=1, keys=pairs).sum()

# make a dataframe with distance for each pair
pw_dist_df = (
    df
    .groupby('start')
    .apply(pw_dist)
    .reset_index()
    .melt(id_vars=['start'], var_name=['indiv_1', 'indiv_2'], value_name='dist')
    )

print('Computed distances:', process.memory_info().rss/1024**3)
sys.stdout.flush()

# compute proper distance as number of diffs divided by window size
pw_dist_df['dist'] /= window_size

# add end column
pw_dist_df.insert(loc=1, column='end', value=pw_dist_df.start + window_size)

# convert indiv labels from object to int and and write hdf
pw_dist_df['indiv_1'] = pw_dist_df['indiv_1'].astype('uint16')
pw_dist_df['indiv_2'] = pw_dist_df['indiv_2'].astype('uint16')
pw_dist_df.to_hdf(args.dist_file, 'df', format='table', mode='w')

print('Dumped distances to hdf:', process.memory_info().rss/1024**3)
sys.stdout.flush()
