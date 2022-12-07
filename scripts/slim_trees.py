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
import time

import psutil
process = psutil.Process(os.getpid())

initialization_with_rec_map = r'''
initialize() {
	initializeTreeSeq();
	initializeMutationRate(0);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeMutationType("m2", 1.0, "f", SEL_COEF);        // introduced
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, 10e6-1);

	lines = readFile("REC_MAP_FILE");
	rates = NULL;
	ends = NULL;
	for (line in lines)
	{
		components = strsplit(line, "\t");
		ends = c(ends, asInteger(components[0]));
		rates = c(rates, asFloat(components[1]));
	}
	ends = c(ends[1:(size(ends)-1)] - 2, 10e6-1);
	rates = rates * 1e-8;
	initializeRecombinationRate(rates, ends);

	initializeSex("CHROM");
}
'''


initialization_with_constant_rec_rate = r'''
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
	sim.treeSeqOutput("TMPDIR/slim_" + simID + ".trees");

	target = sample(p1.genomes, 1);	
	while (target.isNullGenome)
	{
		target = sample(p1.genomes, 1);	
	} 
	target.addNewDrawnMutation(m2, 5e6);
}
SWEEP_START: late() {
	if (sim.countOfMutationsOfType(m2) == 0) {
		if (sum(sim.substitutions.mutationType == m2) == 1) {
			cat(simID + ": FIXED\n");
			sim.deregisterScriptBlock(self);
		} else {
			cat(simID + ": LOST - RESTARTING\n");
			sim.readFromPopulationFile("TMPDIR/slim_" + simID + ".trees");
			setSeed(rdunif(1, 0, asInteger(2^62) - 1));

			target = sample(p1.genomes, 1);	
			while (target.isNullGenome)
			{
				target = sample(p1.genomes, 1);	
			} 
			target.addNewDrawnMutation(m2, 5e6);
		}
	}
}
'''

# partial_sweep = r'''
# SWEEP_START late() {
# 	sim.treeSeqOutput("TMPDIR/slim_" + simID + ".trees");


	# target = sample(p1.genomes, 1);	
	# while (target.isNullGenome)
	# {
	# 	target = sample(p1.genomes, 1);	
	# } 
# 	target.addNewDrawnMutation(m2, 5e6);
# }
# SWEEP_START:SWEEP_ESTABLISHED late() {
# 	mut = sim.mutationsOfType(m2);
# 	if (size(mut) == 1)
# 	{
# 		if (sim.mutationFrequencies(NULL, mut) > 0.5)
# 		{
# 			cat(simID + ": ESTABLISHED\n");
# 			sim.deregisterScriptBlock(self);
# 		}
# 	}
# 	else
# 	{
# 		cat(simID + ": LOST – RESTARTING\n");
# 		sim.readFromPopulationFile("TMPDIR/slim_" + simID + ".trees");
#  		setSeed(rdunif(1, 0, asInteger(2^62) - 1)); 	

	# target = sample(p1.genomes, 1);	
	# while (target.isNullGenome)
	# {
	# 	target = sample(p1.genomes, 1);	
	# } 
    # target.addNewDrawnMutation(m2, 5e6);	
# 	}
# }
# SWEEP_ESTABLISHED: late() {
# 	mut = sim.mutationsOfType(m2);
# 	if (size(mut) == 1)
# 	{
# 		if (sim.mutationFrequencies(NULL, mut) <= 0.5)
# 		{
# 			cat(simID + ": LOST – RESTARTING\n");
# 			sim.readFromPopulationFile("TMPDIR/slim_" + simID + ".trees");
# 			setSeed(rdunif(1, 0, asInteger(2^62) - 1)); 	

# 			target = sample(p1.genomes, 1);
	# while (target.isNullGenome)
	# {
	# 	target = sample(p1.genomes, 1);	
	# } 
# 			target.addNewDrawnMutation(m2, 5e6);	
# 		}
# 	}
# }
# '''

partial_sweep = r'''
SWEEP_START late() {
	sim.treeSeqOutput("TMPDIR/slim_" + simID + ".trees");

	target = sample(p1.genomes, 1);
	while (target.isNullGenome)
	{
		target = sample(p1.genomes, 1);	
	} 	
	target.addNewDrawnMutation(m2, 5e6);
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
 		setSeed(rdunif(1, 0, asInteger(2^62) - 1)); 	

		target = sample(p1.genomes, 1);	
		while (target.isNullGenome)
		{
			target = sample(p1.genomes, 1);	
		} 		
		target.addNewDrawnMutation(m2, 5e6);	
	}
}
'''

selection_episode = r'''
SWEEP_START late() {
	target = sample(p1.genomes, 1);	
	while (target.isNullGenome)
	{
		target = sample(p1.genomes, 1);	
	} 	
	target.addNewDrawnMutation(m2, 5e6);
}
SELECTION_END late() {
	m2muts = sim.mutationsOfType(m2);
	for (index in seqAlong(m2muts))
		m2muts[index].setSelectionCoeff(0.0);
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

#random.seed(7)

parser = argparse.ArgumentParser()
parser.add_argument("--selcoef", type=float)
parser.add_argument("--window", type=int)
parser.add_argument("--samples", type=int)
parser.add_argument("--mutationrate", type=float)
parser.add_argument("--recrate", type=float)
parser.add_argument("--recmap", type=str)
parser.add_argument("--generationtime", type=int)
parser.add_argument("--sweep", type=str, choices=['partial', 'complete', 'nosweep', 'episode'])
parser.add_argument("--sweepstart", type=int)
parser.add_argument("--selectionend", type=int)
parser.add_argument("--sweepestablished", type=int)
parser.add_argument("--xdrive")
parser.add_argument('--demography', nargs='+', type=str)
parser.add_argument("--popsize", type=int)
parser.add_argument("--chrom", type=str, choices=['X', 'A'])
parser.add_argument("--size-reduction", dest='size_reduction', type=float)
parser.add_argument("--totalgenerations", type=int)
parser.add_argument("--dumpscript", action='store_true')
parser.add_argument("trees_file", type=str)
parser.add_argument("dist_file", type=str)
parser.add_argument("sites_file", type=str)
parser.add_argument("vcf_file", type=str)
parser.add_argument("vcf_geno_file", type=str)
args = parser.parse_args()


if args.sweep == 'episode' and args.selectionend is None:
    # if selectionend is not specified we assume episide is 10,000 years:
    args.selectionend = args.sweepstart + int(10000/args.generationtime)

    
window_size = args.window

# slim needs output file to be absolute
if not os.path.isabs(args.trees_file):
    args.trees_file = os.path.abspath(args.trees_file)

assert bool(args.recrate) != bool(args.recmap), "Give either recmap or recrate"
# initialization
if args.recmap:
	slurm_script = (initialization_with_rec_map.replace('SEL_COEF', str(args.selcoef))
								.replace('CHROM', args.chrom)
								.replace('REC_MAP_FILE', str(args.recmap))
					)
else:
	slurm_script = (initialization_with_constant_rec_rate.replace('SEL_COEF', str(args.selcoef))
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
    # slurm_script += partial_sweep.replace('SWEEP_START', str(args.sweepstart)).replace('SWEEP_ESTABLISHED', str(args.sweepestablished))
    slurm_script += partial_sweep.replace('SWEEP_START', str(args.sweepstart))
elif args.sweep == 'complete':
    # slurm_script += complete_sweep.replace('SWEEP_START', str(args.sweepstart)).replace('SWEEP_ESTABLISHED', str(args.sweepestablished))
    slurm_script += complete_sweep.replace('SWEEP_START', str(args.sweepstart))
elif args.sweep == 'episode':
    # slurm_script += complete_sweep.replace('SWEEP_START', str(args.sweepstart)).replace('SWEEP_ESTABLISHED', str(args.sweepestablished))
    slurm_script += selection_episode.replace('SWEEP_START', str(args.sweepstart)).replace('SELECTION_END', str(args.selectionend))

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

# run slim
cmd = './slim3.5/build/slim {}'.format(slurm_script_file.name)
p = subprocess.Popen(cmd.split(), 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
print(stdout)
print(stderr)

# load trees from slim (simplify removes founder invididual in each sub population and thus support for fixed mutations)
ts = pyslim.load(args.trees_file).simplify()

# get nodes/chromosomes for female individuals:
female_nodes = list()
for ind in ts.individuals():
    # if pyslim.decode_individual(ind.metadata).sex == 0:
    if ind.metadata['sex'] == 0:
        female_nodes.extend(ind.nodes)

# get the asmple ids among females chromosomes:
sample_nodes = random.sample(female_nodes, k=args.samples)

# overlay mutations
mutated_ts = msprime.mutate(ts, rate=args.mutationrate*args.generationtime)

# get the positions of each segregating site
positions = [site.position for site in mutated_ts.sites()]  

# get genotypes for sample at variant sites in population:
variants = mutated_ts.variants(samples=sample_nodes, as_bytes=False, 
	impute_missing_data=False)
	# isolated_as_missing=False)
table = np.array([var.genotypes for var in variants])

# turn table into dataframe with positions
df = DataFrame(table, dtype='int8')
df['pos'] = positions
# write sites to hdf
df.to_hdf(args.sites_file, 'df', format='table', mode='w')

# write a VCF
vcf_df = DataFrame(table, dtype='int8')
# remove rows without derived variants
polymorphic = (vcf_df != 0).any(axis=1)
vcf_df = vcf_df.loc[polymorphic]
positions = pd.Series([round(p) for p in positions])
positions = positions[polymorphic]
vcf_df.insert(0, 'FORMAT', 'GT')
vcf_df.insert(0, 'INFO', 'AA=A')
vcf_df.insert(0, 'FILTER', '')
vcf_df.insert(0, 'QUAL', '')
vcf_df.insert(0, 'ALT', 'T')
vcf_df.insert(0, 'REF', 'A')
vcf_df.insert(0, 'ID', list(range(len(positions))))
vcf_df.insert(0, 'POS', positions.values)
vcf_df.insert(0, '#CHROM', '1')
vcf_df.to_csv(args.vcf_file, sep='\t', index=False)

# write a VCF where the haplotypes are artificially treated as unphased diploid genotypes
geno_vcf_df = vcf_df.iloc[:, 0:9]
samples = vcf_df.columns.values[9:].tolist()
for i, (a, b) in enumerate(zip(samples[0::2], samples[1::2])):
    geno_vcf_df[i] = vcf_df[a].astype('str') + '/' +  vcf_df[b].astype('str')
geno_vcf_df.to_csv(args.vcf_geno_file, sep='\t', index=False)


# add a row with zeros for the start of each window so there is at least
# one row in each window
zeros = dict((x, 0) for x in range(args.samples))
extra_df = pd.DataFrame({'pos': range(0, int(mutated_ts.sequence_length), window_size), **zeros})
df = df.append(extra_df, sort=True)

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

# compute proper distance as number of diffs divided by window size
pw_dist_df['dist'] /= window_size

# add end column
pw_dist_df.insert(loc=1, column='end', value=pw_dist_df.start + window_size)

# convert indiv labels from object to int and and write hdf
pw_dist_df['indiv_1'] = pw_dist_df['indiv_1'].astype('uint16')
pw_dist_df['indiv_2'] = pw_dist_df['indiv_2'].astype('uint16')
pw_dist_df.to_hdf(args.dist_file, 'df', format='table', mode='w')

