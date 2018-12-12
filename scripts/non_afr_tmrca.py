import sys, os
from pathlib import Path
import numpy
import pandas
from pandas import DataFrame, Series

from ete3 import Tree

numpy.random.seed(7)

import simons_meta_data


def clade_tmrca(tree, clade_leaves, discrete_time_intervals=None):
    """
    Compute the TMRCA of a a set of leaves relative to the TMRCA 
    of all leaves. Used to compute the relative TMRCA of non-Africans.
    """
    add_node_heights(tree, discrete_time_intervals=discrete_time_intervals)
    clade_root = tree.get_common_ancestor(*clade_leaves)
    root = tree.get_tree_root()
    return clade_root.height, root.height

def add_node_heights(tree, discrete_time_intervals):
    """
    Add height of each internal node as node attribute.
    Assumes an (ultrametric) tree.
    """
    if tree.is_leaf():
        tree.height = 0
        return tree.dist, []
    else:
        child_node_heights = list()
        this_node_heights = list()
        for c in tree.children:
            this_height, child_node_hts = add_node_heights(c, discrete_time_intervals)
            this_node_heights.append(this_height)
            child_node_heights.extend(child_node_hts)

        # take care of round off errors:
        this_node_height = numpy.mean(this_node_heights)
        # and check they are no larger than 0.5
        assert max(this_node_heights) - min(this_node_heights) <= 0.5, this_node_heights

        # take care of round off errors to coalescences in same interval get same height
        if discrete_time_intervals:
            this_node_height = round_to_closest(this_node_height, discrete_time_intervals)

        tree.height = this_node_height

        parent_node_height = this_node_height + tree.dist
        return parent_node_height, [this_node_height] + child_node_heights

        
def tree_stats(df):
    stats_list = list()
    for row in df.itertuples():
        start, end = int(row.start), int(row.end)
        tree = Tree(row.tree)  
        stats_list.append((start, end, *clade_tmrca(tree, non_afr_leaves)))
    return DataFrame.from_records(stats_list, 
                                  columns=['start', 'end', 
                                           'non_afr_tmrca', 'tmrca']
                                 )


root_dir = Path(os.environ['HOME'], 'simons/faststorage/people/kmt')
meta_data_dir = Path(os.environ['HOME'], 'simons/faststorage/data/metadata')


scripts_dir = root_dir / 'scripts'
if str(scripts_dir) not in sys.path:
    sys.path.append(str(scripts_dir))

individuals, populations, regions = simons_meta_data.get_meta_data(meta_data_dir=meta_data_dir)
    
non_afr_leaves = list()
for indiv, d in individuals.items():
    if d['Region'] != 'Africa' and d['Genetic sex assignment'] == 'XY':
        non_afr_leaves.append(indiv + '-A')


_, input_table_file, output_hdf_file = sys.argv
#input_table_file = '/home/kmt/simons/faststorage/people/kmt/steps/argweaver/output/World/X-037500000-037600000.tsv.gz'
input_table_df = pandas.read_table(str(input_table_file))

df = input_table_df.groupby(['chain', 'sample']).apply(tree_stats)
df['start'] = input_table_df.start.min()
df['end'] = input_table_df.end.max()

df.to_hdf(output_hdf_file, 'df', format="table", mode='w')    
