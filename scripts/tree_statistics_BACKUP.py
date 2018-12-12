
import math, re, collections, random, bisect
from pathlib import Path
from ete3 import Tree
import numpy
import scipy
import pandas
import argparse

# def node_heights(tree):
#     """
#     Get height of each internal node in an ultramemtric 
#     tree from a postorder traversal. This traversal uses
#     the order of children as in the t.children attribute 
#     (not the same order as t.traverse('postorder') ).
#     """
#     if tree.is_leaf():
#         return tree.dist, []
#     else:
#         child_node_heights = list()
#         this_node_heights = list()
#         for c in tree.children:
#             this_height, child_node_hts = node_heights(c)
#             this_node_heights.append(this_height)
#             child_node_heights.extend(child_node_hts)

#         # take care of round off errors:
#         this_node_height = np.mean(this_node_heights)
#         # and check they are no larger than 0.5
#         assert max(this_node_heights) - min(this_node_heights) <= 0.5, this_node_heights

#         parent_node_height = this_node_height + tree.dist
#         return parent_node_height, [this_node_height] + child_node_heights

    
# def node_descendant_counts(tree):
#     """
#     Get number of leaves below each internal node in an 
#     ultramemtric tree from a postorder traversal. The order of 
#     nodes returned is the same as for node_heights. This 
#     traversal uses the order of children as in the t.children 
#     attribute (not the same order as t.traverse('postorder') ).
#     """
#     if tree.is_leaf():
#         return 1, []
#     else:
#         child_leaf_counts = list()
#         total_leaves = 0
#         for c in tree.children:
#             child_leaf_count, child_leaf_list = node_descendant_counts(c)
#             child_leaf_counts.extend(child_leaf_list)
#             total_leaves += child_leaf_count
#         return total_leaves, [total_leaves] + child_leaf_counts

    
# def tmrca_stats(tree):
#     """
#     Compute:
#     TMRCA: time to TMRCA of tree.
#     TMRCA_half: time until at a clade consist of at least 50% of sample.
#     Coal_half: time until at least half the sample is part of any subtree (have coalesced).
#                Same as TMRCA_half but not requireing that they are in a single subtree.
#     """
#     total_height, heights = node_heights(tree)
#     total_leaves, leaf_counts = node_descendant_counts(tree)
#     min_half = math.ceil(len(tree)/2)
    
#     coal_half = sorted(heights)[min_half - 2] # -2 because index 0 is the height where two has coalesced
#     tmrca = max(heights)
#     for lc, h in sorted(zip(leaf_counts, heights)):
#         if lc >= min_half:
#             tmrca_half = h
#             break
#     return tmrca, tmrca_half, coal_half

    
# def count_sweep_nodes(tree, max_branch_len):
#     """
#     Compute sweep node count (SwNC):
#     Number of nodes with a short dist up and at least one sort dist down.
#     Top node with no parent is not included.
#     """
#     count = 0
#     for n in tree.traverse('postorder'):
#         has_short_child_branch = any(c.dist <= max_branch_len for c in n.children)
#         if n.up and n.dist <= max_branch_len and has_short_child_branch:
#             count += 1
#     return count


# def count_shortest_branches_within_tree_fraction(tree, tree_fraction):
#     """
#     Compute short branch enrichment (ShBE):
#     Count how many of the shortest branches sort brances
#     that together make up no more than e.g. 10% of total tree length.
#     """
#     all_branch_lengths = sorted(n.dist for n in tree.traverse('postorder') if n.up)
#     total_tree = sum(all_branch_lengths)
#     short_branches = list()
#     tot = 0
#     for b in all_branch_lengths:
#         if tot <= tree_fraction * total_tree:
#             short_branches.append(b)
#         else:
#             break
#         tot += b
#     return len(short_branches)
 

# def find_max_connected_short_branches(tree, max_branch_len):
#     """
#     Compute largest short branch tree (LShBT):
#     Find the largest set of connected branches that are all
#     smaller or equal to some maximum length.
#     """
#     def get_connected_short_branches(tree, max_branch_len):
#         if tree.is_leaf():
#             return 1
#         else:
#             short_count = 0
#             if tree.up and tree.dist <= max_branch_len:
#                 short_count += 1
#             for c in tree.children:       
#                 if c.dist <= max_branch_len:
#                     short_count += get_connected_short_branches(c, max_branch_len)
#             return short_count

#     max_count = 0
#     for n in tree.traverse():
#         count = get_connected_short_branches(n, max_branch_len)
#         max_count = max(count, max_count)
#     return max_count

#################################################

##################################################
# Helper functions
##################################################

def node_heights(tree):
    """
    Get height of each internal node in an ultramemtric 
    tree from a postorder traversal. This traversal uses
    the order of children as in the t.children attribute 
    (not the same order as t.traverse('postorder') ).
    """
    if tree.is_leaf():
        return tree.dist, []
    else:
        child_node_heights = list()
        this_node_heights = list()
        for c in tree.children:
            this_height, child_node_hts = node_heights(c)
            this_node_heights.append(this_height)
            child_node_heights.extend(child_node_hts)

        # take care of round off errors:
        this_node_height = numpy.mean(this_node_heights)
        # and check they are no larger than 0.5
        assert max(this_node_heights) - min(this_node_heights) <= 0.5, this_node_heights

        parent_node_height = this_node_height + tree.dist
        return parent_node_height, [this_node_height] + child_node_heights

    
def node_descendant_counts(tree):
    """
    Get number of leaves below each internal node in an 
    ultramemtric tree from a postorder traversal. The order of 
    nodes returned is the same as for node_heights. This 
    traversal uses the order of children as in the t.children 
    attribute (not the same order as t.traverse('postorder') ).
    """
    if tree.is_leaf():
        return 1, []
    else:
        child_leaf_counts = list()
        total_leaves = 0
        for c in tree.children:
            child_leaf_count, child_leaf_list = node_descendant_counts(c)
            child_leaf_counts.extend(child_leaf_list)
            total_leaves += child_leaf_count
        return total_leaves, [total_leaves] + child_leaf_counts

    
def round_to_closest(n, vals):
    """
    Round to closest discrete time interval.
    """
    pos = bisect.bisect_left(vals, n)
    if pos == 0:
        return vals[0]
    if pos == len(vals):
        return vals[-1]
    before = vals[pos - 1]
    after = vals[pos]
    if after - n < n - before:
        return after
    else:
        return before


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

    
def traverse_tree_by_height(t):
    return sorted(t.traverse(), key=lambda x: x.height, reverse=True)


def traverse_descendants_by_height(t):
    return [n for n in traverse_tree_by_height(t) if n != t]


def add_node_levels(tree):
    """
    Add the coalescence level of each node as attribute.
    Root node is level 1.
    """
    level = 0 # level in time (may represent more than one
              #simultaneous coalscence as happens with discretization)
    pseudolevel = 0 # level randomly separating any simultaneous 
                    # coalescences into seperate levels
    prev_height = float('inf')
    for node in traverse_tree_by_height(tree):
        pseudolevel += 1
        if node.height < prev_height:
            level += 1
            prev_height = node.height
        node.level = level
        node.pseudolevel = pseudolevel

        
##################################################
# Misc tree stats
##################################################
  
def tmrca_stats(tree):
    """
    Compute:
    TMRCA: time to TMRCA of tree.
    TMRCA_half: time until at a clade consist of at least 50% of sample.
    Coal_half: time until at least half the sample is part of any subtree (have coalesced).
               Same as TMRCA_half but not requireing that they are in a single subtree.
    """
    total_height, heights = node_heights(tree)
    total_leaves, leaf_counts = node_descendant_counts(tree)
    min_half = math.ceil(len(tree)/2)
    
    coal_half = sorted(heights)[min_half - 2] # -2 because index 0 is the height where two has coalesced
    tmrca = max(heights)
    for lc, h in sorted(zip(leaf_counts, heights)):
        if lc >= min_half:
            tmrca_half = h
            break
    return tmrca, tmrca_half, coal_half

    
def count_sweep_nodes(tree, max_branch_len):
    """
    Compute sweep node count (SwNC):
    Number of nodes with a short dist up and at least one sort dist down.
    Top node with no parent is not included.
    """
    count = 0
    for n in tree.traverse('postorder'):
        has_short_child_branch = any(c.dist <= max_branch_len for c in n.children)
        if n.up and n.dist <= max_branch_len and has_short_child_branch:
            count += 1
    return count


def count_shortest_branches_within_tree_fraction(tree, tree_fraction):
    """
    Compute short branch enrichment (ShBE):
    Count how many of the shortest branches sort brances
    that together make up no more than e.g. 10% of total tree length.
    """
    all_branch_lengths = sorted(n.dist for n in tree.traverse('postorder') if n.up)
    total_tree = sum(all_branch_lengths)
    short_branches = list()
    tot = 0
    for b in all_branch_lengths:
        if tot <= tree_fraction * total_tree:
            short_branches.append(b)
        else:
            break
        tot += b
    return len(short_branches)
 

def find_max_connected_short_branches(tree, max_branch_len):
    """
    Compute largest short branch tree (LShBT):
    Find the largest set of connected branches that are all
    smaller or equal to some maximum length.
    """
    def get_connected_short_branches(tree, max_branch_len):
        if tree.is_leaf():
            return 1
        else:
            short_count = 0
            if tree.up and tree.dist <= max_branch_len:
                short_count += 1
            for c in tree.children:       
                if c.dist <= max_branch_len:
                    short_count += get_connected_short_branches(c, max_branch_len)
            return short_count

    max_count = 0
    for n in tree.traverse():
        count = get_connected_short_branches(n, max_branch_len)
        max_count = max(count, max_count)
    return max_count


def rel_clade_tmrca(tree, clade_leaves, discrete_time_intervals=[]):
    """
    Compute the TMRCA of a a set of leaves relative to the TMRCA 
    of all leaves. Meant to compute the relative TMRCA of non-Africans.
    """
    add_node_heights(tree, discrete_time_intervals=discrete_time_intervals)
    clade_root = tree.get_common_ancestor(*clade_leaves)
    root = tree.get_tree_root()
    return clade_root.height / root.height


def rel_cross_clade_tmrca(tree, clade_leaves, discrete_time_intervals=[], cross_to=None):
    """
    Compute time of the coalescence event between the common ancestor of 
    the given set of leaves and a leaf not in the specified 
    set of leaves. Normalized with TMRCA of full tree. 
    (Meant to compute the first coalesce event connecting *all* non-european 
    lineages with an African lineage.)
    If cross_to is a list of nodes then we only consider the tree that 
    connects the set of clade_leaves + cross_to. That way the final "cross" 
    coalescence can only be to a lineage in the cross_to set.
    (Meant to be used to only consider cross to a particular set of africans)
    """

    # only consider the subset of clade_leaves and cross_to that are actually in the tree:
    tree_leaves = set(tree.get_leaf_names())
    clade_leaves = [x for x in clade_leaves if x in tree_leaves]
    if cross_to is not None:
        cross_to = [x for x in cross_to if x in tree_leaves]

    if cross_to is None:
        add_node_heights(tree, discrete_time_intervals=discrete_time_intervals)
        clade_root = tree.get_common_ancestor(*clade_leaves)
        root = tree.get_tree_root()
    else:
        # We only consider the part of the tree made up of clade_leaves and cross_to leaves
        subtree_leaves = clade_leaves + cross_to
        tree = tree.copy()
        add_node_heights(tree, discrete_time_intervals=discrete_time_intervals)
        tree.prune(subtree_leaves, preserve_branch_length=True)
        clade_root = tree.get_common_ancestor(*clade_leaves)
        # we cannot use the root method here because the root
        # of the pruned tree is the original root:
        root = tree.get_common_ancestor(*subtree_leaves)

    if clade_root == root:
        return 1
    else:
        return clade_root.up.height / root.height


##################################################
# Sweep statistic
##################################################

prop_coalescences_n_to_l_cache = dict()
prop_ancestor_to_i_cache = dict()

def sweep_stat(l, i, n, t, pop_size):
    """
    Compute probability that i among n lingeages 
    coalesce into one among l lineages in time t.
    """

    tup = (n, l, t)
    if tup in prop_coalescences_n_to_l_cache:
        prop_coalescences_n_to_l = prop_coalescences_n_to_l_cache[tup]
    else:
        # compute prob of coalescing from n to l lineages
        Q = numpy.zeros(shape=(n-l+1, n-l+1))
        for x in range(n-l):
            rate = scipy.misc.comb(n-x, 2) / (2 * pop_size)
            Q[x][x] = -rate
            Q[x][x+1] = rate
        prop_coalescences_n_to_l = scipy.linalg.expm(Q*t)[0][-1]
        prop_coalescences_n_to_l_cache[tup] = prop_coalescences_n_to_l

    # prob that one among l lineages is ancestral to i among n leaves
    tup = (n, i, l)
    if i == 1 or l == 1:
        prop_ancestor_to_i = 1 # CHECK THAT THIS IS OK
    else:
        if tup in prop_ancestor_to_i_cache:
            prop_ancestor_to_i = prop_ancestor_to_i_cache[tup]
        else:
            prop_ancestor_to_i = scipy.misc.comb(n-i-1, l-2) / scipy.misc.comb(n-1, l-1)
            prop_ancestor_to_i_cache[tup] = prop_ancestor_to_i

    assert prop_ancestor_to_i <= 1, prop_ancestor_to_i
    assert prop_coalescences_n_to_l <= 1, prop_coalescences_n_to_l
        
    prob = prop_coalescences_n_to_l * prop_ancestor_to_i
    p_score = -numpy.log10(prob)

    return p_score


def sweep_score(tree, pop_size, min_score, discrete_time_intervals=None):
    
    # Compute the time ranges that the discrete time represents
    if discrete_time_intervals is not None:
        a = numpy.array([0] + list((numpy.diff(discrete_time_intervals) / 2)) + [numpy.inf])
        b = a.cumsum()
        time_intervals = list(zip(b[:-1], b[1:]))
        time_ranges_for_discrete_times = dict(zip(discrete_time_intervals, time_intervals))

    # annotate tree with node heights
    add_node_heights(tree, discrete_time_intervals)
    # annotate tree with node levels (assumes a tree with node heights)
    add_node_levels(tree)

    # compute the statistic for each clade and time interval in the tree (n**2)
    sscore = 0
    sweep_node = None
    sweep_time = None
    sweep_scores = list()
    sweep_times = list()
    sweep_leaf_sets = list()

    sweep_descendants = set()

    for node in traverse_tree_by_height(tree):

        # if node is a descendant of a previously identified sweep node, then continue
        # use node.iter_descendants() to maintain a set() of all nodes under previously identified sweep nodes

        # then we might not need to do the checks for nested leaf_sets at the end...

        if node in sweep_descendants:
            continue

        # leaves are visited last and cannot be sweep nodes
        if node.is_leaf():
            break

        l = node.level
        i = 1
        smax = 0
        smax_t = None

#        for descendant in traverse_descendants_by_height(node):
        for descendant in traverse_tree_by_height(node):
            if descendant is node:
                # do not evaluate node
                continue

            i += 1
            if discrete_time_intervals is not None:
                # find interval that node is in
                node_time_interval = time_ranges_for_discrete_times[node.height]
                # find interval that descendant is in
                descendant_time_interval = time_ranges_for_discrete_times[descendant.height]
                t = max(node_time_interval) - min(descendant_time_interval)
            else:
                t = node.height - descendant.height

            n = descendant.pseudolevel # we use pseudolevel to ensure that n is incrementally
                                       # increased as we include more nodes in the clade.            
            assert l < n, (l, i, n)
            assert n - l >= i-1, (l, i, n) # nr coalescences must be at least as large as clade size
            assert t > 0, t

            s = sweep_stat(l, i, n, t, pop_size=pop_size)
            
            if s > smax:
                smax = s
                smax_t = t
            else:
                break

            # we break when we hit the first leaf, in which case
            # we have reached the largest possible clade size
            if descendant.is_leaf():
                break


        if sscore is None or smax > sscore:
            sscore = smax
            sweep_time = smax_t
            sweep_node = node
            
        # if the score pass threshhold and the (smaller) leaf set is part of a previous
        # larger one, we add the set of leafs as a sweep
        if sscore >= min_score:
            sweep_leaves = set(node.get_leaf_names())
            if not any(sweep_leaves.issubset(x) for x in sweep_leaf_sets):
                sweep_scores.append(sscore)
                sweep_times.append(sweep_time)
                sweep_leaf_sets.append(sweep_leaves)

#            print(len(sweep_leaf_sets))
            

    return sweep_scores, sweep_times, sweep_leaf_sets


def collect_components(graph, threshold):
    components = []
    processed = set()

    def dfs(v, current_component):
        current_component.add(v)
        processed.add(v)
        neighbours = numpy.where(graph[v,:] >= threshold)[0]
        for w in neighbours:
            if w not in processed:
                dfs(w, current_component)

    for v in range(graph.shape[0]):
        if v not in processed:
            current_component = set()
            dfs(v, current_component)
            components.append(current_component)

    return components


def component_stats(component, graph, threshold):
    edge_weights = list()
    for i, n1 in enumerate(component):
        for j, n2 in enumerate(component):
            if i != j and graph[i][j] >= threshold:
                edge_weights.append(graph[i][j])
                
    nr_nodes = len(component)
    connectedness = len(edge_weights) / (nr_nodes**2 - nr_nodes)
    min_weight = edge_weights
    return connectedness, min(edge_weights), max(edge_weights), sum(edge_weights) / len(edge_weights)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--min-clade-score", dest='min_score', type=float)
    parser.add_argument("--min-post-prob", dest='min_post_prob', type=float)
    parser.add_argument("--pop-size", dest='pop_size', type=int)
    parser.add_argument("--discretization", dest='discretization', type=str)
    parser.add_argument("input_table_file", type=Path)
    parser.add_argument("stats_hdf_file", type=Path)
    parser.add_argument("components_hdf_file", type=Path)
    parser.add_argument("component_stats_hdf_file", type=Path)
    args = parser.parse_args()

    discretization = list(map(float, args.discretization.split(',')))

    import simons_meta_data
    individuals, populations, regions = simons_meta_data.get_meta_data()

    # args are used inside this function:
    def add_stats(df):

        graph = None        
        one_int16 = numpy.int16(1)
        zero_int16 = numpy.int16(0)

        leaf_idx = None

        chrom = str(df.chrom.unique()[0])
        window_start = df.start.min()
        window_end = df.end.max()
        window_size = window_end - window_start

        # we need to loop over the mcmc samples too
        groups = df.groupby('sample')

        # number of mcmc samples
        nr_samples = len(groups)

        # list to collect data frames with per-tree stats for each mcmc sample
        group_df_list = list()

        for sample, group in groups:

            group_df = group.copy()

            stats = collections.defaultdict(list)

            for row in group_df.itertuples():

                #start, end coordinates for tree
                start, end = int(row.start), int(row.end)

                # tree parse is slow so we only do this once for all stats
                tree = Tree(row.tree) 
                if graph is None:
                    # we need the first tree to read the number of individuals in the trees

                    # leaf names
                    leaf_names = [x.name for x in tree.get_leaves()]

                    # mapping from names to an index
                    leaf_idx = dict(zip(sorted(leaf_names), range(len(leaf_names))))

                    # intialize the graph to collect coponents from at the end
                    graph = numpy.zeros(shape=(window_size, len(leaf_idx), len(leaf_idx)), dtype=numpy.int16)

                    # get non-african individual names
                    african_leaves = list()
                    nonafrican_leaves = list()
                    for n in leaf_names:
                        if any(n.startswith(x) for x in regions['Africa']):
                            african_leaves.append(n)
                        else:
                            nonafrican_leaves.append(n)

                # # collect various stats for tree:
                # stats['sgm_len'].append(row.end - row.start)
                # tmrca, tmrca_half, coal_half = tmrca_stats(tree)
                # stats['my_tmrca'].append(tmrca)
                # stats['my_tmrca_half'].append(tmrca_half)
                # stats['coal_half'].append(coal_half)
                # stats['shbe'].append(count_shortest_branches_within_tree_fraction(tree, 0.1))
                # stats['lshbt'].append(find_max_connected_short_branches(tree, 0))
                # stats['swnc'].append(count_sweep_nodes(tree, 0))

                # if african_leaves and nonafrican_leaves:
                #     # this is hte "World" analysis
                #     stats['wavestat1'] = rel_clade_tmrca(tree, nonafrican_leaves, discretization)
                #     stats['wavestat2'] = rel_cross_clade_tmrca(tree, nonafrican_leaves, discretization)
                #     stats['wavestat3'] = rel_cross_clade_tmrca(tree, nonafrican_leaves, discretization, 
                #         cross_to=[x+'-A' for x in simons_meta_data.pure_africans])
                #     stats['wavestat4'] = rel_cross_clade_tmrca(tree, nonafrican_leaves, discretization, 
                #         cross_to=[x+'-A' for x in simons_meta_data.group1_pure_africans])
                #     stats['wavestat5'] = rel_cross_clade_tmrca(tree, nonafrican_leaves, discretization, 
                #         cross_to=[x+'-A' for x in simons_meta_data.group2_pure_africans])
                # else:
                #     stats['wavestat1'] = numpy.nan
                #     stats['wavestat2'] = numpy.nan
                #     stats['wavestat3'] = numpy.nan
                #     stats['wavestat4'] = numpy.nan
                #     stats['wavestat5'] = numpy.nan

                # compute the sweep stat and get sweep leaf sets for this tree
                scores, times, leaf_sets = sweep_score(tree, 
                    pop_size=args.pop_size, 
                    min_score=args.min_score, 
                    discrete_time_intervals=discretization)

                if scores:
                    max_sweep_score = max(scores)
                else:
                    max_sweep_score = numpy.nan
                stats['max_sweep_score'].append(max_sweep_score)

                # print("#", start, end, len(leaf_sets))
                if leaf_sets:
                    # add sweep clades to graph as connections between nodes in same sweep
                    for leaf_set in leaf_sets:
                        # print(start, end, len(leaf_set))
                        # idxs = [leaf_idx[x] for x in leaf_set]
                        # for i in range(start, end):
                        #     for j in idxs:
                        #         for k in idxs:
                        #             if j != k:
                        #                 graph[i, j, k] += one_int16
                        idxs = [leaf_idx[x] for x in leaf_set]
                        # vec = numpy.zeros(shape=len(leaf_idx), dtype=numpy.int16)
                        # for x in leaf_set:
                        #     vec[leaf_idx[x]] = one_int16
                        # for i in range(start, end):                            
                        #     # for j in idxs:
                        #     #     graph[i, j, :] += vec
                        #     #     graph[i, :, j] += vec
                        #     graph[i, idxs, :] += vec
                        #     graph[i, :, idxs] += vec

                        graph[start:end, idxs, idxs] += one_int16

                        # print('==')
            # add the tree stats to the row for each tree:
            # group_df['my_tmrca'] = stats['my_tmrca']
            # group_df['my_tmrca_half'] = stats['my_tmrca_half']
            # group_df['coal_half'] = stats['coal_half']
            # group_df['shbe'] = stats['shbe']
            # group_df['lshbt'] = stats['lshbt']
            # group_df['swnc'] = stats['swnc']
            # group_df['sgm_len'] = stats['sgm_len']
            # group_df['wavestat1'] = stats['wavestat1']
            # group_df['wavestat2'] = stats['wavestat2']
            # group_df['wavestat3'] = stats['wavestat3']
            # group_df['wavestat4'] = stats['wavestat4']
            # group_df['wavestat5'] = stats['wavestat5']

            group_df['max_sweep_score'] = stats['max_sweep_score']

            group_df_list.append(group_df.reset_index(level=['sample']))

            # print("** REMOVE THIS BREAK ** REMOVE THIS BREAK ** REMOVE THIS BREAK **")
            # break


        threshold = args.min_post_prob * nr_samples

        component_df_rows = list()
        component_stats_df_rows = list()
        for i in range(window_size):
            numpy.fill_diagonal(graph[i], zero_int16) # set diagonal to zero
            site_components = collect_components(graph[i], threshold)
            # only non-trivial components larger than one
            site_components = [x for x in site_components if len(x) > 1]

            # loop sweeps
            for j, comp in enumerate(site_components):
                # add a row for each individual in the sweep
                for idx in comp:
                    component_df_rows.append((chrom, i, j, idx))
                # add a row with stats for the sweep
                connectedness, min_edge_weight, max_edge_weight, mean_edge_weight = component_stats(comp, graph[i], threshold)
                component_stats_df_rows.append((chrom, i, j, connectedness, min_edge_weight, max_edge_weight, mean_edge_weight))

        components_df = pandas.DataFrame.from_records(component_df_rows, columns=['chrom', 'pos', 'sweep', 'pop_idx'])
        component_stats_df = pandas.DataFrame.from_records(component_stats_df_rows, columns=['chrom', 'pos', 'sweep', 
            'connectedness', 'min_edge_weight', 'max_edge_weight', 'mean_edge_weight'])
        component_stats_df['threshold'] = threshold

        return pandas.concat(group_df_list), components_df, component_stats_df



    input_table_df = pandas.read_table(str(args.input_table_file))

    stats_df_list = list()
    components_df_list = list()
    component_stats_df_list = list()
    for chain, group in input_table_df.groupby('chain'):
        stats_df, components_df, component_stats_df = add_stats(group)
        components_df['chain'] = chain
        stats_df_list.append(stats_df)
        components_df_list.append(components_df)
        component_stats_df_list.append(component_stats_df)
    pandas.concat(stats_df_list).to_hdf(str(args.stats_hdf_file), 'df', mode='w')
    pandas.concat(components_df_list).to_hdf(str(args.components_hdf_file), 'df', mode='w')
    pandas.concat(component_stats_df_list).to_hdf(str(args.component_stats_hdf_file), 'df', mode='w')

#     (pandas.read_table(str(args.input_table_file))
# #        .groupby(['chrom', 'sample', 'start', 'chain'])
#         .groupby(['chain'])
#         .apply(add_stats)
#         .reset_index(drop=True)
#         .to_hdf(str(args.output_hdf_file), 'df', mode='w')
#     )

# python scripts/tree_statistics.py --pop-size 10000 --min-post-prob 0.2 --max-clade-prob 1e-3 --discretization 0.000000,49.193481,122.586947,232.085215,395.449492,639.178343,1002.805899,1545.314509,2354.701987,3562.255340,5363.846221,8051.702368,12061.808515,18044.625462,26970.598323,40287.567936,60155.618452,89797.454603,134021.141756,200000.000000 steps/argweaver/output/Africa/X-000000000-000100000.tsv.gz stats.hdf comp.hdf



