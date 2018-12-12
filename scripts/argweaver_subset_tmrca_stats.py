
from tree_statistics import tmrca_stats
from ete3 import Tree
import gzip, sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument("--exclpops", type=str, default=None)
parser.add_argument("--exclindivs", type=str, default=None)
parser.add_argument("input_file_name", type=str)
parser.add_argument("output_file_name", type=str)
args = parser.parse_args()


def prune(node, included):
    if node.is_leaf():
        if node.name not in included:
            node.delete(preserve_branch_length=True)
    else:
        for child in list(node.children):
            prune(child , included)

#@profile
def main():
    if args.exclpops is not None:
        excluded_populations = args.exclpops.split(',')
    else:
        excluded_populations = []
    if args.exclindivs is not None:
        excluded_individuals = args.exclindivs.split(',')
    else:
        excluded_individuals = []

    # i = 0

    with gzip.open(args.input_file_name, 'rb') as input_file:
        with gzip.open(args.output_file_name, 'wb') as output_file:
            header = True
            for line in input_file:

                fields = [x.decode() for x in line.split()]

                if header:
                    header = False
                    fields += ['pruned_tmrca', 'pruned_tmrca_half', 'pruned_coal_half']
                    s = '\t'.join(fields) + '\n'
                    output_file.write(s.encode())
                    continue

                tree = Tree(fields[32])#.decode())

                included_leaves = list()
                for leaf in tree.get_leaves():
                    if not any(pop in leaf.name for pop in excluded_populations) \
                        and not any(indiv in leaf.name for indiv in excluded_individuals):
                        # included_leaves.append(leaf)
                        included_leaves.append(leaf.name)

                #tree.prune(included_leaves, preserve_branch_length=True)
                prune(tree, included_leaves)


                # hack to ensure there is no nondicotomic node under the root:
                if len(tree.children) == 1 and not tree.children[0].is_leaf():
                    tree.children[0].delete(preserve_branch_length=True)

                assert set(tree.get_leaf_names()) == set(included_leaves)


                # if not node.is_leaf() and len(node.children) == 1 and not node.children[0].is_leaf():
                #     node.children[0].delete(preserve_branch_length=True)


                tmrca, tmrca_half, coal_half = tmrca_stats(tree)

                fields += [str(tmrca), str(tmrca_half), str(coal_half)]
                s = '\t'.join(fields) + '\n'
                output_file.write(s.encode())

                # i += 1
                # if i > 100:
                    # sys.exit()

if __name__ == "__main__":
    main()
