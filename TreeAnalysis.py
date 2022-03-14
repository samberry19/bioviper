"""
TREEANALYSIS.PY
Sam Berry, last edited March 2022

Functions for analyzing posterior distributions of trees as generated from MrBayes
(may someday have support for other kinds of trees). Primarily uses ete3 to analyze
trees, which has significant functionality above and beyond what Bio.Phylo can do,
but trees can also be coverted to and from these two classes in order to do different
kinds of analysis.

(e.g. I use Bio.Phylo for attaching trees to alignments in MultipleSequenceAlignment.py)
"""

import numpy as np
from Bio import Phylo
from Bio.Phylo import Consensus
from ete3 import Tree
from io import StringIO

def bin_trees(forest, rooted=False):

    """Given a set of ETE3 trees, some of which are identical, returns the:
        -- unique trees from the set
        -- counts of each unique tree
        -- indices in the original list corresponding to each unique trees.

    It does so while avoiding calculating the full matrix of robinson-fould distances
    between every pair for efficiency, under the assumption that this matrix will be full of zeros."""

    # calculate the robinson-fould distances of each tree to the first one
    rfdists = rf_distances(forest[0], forest)

    # many of those may be the same; initialize their indices and add them to a growing list
    tree_bins = [np.where(rfdists==0)[0]]

    # these are the indices of trees that are NOT identical, that we now want to consider.
    inds = np.where(rfdists > 0)[0]

    # we will not pay any more attention to the trees matching architecture #1, and only focus on the others now
    new_forest = [forest[i] for i in inds]

    # iterate until we've identified all of the unique architectures
    while len(new_forest) > 0:

        # do basically the same thing
        rfdists = rf_distances(new_forest[0], new_forest)
        tree_bins.append(inds[np.where(rfdists==0)[0]])
        inds = inds[np.where(rfdists>0)[0]]
        new_forest = [forest[i] for i in inds]

    # pick th first tree of each set as its representative - they'll all be identical
    representative_trees = [forest[tree_bins[k][0]] for k in range(len(tree_bins))]
    counts = np.array([len(tree_bins[k]) for k in range(len(tree_bins))])

    return representative_trees, counts, tree_bins


def read_mrbayes_trprobs(filename):

    file = ''

    with open(filename) as f:
        for line in f:
            file = file+line

    ids = {}
    for i in file.split(';')[1].split('\n'):
        l = i.split(' ')
        if len(l) > 4:
            l = [i for i in l if len(i) > 0]
            ids[l[0]] = l[1].replace(',','')

    trees = []; ps = []; Ps = []

    for line in file.split(';')[2:]:
        x = line[line.find('['):]
        if 'p = ' in x and ',' in x:
            p = float(x.split("&W ")[1].split(']')[0])
            P = float(x.split("P = ")[1].split('] = [')[0])
            tree = x.split(' ')[-1]+';'
            tree = Tree(tree)
            for k,leaf in enumerate(tree.get_leaves()):
                leaf.name = ids[leaf.name]
            trees.append(tree)
            ps.append(p)
            Ps.append(P)

    return Ete3Forest(trees, np.array(ps), np.array(Ps), list(ids.values()))

class Ete3Forest:

    def __init__(self, trees, probs, Probs, names):

        self.trees = trees
        self.p = probs
        self.P = Probs
        self.names = names
        self.N = len(self.trees)
        #self.consensus_trees = None

    def __getitem__(self, index):
        return self.trees[index]

    def __iter__(self):
        return __iter__(self.trees)

    def root_all(self, outgroup):

        for tr in self.trees:
            if 'mid' in outgroup:
                root_point = tr.get_midpoint_outgroup()
                tr.set_outgroup(root_point)
            else:
                tr.set_outgroup(outgroup)

    def marginal(self, leaves, outgroup=None, v=True):

        """Calculate the marginal probabilities of tree topologies considering only a subset of leaves.
           Takes ETE3 trees collected into a Forest object, presumably loaded from a MrBayes .trprobs file.
           Relies on trees representing a converged MCMC sample with burnout phase removed such that values
           represent true posterior marginal probabilities.

           Parameters:
            --leaves: the names of the leaves of interest, as a list. should be >3 (>2 if rooted) or else
                       there is only one meaningful topology. all leaves must be exactly as in the tree.
            optional:
                --outgroup: identify one of your leaves as an outgroup to do rooted RF distances. can also
                            pass "midpoint", in which case each tree will be midpoint rooted.
                --v: verbose, true displays results in ascii, false only returns probs and topologies"""

            # will change later if you passed an outgroup
        rooted=False

        # Prune trees to include only leaves of interest
        pruned_trees = []
        for tree in self.trees:
            tr = tree.copy()  # copy to avoid modifying original, and then
            tr.prune(leaves)   # prune - ete3 is very good at this

            # if you gave an outgroup
            if type(outgroup) != type(None):
                rooted=True
                if outgroup=='midpoint':
                    root_point = tr.get_midpoint_outgroup()
                    tr.set_outgroup(root_point)
                else:
                    tr.set_outgroup(outgroup)
            pruned_trees.append(tr)

        unique_trees, counts, x = bin_trees(pruned_trees, rooted=rooted)

        P = []

        for i in x:

            P.append(np.sum(self.p[np.isin(np.arange(self.N), i)]))

        if v:
            for k,tree in enumerate(unique_trees):
                print("Topology",k,": p = "+str(P[k])[:8], "with ", counts[k], "architectures")
                print(tree.get_ascii())
                print('')
                print('')

        return P, unique_trees

    def to_bio(self):
        return BioForest([tree.write() for tree in self.trees], self.p, self.P, self.names)

class BioForest:

    def __init__(self, trees, probs, Probs, names):

        if isinstance(trees[0], Phylo.Newick.Tree):
            self.trees = trees

        elif isinstance(trees[0], str):
            self.trees = [Phylo.read(StringIO(tree), "newick") for tree in trees]

        else:
            raise TypeError("Reading trees from type", trees[0].type, "not yet implemented")

        self.p = probs
        self.P = Probs
        self.names = names
        self.N = len(self.trees)
        #self.consensus_trees = None

    def __getitem__(self, index):
        return self.trees[index]

    def __iter__(self):
        return __iter__(self.trees)

    def find_consensus(self, method, p_thresh=0.95):

        if "maj" in method.lower():
            return Consensus.majority_consensus(trprobs.trees[:np.where(self.P > p_thresh)[0][0]])

        if "strict" in method.lower():
            return Consensus.strict_consensus(trprobs.trees[:np.where(self.P > p_thresh)[0][0]])

    def root_all(self, clade_names):

        for tree in self.trees:
            ca = tree.common_ancestor(next(tree.find_clades(clade_names[0])), next(tree.find_clades(clade_names[1])))
            tree.root_with_outgroup(ca.clades[1])

def rf_distances(tree_1, tree_list):
    return np.array([tree_1.robinson_foulds(tree_2)[0] for tree_2 in tree_list])
