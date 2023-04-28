import numpy as np
import pandas as pd
from bioviper import msa
from Bio import Phylo
from io import StringIO
import ete3
import os
from copy import deepcopy
from matplotlib.colors import to_rgba, to_rgb, ListedColormap
from matplotlib.cm import get_cmap

from tqdm import tqdm

class Tree:

    def __init__(self, biopython_tree, tempfile=".temp.nwk", rooted=False):

        if isinstance(biopython_tree, str):
            biopython_tree = Phylo.read(StringIO(biopython_tree), format='newick')

        elif isinstance(biopython_tree, Phylo.Newick.Clade):
            self.branch_length = biopython_tree.branch_length
            biopython_tree = Phylo.Newick.Tree(biopython_tree)

        self._biopython = biopython_tree

        # Biopython's tree object is missing some easy way of just storing all of the terminals
        #  the iterator is honestly so annoying
        self.leaves = [term for term in self._biopython.get_terminals()]
        self.branches = [nonterm for nonterm in self._biopython.get_nonterminals()]

        # But OK you can also call them terminals and nonterminals if you want, I won't stop you
        self.nonterminals = self.branches
        self.terminals = self.leaves

        self.leaf_names = np.array([term.name for term in self._biopython.get_terminals()])
        self.N = len(self.leaf_names)
        self.ids = self.leaf_names
        self._tempfile = tempfile
        self.depths = self._biopython.depths

        # For compatibility with biopython workflows
        self.clade = self._biopython.clade

        if np.all([type(leaf.color)==type(None) for leaf in self.leaves]):
            self.is_colored = False
        else:
            self.is_colored = True
            self.colors = [leaf.color for leaf in self.leaves]

        #Phylo.write(self._biopython, self._tempfile, "newick")
        self.ete3 = ete3.Tree(self._biopython.__format__("newick"))

        self.rooted = rooted
        self.root = None

    def __getitem__(self, index):
        
        return Tree(self._biopython.clade[index])
    
    def __iter__(self):
        
        self.first_level_clades = [Tree(clade) for clade in self.clade]
        return iter(self.first_level_clades)
    
    def __repr__(self):
        
        return "bioviper.phylo.Tree with "+str(self.N)+" leaves"
    
    def __str__(self):

        return self._biopython.__format__("newick")

    def _update_ete3_tree(self):

        #Phylo.write(self._biopython, self._tempfile, "newick")
        #self.ete3 = ete3.Tree(self._tempfile)

        self.ete3 = ete3.Tree(self._biopython.__format__("newick"))

    def root_at_midpoint(self):

        self._biopython.root_at_midpoint()
        self.rooted = True
        self._update_ete3_tree()

    def set_outgroup(self, outgroup):

        self._biopython.root_with_outgroup(outgroup)
        self.rooted = True
        self.root = self._biopython.root
        self._update_ete3_tree()

    def get_terminals(self):
        return self.leaves

    def get_leaf_index(self, name):

        x = np.where(self.leaf_names==name)[0]
        if len(x)==1:
            return x[0]
        elif len(x)==0:
            return []
        else:
            return x

    def get_leaf(self, name):

        return self.leaves[self.get_leaf_index(name)]

    def find_clades(self, target=None, terminal=None, order="preorder"):

        return self._biopython.find_clades(target=target, terminal=terminal,
            order=order)

    def search_leaf(self, name):

        if name in self.leaf_names:
            return self.get_leaf(name), self.get_leaf_index(name)

        else:

            leaves = []; leaf_indices = []

            for nl, leaf in enumerate(self.leaves):
                if name in leaf.name:
                    leaves.append(leaf)
                    leaf_indices.append(nl)

            if len(leaves)==1:
                return leaves[0], leaf_indices[0]

            else:
                return leaves, leaf_indices


    def prune(self, terms):

        pruned_ete3_tree = deepcopy(self.ete3)
        pruned_ete3_tree.prune([str(i) for i in terms])

        tr = Tree(Phylo.read(StringIO(pruned_ete3_tree.write("newick")), "newick"))

        if self.is_colored:

            tr.color([self.colors[self.get_leaf_index(i)] for i in tr.leaf_names])

        return tr

    def common_ancestor(self, leaves):

        try:
            return Tree(self._biopython.common_ancestor(leaves))

        except ValueError:
            return Tree(self._biopython.common_ancestor([self.search_leaf(leaf)[0].name for leaf in leaves]))

    def color(self, colors, cmap=None, color_branches=True, in_place=True, df_col="color",
                 branch_coloring_method="average", vmin=None, vmax=None, override_cmap=False):

        '''
        Color the tree based on a given input. There are a variety of things that you can pass for the "colors"
            argument:

            a list of array of colors as strings or
        '''

        self.is_colored = True

        # If you ask for it not in-place, first copy the tree, color that new one in place and return it
        if not in_place:
            NewTree = Tree(deepcopy(self._biopython))
            NewTree.color(colors)
            return NewTree

        else:

            # If you pass a colormap, interpret values as being among a color spectrum
            #   and run self.color_spectrum()
            #   override_cmap is used when color_spectrum() calls color again with the floats turned to RGBA
            if type(cmap) != type(None) and override_cmap==False:
                self.color_spectrum(colors, cmap=cmap, color_branches=color_branches, in_place=in_place, df_col=df_col,
                                   branch_coloring_method=branch_coloring_method, vmin=vmin, vmax=vmax)

            else:

                # If the colors argument is a list or array
                if isinstance(colors, (list, tuple, np.ndarray)):

                    # If it's the same length as the number of leaves in the tree, interpret
                    #  as corresponding to the leaves in the same order that they come in
                    if len(colors) == len(self.ids):

                        # Loop through and color each leaf
                        for leaf,color in zip(self.leaves, colors):
                            ColorLeaf(leaf, color)

                    # If colors is length 4, interpret as a constant RGBA value for the whole tree
                    elif len(colors)==4:
                        for leaf in self.leaves:
                            leaf.color = colors

                    # Otherwise it's uninterpretable, raise an error
                    else:
                        raise InputError("Cannot color tree, length of colors array must match number of terminals in tree")

                # If it's a dictionary, we can treat it similarly but we'll use the keys
                elif isinstance(colors, dict):
                    for key in dict:
                        ColorLeaf(get_leaf(key))

            # Pass a pandas dataframe? Not implemented yet
            #elif isinstance(colors, pd.DataFrame):

                #colors_filtered = colors.loc[self.leaf_names]

                #for key,color in colors_filtered[df_col]:
                    #ColorLeaf(self.get_leaf())

                elif isinstance(colors, (str, np.str_)):
                    for leaf in self.leaves:
                        leaf.color = to_branch_color(to_rgb(colors))

                if color_branches:
                    self.recolor_branches(method=branch_coloring_method)

        self.colors = [leaf.color for leaf in self.leaves]

    def color_spectrum(self, color_arr, cmap='viridis', vmin=None, vmax=None, **kwargs):


        '''Color the tree using numeric values using a colormap. First generates colors from the colormap
            and then runs self.color() with those colors.

            Inputs:
                color_arr: an array of the same length - and same order - as self.leaves
                    containining floats or integers
                    (can also be a dictionary where keys are leaves in the tree)

                cmap: a valid matplotlib colormap, as a string or the actual colormap object
                vmin: the minimum value to set to 0 on the colormap. defaults to the minimum of the data
                vmax: the maximum value to set to 1 on the colormap. defaults to the maximum of the data
                other arguments pass through to self.color(), see above for those'''

        if isinstance(color_arr, (list, tuple)):
            color_arr = np.array(color_arr)

        if isinstance(cmap, (str, np.str_)):
            cmap = get_cmap(cmap)

        if type(vmin)==type(None):
            vmin = np.min(color_arr)

        if type(vmax)==type(None):
            vmax = np.max(color_arr)

        if isinstance(color_arr, dict):

            colors = {}

            for key, val in color_arr.items():
                colors[key] = cmap((val - vmin)/(vmax - vmin))
        else:

            normed_color_arr = (color_arr - vmin)/(vmax-vmin)
            colors = cmap(normed_color_arr)

        return self.color(colors, override_cmap=True, **kwargs)


    def recolor_branches(self, method="average", default_color='black'):

        '''
        Recolor internal branches of the tree based on the colors of the leaves.

        Has two methods:
            method="average" colors based on the average color of all leaves contained within the tree
            method="complete" colors only if every leaf is of that color, otherwise colors the default_color

        default_color is set to black unless you change it, pass a string or an RGBA tuple
        '''

        for nonterm in self.branches:

            x = []
            for term in nonterm.get_terminals():

                try:
                    x.append(from_branch_color(term.color))
                except:
                    None

            if len(x) > 0:

                if method=='average':
                    color = np.nanmean(x, axis=0)

                elif method=='complete':
                    if np.all(x)==x[0]:
                        color = x[0]
                    else:
                        color = to_rgb(default_color)
            else:
                color = to_rgb(default_color)

            term.color = to_branch_color(color)

    def copy(self):

        '''Return an identical copy of the tree.'''

        return Tree(deepcopy(self._biopython))


    def set_ids(self, new_ids, in_place=True):

        '''Change the ids of the tree, modifying both the leaves' .name attribute and the
            .ids and .leaf_names attributes of the tree.

            One optional modifier: in_place decides if it modifies the tree in place or returns a copy.'''

        if in_place:

            if isinstance(new_ids, (list, np.ndarray)):
                if len(new_ids) == self.N:

                    for leaf,name in zip(self.leaves, new_ids):
                        leaf.name = name

                    self.ids = new_ids
                    self.leaf_names = new_ids

            elif isinstance(new_ids, dict):

                for leaf_id, new_name in new_ids.items():
                    leaf_ind = self.get_leaf_index(leaf_id)
                    self.leaves[leaf_ind].name = new_name
                    self.ids[leaf_ind] = new_name

                self.leaf_names = self.ids

            self._update_ete3_tree()

        else:
            newtree = set_ids(self.copy())

    def draw_ascii(self):

        print(self.ete3.get_ascii())

    def get_ascii(self):

        return self.ete3.get_ascii()

    def save(self, filename, fmt='default'):

        if fmt=='default':

            if self.is_colored:
                fmt="phyloXML"

            else:
                fmt="newick"

        if self.is_colored and fmt=="newick":
            print("Warning: newick format does not store color information, so colors will be lost!")

        Phylo.write(self._biopython, filename, fmt)


def ColorLeaf(leaf, color):

    '''Color a leaf of a tree. Automatically handles a variety of formats:
        - if the input is already a BranchColor, just assign it
        - if it's a string, convert to RGB and then to a BranchColor
        - if it's an RGB or RGBA tuple, convert to BranchColor'''

    if isinstance(color, Phylo.BaseTree.BranchColor):
        leaf.color = color

    elif isinstance(color, (str, np.str_)):
        leaf.color = to_branch_color(to_rgba(color))

    elif len(color)==4 or len(color)==3:
        leaf.color = to_branch_color(color)

    else:
        raise InputError("Must pass either a named color or RGBA values as a list")


def RFdistance(Tree1, Tree2, normalized=False):

    '''Calculate the Robinson-Foulds distance between two trees'''

    rf = Tree1.ete3.robinson_foulds(Tree2.ete3)

    if normalized:
        return rf[0] / rf[1]

    else:
        return rf[0]

def to_branch_color(rgb_color):

    '''Converts an RGB color to a Bio.Phylo BranchColor object'''

    return [int(255*i) for i in rgb_color[:3]]

def from_branch_color(branch_color):

    '''Converts a Bio.Phylo BranchColor object back to an RGB value'''

    return np.array([branch_color.red/255, branch_color.green/255, branch_color.blue/255])

def readTree(filename, fmt="detect", **kwargs):

    '''
    Read a Tree object in from a file.
    '''

    # Autodetect format from filename, or else try newick and phyloXML
    if fmt=="detect":
        try:
            suff = filename.split('.')[-1]

            if suff in ("nwk", "newick"):
                fmt = "newick"

            elif suff in ("xml", "phyloxml"):
                fmt = "phyloxml"

        except:
            print("Unable to autodetect format, guessing newick...")

            try:
                return readTree(filename, fmt="newick")
            except:
                print("Trying phyloXML...")
                return readTree(filename, fmt="phyloXML")

    return Tree(Phylo.read(filename, fmt), **kwargs)

def readTrees(filename, fmt="detect", **kwargs):

    '''
    Read a Tree object in from a file.
    '''

    # Autodetect format from filename, or else try newick and phyloXML
    if fmt=="detect":
        try:
            suff = filename.split('.')[-1]

            if suff in ("nwk", "newick"):
                fmt = "newick"

            elif suff in ("xml", "phyloxml"):
                fmt = "phyloxml"

        except:
            print("Unable to autodetect format, guessing newick...")

            try:
                return readTrees(filename, fmt="newick")
            except:
                print("Trying phyloXML...")
                return readTrees(filename, fmt="phyloXML")

    return Forest([Tree(tree) for tree in Phylo.parse(filename, fmt, **kwargs)])

def read_mrbayes_trprobs(filename):

    '''Read in an ensemble of trees (a "forest") from a MrBayes .trprobs output file.'''

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

    for line in tqdm(file.split(';')[2:]):
        x = line[line.find('['):]
        if 'p = ' in x and ',' in x:
            p = float(x.split("&W ")[1].split(']')[0])
            #P = float(x.split("P = ")[1].split('] = [')[0])
            tree = x.split(' ')[-1]+';'
            tree = Tree(tree)
            tree.set_ids(ids)
            trees.append(tree)
            ps.append(p)
            #Ps.append(P)

    return Forest(trees, np.array(ps), names=list(ids.values()))


class Forest:

    def __init__(self, trees, probs=None, names=None, renormalize=True):

        self.trees = trees

        if type(probs)==type(None):
            self.p = np.ones(len(self.trees))
        else:
            self.p = probs

        if renormalize and np.sum(self.p) != 1:
            self.p = self.p / np.sum(self.p)

        #self.P = Probs
        self.names = self.trees[0].ids
        self.ids = self.names; self.term_names = self.names
        self.N = len(self.trees)
        #self.consensus_trees = None

    def __getitem__(self, index):
        return self.trees[index]

    def __iter__(self):
        return iter(self.trees)

    def root_all(self, outgroup):

        for tr in self.trees:
            if 'mid' in outgroup:
                root_point = tr.get_midpoint_outgroup()
                tr.set_outgroup(root_point)
            else:
                tr.set_outgroup(outgroup)

    # def find_consensus(self, method="majority"):

    #     if method=="majority":
    #         return Consensus.strict_consensus(self.trees)



    def prune(self, leaves, outgroup=False, rebin=True):


        """Prune the trees to contain only a subset of leaves."""

        rooted=False
        # Prune trees to include only leaves of interest
        pruned_trees = []

        for tree in tqdm(self.trees):
            tr = tree.prune(leaves)

            # if you gave an outgroup
            if type(outgroup) != type(None):
                rooted=True
                if outgroup=='midpoint':
                    root_point = tr.get_midpoint_outgroup()
                    tr.set_outgroup(root_point)
                else:
                    tr.set_outgroup(outgroup)
            pruned_trees.append(tr)

        if rebin:

            unique_trees, counts, x = bin_trees(pruned_trees, rooted=rooted)

            P = []

            for i in x:
                P.append(np.sum(self.p[np.isin(np.arange(self.N), i)]))

            f = Forest(unique_trees, P, np.arange(len(P)))
            f.counts = counts

            return f

        else:
            return Forest(pruned_trees, p=self.p)


    def calc_marginals(self, leaves, outgroup=None):

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

        pruned_forest = self.prune(leaves, outgroup=outgroup, rebin=True)

        for k,tree in enumerate(pruned_forest):
            print("Topology",k+1,": p = "+str(pruned_forest.p[k])[:8], "with ", pruned_forest.counts[k], "architectures")
            print(tree.get_ascii())
            print('')
            print('')

        return pruned_forest

def rf_distances(tree_1, tree_list):
    return np.array([RFdistance(tree_1, tree_2) for tree_2 in tree_list])

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
