import numpy as np
import pandas as pd
from Bio import AlignIO, Phylo, SeqIO, PDB, pairwise2

from Bio.Align import _aligners
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord,_RestrictedDict
from Bio.Align import MultipleSeqAlignment

from .pdb import readPDB
from .phylo import readTree, Tree

# Otherwise we get some annoying pandas warnings
import warnings
warnings.simplefilter("ignore", UserWarning)

from copy import copy,deepcopy

alphabet = '-ACDEFGHIKLMNPQRSTVWY'

def readAlignment(alignment_file, format="fasta", calc_frequencies=False):

    '''
    Read in an alignment in a format accepted by biopython. Arguments:
        -alignment_file: the path to the alignment you want to load
        -format: 
    '''

    msa = MultipleSequenceAlignment(AlignIO.read(alignment_file, format=format))

    if calc_frequencies==True:
        msa.calc_frequencies()
        msa.calc_coverage()

    return msa

class MultipleSequenceAlignment(MultipleSeqAlignment):

    """Represents a classical multiple sequence alignment (MSA).

    By this we mean a collection of sequences (usually shown as rows) which
    are all the same length (usually with gap characters for insertions or
    padding). The data can then be regarded as a matrix of letters, with well
    defined columns."""

    def __init__(
        self, records, ids=None, names=None, descriptions=None,
        alphabet=None, annotations=None, column_annotations=None,
        tree=None, _tree_branches=None
    ):
        """Initialize a new MultipleSeqAlignment object.

        Arguments:
         - records - A list (or iterator) of SeqRecord objects, whose
                     sequences are all the same length.  This may be an be an
                     empty list.
         - alphabet - For backward compatibility only; its value should always
                      be None.
         - annotations - Information about the whole alignment (dictionary).
         - column_annotations - Per column annotation (restricted dictionary).
                      This holds Python sequences (lists, strings, tuples)
                      whose length matches the number of columns. A typical
                      use would be a secondary structure consensus string."""

        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")

        # IF you generate from a matrix rather than a list of SeqRecords
        if isinstance(records, np.ndarray):
            self._records = []
            self.matrix = records

            record_list = []

            for k,rec in enumerate(records):

                record_info = []

                seq = SeqRecord(Seq(''.join(rec)))

                if type(ids) != type(None):
                    seq.id = ids[k]
                else:
                    seq.id = None

                if type(names) != type(None):
                    seq.name = str(names[k])

                if type(descriptions) != type(None):
                    seq.description = str(descriptions[k])

                elif type(_tree_branches) != type(None):
                    seq.branch = _tree_branches[k]

                record_list.append(seq)

            super().__init__(record_list)

        else:
            '''If you pass it a list of record objects, it will initialize it as a
            normal MultipleSeqAlignment object.'''

            super().__init__(records)

        # Annotations about the whole alignment
        if annotations is None:
            annotations = {}
        elif not isinstance(annotations, dict):
            raise TypeError("annotations argument should be a dict")
        self.annotations = annotations

        # Annotations about each column of the alignment
        if column_annotations is None:
            column_annotations = {}
        # Handle this via the property set function which will validate it
        self.column_annotations = column_annotations

        #Useful for indexing
        id_to_index = pd.DataFrame(zip(self.ids, np.arange(self.N)), columns=("id", "index")).set_index("id")
        self.id_to_index = id_to_index

        self.N, self.L = self.matrix.shape

        self.coverage = None
        self.frequencies = None
        self.tree = None
        self.numeric_matrix = None
        self.fij = None
        self._pseudocount = 0
        self.sequence_weights = np.ones(self.N)

    def save(self, filename, format='Auto'):

        if format=='Auto':
            try:
                ending = filename.split('.')[-1]
                if "fa" in ending:
                    format = "fasta"
                elif "sto" in ending:
                    format = "stockholm"
                elif "clust" in ending:
                    format = "clustal"
            except:
                return TypeError("Please pass a format or a filename with an ending containing fa, sto or clust")

        AlignIO.write(MultipleSeqAlignment(self._records), filename, format)

    def update(self):
        if self.matrix.shape[0] != len(self._records):
            self.matrix = np.array(self.records)
            self.N, self.L = self.matrix.shape
            self.ids = np.array([rec.id for rec in self._records])
            self.names = np.array([rec.name for rec in self._records])
            self.descriptions = [rec.description for rec in self._records]

    def dealign(self):

        """Return raw sequences without any gap characters (or lowercase)"""

        ds = []

        for rec in self._records:
            s = Seq(str(rec.seq).replace('-','').replace('.','').upper())
            ds.append(SeqRecord(s, id = rec.id, name=rec.name, description=rec.description))

        return SequenceArray(ds)


    def sequence_lengths(self):

        W = (self.matrix!='-')&(self.matrix!='.')
        return np.sum(W,axis=1)

    def calc_coverage(self):
        if type(self.coverage)==type(None):
            self.coverage = np.mean((self.matrix!='.')&(self.matrix!='-'),axis=0)
            return self.coverage
        else:
            return self.coverage

    def calc_frequencies(self, pseudocount=0, first_pos=0, pos_indices = None, store=True):

        self._pseudocount = pseudocount

        M_upper = self.upper().replace('.','-')
        if pseudocount==0:
            frequencies = np.array([ \
                            [np.average(M_upper.matrix[:,k]==i, weights=self.sequence_weights) for i in alphabet] \
                          for k in range(M_upper.matrix.shape[1])])

        else:
            frequencies = np.array([ \
                             [np.average(M_upper.matrix[:,k]==i, weights=self.sequence_weights) + \
                                pseudocount/np.sum(self.sequence_weights) for i in alphabet] \
                            for k in range(M_upper.matrix.shape[1])])
            frequencies = frequencies / np.sum(frequencies, axis=1)[:,None]

        frequencies = pd.DataFrame(frequencies, columns=list(alphabet))

        frequencies.sequence_weights = self.sequence_weights

        if type(pos_indices) != type(None):
            frequencies.index = pos_indices
        elif first_pos != 0:
            frequencies.index = np.arange(first_index, self.L + first_index)

        if store:
            self.frequencies = frequencies

        return frequencies

    def as_numeric(self, recalculate=False):

        if type(self.numeric_matrix) == type(None) or recalculate:
            m = copy(self.upper().replace('.','-').matrix)
            mnum = np.zeros(m.shape)
            for n,letter in enumerate(alphabet):
                mnum[m==letter]=n
            self.numeric_matrix = mnum.astype('int')

        return self.numeric_matrix

    def as_df(self, numeric=False):

        '''Return alignment as a pandas dataframe.'''

        if numeric:
            df = pd.DataFrame.from_records(self.as_numeric(), index=self.ids)
        else:
            df = pd.DataFrame.from_records(self.matrix, index=self.ids)

        return df

    def one_hot_encoding(self, full=False, flat=True):

        '''Generate a one-hot encoded pandas dataframe, e.g. for dimensionality reduction. Two arguments:

            full: whether to use all amino acid possibilities (TRUE) or only those seen in the alignment (FALSE)
            flat: whether to flatten the one hot encoding to 2D (shape L x 21*N) or keep it 3D (final array L x N x 21)
                (flat=False is only an option for full=True, otherwise it won't be a full matrix)

            If full=False, returns a pandas dataframe
            If full=True, returns a numpy array.'''

        if full:

            # Initialize an array of all zeros of the correct shape
            x = np.zeros((self.N, self.L, 21))

            for n, seq_num in enumerate(self.as_numeric()):
                x[n, np.arange(self.L), seq_num] = 1

                if flat:
                    return x.reshape((self.N, self.L*21))
                else:
                    return x

        else:

            return pd.get_dummies(self.upper().as_df())

    def calc_pair_frequencies(self, pseudocount=0):

        if type(self.numeric_matrix) == type(None):
            self.as_numeric()

        X = [[[] for i in range(self.L)] for i in range(self.L)]

        for i in range(self.L):
            for j in range(i, self.L):
                A = np.histogram2d(self.numeric_matrix[:,i], self.numeric_matrix[:,j], bins=[np.arange(22), np.arange(22)])[0] + pseudocount
                A = A/np.sum(A)
                X[i][j] = A
                X[j][i] = A

        self.fij = np.array(X)

        return self.fij

    def _set_per_column_annotations(self, value):
        if not isinstance(value, dict):
            raise TypeError(
                "The per-column-annotations should be a (restricted) dictionary."
            )
        # Turn this into a restricted-dictionary (and check the entries)
        if len(self):
            # Use the standard method to get the length
            expected_length = self.get_alignment_length()
            self._per_col_annotations = _RestrictedDict(length=expected_length)
            self._per_col_annotations.update(value)
        else:
            # Bit of a problem case... number of columns is undefined
            self._per_col_annotations = None
            if value:
                raise ValueError(
                    "Can't set per-column-annotations without an alignment"
                )

    def _get_per_column_annotations(self):
        # if self._per_col_annotations is None:
        #     # This happens if empty at initialisation
        #     if len(self):
        #         # Use the standard method to get the length
        #         expected_length = self.get_alignment_length()
        #     else:
        #         # Should this raise an exception? Compare SeqRecord behaviour...
        #         expected_length = 0
        #     self._per_col_annotations = _RestrictedDict(length=expected_length)
        #return self._per_col_annotations

        return None

    column_annotations = property(
        fget=_get_per_column_annotations,
        fset=_set_per_column_annotations,
        doc="""Dictionary of per-letter-annotation for the sequence.""",
    )

    def extend(self, records):
        """Add more SeqRecord objects to the alignment as rows.

        They must all have the same length as the original alignment. For
        example,

        Because the alignment object allows iteration over the rows as
        SeqRecords, you can use the extend method with a second alignment
        (provided its sequences have the same length as the original alignment).
        """
        if len(self):
            # Use the standard method to get the length
            expected_length = self.get_alignment_length()
        else:
            # Take the first record's length
            records = iter(records)  # records arg could be list or iterator
            try:
                rec = next(records)
            except StopIteration:
                # Special case, no records
                return
            expected_length = len(rec)
            self._append(rec, expected_length)
            # Can now setup the per-column-annotations as well, set to None
            # while missing the length:
            self.column_annotations = {}
            # Now continue to the rest of the records as usual

        for rec in records:
            self._append(rec, expected_length)

        try:
            self.matrix = np.concatenate([self.matrix, MultipleSequenceAlignment(records).matrix])
        except:
            self.matrix = np.array(self._records)

        self.ids = np.array([rec.id for rec in self._records])
        self.names = np.array([rec.name for rec in self._records])
        self.descriptions = np.array([rec.description for rec in self._records])

    def upper(self):

        """Make all letters in an alignment uppercase."""

        return MultipleSequenceAlignment(np.char.upper(self.matrix), ids=self.ids, names=self.names, descriptions=self.descriptions)


    def replace(self, a, b):

        """Replace one letter in an alignment with another."""

        return MultipleSequenceAlignment(np.char.replace(self.matrix, a, b), ids=self.ids, names=self.names, descriptions=self.descriptions)

    def standardize(self):

        ''' Make all amino acids upper case and replace all noncanonical characters with gaps ('-') '''

        matrix_upper = np.char.upper(self.matrix)
        noncanonical_characters = [letter for letter in np.unique(matrix_upper) if letter not in alphabet]

        for character in noncanonical_characters:
            matrix_upper = np.char.replace(matrix_upper, character, '-')

        return MultipleSequenceAlignment(matrix_upper, ids=self.ids,
                                names=self.names, descriptions=self.descriptions)

    def set_ids(self, ids):

        """Reset all IDs in place to a new list. Will both change self.ids and the .id attribute
           of each individual sequence record."""

        if len(ids) == self.N:
            self.ids = ids
            for n,seq in enumerate(self._records):
                seq.id = ids[n]
        else:
            raise ValueError('ERROR: cannot use id list of length', len(ids), 'for an alignment with', self.N, 'sequences')


    def modify_ids(melf, func):

        """Modify all IDs in place with a function. Either pre-define the function or use lambda, e.g.
            alignment.modify_ids(lambda x:x.split('|')[2])"""

        ids = []

        for seq in self._records:
            seq.id = func(seq.id)
            ids.append(seq.id)

        self.ids = np.array(ids)

    def pos_covariance(self):

        '''Get the 4D position by position covariance matrix, considering one-hot encoded sequences'''

        return np.cov(self.one_hot_encoding(full=True).T).reshape((sampler.alignment.L, sampler.alignment.L, 21, 21))

    def seq_covariance(self):

        '''Get the covariance matrix across sequences'''

        return np.cov(np.array(self.one_hot_encoding()))

    def __getitem__(self, index):
        """Access part of the alignment. Indexes like a 2D numpy array and returns
        another MultipleSequenceAlignment object unless:
            1) you ask for a single row, when you get a SeqRecord object
            2) you ask for a single column, when you get a numpy array
            3) you ask for a single row AND a single column, when you get a single letter as a string."""

        if isinstance(index, int):
            return self._records[index]

        elif isinstance(index, slice):
            return self.__getitem__(np.arange(self.L)[index])

        elif isinstance(index, np.ndarray):

            if isinstance(index[0], (int, np.int64, float, np.float64)):

                ids = self.ids[index]; names = self.names[index]
                descs = list(np.array(self.descriptions)[index])

                if type(self.tree) != type(None):

                    if type(index[0])==np.bool_:
                        branches = [self._records[k] for k in np.where(index)[0]]
                    else:
                        branches = [self._records[k] for k in index]

                    new_msa = MultipleSequenceAlignment(self.matrix[index], ids=ids,
                        names=names, descriptions=descs,_tree_branches=branches)
                    new_msa.tree = self.tree

                    return new_msa

                else:
                    return MultipleSequenceAlignment(self.matrix[index], ids=ids, names=names, descriptions=descs)

            else:

                return self.__getitem__(np.array(self.id_to_index.loc[index]))

        elif isinstance(index, list):
            return self.__getitem__(np.array(id))

        elif isinstance(index, str):
            return self.__getitem__(np.where(self.ids==index)[0])[0]

        elif len(index) != 2:
             raise TypeError("Invalid index type.")

        else:
            # Handle double indexing
            row_index, col_index = index

            # If you ask for a single row...
            if isinstance(row_index, int):

                # and a single, column, it will return a string (same as indexing self.matrix)
                if isinstance(col_index, int):
                    return self.matrix[row_index, col_index]

                # and a slice, return a single record only at that particular set of positions
                elif isinstance(col_index, (slice, np.ndarray)):
                    return self.__getitem__((np.arange(self.N), col_index))._records[row_index]

            # If you pass a single column, it will return an array, not an MSA object with a single column because that's dumb
            #   here is one difference from a MultipleSeqAlignment in Biopython, where a string is returned
            #   it seems to me that an array is a much more logical data structure for this

            elif isinstance(col_index, int):
                return self.matrix[:,col_index]

            # When you pass complicated indexing
            else:

                ids = self.ids[row_index]
                names = self.names[row_index]
                descs = list(np.array(self.descriptions)[row_index])
                #new_msa = MultipleSequenceAlignment(self.matrix[row_index, col_index], ids=ids)

                # if there's an attached tree, copy over the tree without re-attaching
                if type(self.tree) != type(None):
                    new_msa = MultipleSequenceAlignment(self.matrix[row_index, col_index], ids=ids,
                        names=names, descriptions=descs,
                        _tree_branches=[rec.branch for rec in self._records])
                    new_msa.tree = self.tree

                    return new_msa

                else:
                    return MultipleSequenceAlignment(self.matrix[row_index, col_index],
                        ids=ids, names=names, descriptions=descs)


    def search_id(self, id_str):

        """Search for an incomplete ID in the alignment"""

        in_id = np.array([str(id_str) in id_val for id_val in self.ids])
        matched = np.where(in_id)[0].astype('int')

        if len(matched) == 0:
            return None

        elif len(matched) == 1:
            return self.__getitem__(int(matched[0]))

        else:
            return self.__getitem__(matched)

    def search_ids(self, ids):

        seqs = []

        for id_val in ids:
            matched = self.search_id(id_val)

            if not isinstance(matched, type(None)):

                if isinstance(matched, MultipleSequenceAlignment):
                    for match in matched:
                        seqs.append(match)

                elif isinstance(matched, SeqRecord):
                    seqs.append(matched)


        return MultipleSequenceAlignment(seqs)


    def search_sequence(self, sequence, include_gaps=False, case_sensitive=False):


        """Search for a particular piece of sequence in the alignment. Only matches the sequence exactly."""

        x = []
        for nr,record in enumerate(self._records):
            if include_gaps==True:
                if case_sensitive:
                    seq = record.seq
                else:
                    seq = record.seq.upper()
                    sequence=sequence.upper()

                if sequence in seq:
                    x.append(nr, seq.index(sequence))

            else:
                if case_sensitive:
                    seq = record.seq.replace('-','').replace('.').upper()
                else:
                    seq = record.seq.replace('-','').replace('.','').upper()
                    sequence=sequence.upper()

                if sequence in seq:
                    where_ungapped = np.where((self.matrix[nr]!='-')&(self.matrix[nr]!='.'))[0]
                    x.append(nr, where_ungapped[seq.index(sequence)])

        if len(x)==0:
            return None

        if len(x)==1:
            return x[0]

        else:
            return x

    ### TREE FUNCTIONALITY ##

    def attach_tree(self, tree, format='newick', prune_unmatched=False):

        """Attach a precalculated tree to the alignment. See documentation for how usage.
            As of now, IDs in the tree must exactly match IDs in the alignment."""

        if type(tree)==str:
            tree = readTree(tree, format)

        if isinstance(tree, (Phylo.Newick.Tree, Phylo.PhyloXML.Phylogeny)):
            tree = Tree(tree)

        self.tree = deepcopy(tree)

        for leaf in self.tree.leaves:

            try:
                n = list(self.ids).index(leaf.name)

            except:

                ns = [n for n,i in enumerate(self.ids) if leaf.name in i]
                if len(ns) == 0:
                    print("Terminal name", leaf.name, "not in alignment - will not be attached!")
                    n = -1

                elif len(ns) == 1:
                    n = n[0]

                else:
                    print("Ambiguity with terminal time", leaf.name, "in alignment", len(ns), "times!")
                    n = -1

            if n > -1:
                self._records[n].branch = leaf
                leaf.seq = self._records[n].seq
                leaf.nseq = n

            elif prune_unmatched:
                self.tree.prune(leaf)

    def set_sequence_weights(self, weights):

        self.sequence_weights = np.array(weights)
        for n,record in enumerate(self._records):
            record.weight = weights[n]

    def get_index(self, id):

        """Return the index of a particular sequence id. Should work for partial IDs;
            for indices that appear multiple times it will return an array, and for missing
            indices it will return None."""


        try:
            # pandas indexing here is by far the fastest way to do this
            return self.id_to_index.loc[id].iloc[0]

        except:

            a = self.search_id(id)
            if type(a)==type(None):
                return None

            elif type(a) == SeqRecord:
                return self.get_index(a.id)

            else:
                return np.array([self.get_index(seq.id) for seq in a])

    def fix_tree(self):

        """This updates the current attached tree to (1) make sure that the number matching is correct
        and (2) remove any nodes that are not present in the alignment.

        Technically, it makes a copy of the original tree and re-attaches it before modifying.

        You'll need to do this if you subset particular sequences in an alignment with an attached tree
        and now want to use the new tree. I don't do it by default because for very large alignments
        I suspect it may take a little while and maybe you don't care."""

        if type(self.tree) != type(None):
            new_tree = deepcopy(self.tree)
            self.attach_tree(new_tree)
            branches = [record.branch for record in self._records]

            for term in self.tree.get_terminals():
                if tree not in branches:
                    self.tree.prune(branch)


    def subset_by_ids(self, ids, sort=False, match_order=True):

        """
        Subset the MSA using a set of IDs. Has several options related to the order they are returned:

            match_order=TRUE (default): returns an MSA where sequences are ordered according to the list you passed.
                This is the slowest option. If the IDs in the passed array are only partial, e.g. don't fully match
                the ids in the MSA, this will still work, but it will be even slower. At the moment, it does NOT work
                if some of the ids are not present in the alignment - if you think this is the case, I would first use
                np.intersect1d to get only the overlapping sequence IDs.

            match_order=FALSE, sort = FALSE: preserve the original order of the MSA, only excluding sequences that
                    aren't in IDs. This is an intermediate speed. Only works for complete ids.

            sort = TRUE (overrides match_order): first sort each list and then intersect. This is much faster than either
                    of the other options for very large alignments (e.g. 100,000+ sequences), but returns sequences alphabetically
                    rather than in any of the original orders. Only works for complete ids.
        """


        self.__getitem__[np.array(self.id_to_index.loc[ids]["index"])]


    def subset_by_clade(self, clade):

        """For a clade in the attached tree, subset the alignment to get only those sequences."""

        try:
            return self.__getitem__(np.array([term.nseq for term in clade.get_terminals()]))

        # This will happen if some of the terminals don't have a .nseq attribute
        except AttributeError:

            # Define the terminals and their names
            terminals = clade.get_terminals()
            terminal_names = [term.name for term in terminals]

            # Find which ones are actually in
            isin = np.isin(terminal_names, self.ids)

            # If any were found, return them
            if np.sum(isin) > 0:
                return self.__getitem__(np.array([terminals[k].nseq for k in np.where(isin)[0]]))

            # Else nothing is returned
            else:
                return None

    def build_tree(self, method='nj', model='autodetect', root=True, ladderize=True, filename="fasttree.fa"):

        "Build a tree and attach it. Calls functions from the TreeBuilders.py script"""

        if method=='nj':
            self.attach_tree(NJTree(MultipleSeqAlignment(self.upper()._records), matrix=model, root=root, ladderize=ladderize))

        elif method=='upgma':
            self.attach_tree(UPGMATree(self, matrix=model, root=root, ladderize=ladderize))

        elif method=='FastTree':
            if model=='autodetect':
                model = "LG"
            format_id = filename.split(".")[-1]
            if "fa" in format_id:
                format = "fasta"
            elif "sto" in format_id:
                format = "stockholm"
            elif "clust" in format_id:
                format = "clustal"
            else:
                print("Cannot interpret file format for", filename, "defaulting to fasta")
                format = "fasta"

            self.attach_tree(FastTree(self, model=model, root=root, ladderize=ladderize, cat=cat,
                                temp_aln=filename))

    def attach_structure(self, structure, nseq=None, name=None, chains='A', model=0):

        if isinstance(structure, str):
            if type(name) == type(None):
                if type(nseq) != type(None):
                    name = self.ids[nseq]
                else:
                    name = ''

            structure = readPDB(structure, name, chains, model)

        if type(nseq) == None:

            x = self.search_sequence(structure.sequence)

            if isinstance(x, NoneType):
                print('Could not find sequence')
                return None

            elif isinstance(x, int):
                nseq = x

            elif isinstance(x, list):
                print(len(x), 'matching sequences found! Appending to first...')
                nseq = x[0]

        elif type(nseq) == str:

            a = np.where(self.ids==nseq)[0]
            if len(a)==0:
                a = self.get_index(nseq)
                if type(a) == type(None):
                    print("string", nseq, "is not located anywhere in alignment IDs!")
            else:
                a = a[0]

            nseq = a

        record = self._records[nseq]

        # Align the structure's sequence to the record's sequence
        #  replace gaps in the record's sequence with "x" to not confuse these
        #  gaps (which we want to remember) with possible new gaps introduced
        #  in the alignment
        A = PairwiseAlign(str(record.seq).replace('-','x'), structure.ordered_sequence,
                            alignment_method='globalms', alignment_params=(2,-1,-0.75,-.1))

        A1 = A[:,np.where(A.matrix[0]!='-')[0]]
        A2 = A[:,np.where(A.matrix[1]!='-')[0]]

        structure._alipos = np.where(A1.matrix[1]!='-')[0]
        structure._pdbpos = np.where(A2.matrix[0]!='-')[0]

        try:
            record.structures
        except:
            record.structures = {}

        if type(name) != None:

            if 0 in record.structures:
                record.structures[np.max([key for key in record.structures.keys() if type(key)==int])+1] = structure
            else:
                record.structures[0] = structure

        else:
            record.structures[name] = structure

        try:
            self.structures.append((nseq, structure))
        except:
            self.structures = [(nseq, structure)]

def readSequences(filename, format="Detect"):

    if format=='Detect':
        ending = filename.split('.')[-1]
        if 'fastq' in ending:
            format="fastq"
        elif 'fa' in ending:
            format="fasta"
        elif 'sto' in ending:
            format = "stockholm"
        elif 'nex' in ending:
            format = "nexus"
        elif "a2m" in ending:
            format = "fasta"

    return SequenceArray(list(SeqIO.parse(filename, format)), fmt=format)

class SequenceArray:

    def __init__(self, records, fmt="fasta"):

        self._records = records
        self.ids = np.array([record.id for record in self._records])
        self.names = np.array([record.id for record in self._records])
        self.descriptions = [record.id for record in self._records]
        self.N = len(self._records)
        self.L = np.array([len(rec.seq) for rec in self._records])
        self.format = fmt

        # For fastq sequencing reads with phred quality scores, save those phred quality scores in a centralized place
        if fmt == "fastq":
            self.phred_quality = [np.array(rec.letter_annotations["phred_quality"]) for rec in self._records]
            self.mean_phred_quality = np.array([np.mean(qual_arr) for qual_arr in self.phred_quality])

    def __len__(self):
        return self.N

    def __getitem__(self, index):

        if isinstance(index, int):
            return self._records[index]

        elif isinstance(index, slice):
            return SequenceArray(self._records[index])

        elif isinstance(index, (list, np.ndarray)):

            if type(index[0])==np.bool_:
                return SequenceArray([self._records[k] for k in np.where(index)[0]])

            else:
                return SequenceArray([self._records[k] for k in index])

        elif isinstance(index, str):

            if index in self.ids:
                return self._records[self.ids.index(index)]
            else:
                raise IndexError("No such ID in SequenceArray. If ID is incomplete, try seqs.search_id(name)")

    def from_file(self, filename, format="Detect"):

        if format=='Detect':
            ending = filename.split('.')[-1]
            if 'fa' in ending:
                format="fasta"
            elif 'sto' in ending:
                format = "stockholm"
            elif 'nex' in ending:
                format = "nexus"
            elif "a2m" in ending:
                format = "fasta"

        self._records = list(SeqIO.parse(filename, format))
        self.ids = np.array([record.id for record in self._records])
        self.names = np.array([record.id for record in self._records])
        self.descriptions = [record.id for record in self._records]

    def search_id(self, id_str):

        """Search for an incomplete ID in the alignment"""

        in_id = np.array([str(id_str) in id_val for id_val in self.ids])

        if np.sum(in_id) > 1:
            return self.__getitem__(np.where(in_id)[0].astype('int'))
        elif np.sum(in_id)==1:
            return self.__getitem__(np.where(in_id)[0].astype('int'))[0]
        else:
            return []

    def __iter__(self):
        return __iter__(self._records)

    def search_ids(self, ids):

        return SequenceArray([self.search_id(id_val) for id_val in ids])

    def search_sequence(self, sequence):

        x = []
        for nr,record in enumerate(self._records):

            if sequence in record.seq:
                x.append(nr, seq.index(sequence))

        if len(x)==0:
            return None

        if len(x)==1:
            return x[0]

        else:
            return x

    def attach_structure(self, structure, nseq=None, name=None, chains='A', model=0):

        if isinstance(structure, str):
            structure = readPDB(filename, name, chains, model)

        if type(nseq) == None:

            x = self.search_sequence(structure.sequence)

            if isinstance(x, NoneType):
                print('Could not find sequence')
                return None

            elif isinstance(x, int):
                nseq = x

            elif isinstance(x, list):
                print(len(x), 'matching sequences found! Appending to first...')
                nseq = x[0]

        elif type(nseq) == str:

            a = np.where(self.ids==nseq)[0]
            if len(a)==0:
                a = self.search_id(nseq)[0]
                if len(a)==0:
                    print("string", nseq, "is not located anywhere in alignment IDs!")
                else:
                    a = a[0]
            else:
                a = a[0]

            nseq = a

        record = self._records[nseq]

        # Align the structure's sequence to the record's sequence
        #  replace gaps in the record's sequence with "x" to not confuse these
        #  gaps (which we want to remember) with possible new gaps introduced
        #  in the alignment
        A = PairwiseAlign(str(record.seq).replace('-','x'), structure.ordered_sequence,
                            alignment_method='globalms', alignment_params=(2,-1,-0.75,-.1))

        A1 = A[:,np.where(A.matrix[0]!='-')[0]]
        A2 = A[:,np.where(A.matrix[1]!='-')[0]]

        structure._seqpos = np.where(A1.matrix[1]!='-')[0]
        structure._pdbpos = np.where(A2.matrix[1]!='-')[0]

        try:
            record.structures
        except:
            record.structures = {}

        if type(name) != None:

            if 0 in record.structures:
                record.structures[np.max([key for key in record.structures.keys() if type(key)==int])+1] = structure
            else:
                record.structures[0] = structure

        else:
            record.structures[name] = structure

        try:
            self.structures.append((nseq, structure))
        except:
            self.structures = [(nseq, structure)]

    def save(self, filename, format='Auto'):

        if format=='Auto':
            try:
                ending = filename.split('.')[-1]
                if "fa" in ending:
                    format = "fasta"
                elif "sto" in ending:
                    format = "stockholm"
                elif "clust" in ending:
                    format = "clustal"
            except:
                return TypeError("Please pass a format or a filename with an ending containing fa, sto or clust")

        SeqIO.write(self._records, filename, format)


def PairwiseAlign(sequence1, sequence2, alignment_method, alignment_params):

    '''
    Pairwise alignment of two sequences using simple alignment tools from pairwise2.align. A standalone function
    to make it easier to incorporate choices of different methods or parameters into a more complex workflow
    '''

    if alignment_method == 'globalms':
        alignment = pairwise2.align.globalms(sequence1, sequence2, *alignment_params)[0]
    elif alignment_method == 'globalxx':
        alignment = pairwise2.align.globalxx(sequence1, sequence2, *alignment_params)[0]
    elif alignment_method == 'globalmx':
        alignment = pairwise2.align.globalmx(sequence1, sequence2, *alignment_params)[0]
    elif alignment_method == 'localxx':
        alignment = pairwise2.align.localxx(sequence1, sequence2, *alignment_params)[0]
    elif alignment_method == 'muscle':
        alignment = MuscleAlign(sequence1, sequence2, muscle_path)
    else:
        print('Not a valid option for alignment! Currently supported options are globalxx, globalmx, globalms, and localxx')

    ali = np.array([np.array(list(alignment[0])), np.array(list(alignment[1]))])
    return MultipleSequenceAlignment(ali)


def cluster_sizes(ali, thresh):

    try:
        imat = ali.identity_matrix
    except:
        imat = CalcIdentityMatrix(ali)

    return np.array([np.sum(imat[i] > thresh) for i in range(ali.N)])

def CalcIdentityMatrix(msa, gap_handling=1):

    '''Code to calculate a percent identity matrix between every pair of sequences
    in a multiple sequence alignment. Uses numpy array broadcasting and so is much faster
    than Biopython's built-in method.

    Comes with three options for handling gap characters when calculating identity:
        gap_handling = 0 ("keep_all") :
                            treats gaps as a character, so two gaps are "equal"

        gap_handling = 1 ("ignore_doubles") :
                            doesn't include positions where both sequences are gapped but
                            if only one is gapped this counts as a difference
                            (DEFAULT, this one makes the most sense to me)

        gap_handling = 2 ("ignore_all"):
                            ignore positions that are gapped in either sequence
    '''

    # Initialize the identity matrix
    IM = []
    for nseq, seq in enumerate(msa.matrix):

        # If we're treating gaps as a character, this is easy
        #  loop through sequences first rather than doing it all at once with numpy to lower memory usage
        if str(gap_handling) in ('0', 'keep_all', 'False'):

            # Take the mean identity of the sequence with every other sequence
            IM.append(np.mean(seq==msa.matrix, axis=1))

        # We need to do a little bit more to ignore gaps
        if str(gap_handling) in ('1', 'ignore_doubles'):

            # We're going to loop through sequences again
            ids = []

            # Don't need to precalculate the ._coverage_mask if it's already defined for that MSA
            try:
                msa._coverage_mask[nseq]
            except:
                msa._coverage_mask = np.array([msa.matrix[n]!='-' for n in range(msa.N)])

            for n in range(msa.N):

                # If either is TRUE, return TRUE
                use_pos = np.add(msa._coverage_mask[n], msa._coverage_mask[nseq])

                #
                ids.append(np.mean(seq[use_pos]==msa.matrix[n][use_pos]))

            IM.append(ids)

        if str(gap_handling) in ('2', 'ignore_all'):

            ids = []

            try:
                msa._coverage_mask[nseq]
            except:
                msa._coverage_mask = np.array([msa.matrix[n]!='-' for n in range(msa.N)])

            for n in range(msa.N):

                # If either is FALSE, return FALSE
                use_pos = msa._coverage_mask[n] & msa._coverage_mask[nseq]

                ids.append(np.mean(seq[use_pos]==msa.matrix[n][use_pos]))

            IM.append(ids)

    return np.array(IM)
