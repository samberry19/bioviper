
import numpy as np
from Bio import AlignIO, Phylo, PDB, pairwise2

from Bio.Align import _aligners
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord,_RestrictedDict
from Bio.PDB.DSSP import DSSP
from Bio.Align import MultipleSeqAlignment

from TreeBuilders import *

from copy import copy,deepcopy

alphabet = '-ACDEFHIKLMNPQRSTVWXY'

def readAlignment(alignment_file, format="fasta", calc_frequencies=True):

    msa = MultipleSequenceAlignment(AlignIO.read(alignment_file, format=format))

    if calc_frequencies==True:
        msa.calc_frequencies()
        msa.calc_coverage()

    return msa


class MultipleSequenceAlignment:

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

        self._records = []

        # IF you generate from a matrix rather than a list of SeqRecords
        if isinstance(records, np.ndarray):
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

            self.extend(record_list)

        else:
            #if records:
            self.extend(records)

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

        self.N, self.L = self.matrix.shape

        self.coverage = None
        self.frequencies = None
        self.tree = None

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

        return ds


    def sequence_lengths(self):

        W = (self.matrix!='-')&(self.matrix!='.')
        return np.sum(W,axis=1)

    def calc_coverage(self):
        if type(self.coverage)==type(None):
            self.coverage = np.mean((self.matrix!='.')&(self.matrix!='-'),axis=0)
            return self.coverage
        else:
            return self.coverage

    def calc_frequencies(self):
        if type(self.frequencies)==type(None):
            M_upper = self.upper().replace('.','-')
            self.frequencies = np.array([[np.sum(M_upper.matrix[:,k]==i) for i in alphabet] for k in range(M_upper.matrix.shape[1])]) / self.N
            return self.frequencies
        else:
            return self.frequencies

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
        if self._per_col_annotations is None:
            # This happens if empty at initialisation
            if len(self):
                # Use the standard method to get the length
                expected_length = self.get_alignment_length()
            else:
                # Should this raise an exception? Compare SeqRecord behaviour...
                expected_length = 0
            self._per_col_annotations = _RestrictedDict(length=expected_length)
        return self._per_col_annotations

    column_annotations = property(
        fget=_get_per_column_annotations,
        fset=_set_per_column_annotations,
        doc="""Dictionary of per-letter-annotation for the sequence.""",
    )

    def _str_line(self, record, length=50):
        """Return a truncated string representation of a SeqRecord (PRIVATE).

        This is a PRIVATE function used by the __str__ method.
        """
        if record.seq.__class__.__name__ == "CodonSeq":
            if len(record.seq) <= length:
                return "%s %s" % (record.seq, record.id)
            else:
                return "%s...%s %s" % (
                    record.seq[: length - 3],
                    record.seq[-3:],
                    record.id,
                )
        else:
            if len(record.seq) <= length:
                return "%s %s" % (record.seq, record.id)
            else:
                return "%s...%s %s" % (
                    record.seq[: length - 6],
                    record.seq[-3:],
                    record.id,
                )

    def __str__(self):
        """Return a multi-line string summary of the alignment.

        See also the alignment's format method.
        """
        rows = len(self._records)
        lines = [
            "Alignment with %i rows and %i columns"
            % (rows, self.get_alignment_length())
        ]
        if rows <= 20:
            lines.extend(self._str_line(rec) for rec in self._records)
        else:
            lines.extend(self._str_line(rec) for rec in self._records[:18])
            lines.append("...")
            lines.append(self._str_line(self._records[-1]))
        return "\n".join(lines)

    def __repr__(self):
        """Return a representation of the object for debugging.
        """
        # A doctest for __repr__ would be nice, but __class__ comes out differently
        # if run via the __main__ trick.
        return "<%s instance (%i records of length %i) at %x>" % (
            self.__class__,
            len(self._records),
            self.get_alignment_length(),
            id(self),
        )
        # This version is useful for doing eval(repr(alignment)),
        # but it can be VERY long:
        # return "%s(%r)" \
        #       % (self.__class__, self._records)

    def __format__(self, format_spec):
        """Return the alignment as a string in the specified file format.
        """
        if format_spec:
            from io import StringIO
            from Bio import AlignIO

            handle = StringIO()
            AlignIO.write([self], handle, format_spec)
            return handle.getvalue()
        else:
            # Follow python convention and default to using __str__
            return str(self)

    def __iter__(self):
        """Iterate over alignment rows as SeqRecord objects.
        """
        return iter(self._records)

    def __len__(self):
        """Return the number of sequences in the alignment.

        Use len(alignment) to get the number of sequences (i.e. the number of
        rows), and alignment.get_alignment_length() to get the length of the
        longest sequence (i.e. the number of columns).

        This is easy to remember if you think of the alignment as being like a
        list of SeqRecord objects.
        """
        return len(self._records)

    def get_alignment_length(self):
        """Return the maximum length of the alignment.
        """
        max_length = 0

        for record in self._records:
            if len(record.seq) > max_length:
                max_length = len(record.seq)

        return max_length

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

    def append(self, record):
        """Add one more SeqRecord object to the alignment as a new row."""

        if self._records:
            self._append(record, self.get_alignment_length())
        else:
            self._append(record)

    def _append(self, record, expected_length=None):
        """Validate and append a record (PRIVATE)."""
        if not isinstance(record, SeqRecord):
            raise TypeError("New sequence is not a SeqRecord object")

        # Currently the get_alignment_length() call is expensive, so we need
        # to avoid calling it repeatedly for __init__ and extend, hence this
        # private _append method
        if expected_length is not None and len(record) != expected_length:
            # TODO - Use the following more helpful error, but update unit tests
            # raise ValueError("New sequence is not of length %i"
            #                  % self.get_alignment_length())
            raise ValueError("Sequences must all be the same length")

        self._records.append(record)

    def __add__(self, other):
        """Combine two alignments with the same number of rows by adding them.

        If you have two multiple sequence alignments (MSAs), there are two ways to think
        about adding them - by row or by column. Using the extend method adds by row.
        Using the addition operator adds by column.

        """
        if not isinstance(other, MultipleSequenceAlignment):
            raise NotImplementedError
        if len(self) != len(other):
            raise ValueError(
                "When adding two alignments they must have the same length"
                " (i.e. same number of rows)"
            )
        merged = (left + right for left, right in zip(self, other))
        # Take any common annotation:
        annotations = {}
        for k, v in self.annotations.items():
            if k in other.annotations and other.annotations[k] == v:
                annotations[k] = v
        column_annotations = {}
        for k, v in self.column_annotations.items():
            if k in other.column_annotations:
                column_annotations[k] = v + other.column_annotations[k]
        return MultipleSequenceAlignment(
            merged, annotations=annotations, column_annotations=column_annotations
        )

    def upper(self):

        """Make all letters in an alignment uppercase."""

        return MultipleSequenceAlignment(np.char.upper(self.matrix), ids=self.ids)

    def replace(self, a, b):

        """Replace one letter in an alignment with another."""

        return MultipleSequenceAlignment(np.char.replace(self.matrix, a, b), ids=self.ids)

    def __getitem__(self, index):
        """Access part of the alignment. Indexes like a 2D numpy array and returns
        another MultipleSequenceAlignment object unless:
            1) you ask for a single row, when you get a SeqRecord object
            2) you ask for a single column, when you get a numpy array
            3) you ask for a single row AND a single column, when you get a single letter as a string."""

        if isinstance(index, int):
            return self._records[index]

        elif isinstance(index, slice):
            return self.__getitem__(index, np.arange(self.L))

        elif isinstance(index, np.ndarray):
            ids = self.ids[index]; names = self.names[index]
            descs = list(np.array(self.descriptions)[row_index])

            if type(self.tree) != type(None):

                if type(x[0])==np.bool_:
                    branches = [self._records[k] for k in np.where(index)[0]]
                else:
                    branches = [self._records[k] for k in index]

                new_msa = MultipleSequenceAlignment(self.matrix[index], ids=ids,
                    names=names, descriptions=descs,_tree_branches=branches)
                new_msa.tree = self.tree

                return new_msa

            else:
                return MultipleSequenceAlignment(self.matrix[index], ids=ids, names=names, descriptions=descs)

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

        #return MultipleSequenceAlignment(self.matrix[index], ids=self.ids)

        # if isinstance(index, int):
        #     # e.g. result = align[x]
        #     # Return a SeqRecord
        #     return self._records[index]
        #
        # elif isinstance(index, slice):
        #     # e.g. sub_align = align[i:j:k]
        #     new = MultipleSequenceAlignment(self._records[index])
        #     if self.column_annotations and len(new) == len(self):
        #         # All rows kept (although could have been reversed)
        #         # Preserve the column annotations too,
        #         for k, v in self.column_annotations.items():
        #             new.column_annotations[k] = v
        #     return new
        # elif len(index) != 2:
        #     raise TypeError("Invalid index type.")
        #
        # # Handle double indexing
        # row_index, col_index = index
        # if isinstance(row_index, int):
        #     # e.g. row_or_part_row = align[6, 1:4], gives a SeqRecord
        #     return self._records[row_index][col_index]
        #
        # elif isinstance(col_index, int):
        #     # e.g. col_or_part_col = align[1:5, 6], gives a string
        #     return np.array([rec[col_index] for rec in self._records[row_index]])
        #
        # else:
        #     # e.g. sub_align = align[1:4, 5:7], gives another alignment
        #     new = MultipleSequenceAlignment(
        #         rec[col_index] for rec in self._records[row_index]
        #     )
        #
        #     if self.column_annotations and len(new) == len(self):
        #         # All rows kept (although could have been reversed)
        #         # Preserve the column annotations too,
        #         for k, v in self.column_annotations.items():
        #             new.column_annotations[k] = v[col_index]
        #     return new

    def sort(self, key=None, reverse=False):
        """Sort the rows (SeqRecord objects) of the alignment in place.

        This sorts the rows alphabetically using the SeqRecord object id by
        default. The sorting can be controlled by supplying a key function
        which must map each SeqRecord to a sort value.

        This is useful if you want to add two alignments which use the same
        record identifiers, but in a different order.

        """
        if key is None:
            self._records.sort(key=lambda r: r.id, reverse=reverse)
        else:
            self._records.sort(key=key, reverse=reverse)

    @property
    def substitutions(self):
        """Return an Array with the number of substitutions of letters in the alignment.

        As an example, consider a multiple sequence alignment of three DNA sequences:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from Bio.Align import MultipleSeqAlignment
        >>> seq1 = SeqRecord(Seq("ACGT"), id="seq1")
        >>> seq2 = SeqRecord(Seq("A--A"), id="seq2")
        >>> seq3 = SeqRecord(Seq("ACGT"), id="seq3")
        >>> seq4 = SeqRecord(Seq("TTTC"), id="seq4")
        >>> alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])
        >>> print(alignment)
        Alignment with 4 rows and 4 columns
        ACGT seq1
        A--A seq2
        ACGT seq3
        TTTC seq4

        >>> m = alignment.substitutions
        >>> print(m)
            A   C   G   T
        A 3.0 0.5 0.0 2.5
        C 0.5 1.0 0.0 2.0
        G 0.0 0.0 1.0 1.0
        T 2.5 2.0 1.0 1.0
        <BLANKLINE>

        Note that the matrix is symmetric, with counts divided equally on both
        sides of the diagonal. For example, the total number of substitutions
        between A and T in the alignment is 3.5 + 3.5 = 7.

        Any weights associated with the sequences are taken into account when
        calculating the substitution matrix.  For example, given the following
        multiple sequence alignment::

            GTATC  0.5
            AT--C  0.8
            CTGTC  1.0

        For the first column we have::

            ('A', 'G') : 0.5 * 0.8 = 0.4
            ('C', 'G') : 0.5 * 1.0 = 0.5
            ('A', 'C') : 0.8 * 1.0 = 0.8

        """
        letters = set.union(*[set(record.seq) for record in self])
        try:
            letters.remove("-")
        except KeyError:
            pass
        letters = "".join(sorted(letters))
        m = substitution_matrices.Array(letters, dims=2)
        for rec_num1, alignment1 in enumerate(self):
            seq1 = alignment1.seq
            weight1 = alignment1.annotations.get("weight", 1.0)
            for rec_num2, alignment2 in enumerate(self):
                if rec_num1 == rec_num2:
                    break
                seq2 = alignment2.seq
                weight2 = alignment2.annotations.get("weight", 1.0)
                for residue1, residue2 in zip(seq1, seq2):
                    if residue1 == "-":
                        continue
                    if residue2 == "-":
                        continue
                    m[(residue1, residue2)] += weight1 * weight2

        m += m.transpose()
        m /= 2.0

        return m


    def search_id(self, id_str):

        """Search for an incomplete ID in the alignment"""

        in_id = np.array([str(id_str) in id for id in self.ids])
        return self.__getitem__(np.where(in_id)[0].astype('int'))

    def search_ids(self, ids):

        return MultipleSequenceAlignment([search_id(id)[0] for id in ids])

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
            tree = Phylo.read(tree, format)

        self.tree = deepcopy(tree)

        for term in self.tree.get_terminals():
            try:
                n = list(self.ids).index(term.name)
            except:
                ns = [n for n,i in enumerate(self.ids) if term.name in i]
                if len(ns) == 0:
                    print("Terminal name", term.name, "not in alignment - will not be attached!")
                    n = -1
                elif len(ns) == 1:
                    n = n[0]
                else:
                    print("Ambiguity with terminal time", term.name, "in alignment", len(ns), "times!")
                    n = -1

            if n > -1:
                self._records[n].branch = term
                term.seq = self._records[n].seq
                term.nseq = n
            elif prune_unmatched:
                self.tree.prune(term)


    def get_index(self, id):

        """Return the index of a particular sequence id. Should work for partial IDs;
            for indices that appear multiple times it will return an array, and for missing
            indices it will return None."""


        if id in self.ids:
            return list(self.ids).index(id)

        else:

            a = self.search_id(id)
            if a.N == 0:
                return None

            elif a.N == 1:
                return self.get_index(a[0].id)

            else:
                return np.array([self.get_index(record.name) for record in a])

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


    def subset_by_clade(self, clade):

        """For a clade in the attached tree, subset the alignment to get only those sequences."""

        return self.__getitem__(np.array([term.nseq for term in clade.get_terminals()]))

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

        A = A[:,np.where(A.matrix[0]!='-')[0]]

        if '-' in A[0].seq:
            print(str(A[0].seq).count('-'), "positions not aligned - removing...")

            x = []
            nres = 0
            for n,res in enumerate(A.matrix[1]):
                if res!='-':
                    nres+=1
                    x.append(A.matrix[0][n] != '-')

            new_structure = structure[np.array(x)]
            self.attach_structure(structure, name=name)

        else:
            structure._alipos = np.where(A.matrix[1]!='-')[0]

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


def readPDB(filename, name, chains='A', model=0):

    '''
    Given a PDB code and a directory of where to find the structure, returns
    the distance matrix and residue ids.
    '''

    # Load the structure using Bio.PDB
    structure = PDB.PDBParser(QUIET=True).get_structure(name, filename)

    # The number model to use in the structure; usually there is only one, and this is 0 (default)
    m = structure[model]

    # Which chains to extract - either all chains, or a single chain
    if isinstance(chains, str):
        if chains.lower()=='all':    # Just extract all chains
            structure_list = []
            for ch in m:
                structure_list.append(ProteinStructure(ch, name=name+'-'+ch.id))
            return structure_list

        else:    # If it's a string and not 'All', we assume it's a single chain index
            try:    # A try... except here so that it can tell you if you put in an incorrect chain ID
                ch = m[chains]
            except:
                print(chains, 'is not a valid chain ID in', filename)
                return None

            return ProteinStructure(ch, name=name)

    else:
        # If you pass it a list, it will extract those chains
        structure_list = []
        for ch in m:
            structure_list.append(ProteinStructure(ch, name=name+'-'+ch.id))
        return structure_list


def readAlphaFold(filename, name, pae_file=None):

    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(name, filename)
    prot = ProteinStructure(structure[0]["A"], name=name)

    if type(pae_file)!=None:
        pae_matrix = np.genfromtxt(pae_file)
        prot.pae = pae_matrix

    return prot


letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}

class ProteinStructure:

    '''
    A class for protein structures, which uses Biopython residue and atom objects but is intended to
    be much easier to use that a Biopython <structure> object, which are just a complete mess.
    '''

    def __init__(self, chain, annotation_type='bfactor', annotate_by='residue', resolution=None, name=None):

        self.chain = chain
        self.residues = [residue for residue in chain]
        self.residue_ids = np.array([residue.id[1] for residue in chain])
        self.sequence = ''.join([letters[residue.resname] for residue in chain if residue.resname in letters])
        self.ordered_sequence = ''.join([letters[residue.resname] for residue in chain if residue.resname in letters and not residue.is_disordered()])

        self.n_residues = len(self.residues)
        self.atoms = np.concatenate([[atom for atom in residue] for residue in self.residues])
        self.n_atoms = len(self.atoms)
        self.xyz = np.array([atom.coord for atom in self.atoms])
        self.dmatrix = None
        self._annotation_type = annotation_type
        self._annotate_by = annotate_by
        self.resolution = resolution
        self.name = name
        if len(self.name) == 4:
            self.pdb = self.name.upper()

        if annotation_type.lower()=='bfactor':

            if 'res' in annotate_by.lower():
                bfactors = []
                for residue in self.residues:
                    bfactors.append(np.mean([atom.bfactor for atom in residue]))
                self.bfactors = np.array(bfactors)
            else:
                self.bfactors = np.array([atom.bfactor for atom in self.atoms])

        elif annotation_type.lower()=='plddm':

            plddms = []
            for residue in self.residues:
                plddms.append(np.max([atom.bfactor for atom in residue]))

            self.plddm = np.array(plddms)

    def __len__(self):
        return self.n_residues

    def __str__(self):
        return "Protein structure "+self.name+" with "+str(self.n_residues)+" residues and "+str(self.n_atoms)+" atoms"

    def __getitem__(self, index):

        if isinstance(index, int):
            try:
                ind = int(np.where(self.residue_ids==index)[0][0])
                return self.residues[ind]
            except:
                return ValueError("No such residue in protein!")

        elif isinstance(index, slice):
            sel = np.array([np.where(self.residue_ids==i)[0][0] for i in index])
            return ProteinStructure(self.residues[sel])

        elif isinstance(index, (list, np.ndarray)):
            return ProteinStructure([self.residues[k] for k in index])

    def __iter__(self):
        return iter(self.residues)

    def select_atoms(self, key):

        if isinstance(key, str):

            # if ' or ' in str:
            #     as = []
            #     for x in key.split(' or '):
            #         as.append([n for n,atom in enumerate(self.atoms) if atom.name==x])
            #     return self.select_atoms(np.concatenate(as))

            atom_sel = [n for n,atom in enumerate(self.atoms) if atom.name==key]
            return self.select_atoms(atom_sel)

        elif isinstance(key, int):
            return self.atoms[key]

        else:
            residues_in_selection = [self.atoms[k].get_parent().id[1] for k in key]
            unique_residues = np.unique(residues_in_selection)
            where_res = [np.where(np.array(residues_in_selection)==res)[0] for res in unique_residues]

            final_residue_list = []

            for n,resn in enumerate(unique_residues):
                res = self.residues[list(self.residue_ids).index(resn)]
                newres = PDB.Residue.Residue(res.id, res.resname, res.segid)

                for atom_index in where_res[n]:
                    newres.add(copy(self.atoms[key[atom_index]]))

                final_residue_list.append(newres)

            return ProteinStructure(final_residue_list, annotation_type=self._annotation_type,
                                    annotate_by=self._annotate_by, resolution=self.resolution, name=self.name)


    def distance_matrix(self, method="CA", recalculate=False):

        if type(self.dmatrix)==type(None) or recalculate==True:

            if method=="CA":
                self.dmatrix= CalcDistanceMatrix(self.select_atoms("CA"))
                self.distance_method = "CA"
            elif method.lower() in ('full', 'all'):
                self.dmatrix = CalcDistanceMatrix(self)
                self.distance_method = "full"

        elif self.distance_method != method:
            if method=="CA":
                self.dmatrix= CalcDistanceMatrix(self.select_atoms("CA"))
                self.distance_method = "CA"
            elif method.lower() in ('full', 'all'):
                self.dmatrix = CalcDistanceMatrix(self)
                self.distance_method = "full"

        return self.dmatrix

    def contacts(self, threshold, dist_threshold=None):

        if type(self.dmatrix)==type(None):
            self.dmatrix = CalcDistanceMatrix(self)

        if type(dist_threshold)==type(None):
            return np.where(self.dmatrix < threshold)
        else:
            s = self.dmatrix.shape[0]
            off_diag = np.abs(np.arange(s)[None,:] - np.arange(s)[:,None]) > dist_threshold
            return np.where((self.dmatrix < threshold)&(off_diag))


    def annotate_secondary_structure(dssp_path='mkdssp', filename='temp.pdb'):

        self.save(self, filename)
        self.ss = PDB.DSSP.DSSP(chain, filename, dssp_path)

    def save(self, pdb_file):
        io = PDB.PDBIO()
        io.set_structure(self.chain)
        io.save(pdb_file)

def CalcDistanceMatrix(Structure):

    D = []

    for i in range(len(Structure.xyz)):
        D.append(np.sqrt(np.sum((Structure.xyz[i]-Structure.xyz)**2,axis=1)))

    return np.array(D)

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
