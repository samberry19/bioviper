# bioviper

**Warning**: this code is still under development - while most parts are functional, there are still likely bugs and features that are yet to be implemented. Make a pull request if you have an issue!

A set of classes and functions built on top of Biopython to make working with biological data - primarily sequence alignments, protein structures and phylogenetic trees - simpler and more straightforward.

The core of this code is a set of three new classes: a MultipleSequenceAlignment class, a SequenceArray class (for unaligned sequences) and a ProteinStructure class.

I am working on a set of Jupyter notebooks that should explain the basic usage of these objects.

## Installation

On Mac/Linux, I would recommend cloning or downloading this repository and then installing it into a conda environment. You can create a conda environment, let's called it "evo," assuming you have conda installed with:

```
conda create -n evo pip numpy pandas Biopython
```

Then navigate to the folder alignment_tools and type:

```
conda install .
```

I have no idea how to do this on Windows, but let me know if you install it and how and I can add it here :)

## Basics

To load an alignment from a file:

```
msa = readAlignment("alignment.fa", format="fasta")
```

Much like a traditional Biopython MSA, this object is essentially a list of Bio.SeqRecord.SeqRecord objects, but it can do a few nice new things:

```
msa.ids                      # gives you all of the ids in the sequence alignment
msa.matrix                   # gives you the alignment as a 2D numpy array
msa[np.array([3,5,6,9])]     # gives a new MSA with the 3rd, 5th, 6th and 9th sequences
msa[:,np.array([3,5,6,9])]   # gives a new MSA with the 3rd, 5th, 6th and 9th columns
msa.calc_coverage()          # calculates the alignment coverage and stores it as msa.coverage
msa.calc_frequencies()       # calculates the frequency of each amino acid and stores it as msa.frequencies
msa.search_id("HUMAN")       # returns all sequences with "HUMAN" in the id
msa.search_sequence("ACYWL") # searches for sequence record(s) with the following subsequence
msa.sequence_lengths()       # returns the length of each sequence in the alignment, not countin gaps
msa.dealign()                # gives you all dealigned sequences (gaps removed)
msa.save(filename)           # saves to a file
```

To load a structure from a file:

```
structure = readPDB("6bu5.pdb", name="Nramp outward-open", model=0, chain="A")
```

There is no easy-to-use protein structure object in Biopython, so this one has mostly entirely new functionality (while incorporating Biopython
residue and atom objects). Simple usage includes:

```
structure.sequence              # the full sequence from the PDB file
structure.ordered_sequence      # the part of the sequence for which there is structural data (the ordered residues)
structure.xyz                   # the xyz coordinates of the structure as a 2D numpy array
structure.residues              # a list of all of the residue objects
structure.atoms                 # a list of all of the atom objects
structure[np.array([3,5,6,9])   # returns a ProteinStructure with only the 34d, 5th, 6th and 9th residues 
structure.select_atoms("CA")    # returns a ProteinStructure with only the C-alpha atoms
structure.distance_matrix()     # calculates (if new) or returns (if already calculated) a distance matrix for the structure
                                   # by default, does so for only the C-alphas, but can also accept other arguments (see notebooks)
structure.contacts(7)           # returns all structural contacts within 7 Å (based on the distance matrix)
structure.assign_ss()           # assign secondary structures using DSSP (requires DSSP to be installed)
structure.save(filename.pdb)    # saves the structure as a pdb file
```

## Attaching structures to a sequence alignment

One new feature of this code is that structures and phylogenetic trees can be "attached" to a multiple sequence alignment in order to facilitate
analysis that integrates the different pieces of information.

A structure can be attached to first sequence of the MSA with ``msa.attach_structure(structure, 0)``. If you don't provide a sequence number, the program will attempt to find a 
matching sequence using ``msa.search_sequence()`` - this requires that the structure not have any sequence mutations. 
Multiple structures can be attached to the same sequence and can be accessed with ``msa[0].structures``; if a single structure is attached, this will be found at ``msa[0].structures[0]``.
The structure's sequence will be automatically aligned to its sequence in the MSA and the corresponding 
MSA positions to each position in the structure can be accessed with ``structure._alipos``, which is used in functions such as ``distance_difference``.

## Attaching phylogenetic trees to a sequence alignment

A phylogenetic tree with terminal ids matching the ids in the sequence alignment can be loaded as follows:
```
tree = Phylo.read("example_tree.nwk", "newick")
msa.attach_tree(tree)
```

Both the tree and the MSA's sequence records will now gain new attributes. For each sequence record that corresponded to a branch in the tree, that branch can
now be accessed as ``record.branch``. For each branch in the tree, the index of its corresponding sequence is stored in ``terminal_node.nseq``. In order to select a 
portion of the MSA based on a clade in the tree, you can call e.g. ``msa.subset_by_clade(msa.tree.clade[0][0][1])``. 

One note to keep in mind is that if you subset the sequences in the MSA, it will pass along the branch attributes but not automically correct the ``nseq`` values in
the tree or prune the tree to contain only the sequences in the new MSA. This can be done, however, by calling ``msa.fix_tree()`` on the new alignment; it is only not
implemented by default because for larger alignments this may be relatively slow.
