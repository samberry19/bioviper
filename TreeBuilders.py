import numpy as np
from Bio import Phylo,AlignIO,SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
from Bio.Phylo.Applications import FastTreeCommandline,PhymlCommandline
#from MSA import MultipleSequenceAlignment
import subprocess

blosum_options = np.array([45, 50, 55, 60, 62, 65, 70, 75, 80, 85, 90, 95])

def CalcIdentityMatrix(msa):

    '''Code to calculate a percent identity matrix between every pair of sequences
    in a multiple sequence alignment. Uses numpy array broadcasting and so is relatively
    memory-intensive but much faster than BioPython's build in method'''

    IM = []
    for seq in msa.matrix:

        # Take the mean identity of the sequence with every other sequence
        IM.append(np.mean(seq==msa.matrix, axis=1))

    return IM

def autodetect_blosum_matrix(msa):

    '''Construct an identity matrix to determine what the best BLOSUM matrix to use is.'''

    ssm = CalcIdentityMatrix(msa)
    best_blosum = blosum_options[np.argmin(np.abs(blosum_options - np.mean(ssm)))]
    return 'blosum'+str(best_blosum)

def NJTree(msa, matrix='autodetect', root=True, ladderize=True):

    """Construct a neighbor-joining distance tree from a sequence alignment."""

    if matrix == 'autodetect':
        matrix = autodetect_blosum_matrix(msa)

    calculator = DistanceCalculator(matrix)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(msa)

    if root:
        tree.root_at_midpoint()
    if ladderize:
        tree.ladderize()

    return tree


def UPGMATree(msa, matrix='autodetect', root=True, ladderize=True):

    """Construct a UPGMA distance tree from a sequence alignment."""

    if matrix=='autodetect':
        matrix = autodetect_blosum_matrix(msa)

    calculator = DistanceCalculator(matrix)
    constructor = DistanceTreeConstructor(calculator, 'upgma')
    tree = constructor.build_tree(ali)

    if root:
        tree.root_at_midpoint()
    if ladderize:
        tree.ladderize()

    return tree


def MuscleAlign(seq1, seq2, muscle_path, seq_filename='temp.fa', aln_filename='temp_aln.fa'):

    '''
    Align two sequences using command-line MUSCLE from Python
    '''

    #Convert sequences as trings to Biopython sequence objects
    s1=SeqIO.SeqRecord(Seq(seq1), id='s1', description='')
    s2=SeqIO.SeqRecord(Seq(seq2), id='s2', description='')

    # Write them to a file
    SeqIO.write([s1,s2], seq_filename, 'fasta')

    # Run muscle using subprocess.call()
    subprocess.call(muscle_path + ' -in ' + seq_filename + ' -out ' + aln_filename, shell=True)

    # Load back in and return the raw sequences
    return [str(i.seq) for i in AlignIO.read(aln_filename, 'fasta')]

def FastTree(aln, fasttree_path, gamma=True, model='LG', cat=20, intree=None, fixed_topology=False,
             out_dir='./', outtree = 'fasttree.nwk', outfile='fasttree_out.txt', temp_aln='temp_aln.fa', aln_format='fasta'):

    '''
    Run command-line FastTree from within your Python script, and get back the tree and the rate parameters fit to it.

    Inputs:
        aln: your alignment, as a MultipleSequenceAlignment object
        fasttree_path: full path to your local download of fasttree
        optional:
            gamma: fit gamma rate parameters to the data. Defaults to true.
            model: evolutionary model for the data. Defaults to 'LG'
            cat: number of rate categories. Defaults to 20. Currently may not work with any other number...
            intree: pass a starting tree for the FastTree simulation. Defaults to None
            fixed_topology: don't infer trees, simply infer parameters based on the intree.
                This currently doesn't seem to work right; I am troubleshooting it.
            out_dir: output directory to put files in; defaults to the current directory
            outtree: name of the tree file that will be created; defaults to fasttree.nwk
            outfile: name of the output file that will be created; defaults to fasttree_out.txt
            temp_aln: name of the file to house the temporary alignment file made from the saved alignment.
            aln_format: format of the alignment file. Defaults to fasta. I have no idea why you'd change this but it's here

    Returns:
        tree: your tree as a Bio.Phylo tree
        (if gamma==True): rp, a dictionary of rate parameters where:
            rp[alpha] gives the alpha parameter for the gamma distribution
            rp[rates] gives the rates across each site
            RateCategoryProbs(rp, L) gives the probabilities of each category
    '''

    # save your alignment as a temporary fasta file
    AlignIO.write(aln, temp_aln, format=aln_format)

    # we want to make sure that the model is lowercase, and not fail if it started uppercase
    m = model.lower()

    if gamma==True:
        g= '-gamma '
    else:
        g=''

    if type(intree)==type(None):
        tree_args = ' '
    else:
        # If you passed it a tree, also have to write a treefile
        Phylo.write(intree, out_dir+'temp_tree.nwk', 'newick')
        if fixed_topology==True:
            tree_args = ' -intree '+out_dir+'temp_tree.nwk -nome -mllen '
        else:
            tree_args = ' -intree '+out_dir+'temp_tree.nwk '

    # this is where fasttree is called with all of the parameters set above
    subprocess.call(fasttree_path+tree_args+g+'-'+m+' -cat '+str(cat)+' -log '+out_dir+outfile+' '+out_dir+temp_aln+' > '+out_dir+outtree, shell=True)

    # read in the tree from the tree file
    tree = Phylo.read(out_dir+outtree, 'newick')

    # read back in the rate parameters from the outfile
    if gamma==True:
        rp = ReadRateParams(out_dir+outfile)

        return tree,rp

    else:
        return tree

    if gamma==True:
        g= '-gamma '
    else:
        g=''

    if type(intree)==type(None):
        tree_args = ' '
    else:
        # If you passed it a tree, also have to write a treefile
        Phylo.write(intree, out_dir+'temp_tree.nwk', 'newick')
        if fixed_topology==True:
            tree_args = ' -intree '+out_dir+'temp_tree.nwk -nome -mllen '
        else:
            tree_args = ' -intree '+out_dir+'temp_tree.nwk '

    # this is where fasttree is called with all of the parameters set above
    subprocess.call(fasttree_path+tree_args+g+'-'+m+' -cat '+str(cat)+' -log '+out_dir+outfile+' '+out_dir+temp_aln+' > '+out_dir+outtree, shell=True)

    # read in the tree from the tree file
    tree = Phylo.read(out_dir+outtree, 'newick')

    # read back in the rate parameters from the outfile
    if gamma==True:
        rp = ReadRateParams(out_dir+outfile)

        return tree,rp

    else:
        return tree
