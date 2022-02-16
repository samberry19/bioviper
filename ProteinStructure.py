import numpy as np
import pandas as pd
from Bio import PDB, pairwise2

from Bio.Align import _aligners
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord,_RestrictedDict
from Bio.PDB.DSSP import DSSP

from copy import copy

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

    def ordered_residues(self):
        return [res for res in self.residues if not res.is_disordered()]

    def save(self, pdb_file):
        io = PDB.PDBIO()
        io.set_structure(self.chain)
        io.save(pdb_file)

def CalcDistanceMatrix(Structure):

    D = []

    for i in range(len(Structure.xyz)):
        D.append(np.sqrt(np.sum((Structure.xyz[i]-Structure.xyz)**2,axis=1)))

    return np.array(D)
