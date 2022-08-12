import numpy as np
import pandas as pd
from Bio import PDB, pairwise2
from Bio.PDB.DSSP import DSSP

from copy import copy

def readPDB(filename, name=None, chains='A', model=0):

    '''
    Given a PDB code and a directory of where to find the structure, returns
    the distance matrix and residue ids.

    Inputs:
        filename: the filepath to your .pdb file
        name: an optional name to assign the structure, can be the PDB code
            (generated automatically from the filename if you don't pass one)
        chains: which chains of the structure to load. Defaults to "A."
            Can also pass 'all', which will load all chains.
        model: which model to use from multi-model PDBs. Defaults to 0."

    Outputs:
        sturcture: a ProteinStructure object
    '''

    # If you don't pass a name, automatically fill in a name from the pdb file
    #  (for some reason Biopython's PDBParser requires a "PDB code" which can really
    #    just be any name, doesn't have to be the real PDB code)
    if type(name) == type(None):
        if '/' in filename:
            end_path = filename.split('/')[-1]
            if '.' in end_path:
                name = end_path.split('.')[0]
            else:
                name = end_path
        else:
            if '.' in filename:
                name = filename.split(".")[-2]
            name = filename

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
        self.dmatrix = None; self.distance_method = None
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
            return ProteinStructure(residues[index], name=self.name,
                    annotation_type=self._annotation_type, annotate_by=self._annotate_by,
                    resolution=self.resolution)

        elif isinstance(index, (list, np.ndarray)):
            return ProteinStructure([self.residues[k] for k in index], name=self.name,
                    annotation_type=self._annotation_type, annotate_by=self._annotate_by,
                    resolution=self.resolution)

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

    def distance(self, i, j, method="CA", atoms=None):

        '''
        Get the distance between two residues in a structure by passing the numbers of two residues.

        Has three available 'modes' of calculating the distance:

            "C-alphas" : gives the distance between the alpha carbons of each residue (DEFAULT)
                           set method="CA", atoms=None (Default)
            "min" :      minimum distance between any atoms in the residues
                           set method="min", atoms=None (Default)
            "atoms" :      pick user-defined atoms from each residue to get a very specific distance
                            atoms = (ID_1, ID_2), e.g. atoms=("O1", "O2") gives distance between O1 in the
                            first residue and O2 in the second residue.
                            the method keyword can be set to anything, it is ignored

        '''

        if isinstance(atoms, (list, tuple)):

            if len(atoms)==2:
                a1 = ProteinStructure([structure_1.residues[i]], name='p1').select(atoms(atoms[0])).xyz[0]
                a2 = ProteinStructure([structure_1.residues[j]], name='p2').select(atoms(atoms[1])).xyz[0]

                return np.sqrt(np.sum((a1 - a2)**2))

            else:
                return IndexError(
                    "To get distances between user-selected atoms, you must pass a list of 2 valid PDB atom indices"
                )

        elif method.upper() == 'CA':
            calphas = self.select_atoms("CA")
            return np.sqrt(np.sum((calphas.xyz[i]-calphas.xyz[j])**2))

        elif method.upper() == 'MIN':
            a1 = ProteinStructure([structure_1.residues[i]], name='p1').xyz
            a2 = ProteinStructure([structure_1.residues[j]], name='p2').xyz
            return np.min(np.sqrt(np.sum((a1[:,None,:] - a2[None,:,:])**2, axis=2)))

    def contacts(self, threshold, dist_threshold=None):

        return get_contacts(self, threshold, dist_threshold)

    def annotate_secondary_structure(dssp_path='mkdssp', filename='temp.pdb'):

        self.save(self, filename)
        self.ss = PDB.DSSP.DSSP(chain, filename, dssp_path)

    def ordered_residues(self):
        return [res for res in self.residues if not res.is_disordered()]

    def save(self, pdb_file):
        io = PDB.PDBIO()
        io.set_structure(self.chain)
        io.save(pdb_file)

#class Complex(Structure):

def get_contacts(structure, thresh, dist_threshold=None):

    '''
    Calculate contacts for which the C-alphas are within some threshold for a ProteinStructure object.
    '''

    if type(self.dmatrix)==type(None):
        self.dmatrix = CalcDistanceMatrix(self)

    if type(dist_threshold)==type(None):

        return structure.residue_ids[np.where(structure.distance_matrix() < thresh)[0]], \
               structure.residue_ids[np.where(structure.distance_matrix() < thresh)[1]]

    else:
        s = self.dmatrix.shape[0]
        off_diag = np.abs(np.arange(s)[None,:] - np.arange(s)[:,None]) > dist_threshold

        return structure.residue_ids[np.where((structure.distance_matrix() < thresh)&(off_diag))[0]], \
               structure.residue_ids[np.where((structure.distance_matrix() < thresh)&(off_diag))[1]]



def CalcDistanceMatrix(Structure):

    D = []

    for i in range(len(Structure.xyz)):
        D.append(np.sqrt(np.sum((Structure.xyz[i]-Structure.xyz)**2,axis=1)))

    return np.array(D)
