# TODO:
# 1) [DONE] figure out how to read a mol2 file with multiple molecules into rdkit
# 2) [SEMI-DONE] solve or accound for all errors that occur while reading in files
# 3) [DONE] maybe implement a dataframe or something to help organize files/molecules
# 4) [SEMI-DONE: see #2] go through and test that all molecules were read in correctly
# 5) [NOT DONE: no need this week] unpack or write to a compressed gz (learn gzip moldule)file
# 6) MiniDrugBank.sdf: do rdkit and openeye read all these molecules the same?
# 7) - how to compare molecules that are read into rdkit vs. openeye?
# 8) - compare, graphs(stereo, aromaticity, connectivity, atoms), partial charge, isomers(?), conformers etc.
# 9) - better coding pratices, less experimental code, more OO programming, get a final product
# 10) - make the software very readable, easy to follow the logical flow.
# 11) - make functions first, think about classes later 
# 12) - what are the details of a jupyter(?) notebook. 

#RDKit modules
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
#OpenEye modules
from openeye import oechem
from openeye import oedepict
#Other modules
import csv
import os
from pathlib import Path
import pandas as pd
import numpy as np
import warnings
import inspect
IPythonConsole.ipython_useSVG=False
# class and function definitions

class ForwardMol2MolSupplier:
    def __init__(self, file):
        self.file = file
        self.at_end_of_file = False
        self.file_length = sum(1 for line in open(file, "r"))
        self.is_mol2()  # checks for file identity
        self.current_mol = 0

    def is_mol2(self):
        if self.file.suffix != '.mol2':
            raise NameError("mol2 file required for class \"ForwardMol2MolSupplier\"")
        # check to make sure only TRIPOS files are used
        file = open(self.file, "r")
        line = file.readline()
        if line != "@<TRIPOS>MOLECULE\n":
            raise NameError("only TRIPOS format is supported")

    def get_next_mol2_block(self):
        """
        read from the current file index "curr" to the start of the next molecule
        """
        f = open(self.file, "r")

        mol2string = "@<TRIPOS>MOLECULE"
        for i, line in enumerate(f):
            # parse through until the correct line has been reached
            if i < self.current_mol:
                continue
            # if we have reached the end of the current molecule, break
            if line == "@<TRIPOS>MOLECULE\n":
                if i == self.current_mol:   #if at start of the molecule
                    mol2string = "@<TRIPOS>MOLECULE\n"
                    continue
                else:                       #if at the end
                    self.current_mol = i
                    break
            if i == self.current_mol + 1:   #get the identifier of the molecule
                name = line.split("\n")[0]  #remove newline 
            mol2string = "".join((mol2string, line))
            if (i + 1 == self.file_length): # if line is the last line of the file
                self.at_end_of_file = True
                self.current_mol = i + 1
        return name, mol2string

    def __iter__(self):
        self.current_mol = 0
        self.at_end_of_file = False
        return self
    def __next__(self):
        if self.at_end_of_file:
            raise StopIteration
        mol_name, mol2_block = self.get_next_mol2_block()
        return mol_name, Chem.rdmolfiles.MolFromMol2Block(mol2_block, sanitize=False, removeHs=False)

# we should be in '/home/coda3831/anaconda3/envs/openff-dev'
try:
    os.chdir(os.getcwd() + '/openforcefield/openforcefield/Connor')
except Exception as e:
    print(f"cwd = {os.getcwd()} \n check your workspace directory")
    print(e)

drug_bank_mol2 = Path('/home/coda3831/anaconda3/envs/openff-dev/openforcefield/openforcefield/data/molecules/MiniDrugBank_tripos.mol2')
drug_bank_sdf = Path('/home/coda3831/anaconda3/envs/openff-dev/openforcefield/openforcefield/data/molecules/MiniDrugBank.sdf')

mols = pd.DataFrame(columns=['name', 'rdkit_mol2', 'rdkit_sdf', 'openeye_mol2', 'openeye_sdf'])

rdkit_mol2s = ForwardMol2MolSupplier(drug_bank_mol2)  # returns an iterator class    
rdkit_sdfs = Chem.rdmolfiles.ForwardSDMolSupplier(str(drug_bank_sdf), sanitize=True, removeHs=False) # returns an iterator class

openeye_mol2s = oechem.oemolistream(".mol2")
openeye_mol2s.open(str(drug_bank_mol2))
openeye_sdfs = oechem.oemolistream(".sdf")
openeye_sdfs.open(str(drug_bank_sdf))

# test to make sure all iterators have the same length and that the same # of molecuels were read
# sum1 = sum(1 for molecule in rdkit_mol2s)
# sum2 = sum(1 for molecule in rdkit_sdfs)
# sum3 = sum(1 for molecule in openeye_mol2s.GetOEGraphMols())
# sum4 = sum(1 for molecule in openeye_sdfs.GetOEGraphMols())

# assert sum1==sum2==sum3==sum4        # all are equal to 371 (are there 371 molecules in the DrugBank?)

print("\n\n\n_______________________________________________________________________________________")

for rdkit_mol2, rdkit_sdf_mol in zip(rdkit_mol2s, rdkit_sdfs):
    # TODO: Better way to do this? I'm having to rely on the fact that rdkit and my class will
    # read the files in exactly the same way. I will include some tests to make sure there are no reading differences
    name, rdkit_mol2_mol = rdkit_mol2    # this needs to be unpacked since my class returns a tuple 
    # A few tests to make sure all molecules are read in the same
    try:
        mols = mols.append({'name': name,
                     'rdkit_mol2': rdkit_mol2_mol,
                      'rdkit_sdf': rdkit_sdf_mol,
                   'openeye_mol2': None,
                    'openeye_sdf': None},
                     ignore_index=True)
    except Exception:
        print(f"\n---failure with {name} during assertion---\n")
        mols = mols.append({'name': name,
                     'rdkit_mol2': np.NaN,
                      'rdkit_sdf': np.NaN,
                   'openeye_mol2': np.NaN,
                    'openeye_sdf': np.NaN},
                     ignore_index=True)


openeye_mol2_list = []
for openeye_mol2_mol in openeye_mol2s.GetOEGraphMols():
    openeye_mol2_list.append(oechem.OEGraphMol(openeye_mol2_mol))

openeye_sdf_list = []
for openeye_sdf_mol in openeye_sdfs.GetOEGraphMols():
    openeye_sdf_list.append(oechem.OEGraphMol(openeye_sdf_mol))
test = 1
mols["openeye_mol2"] = openeye_mol2_list
mols["openeye_sdf"] = openeye_sdf_list

# openeye reads file from bottom to top. Must reverse these columns to align with rdkit data with openeye data
# for whatever reason, I cannot reverse the order of the openeye columns
# mols['rdkit_mol2'] = mols['rdkit_mol2'].values[::-1]
# mols['rdkit_sdf'] = mols['rdkit_sdf'].values[::-1]
test = 1
# now, we have read in all the molecules in all the combinations that they could be read in
# we can begin testing to see if all of these molecule representations are the same. This will
# be done using a number of functions:
# 1) 

# the point of this class is to get rdkit and openeye information into a format that allows for
# the data to be easily compared. This is also just so I can get some practice with python classes 
class Atom:
    def __init__(self, atomic_num, Idx, is_aromatic, formal_charge):
        self.edges = []
        self.atomic_num = atomic_num 
        self.Idx = Idx
        self.is_aromatic = is_aromatic


def compare(mols):
    """
    params: 
        mols - list of mols to compare. Can be all rdchem.Mol, all OEGraphMol, or a mixture of the two
    returns: 
        (bool) True if all mols have the same graph, False otherwise
    """
    equal = 1
    graphs = []    # list of lists, list of molecule atoms
    for mol in mols:
        molecule_atoms = []
        if isinstance(mol, Chem.rdchem.Mol):
            # def __init__(self, name, Idx, is_aromatic, formal_charge)
            # TODO: add formal charge, stereochem, weight
            for a in mol.GetAtoms():   # "a" for atom
                atomic_num = a.GetAtomicNum()
                Idx = a.GetIdx()
                is_aromatic = a.GetIsAromatic()
                formal_charge = a.GetFormalCharge()
                atom_class = Atom(atomic_num, Idx, is_aromatic, formal_charge)

                neighbor_atoms = [neighbor.GetAtomicNum() for neighbor in a.GetNeighbors()]
                neighbor_atoms.sort()
                atom_class.edges = neighbor_atoms
                molecule_atoms.append(atom_class)

        elif isinstance(mol, oechem.OEGraphMol):
            for a in mol.GetAtoms():
                atomic_num = a.GetAtomicNum()
                Idx = a.GetIdx()
                is_aromatic = a.IsAromatic()
                formal_charge = a.GetFormalCharge()
                atom_class = Atom(atomic_num, Idx, is_aromatic, formal_charge)

                neighbor_atoms = [neighbor.GetAtomicNum() for neighbor in a.GetAtoms()]
                neighbor_atoms.sort()
                atom_class.edges = neighbor_atoms
                molecule_atoms.append(atom_class)

        else:
            print(f"improper class <{type(mol)}> passed to compare_graph \n \
            only classes Chem.rdchem.Mol and oechem.OEGraphMol allowed")
            return -1
        graphs.append(molecule_atoms)

    test=1
    # now we can go through and compare each molecule atom-by-atom
    # print('          rdkit            openeye ')
    # print('     mol2        sdf    mol2     sdf')
    for atoms in zip(*graphs):     # this is just a really cool trick-REMEMBER THIS: "4.7.5 Unpacking Argument Lists"
        prev, *rest = atoms
        # print(f"{prev.Idx:02}", end=":  ")
        for atom in atoms:     # this compares the first element to itself but oh well 
            # compare all attributes to atom before it, if anything is not equal, print a statment and return
            if atom.edges != prev.edges:
                print("atom neighbors not equal")
                equal = 0
            if atom.atomic_num != prev.atomic_num:
                print("atom atomic number not equal")
                equal = 0
            if atom.Idx != prev.Idx:
                print("atom Idx not equal")
                equal = 0
            if atom.is_aromatic != prev.is_aromatic:
                # print("atom aromaticity not equal")
                equal = 0
            
            # print(atom.is_aromatic, end="\t")
            prev = atom
        # print("")
    return equal

num_equal = 0
num_unequal = 0
num_failed = 0
for index, row in mols.iterrows():
    name = row[0]
    m = row[2:]
    comparison = compare(m)
    if comparison == 1:
        num_equal += 1
        # print(f"equal!!!!: {name} ({index})")
    elif comparison == -1:
        num_failed += 1
    elif comparison == 0:
        num_unequal += 1
    else:
        print("something weird happened")
print(f"number equal: {num_equal}")
print(f"number not equal: {num_unequal}")
print(f"number falied: {num_failed}")

test = mols.iloc[0, :]
test = test.to_list()
name = test[0]
test = test[1:]
compare(test)

# RANDOM NOTES:
# [print([atom.atomic_num for atom in atoms]) for atoms in zip(*graphs)]
# 7) - how to compare molecules that are read into rdkit vs. openeye?
# 8) - compare, graphs(stereo, aromaticity, connectivity, atoms), partial charge, isomers(?), conformers etc.




