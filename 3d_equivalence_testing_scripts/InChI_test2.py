#RDKit modules
from rdkit import Chem
#OpenEye modules
from openeye import oechem
from openeye import oedepict
#Other modules
import csv
import os
from pathlib import Path
# we should be in '/home/coda3831/anaconda3/envs/openff-dev'
try:
    os.chdir(os.getcwd() + '/openforcefield/openforcefield/Connor')
except Exception as e:
    print(f"cwd = {os.getcwd()} \n check your workspace directory")
    print(e)

    # first, load in a series of smiles mols from a file
mols = [] # will be a list of dictionaries 
with open('smiles_molecules.txt', 'r', newline='') as file:
    fileReader = csv.DictReader(file, delimiter= ' ', 
                        skipinitialspace=True, fieldnames=['iupac', 'smiles', 'inchi', 'is_aromatic'])
        # iupac: IUPAC name
        # smiles: smiles format, to be read into RDKit and openeye
        # inchi: Inchi name, mostly to compare if it can distinguish benzene and cyclohexane (aromatic vs. not aromatic)
        # is_aromatic: stores whether I can tell if the molecules is aromatic (1) or not (0)
    for row in fileReader:
        mols.append(row)
    # Let's load a series of molecules into RDKit and OpenEye...
for mol in mols:
    mol['rdkit_mol'] = Chem.MolFromSmiles(mol['smiles'])  
    mol['rdkit_mol'].GetAtoms()  #rdkit directly returns a mol
    openEyeMol = oechem.OEGraphMol()                        
    oechem.OESmilesToMol(openEyeMol, mol['smiles'])         #openeye stores the interpreted molecule in an OEGraphMol object, returns 
    mol['openeye_mol'] = openEyeMol                         #true if no errors were recieved
test = 1

# Do OEMols/RDMols store aromaticity? 
# ***Answer***: YES
for mol in mols:
    print(f"molecule: _____{mol['iupac']}_____, is_aromatic: {mol['is_aromatic']}")
    # Find number of aromatic atoms in the molecule
    num_aromatic_atoms = 0
    for atom in mol['rdkit_mol'].GetAromaticAtoms():
        num_aromatic_atoms += 1
    print(f"\t RDKit found {num_aromatic_atoms} aromatic atoms")

    # openeye requires that you check every atom for aromaticity 
    num_aromatic_atoms = 0
    for atom in mol['openeye_mol'].GetAtoms():
        if atom.IsAromatic() == True:
            num_aromatic_atoms += 1
    print(f"\t OpenEye found {num_aromatic_atoms} aromatic atoms")

# Can you derive it from other info? 
# ***Answer***: Maybe? It would be hard and openeye/rdkit percieve this already if given enough info
    # Criteria for aromaticity: 
        # 1) circular: you could search through the atoms and their neighbors recursively until you found a repeat Idx? 
        # 2) conjugated: No idea. Search for alternating single/double bonds? See if everything is sp2 hybridized?
        # 3) 4n+2 rule: for every double bond, add 2 pi electrons. for every negative formal charge on an atom, add 2 pi electrons 
# What does an inchi look like when aromaticity is present in the molecule?
# ***Answer***: reading aromaticity from InchI must be inferred from whatever software is writing a molecule from
    # the given InchI string. In other words, InchI does not give bond order info nor does it specifiy which atoms 
    # are aromatic explicitly. I'm guessing RDKit differentiates between an aromatic and nonaromatic molecule by
    # perceiving bond orders from degrees of unsaturation, general structure, etc. 
# let's compare cyclohexane and benzene (note that only RDKit can read in an InChI string):
for mol in mols:
    if mol['iupac'] == 'cyclohexane':
        cyclohexane = Chem.inchi.MolFromInchi(mol['inchi'])
    if mol['iupac'] == 'benzene':
        benzene = Chem.inchi.MolFromInchi(mol['inchi'])
print("\n\ncyclohexane      benzene \n atoms            atoms")
for c_atom, b_atom in zip(cyclohexane.GetAtoms(), benzene.GetAtoms()):
    print(f"{c_atom.GetSymbol()}:", end='')
    # Testing to see if it could recognize the aromatic atoms
    print(f"({c_atom.GetIsAromatic()})", end="\t")

    print(f"{b_atom.GetSymbol()}:", end='')
    # Testing to see if it could recognize the aromatic atoms
    print(f"({b_atom.GetIsAromatic()})", end="\n")
# so yes, rdkit somehow perceives aromaticity from an InchI string format 


# How does InChI represent resonance states? And formal charge? Eg guandinium or arginine?
# InChI has an optional charge layer that can be used to store OVERALL charge info. In guanidinium, the charge
# layer can be seen with the "/p+1" tag: 
    # guanidinium: InChI=1S/CH5N3/c2-1(3)4/h(H5,2,3,4)/p+1
# compare this to SMILES which specifies charge info for each atom:
    # guanidinium: C(=[NH2+])(N)N   
    #                     ^
    # notice the "+" and specified 2 protons on the given nitrogen 
# InChI does not store formal charge or bond orders, so it is a logical conclusion (aka "I haven't looked this up but it seems obvious") 
# that InChI cannot store resonance states like SMILES 

# Does the same molecule always make the same inchi, or could it be represented by multiple InChIs?
    # Inchi keys have multiple layers that can store information, and the programmer can choose to what detail--how
    # many levels--he/she chooses to represent their molecule with. The only required levels are the chemical formula
    # level, and the connectivity level. Additional levels such as charge, explicit protons, isotopes can be used or 
    # not depending on the level of detail needed. 

    # Excluding isotopes, strange charges, stereochem, etc. an InChI key will only ever represent one molecule. It is
    # impossible to construct two different molecules from the same InChI string (again, not taking stereochemistry, isomers,
    # etc. into account)
# What is canonicalization? 
    # Standardizing the format of some data so that everyone can be on the same page. This is done to ease the transition 
    # of using data from one website, database, software, etc. to another. It would be a pain if you had to manually reformat 
    # your data every time you had to search for a molecule of use a different drug database
# Is there a canonical InChI?
    # Yes: Standard InChI. This is denoted by the "S" in "InChI=1S". For searching/database purposes, you can use the InChIKey 
    # format. 

# Load MiniDrugBank into RDKit and OpenEye from SDF and MOL2 formats:
#  https://github.com/openforcefield/openforcefield/tree/master/openforcefield/data/molecules

print("\n\n")
path2DrugBank = '/home/coda3831/anaconda3/envs/openff-dev/openforcefield/openforcefield/data/molecules'
drug_bank = Path(path2DrugBank)
failed = []
mol2 = []   # list of rdkit mols and openeye OEGraphs
sdf = []    # list of rdkit mols and openeye OEGraphs
for f in drug_bank.iterdir():
    try:
        if (f.suffix == '.mol2'):
            continue
            # read in the mol2 format to rdkit and openeye
            rdkit_list = []
            rdkit_mol = Chem.rdmolfiles.MolFromMol2File(str(f))   
                # The parser expects the atom-typing scheme used by Corina. Atom types from Tripos’ dbtranslate are less supported
                # Could only get this to work with mol files that contained a single molecule 
            rdkit_list.append(rdkit_mol)

            # now try openeye's file reading 
            # in future: use openeye with mol2 
            ifs = oechem.oemolistream(".mol2")
            ifs.open(str(f))
            openeye_list = []
            for mol in ifs.GetOEGraphMols():
                openeye_list.append(mol)
            # now we have two lists full of all the molecules rdkit and openeye could get from the mol file
            new_file = {'file': f.stem ,'rdkit': rdkit_list, 'openeye': openeye_list}
            mol2.append(new_file)
                
        elif(f.suffix == '.sdf'):
            # read in the sdf format to rdkit and openeye
            rdkit_list = []
            suppl = Chem.rdmolfiles.ForwardSDMolSupplier(str(f))

            for mol in suppl:
                rdkit_list.append(mol)

            # now try openeye's file reading 
            ifs = oechem.oemolistream(".sdf")
            ifs.open(str(f))
            openeye_list = []
            for mol in ifs.GetOEGraphMols():
                openeye_list.append(mol)
            # now we have two lists full of all the molecules rdkit and openeye could get from the mol file
            new_file = {'file': f.stem ,'rdkit': rdkit_list, 'openeye': openeye_list}
            sdf.append(new_file)
        else:
            continue 
    except Exception as e:
        failed.append(f)
# TODO:
# 1) figure out how to read a mol2 file with multiple molecules into rdkit
# 2) solve or accound for all errors that occur while reading in files
# 3) maybe implement a dataframe or something to help organize files/molecules
# 4) go through and test that all molecules were read in correctly
# 5) unpack or write to a compressed gz (learn gzip moldule)file
# 6) MiniDrugBank.sdf: do rdkit and openeye read all these molecules the same?
# 7) - how to compare molecules that are read into rdkit vs. openeye?
# 8) - compare, graphs(stereo, aromaticity, connectivity, atoms), partial charge, isomers(?), conformers etc.
# 9) - better coding pratices, less experimental code, more OO programming, get a final product
# 10) - make the software very readable, easy to follow the logical flow.
# 11) - make functions first, think about classes later 
# 12) - what are the details of a jupyter(?) notebook. 


test = 1
for failure in failed:
    print(f"file {failure.name} failed to be added")

