from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
IPythonConsole.ipython_useSVG=False


InChI_String = "InChI=1S/C6H8O6/c7-1-2(8)5-3(9)4(10)6(11)12-5/h2,5,7-8,10-11H,1H2/t2-,5+/m0/s1"
# InChI_String = "InChI=1S/C16H20/c1-11(2)13-5-7-16-10-14(12(3)4)6-8-15(16)9-13/h5-12H,1-4H3"
mol = Chem.inchi.MolFromInchi(InChI_String, sanitize=False, removeHs=False)
test =1
# Connectivity
# Atom names / chemical structure
for atom in mol.GetAtoms():
    print(f"{atom.GetSymbol()}({atom.GetMass()}):", end='')
    for neighbor in atom.GetNeighbors():
        print(f"->{neighbor.GetAtomicNum()} ({neighbor.GetIdx()})", end='')
    print("")
    # Testing to see if it recognizes aromaticies 
    # TODO: test on aromatic molecule
    print(f"    Is Aromatic?: {atom.GetIsAromatic()}")
    # stereochem
    print(f"    chirality: {atom.GetChiralTag().name}")
    
# Draw the molecule
Draw.MolToFile(mol, "openforcefield/openforcefield/Connor/test.png")


# InchI Key
Key = Chem.inchi.MolToInchiKey(mol)
print(Key)

# now for openeye
from openeye import oechem
from openeye import oedepict


smiles = '''\
CCO
c1cnccc1'''

ims = oechem.oemolistream()
ims.SetFormat(oechem.OEFormat_SMI)
ims.openstring(smiles)

mols = []
mol = oechem.OEGraphMol()
for mol in ims.GetOEMols():
    mols.append(oechem.OEGraphMol(mol))
# Since openeye cannot read InChI, I'm just going to mess around with
# the give list of OEMol objects
for atom in mols[0].GetAtoms():
    print(atom.GetAtomicNum(), end="")
    for neighbor in atom.GetAtoms():
        print(f"->{neighbor.GetAtomicNum()}", end="")
    print("")
    print(f"  Formal Charge: {atom.GetFormalCharge()}")
    print(f"  Implicit H Count: {atom.GetImplicitHCount()}")
    # apparently you can set params too...
    if atom.GetAtomicNum() == 6:
        atom.SetIsotope(12)
    print(f"  Mass: {atom.GetIsotope()}")
    print(f"  Partial Charge: {atom.GetPartialCharge()}")
    print(f"  Hybridization: {atom.GetHyb()}")
    print(f"  Name: {atom.GetName()}")
    print(f"  Aromaticitiy: {atom.IsAromatic()}")
print ("Bond Info")
for bond in mols[0].GetBonds():
    print(bond.GetOrder())
# now try to generate and graph 2d coords for the first molecule
mol = mols[0]
oedepict.OEPrepareDepiction(mol)
disp = oedepict.OE2DMolDisplay(mol)

# TODO: Understand this

oedepict.OERenderMolecule("openforcefield/openforcefield/Connor/test2.png", disp)

# now write to InChI format
oms = oechem.oemolostream()
oms.SetFormat(oechem.OEFormat_INCHI)
oms.openstring()

for mol in mols:
    oechem.OEWriteMolecule(oms, mol)

molfile = oms.GetString()
print("MOL string\n", molfile.decode('UTF-8'))
