# Author: Connor Davel
# Date Last modified: 1/9/20
# 
# This script will serve as the beginnings of hopefully a future jupyter notebook about loading openFF mols from 3D.
# Much of this is adapted/expanded from Jeff's notebook made a year ago to answer the same question: 
# https://github.com/openforcefield/openforcefield/blob/master/examples/deprecated/rdkit_openeye_comparison/compare_rdk_oe_mol_loading.ipynb
# Steps: 
# 1) [:D] Load in sets of molecules from 3D using OE and RDKit wrapper with the OFFMol.from_file() function
#    This pandas dataFrame will have 2 columns in the form:
#                 OFF_from_OE  |  OFF_from_RD
#               1|     ...     |     ...
#               2|     ...     |     ...
#              ...     ...     |     ...
# 2) [:D] Function 1 - OFF to SMILES string then string comparison
# 3) [:D] Function 2 - OFF to InChI string then string comparison
# 4) [:D] Function 3 - are_aromatic() from existing toolkit code
# 5) [:D] Function 4 - visual graph of NetworkX graphs with basic visual markers for different bonds/atoms
# 6) [X] calculate and store comparison results in another pandas dataframe (true if comparison function says that the two OFF molecules are
#    the same, false otherwise). Export this dataframe to a 

# STEP #1:
# It is important that this only be tested with an SDF files so that we do not introduce possible errors associated with 
# converting from mol2 to SDF or from SDF to mol2. 

# ????????????????? Since the RDKit wrapper can only read "SDF", "MOL", or "SMI", I'm not sure
# if PDF files will be viable for this project

# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import operator
import uuid
import warnings
from collections import Counter, OrderedDict
from copy import deepcopy
from io import StringIO
import os

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism import GraphMatcher
from simtk import unit
from simtk.openmm.app import Element, element

import openff.toolkit
from openff.toolkit.utils import (
    MessageException,
    check_units_are_compatible,
    deserialize_numpy,
    quantity_to_string,
    serialize_numpy,
    string_to_quantity,
)
from openff.toolkit.utils.serialization import Serializable
from openff.toolkit.utils.toolkits import (
    DEFAULT_AROMATICITY_MODEL,
    GLOBAL_TOOLKIT_REGISTRY,
    InvalidToolkitRegistryError,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
    ToolkitWrapper,
    UndefinedStereochemistryError,
)

from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper, UndefinedStereochemistryError
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.topology.molecule import Molecule, FrozenMolecule
import pandas as pd
import difflib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# following is needed so that the child classes bellow run with all the libraries found in toolkits.py
import logging 
logger = logging.getLogger(__name__)

# for testing purposes, I will be using new classes that inherit from RDKitToolkitWrapper and OpenEyeToolkitWrapper. This way I can
# override class functions and make some minimal modifications. 

class TestRDKitToolkitWrapper(RDKitToolkitWrapper):
    # An exact copy from original except for some minor changes:
    # When a file is not read correctly, instead of continuing to the next molecule, append a "None" type 
    # to the mols list. This ensures that both Toolkits return exactly 371 things, even if some are "None" 
    def from_file(
        self, file_path, file_format, allow_undefined_stereo=False, _cls=None
    ):
        """
        Create an openff.toolkit.topology.Molecule from a file using this toolkit.



        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.
        _cls : class
            Molecule constructor
        Returns
        -------
        molecules : iterable of Molecules
            a list of Molecule objects is returned.

        """
        from rdkit import Chem

        file_format = file_format.upper()

        mols = list()
        if (file_format == "MOL") or (file_format == "SDF"):
            for rdmol in Chem.SupplierFromFilename(
                file_path, removeHs=False, sanitize=False, strictParsing=True
            ):
                if rdmol is None:
                    mols.append(np.NaN)     # --------------------------------------------> added from original
                    continue

                # Sanitize the molecules (fails on nitro groups)
                try:
                    Chem.SanitizeMol(
                        rdmol,
                        Chem.SANITIZE_ALL
                        ^ Chem.SANITIZE_SETAROMATICITY
                        ^ Chem.SANITIZE_ADJUSTHS,
                    )
                    Chem.AssignStereochemistryFrom3D(rdmol)
                except ValueError as e:
                    logger.warning(rdmol.GetProp("_Name") + " " + str(e))
                    mols.append(np.NaN)   #-------------------------------------------------> added from original
                    continue
                Chem.SetAromaticity(rdmol, Chem.AromaticityModel.AROMATICITY_MDL)
                mol = self.from_rdkit(
                    rdmol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
                )
                mols.append(mol)

        elif file_format == "SMI":
            # TODO: We have to do some special stuff when we import SMILES (currently
            # just adding H's, but could get fancier in the future). It might be
            # worthwhile to parse the SMILES file ourselves and pass each SMILES
            # through the from_smiles function instead
            for rdmol in Chem.SmilesMolSupplier(file_path, titleLine=False):
                rdmol = Chem.AddHs(rdmol)
                mol = self.from_rdkit(
                    rdmol, allow_undefined_stereo=allow_undefined_stereo, _cls=_cls
                )
                mols.append(mol)

        elif file_format == "PDB":
            raise Exception(
                "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                "is likely to be lost. To read a PDB using RDKit use Molecule.from_pdb_and_smiles()"
            )
            # TODO: See if we can implement PDB+mol/smi combinations to get complete bond information.
            #  testing to see if we can make a molecule from smiles and then use the PDB conformer as the geometry
            #  and just reorder the molecule
            # https://github.com/openforcefield/openff-toolkit/issues/121
            # rdmol = Chem.MolFromPDBFile(file_path, removeHs=False)
            # mol = Molecule.from_rdkit(rdmol, _cls=_cls)
            # mols.append(mol)
            # TODO: Add SMI, TDT(?) support

        return mols
class TestOpenEyeToolkitWrapper(OpenEyeToolkitWrapper):
    def from_file(
        self, file_path, file_format, allow_undefined_stereo=False, _cls=None
    ):
        """
        Return an openff.toolkit.topology.Molecule from a file using this toolkit.

        Parameters
        ----------
        file_path : str
            The file to read the molecule from
        file_format : str
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for details.
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.
        _cls : class
            Molecule constructor

        Returns
        -------
        molecules : List[Molecule]
            The list of ``Molecule`` objects in the file.

        Raises
        ------
        GAFFAtomTypeWarning
            If the loaded mol2 file possibly uses GAFF atom types, which
            are not supported.

        Examples
        --------

        Load a mol2 file into an OpenFF ``Molecule`` object.

        >>> from openff.toolkit.utils import get_data_file_path
        >>> mol2_file_path = get_data_file_path('molecules/cyclohexane.mol2')
        >>> toolkit = OpenEyeToolkitWrapper()
        >>> molecule = toolkit.from_file(mol2_file_path, file_format='mol2')

        """
        from openeye import oechem

        ifs = oechem.oemolistream(file_path)
        return self._read_oemolistream_molecules(
            ifs, allow_undefined_stereo, file_path=file_path, _cls=_cls
        )
class TestMolecule(Molecule):
    # this is the exact same function but wrapped in a try statement to help debug 
    @staticmethod
    def are_isomorphic(
        mol1,
        mol2,
        return_atom_map=False,
        aromatic_matching=True,
        formal_charge_matching=True,
        bond_order_matching=True,
        atom_stereochemistry_matching=True,
        bond_stereochemistry_matching=True,
        strip_pyrimidal_n_atom_stereo=True,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Determines whether the two molecules are isomorphic by comparing their graph representations and the chosen
        node/edge attributes. Minimally connections and atomic_number are checked.

        If nx.Graphs() are given they must at least have atomic_number attributes on nodes.
        other optional attributes for nodes are: is_aromatic, formal_charge and stereochemistry.
        optional attributes for edges are: is_aromatic, bond_order and stereochemistry.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mol1 : an openforcefield.topology.molecule.FrozenMolecule or TopologyMolecule or nx.Graph()
        mol2 : an openforcefield.topology.molecule.FrozenMolecule or TopologyMolecule or nx.Graph()
            The molecule to test for isomorphism.

        return_atom_map: bool, default=False, optional
            will return an optional dict containing the atomic mapping.

        aromatic_matching: bool, default=True, optional
            compare the aromatic attributes of bonds and atoms.

        formal_charge_matching: bool, default=True, optional
            compare the formal charges attributes of the atoms.

        bond_order_matching: bool, deafult=True, optional
            compare the bond order on attributes of the bonds.

        atom_stereochemistry_matching : bool, default=True, optional
            If ``False``, atoms' stereochemistry is ignored for the
            purpose of determining equality.

        bond_stereochemistry_matching : bool, default=True, optional
            If ``False``, bonds' stereochemistry is ignored for the
            purpose of determining equality.

        strip_pyrimidal_n_atom_stereo: bool, default=True, optional
            If ``True``, any stereochemistry defined around pyrimidal
            nitrogen stereocenters will be disregarded in the isomorphism
            check.

        toolkit_registry : openforcefield.utils.toolkits.ToolkitRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for
            removing stereochemistry from pyrimidal nitrogens.

        Returns
        -------
        molecules_are_isomorphic : bool

        atom_map : default=None, Optional,
            [Dict[int,int]] ordered by mol1 indexing {mol1_index: mol2_index}
            If molecules are not isomorphic given input arguments, will return None instead of dict.
        """

        try:
            # Do a quick hill formula check first
            if Molecule.to_hill_formula(mol1) != Molecule.to_hill_formula(mol2):
                return False, None

            # Build the user defined matching functions
            def node_match_func(x, y):
                # always match by atleast atomic number
                is_equal = x["atomic_number"] == y["atomic_number"]
                if aromatic_matching:
                    is_equal &= x["is_aromatic"] == y["is_aromatic"]
                if formal_charge_matching:
                    is_equal &= x["formal_charge"] == y["formal_charge"]
                if atom_stereochemistry_matching:
                    is_equal &= x["stereochemistry"] == y["stereochemistry"]
                return is_equal

            # check if we want to do any bond matching if not the function is None
            if aromatic_matching or bond_order_matching or bond_stereochemistry_matching:

                def edge_match_func(x, y):
                    # We don't need to check the exact bond order (which is 1 or 2)
                    # if the bond is aromatic. This way we avoid missing a match only
                    # if the alternate bond orders 1 and 2 are assigned differently.
                    if aromatic_matching and bond_order_matching:
                        is_equal = (x["is_aromatic"] == y["is_aromatic"]) or (
                            x["bond_order"] == y["bond_order"]
                        )
                    elif aromatic_matching:
                        is_equal = x["is_aromatic"] == y["is_aromatic"]
                    elif bond_order_matching:
                        is_equal = x["bond_order"] == y["bond_order"]
                    else:
                        is_equal = None
                    if bond_stereochemistry_matching:
                        if is_equal is None:
                            is_equal = x["stereochemistry"] == y["stereochemistry"]
                        else:
                            is_equal &= x["stereochemistry"] == y["stereochemistry"]

                    return is_equal

            else:
                edge_match_func = None

            # Here we should work out what data type we have, also deal with lists?
            def to_networkx(data):
                """For the given data type, return the networkx graph"""
                from openff.toolkit.topology import TopologyMolecule

                if strip_pyrimidal_n_atom_stereo:
                    SMARTS = "[N+0X3:1](-[*])(-[*])(-[*])"

                if isinstance(data, FrozenMolecule):
                    # Molecule class instance
                    if strip_pyrimidal_n_atom_stereo:
                        # Make a copy of the molecule so we don't modify the original
                        data = deepcopy(data)
                        data.strip_atom_stereochemistry(
                            SMARTS, toolkit_registry=toolkit_registry
                        )
                    return data.to_networkx()
                elif isinstance(data, TopologyMolecule):
                    # TopologyMolecule class instance
                    if strip_pyrimidal_n_atom_stereo:
                        # Make a copy of the molecule so we don't modify the original
                        ref_mol = deepcopy(data.reference_molecule)
                        ref_mol.strip_atom_stereochemistry(
                            SMARTS, toolkit_registry=toolkit_registry
                        )
                    return ref_mol.to_networkx()
                elif isinstance(data, nx.Graph):
                    return data

                else:
                    raise NotImplementedError(
                        f"The input type {type(data)} is not supported,"
                        f"please supply an openff.toolkit.topology.molecule.Molecule,"
                        f"openff.toolkit.topology.topology.TopologyMolecule or networkx "
                        f"representation of the molecule."
                    )

            mol1_netx = to_networkx(mol1)
            mol2_netx = to_networkx(mol2)
            from networkx.algorithms.isomorphism import GraphMatcher

            GM = GraphMatcher(
                mol1_netx, mol2_netx, node_match=node_match_func, edge_match=edge_match_func
            )
            isomorphic = GM.is_isomorphic()

            if isomorphic and return_atom_map:
                topology_atom_map = GM.mapping

                # reorder the mapping by keys
                sorted_mapping = {}
                for key in sorted(topology_atom_map.keys()):
                    sorted_mapping[key] = topology_atom_map[key]

                return isomorphic, sorted_mapping

            else:
                return isomorphic, None
        except Exception:
            return -99


# 2) Function 1 - OFF to SMILES string then string comparison
def get_str_symbols(string):
    """
    returns a dictionary with the symbols as keys and the number of occurences of each symbol as the values from str. 
    This function does not work well for formulas with Ni, Ts, Nb, and some other elements that will not ever be used
    """
    letters = dict()
    ignore = "[]()123456789-+=\?#$:.Hh@"
    for i in range(len(string)):
        f = string[i:i+1]
        n = string[i+1:i+2]
        if f in ignore:
            continue
        if f.isupper():
            if n in 'ltre':
                i_str = f + n
            else:
                i_str = f
        elif f in 'bcnops':
            i_str = f
        else:
            continue

        if i_str in letters.keys():
            letters[i_str] += 1
        else:
            letters[i_str] = 1

    return letters

def smiles_str_compare(mol1, mol2, toolkit_wrapper, isom=True, explicit_h=True):
    """
    converts the molecules to a SMILES string and then does a string comparison 
    mol1: an OFF Molecule class 
    mol2: an OFF Molecule class
    toolkit_wrapper: either RDKitToolkitWrapper or OpenEyeToolkitWrapper (hopefully testing will show that this distinction does not matter)
    returns: 1 if SMILES are identical, 0 if SMILES have the same number of heavy atoms (everything except H), -1 else 
    """
    try:
        mol1_str = mol1.to_smiles(toolkit_registry=toolkit_wrapper, isomeric=isom, explicit_hydrogens=explicit_h)
        mol2_str = mol2.to_smiles(toolkit_registry=toolkit_wrapper, isomeric=isom, explicit_hydrogens=explicit_h)

        if mol1_str == mol2_str:
            return 1
        # now we start the comparison of all the heavy atoms. Hopefully this will tell us if explicit H's are the cause of the string inequality
        mol1_symbols = get_str_symbols(mol1_str)
        mol2_symbols = get_str_symbols(mol2_str)
        if mol1_symbols == mol2_symbols:
            return 0

        return -1
    except Exception:
        return -999

# 3) Function 2 - OFF to InChI string then string comparison
def parse_std_InChI(string):
    """
    Does a basic parsing of a standard InChI string using what rules I currently know.
    Returns [atom_block, connectivity_block, explicity_hydrogen_block]
    """

    layers = string.split('/')
    atoms = None
    connect = None
    hydro = None
    for layer in layers:
        if 'InChI=' in layer:
            continue    
        if layer[0:1].isupper():
            atoms = layer
        if layer[0:1] == 'c':
            connect = layer
        if layer[0:1] == 'h':
            hydro = layer
    return [atoms, connect, hydro]
    
def inchi_str_compare(mol1, mol2, toolkit_wrapper, fixed_h = False):
    """
        converts the molecules to an InChI string and then does a string comparison 
        mol1: an OFF Molecule class 
        mol2: an OFF Molecule class
        toolkit_wrapper: either RDKitToolkitWrapper or OpenEyeToolkitWrapper (hopefully testing will show that this distinction does not matter)
        returns: 3 if InChI are identical, 
                 2 if InChI have the same atoms and connectivity and explicit hydrogen block, 
                 1 if InChI have the same atoms and connectivity
                 0 if InChI have the same atoms
                 -1 else
    """
    try:
        mol1_str = mol1.to_inchi(fixed_hydrogens=fixed_h, toolkit_registry=toolkit_wrapper)
        mol2_str = mol2.to_inchi(fixed_hydrogens=fixed_h, toolkit_registry=toolkit_wrapper)

        if mol1_str == mol2_str:
            return 3
        # now we start the comparison of all the heavy atoms. Hopefully this will tell us if explicit H's are the cause of the string inequality
        mol1_blocks = parse_std_InChI(mol1_str)
        mol2_blocks = parse_std_InChI(mol2_str)
        if mol1_blocks[0] == mol2_blocks[0]:  #compares the atoms block
            if mol1_blocks[1] == mol2_blocks[1]:  #compares the connectivity block
                if mol1_blocks[2] == mol2_blocks[2]:  #compares the explicity hydrgen block
                    return 2
                return 1
            return 0
        else:
            return -1
    except Exception:
        return -999

# 5) Function 4 - visual graph of NetworkX graphs with basic visual markers for different bonds/atoms
def graph_networkx(mol1, folder, toolkit_wrapper, strip_pyrimidal_n_atom_stereo=False):
    # the following function comes straight from the are_isomorphic function from the FrozenMolecule class (see line 2595 molecule.py)
    def to_networkx(data):
            """For the given data type, return the networkx graph"""

            if strip_pyrimidal_n_atom_stereo:
                SMARTS = "[N+0X3:1](-[*])(-[*])(-[*])"

            if isinstance(data, FrozenMolecule):
                # Molecule class instance
                if strip_pyrimidal_n_atom_stereo:
                    # Make a copy of the molecule so we don't modify the original
                    data = deepcopy(data)
                    data.strip_atom_stereochemistry(
                        SMARTS, toolkit_registry=toolkit_wrapper
                    )
                return data.to_networkx()

            else:
                raise NotImplementedError(
                    f"The input type {type(data)} is not supported,"
                    f"please supply an openff.toolkit.topology.molecule.Molecule,"
                )

    mol1_netx = to_networkx(mol1)

    color_map = {1: 'y',
                6: 'k',
                7: 'b',
                8: 'r',
                15: 'm',
                16: 'g'
                }
    size_map = {1: 20,
                6: 75,
                7: 150,
                8: 150,
                15: 200,
                16: 200
                }
    label_map = {1: "",
                6: "",
                7: "N",
                8: "O",
                15: "P",
                16: "S"
                }
    labeldict = dict()
    for i, val in mol1_netx.nodes(data=True):
        labeldict[i] = label_map.get(val.get('atomic_number'), f"{val.get('atomic_number')}")

    color_values = [color_map.get(val.get('atomic_number'), 'c') for i, val in mol1_netx.nodes(data=True)]
    size_values = [size_map.get(val.get('atomic_number'), 200) for i, val in mol1_netx.nodes(data=True)]
    single_bonds = [edge for edge in mol1_netx.edges(data=True) if (edge[2].get('bond_order') == 1)]
    double_bonds = [edge for edge in mol1_netx.edges(data=True) if (edge[2].get('bond_order') == 2)]
    triple_bonds = [edge for edge in mol1_netx.edges(data=True) if (edge[2].get('bond_order') == 3)]

    # # -----------for layout testing purposes, will eventually only use the best layout --------------
    # top = nx.bipartite.sets(mol1_netx)[0]
    # layouts = {
    #     'bipartite_layout': nx.bipartite_layout(mol1_netx, top),
    #     'circular_layout': nx.circular_layout(mol1_netx),
    #     'kamada_kawai_layout': nx.kamada_kawai_layout(mol1_netx),
    #     'planar_layout': nx.planar_layout(mol1_netx),
    #     'random_layout': nx.random_layout(mol1_netx),
    #     'shell_layout': nx.shell_layout(mol1_netx),
    #     'spring_layout': nx.spring_layout(mol1_netx),
    #     'spectral_layout': nx.spectral_layout(mol1_netx),
    #     'spiral_layout': nx.spiral_layout(mol1_netx)
    # }

    # for key in layouts:
    #     try:
    #         pos = layouts.get(key)
    #         nx.draw_networkx_nodes(mol1_netx, pos, node_color = color_values, node_size = 200)
    #         nx.draw_networkx_labels(mol1_netx, pos)
    #         nx.draw_networkx_edges(mol1_netx, pos, edgelist=single_bonds, edge_color='k')
    #         nx.draw_networkx_edges(mol1_netx, pos, edgelist=double_bonds, edge_color='r')
    #         nx.draw_networkx_edges(mol1_netx, pos, edgelist=triple_bonds, edge_color='c')
    #         plt.savefig(f"openforcefield/openforcefield/Connor/figures/layout_testing/{key}.png")
    #     except Exception as e:
    #         print(f"problem while plotting layout: {key}")
    #     plt.clf()

    # # ------------for kamada kawai layout testing. I want to make the graph will be plotted the same way every time. 
    # # This is way the spring_layout is not an option (It creates random graphs after each iteration) --------------

    # for i in range(1, 10):
    #     pos = nx.kamada_kawai_layout(mol1_netx)
    #     nx.draw_networkx_nodes(mol1_netx, pos, node_color = color_values, node_size = 200)
    #     nx.draw_networkx_labels(mol1_netx, pos)
    #     nx.draw_networkx_edges(mol1_netx, pos, edgelist=single_bonds, edge_color='k')
    #     nx.draw_networkx_edges(mol1_netx, pos, edgelist=double_bonds, edge_color='r')
    #     nx.draw_networkx_edges(mol1_netx, pos, edgelist=triple_bonds, edge_color='c')
    #     plt.savefig(f"openforcefield/openforcefield/Connor/figures/layout_testing/kamada_test_{i}.png")
    #     plt.clf()

    # --------------CONCLUSION-----------------
    # Will only use the kamada kawai layout from now on...
    plt.clf()
    fig, ax = plt.subplots(1, 1, num=1)

    pos = nx.kamada_kawai_layout(mol1_netx)
    nx.draw_networkx_nodes(mol1_netx, pos, node_color = color_values, node_size = size_values, ax=ax)
    nx.draw_networkx_labels(mol1_netx, pos, font_size = 7, alpha = 0.75, ax=ax, labels=labeldict)
    nx.draw_networkx_edges(mol1_netx, pos, edgelist=single_bonds, edge_color='k', ax=ax)
    nx.draw_networkx_edges(mol1_netx, pos, edgelist=double_bonds, edge_color='r', ax=ax)
    nx.draw_networkx_edges(mol1_netx, pos, edgelist=triple_bonds, edge_color='c', ax=ax)
    plt.savefig(f"figures/{folder}/{mol1.name}.png")
    plt.clf()

    # test = 1
    #  "openforcefield/openforcefield/Connor/test2.png"

def graph_mult_networkx(mol_list, name_list, toolkit_wrapper, folder, strip_pyrimidal_n_atom_stereo=False):
    # the following function comes straight from the are_isomorphic function from the FrozenMolecule class (see line 2595 molecule.py)
    def to_networkx(data):
            """For the given data type, return the networkx graph"""

            if strip_pyrimidal_n_atom_stereo:
                SMARTS = "[N+0X3:1](-[*])(-[*])(-[*])"

            if isinstance(data, FrozenMolecule):
                # Molecule class instance
                if strip_pyrimidal_n_atom_stereo:
                    # Make a copy of the molecule so we don't modify the original
                    data = deepcopy(data)
                    data.strip_atom_stereochemistry(
                        SMARTS, toolkit_registry=toolkit_wrapper
                    )
                return data.to_networkx()

            else:
                raise NotImplementedError(
                    f"The input type {type(data)} is not supported,"
                    f"please supply an openff.toolkit.topology.molecule.Molecule,"
                )

    filename = ""
    for mol in mol_list:
        filename = filename + mol.name + "_"

    color_map = {1: 'y',
                6: 'k',
                7: 'b',
                8: 'r',
                15: 'm',
                16: 'g'
                }
    size_map = {1: 20,
                6: 75,
                7: 150,
                8: 150,
                15: 200,
                16: 200
                }
    label_map = {1: "",
                6: "",
                7: "N",
                8: "O",
                15: "P",
                16: "S"
                }
    plt.clf()
    fig, ax = plt.subplots(len(mol_list), 1, num=1)
    plot = 0
    for mol in mol_list:
        mol_netx = to_networkx(mol)

        labeldict = dict()
        for i, val in mol_netx.nodes(data=True):
            labeldict[i] = label_map.get(val.get('atomic_number'), f"{val.get('atomic_number')}")

        color_values = [color_map.get(val.get('atomic_number'), 'c') for i, val in mol_netx.nodes(data=True)]
        size_values = [size_map.get(val.get('atomic_number'), 200) for i, val in mol_netx.nodes(data=True)]
        single_bonds = [edge for edge in mol_netx.edges(data=True) if (edge[2].get('bond_order') == 1)]
        double_bonds = [edge for edge in mol_netx.edges(data=True) if (edge[2].get('bond_order') == 2)]
        triple_bonds = [edge for edge in mol_netx.edges(data=True) if (edge[2].get('bond_order') == 3)]
        
        pos = nx.kamada_kawai_layout(mol_netx)

        nx.draw_networkx_nodes(mol_netx, pos, node_color = color_values, node_size = size_values, ax=ax[plot])
        nx.draw_networkx_labels(mol_netx, pos, font_size = 7, alpha = 0.75, ax=ax[plot], labels=labeldict)
        nx.draw_networkx_edges(mol_netx, pos, edgelist=single_bonds, edge_color='k', ax=ax[plot])
        nx.draw_networkx_edges(mol_netx, pos, edgelist=double_bonds, edge_color='r', ax=ax[plot])
        nx.draw_networkx_edges(mol_netx, pos, edgelist=triple_bonds, edge_color='c', ax=ax[plot])
        ax[plot].title.set_text(name_list[plot])

        plot += 1
    plt.savefig(f"figures/{folder}/{filename}.png")
    plt.clf()

# 6) calculate and store comparison results in another pandas dataframe (true if comparison function says that the two OFF molecules are
#    the same, false otherwise). Export this dataframe to a csv

test_output = False
def test_function(run_function, dframe):
    # trying to find why are_isomorphic returns False at one end of the script and True at the other
    # Will drop this function in a few places throughout the script to figure out what is going on
    if run_function:
        names = ['DrugBank_3817',
                'DrugBank_4032',
                'DrugBank_1971',
                'DrugBank_2140',
                'DrugBank_2563',
                'DrugBank_2585',
                'DrugBank_2687']
        test = dframe[dframe['name'].isin(names)]
        for index, row in test.iterrows():
            m1 = row[0]
            m2 = row[1]
            a, b = TestMolecule.are_isomorphic(m1, m2)
            print(a, end="\t")
        print("")

def run_comparisons():
    """
    this is just a function to wrap the comparison calculations and organize the workspace. 
    """
    def unwrapper_func(are_isomorphic_output):
        """
        This is a helper function that helps unwrap the output from the are_isomorphic function in the TestMolecule class
        """
        if are_isomorphic_output == -99:
            return -99
        else:
            f, l = are_isomorphic_output
            return f
    def compare_func(mol1, mol2):
        """
        another helper function. Not sure if these list comps are even more efficient at this point but whatever, I've dig this deep already... 
        """
        if isinstance(mol1, Molecule) and isinstance(mol2, Molecule):
            return True
        else:
            if isinstance(mol2, float):
                if mol2 != mol2:
                    return False
                elif mol2 < -10:
                    return True
                else:
                    print("unexplected case in compare_func")
                    return False
            elif isinstance(mol1, float):
                if mol1 != mol1:
                    return False
                elif mol1 < -10:
                    return True
                else:
                    print("unexplected case in compare_func")
                    return False

    
    comparison_results = pd.DataFrame(columns=['name', 'smi_rd', 'smi_oe', 'inch_rd', 'inch_oe', 'iso_arom', 'iso_bond_order', 'iso_arom_bond_order'])
    comparison_results['name'] = [mol.name for mol in mols['OFF_from_OE']]

    comparison_results['smi_rd'] = [smiles_str_compare(mol_oe, mol_rd, toolkit_wrapper=rdkit_toolkit_wrapper) if (compare_func(mol_oe, mol_rd)) else "nan" \
             for mol_oe, mol_rd in zip(mols['OFF_from_OE'], mols['OFF_from_RD'])]
    print("1")
    comparison_results['smi_oe'] = [smiles_str_compare(mol_oe, mol_rd, toolkit_wrapper=open_eye_toolkit_wrapper) if (compare_func(mol_oe, mol_rd)) else "nan" \
              for mol_oe, mol_rd in zip(mols['OFF_from_OE'], mols['OFF_from_RD'])]
    print("2")
    comparison_results['inch_rd'] = [inchi_str_compare(mol_oe, mol_rd, toolkit_wrapper=rdkit_toolkit_wrapper) if (compare_func(mol_oe, mol_rd)) else "nan" \
             for mol_oe, mol_rd in zip(mols['OFF_from_OE'], mols['OFF_from_RD'])]
    print("3")
    comparison_results['inch_oe'] = [inchi_str_compare(mol_oe, mol_rd, toolkit_wrapper=open_eye_toolkit_wrapper) if (compare_func(mol_oe, mol_rd)) else "nan" \
              for mol_oe, mol_rd in zip(mols['OFF_from_OE'], mols['OFF_from_RD'])]
    print("4")

    comparison_results['iso_arom'] = [unwrapper_func(TestMolecule.are_isomorphic(mol_oe, mol_rd, bond_order_matching=False)) if (compare_func(mol_oe, mol_rd)) else "nan" \
             for mol_oe, mol_rd in zip(mols['OFF_from_OE'], mols['OFF_from_RD'])]
    print("5")
    comparison_results['iso_bond_order'] = [unwrapper_func(TestMolecule.are_isomorphic(mol_oe, mol_rd, aromatic_matching=False)) if (compare_func(mol_oe, mol_rd)) else "nan" \
                for mol_oe, mol_rd in zip(mols['OFF_from_OE'], mols['OFF_from_RD'])]
    print("6")
    comparison_results['iso_arom_bond_order'] = [unwrapper_func(TestMolecule.are_isomorphic(mol_oe, mol_rd)) if (compare_func(mol_oe, mol_rd)) else "nan" \
              for mol_oe, mol_rd in zip(mols['OFF_from_OE'], mols['OFF_from_RD'])]
    print("7")
    comparison_results.to_csv("comparison_results.csv")
    return comparison_results

class CustomSDFMolSupplier:
    """
    returns each individual SDF block from an SDF file. So far, only tested and working on
    MiniDrugBank.sdf. With the nature of the file parsing, errors are extremely likely. It would
    be nice if I could get ahold of some "official" code that does this. Oh well 
    """
    def __init__(self, file):
        
        self.file = self.check_file(file)
        self.at_end_of_file = False
        self.file_length = sum(1 for line in open(file, "r"))
        self.is_sdf()  # checks for file identity
        self.current_mol = 0

    def check_file(self, file):
        try:
            from pathlib import Path
        except Exception as e:
            print(f"Exception during import in CustomSDFMolSupplier: \n\n {e} \n")
        if isinstance(file, str):
            f = Path(file)
        elif isinstance(file, Path):
            f = file
        else:
            print(f"error in CustomSDFMolSupplier:\n\t type {type(file)} given, which is not accepted. Use a string or Path variable")
        return f    
    
    def is_sdf(self):
        if self.file.suffix != '.sdf':
            raise NameError("sdf file required for class \"CustomSDFMolSupplier\"")

    def get_next_block(self):
        """
        read from the current file index "curr" to the start of the next molecule
        """
        f = open(self.file, "r")

        sdfstring = ""     # this is only defined here to that the top-most molecule does not break
        for i, line in enumerate(f):
            # parse through until the correct line has been reached
            if i < self.current_mol:
                continue
            # if we have reached the end of the current molecule, break
            if line == "$$$$\n":
                if (i + 1 == self.file_length): # if line is the last line of the file
                    self.at_end_of_file = True
                    self.current_mol = i + 1
                    break
                elif i == self.current_mol:   #if at start of the molecule
                    sdfstring = ""
                    continue
                else:                       #if at the end
                    self.current_mol = i
                    break
            if i == self.current_mol + 1 and i != 1:   #get the identifier of the molecule
                name = line.split("\n")[0]  #remove newline 
            elif i == 0:         # special case when reading the first block of the file
                name = line.split("\n")[0]
            sdfstring = "".join((sdfstring, line))
        return name, sdfstring

    def __iter__(self):
        self.current_mol = 0
        self.at_end_of_file = False
        return self
    def __next__(self):
        if self.at_end_of_file:
            raise StopIteration
        mol_name, mol2_block = self.get_next_block()
        return mol_name, mol2_block

class BlockToAtomSupplier(CustomSDFMolSupplier):
    """
    Does the same thing as CustomSDFMolSupplier but returns atom info instead of the raw SDF block
    Same applies here - lots of testing needs to be done, but these classes work for now. I'm not sure
    if I would be comfortable putting this in a Jupyter notebook and airing it to the programming world... 
    """
    def get_str_symbols(self, string):
        """
        returns a dictionary with the symbols as keys and the number of occurences of each symbol as the values from str. 
        This function does not work well for formulas with Ni, Ts, Nb, and some other elements that will not ever be used
        """
        letters = dict()
        ignore = "[]()123456789-+=\?#$:.@"
        for i in range(len(string)):
            f = string[i:i+1]
            n = string[i+1:i+2]
            if f in ignore:
                continue
            if f.isupper():
                if n in 'ltre':
                    i_str = f + n
                else:
                    i_str = f
            elif f in 'bcnops':
                i_str = f
            else:
                continue

            if i_str in letters.keys():
                letters[i_str] += 1
            else:
                letters[i_str] = 1

        return letters

    def parse_atoms(self, b):
        """
        returns a dictionary of atom names (key) and the corresponding number of atoms of that type.
        If you run into issues with this function, it probably means that the sdf block or file you are reading
        this from differs from MiniDrugBank.sdf by a single space or indetation difference. 
        """
        string = ""
        split_array = b.split("\n")
        found = False
        for line in split_array: # only want to get the "big block" of info
            if len(line) > 60:
                found = True
                string += line
            if len(line) < 40 and found:  # stop reading the rest when done 
                break
        return self.get_str_symbols(string)

    def __next__(self):
        if self.at_end_of_file:
            raise StopIteration
        mol_name, mol2_block = self.get_next_block()
        return mol_name, self.parse_atoms(mol2_block)

def get_atom_nums(mols):
    """
    I just discovered function Generators and this is the best thing ever 
    This func takes in mols and returns a generator object that supplies dictionaries with element symbols as keys and number of atoms for each element as the values.
    Mols is a list of OFF Molecules. Can also accept a pandas series. stops iteration and prints error if error occurs 
    """
    # begin with list and Series handling
    try:
        if isinstance(mols, list):
            pass
        elif isinstance(mols, pd.Series):
            mols = mols.to_list()
        else:
            print(f"improper Dtype <{type(mols)}> passed to get_atom_nums(mols)")
            raise StopIteration
    except Exception:
        print(f"improper Dtype <{type(mols)}> passed to get_atom_nums(mols)")
        raise StopIteration
    # define the function that operates on an individual OFFMol and returns a dictionary

    def OFF_dict(mol):
        atom_syms = [atom.element.symbol for atom in mol.atoms]
        atom_dict = dict()
        for a in atom_syms:
            if a in atom_dict.keys():
                atom_dict[a] += 1
            else:
                atom_dict[a] = 1
        return atom_dict

    for mol in mols:
        try:
            yield mol.name, OFF_dict(mol)
        except Exception:
            yield -999, "tears"

def run_manual_comparisons():
    """runs functions to (literally) count the number of atoms from the sdf file. This is supposed
        catch any errors that would be due to adding in protons. 
    """
    comparison_results = pd.DataFrame(columns=['name2', 'manual_reading', 'open_eye', 'rdkit', 'proton_num(man vs oe vs rd)','comparison'])
    for file_results, off_oe, off_rd in zip(BlockToAtomSupplier(drug_bank_file), get_atom_nums(mols['OFF_from_OE']), get_atom_nums(mols['OFF_from_RD'])):
        manual_name, manual_dict = file_results
        oe_name, oe_results = off_oe
        rd_name, rd_results = off_rd
        if manual_name == oe_name == rd_name:
            are_same = (manual_dict == oe_results == rd_results)
            row = pd.Series({'name2': manual_name, 
                   'manual_reading': manual_dict, 
                   'open_eye': oe_results, 
                   'rdkit': rd_results, 
                   'proton_num(man vs oe vs rd)': f"{manual_dict.get('H')} vs {oe_results.get('H')} vs {rd_results.get('H')}",
                   'comparison': are_same})
            comparison_results = comparison_results.append(row, ignore_index=True)
        else:
            row = pd.Series({'name2': manual_name, 
                   'manual_reading': oe_name, 
                   'open_eye': rd_name, 
                   'rdkit': "--", 
                   'proton_num(man vs oe vs rd)': "--",
                   'comparison': "--"})
            comparison_results = comparison_results.append(row, ignore_index=True)
            
    comparison_results.to_csv("manual_comparison_results.csv")
    return comparison_results

def graph_deviants():

    df1 = run_manual_comparisons()
    df2 = run_comparisons()
    combined = pd.concat([mols, df1, df2], axis=1)

    condition = {'deviant_mol_structures/arom_kekule_diff': (combined['comparison'] == True) & (combined['iso_bond_order'] == False) & (combined['iso_arom'] == True),
                  'deviant_mol_structures/large_kekule_diff': (combined['comparison'] == True) & (combined['iso_bond_order'] == False) & (combined['iso_arom'] == False),
                  'deviant_mol_structures/proton_or_atom_diff': (combined['comparison'] == False), 
                  'deviant_mol_structures/failed_to_load': (combined['smi_rd'] == 'nan')}
    for folder in condition.keys():
        error_mols = combined.loc[condition.get(folder), ['OFF_from_OE', 'OFF_from_RD']]
        for mol1, mol2 in zip(error_mols['OFF_from_OE'], error_mols['OFF_from_RD']):
            # try:
            #     if isinstance(mol1, Molecule) and isinstance(mol2, Molecule):
            #         if mol1.name != mol2.name:
            #             print("\n\n\n\nCritical Error in graph_deviants")
            #             return None
            #         else:
            #             graph_mult_networkx([mol1, mol2], ['OpenEye', 'RDKit'], toolkit_wrapper=rdkit_toolkit_wrapper, folder=folder)
            #     else: 
            #         if isinstance(mol1, Molecule):
            #             graph_networkx(mol1, folder=folder, toolkit_wrapper=rdkit_toolkit_wrapper)
            #         elif isinstance(mol2, Molecule):
            #             graph_networkx(mol2, folder=folder, toolkit_wrapper=rdkit_toolkit_wrapper)
            #         else:
            #             continue
            # except Exception:
            #     print(f"unkown exception in loop: {mol1} \t {mol2}")

            if isinstance(mol1, Molecule) and isinstance(mol2, Molecule):
                if mol1.name != mol2.name:
                    print("\n\n\n\nCritical Error in graph_deviants")
                    return None
                else:
                    graph_mult_networkx([mol1, mol2], ['OpenEye', 'RDKit'], toolkit_wrapper=rdkit_toolkit_wrapper, folder=folder)
            else: 
                if isinstance(mol1, Molecule):
                    graph_networkx(mol1, folder=folder, toolkit_wrapper=rdkit_toolkit_wrapper)
                elif isinstance(mol2, Molecule):
                    graph_networkx(mol2, folder=folder, toolkit_wrapper=rdkit_toolkit_wrapper)
                else:
                    continue


# ----------------------Week 3: Investigating stereochemistry, protonization, and loading failures ---------------------
# Overall goal: find out what steps in the toolkits are resulting in errors and differences between the resulting molecules 

def organize_errors():
    # organizes the previous results and adds a column naming the type of error identified.
    # ONLY returns the molecules that exhibit some kind or error. All "good" molecules are discarded. 
    df1 = run_manual_comparisons()
    df2 = run_comparisons()
    combined = pd.concat([mols, df1, df2], axis=1)
    # the columns added here will help organize the data into error types: 
        # arom_kekule_diff: apparaent differences in the bond order around aromatic rings
        # stereo_diff: differences in something other than bond order, stereo, etc. Probably stereochem around thiones
        # proton_diff: differences in the number of atoms, but all of these show up as protonization differences
        # failed_to_load: mols that rdkit failed to load. 
        # all of the above should be mutually exculsive, but we will see
    condition = {'arom_kekule_diff': (combined['comparison'] == True) & (combined['iso_bond_order'] == False) & (combined['iso_arom'] == True),
                  'stereo_diff': (combined['comparison'] == True) & (combined['iso_bond_order'] == False) & (combined['iso_arom'] == False),
                  'proton_diff': (combined['comparison'] == False), 
                  'failed_to_load': (combined['smi_rd'] == 'nan')}    
    
    for error_name in condition.keys():
        combined[error_name] = condition.get(error_name)
    # check one last time to make sure all names match. Not really necessary but oh well... 
    for i, row in combined.iterrows():
        if (not isinstance(row['OFF_from_RD'], float)) and (row['name2'] == row['name'] == row['OFF_from_OE'].name == row['OFF_from_RD'].name):
            continue
        elif (isinstance(row['OFF_from_RD'], float)) and (row['name2'] == row['name']):   # as a reminder, these cases indicate that RDKit failed to load in the molecule 
            continue 
        else:
            print('critical error in organize_errors: names do not match')
            raise Exception

    combined = combined[combined['arom_kekule_diff'] | combined['stereo_diff'] | combined['proton_diff'] | combined['failed_to_load']]

    test_function(test_output, combined)

    return combined


# the following functions group my attempts to investigate each type of error 

def stereo_diff():
    # since the differences in stereochemistry are hard to see from 2D or 3D graphs, the following code organizes the 
    # stereochemistry of atoms from the open_eye and rdkit molecules into dataframes and prints the results 
    # The form of the dataframe will be:
    #    name | oe_N | rd_N | oe_C | rd_C | oe_S | rd_S | oe_etc | rd_etc 
    #         |      |      |      |      |      |      |        |    
    # where: name = name of the compound being loaded into rdkit and openeye
    #        oe_N = list of stereochemistries for each N in the openeye molecule
    #        rd_N = list of stereochemistries for each N in the rdkit molecule 
    #         ...
    #        rd_etc = a list of tuples that catch any stereochemistry options not specified in the first few columns. 
    # Note: a bit of organization will be needed in excel or an interactive window to make sense of this data
    # Returns: the resulting dframe. Also prints the dframe to a csv in the current working directory 
  
    dframe = organize_errors()

    test_function(test_output, dframe)

    dframe = dframe[dframe['stereo_diff'] == True]
    test_function(test_output, dframe)
    # define the dictinary that will store all stereo info
    temp = {'name': None,
            'oe_N': [],
            'rd_N': [],
            'oe_C': [],
            'rd_C': [],
            'oe_S': [],
            'rd_S': [],
            'oe_etc': [],
            'rd_etc': []}
    stereo_dframe = pd.DataFrame(columns=temp.keys())
    test_function(test_output, stereo_dframe)
    # now check for stereochemistry for each molecule, atom-by-atom. Any stereocenters that are not
    # N, C, or S will fall into the "etc" category and stored just in case
    for index, row in dframe.iterrows():
        d = None
        d = deepcopy(temp)
        d['name'] = row['name']
        oe_mol = row['OFF_from_OE']
        rd_mol = row['OFF_from_RD']
        # iter through openeye molecule
        test = 1
        for atom in oe_mol.atoms:
            stereo = atom.stereochemistry
            if bool(stereo):
                el = atom.element.symbol
                if el == 'N':
                    d['oe_N'].append(stereo)
                elif el == 'C':
                    d['oe_C'].append(stereo)
                elif el == 'S':
                    d['oe_S'].append(stereo)
                else:  # this step catches all other cases 
                    d['oe_etc'].append((el, stereo))
        # iter through rdkit molecule
        for atom in rd_mol.atoms:
            stereo = atom.stereochemistry
            if bool(stereo):
                el = atom.element.symbol
                if el == 'N':
                    d['rd_N'].append(stereo)
                elif el == 'C':
                    d['rd_C'].append(stereo)
                elif el == 'S':
                    d['rd_S'].append(stereo)
                else:  # this step catches all other cases 
                    d['rd_etc'].append((el, stereo))
        # sort everything for readability 
        for key in d.keys():
            if isinstance(d[key], list):
                d[key].sort()

        # dear god please work... 
        stereo_dframe = stereo_dframe.append(d, ignore_index=True)
        print(TestMolecule.are_isomorphic(oe_mol, rd_mol))

    test_function(test_output, stereo_dframe)
    print(stereo_dframe)
    big_df = pd.concat([dframe.reset_index(), stereo_dframe], axis = 1)
    big_df.to_csv("stereo_differences.csv")
    test_function(test_output, big_df)
    return dframe

def smi_diff():
    df1 = run_manual_comparisons()
    df2 = run_comparisons()
    dframe = pd.concat([mols, df1, df2], axis=1)
    smiles_df = pd.DataFrame(columns=['index', 'name', 'molecule_type', 'sm_using_oe', "sm_using_rd"])
    # where name = name of the molecule
    #       molecule_type = what toolkit was used to load in the molecule in the first place
    #       sm_using_oe = the smiles string output when openeye is used
    #       sm_using_rd = the smiles string output when rdkit is used
    for index, row in dframe.iterrows():
        name = row['name']
        oe_mol = row['OFF_from_OE']
        rd_mol = row['OFF_from_RD']
        if (row['smi_oe'] == 0 or row['smi_rd'] == 0) and (row['iso_bond_order'] == True) and (row['comparison'] == True):
            # if the very specific set of test cases are satisfied such that we have a smiles error that is not included in any of the other already discovered errors...
            d_from_oe_mol = {'index': index, 'name': name, 'molecule_type': 'OpenEye',
                    'sm_using_oe': oe_mol.to_smiles(toolkit_registry=open_eye_toolkit_wrapper),
                    'sm_using_rd': oe_mol.to_smiles(toolkit_registry=rdkit_toolkit_wrapper)}
            d_from_rd_mol = {'index': index, 'name': name, 'molecule_type': 'RDKit',
                    'sm_using_oe': rd_mol.to_smiles(toolkit_registry=open_eye_toolkit_wrapper),
                    'sm_using_rd': rd_mol.to_smiles(toolkit_registry=rdkit_toolkit_wrapper)}
            smiles_df = smiles_df.append(d_from_oe_mol, ignore_index=True)
            smiles_df = smiles_df.append(d_from_rd_mol, ignore_index=True)
            # add an extra line just for readability
            smiles_df = smiles_df.append({'index': "", 'name': "", 'molecule_type': "", 'sm_using_oe': "", 'sm_using_rd': ""}, ignore_index = True)
    smiles_df.to_csv("smi_differences.csv")
    return smiles_df
            

if __name__ == "__main__":
    if os.getcwd() != '/home/coda3831/anaconda3/envs/openff-dev':
        raise Exception
    os.chdir('openforcefield/openff/toolkit/Connor')

    drug_bank_file = get_data_file_path('molecules/MiniDrugBank.sdf')
    # drug_bank_file = get_data_file_path('molecules/methyl-fluorene.sdf')
    rdkit_toolkit_wrapper = TestRDKitToolkitWrapper()
    open_eye_toolkit_wrapper = TestOpenEyeToolkitWrapper()

    OFF_using_rd = Molecule.from_file(drug_bank_file, 
                                    toolkit_registry=rdkit_toolkit_wrapper, 
                                    file_format='sdf', 
                                    allow_undefined_stereo=True)
    OFF_using_oe = Molecule.from_file(drug_bank_file, 
                                    toolkit_registry=open_eye_toolkit_wrapper, 
                                    file_format='sdf', 
                                    allow_undefined_stereo=True)
    assert(len(OFF_using_rd) == len(OFF_using_oe))

    mols = pd.DataFrame(data={"OFF_from_OE": OFF_using_oe, "OFF_from_RD": OFF_using_rd})
    # perform a basic check to make sure the dataframe lines up. This is just for my sanity

    all_are_the_same = 1

    [(all_are_the_same * int((str(mol1.name)==str(mol2.name)))) for mol1, mol2 in zip(mols['OFF_from_OE'], mols['OFF_from_RD']) if not isinstance(mol2,float)]
    # if there are any differences (besides RDKit molecules being np.NaN due to incorrect loading that I don't know how to fix right now), 
    # all_are_the_same will return as 0. 
    assert(all_are_the_same == 1)

    # dframe = smi_diff()
    # graph_deviants()
    run_manual_comparisons()
    # print(dframe)



# 1/20 TODO
# 1) part of the openeye toolkit calls "add_H's" and this does not happen in RDKit
#       neither rdkit or openeye should be doing thiss automatically 
#       toolkits.py - line 1298-1299. 
#       check if this is being triggered, what does graph look like before and after? 
# 2) large_kekule diff
#      around the sulfer (S=O) group, one of the toolkits considers it to be a stereocenter, one does not. 
# 3) update the openeye and rdkit columns for thione -- see which one is assigning it stereochem and which is not
#        for the "large_kekule_diff", create a separate table with col1 = rd_thione_stereo, col2 = oe_thione_stereo
# 4) for fail_to_load mols, it seems like one of the toolkits are adding protons to correct nitrogen charges. 
#   also, some of the N atoms should have charges, but don't. This is the case when you have opposite formal charges
#   one bond away from eachother. 
#   Try manually modifying mols from the minidrug bank so that the atoms have physically possible formal charges. 
#   Try to "force" RDKit to load in mols by tweeking sanitization (if possible). 
# 5) Create a database of "edge-cases" that test the boundaries between different conventions and toolkits 
#    For each of these edge cases, include info on: what RDKit does, what OpenEye does, and when should happen (ie. identical
#       smarts matching no matter what toolkit or input form you use for all "common cases" > 99% of drug-like molecules). 
# 6) https://github.com/openforcefield/openff-toolkit/issues/146
#    see above for current table on smarts matching. 
# 7) Open an issue on github for openeye erraneously adding protons 
# 8) should be invited to internal -- keep up to date with "RDKit guy" (Greg Landrum)

# for name, atoms in get_atom_nums(OFF_using_rd):
#     print(f"{name}: {atoms}")

# len1 = sum(1 for name, atoms in get_atom_nums(mols['OFF_from_OE']))

# for name, atoms in BlockToAtomSupplier(drug_bank_file):
#     print(f"{name}: {atoms}")


# summer internship through CalTech
#       Apply for comp chem 
# http://announcements.surf.caltech.edu/index.cfm?event=ShowAOPublicList&inFrame=&type=SURF&formType=AO_CIT


# 1/26/21 Notes
# thion stereochem:
#   do they all dissagree? Load in a small molecule and see which toolkit is wrong. 
# failure to load:
#   count atoms in openeye molecules vs the number of atoms in the original file. 
#   RDKit is choosing to not load in. OpenEye is protonizing to fix percieved error. -> + and - charge 
#   on neighboring atoms? Any nitrogen with four bonds? 
# smiles differences:
#   indicates that each toolkit has different methods of assigning stereochem to N. This is not caught
#   in the are_isomorphic test because the strip_pyrimidal_n parameter is set to true. 
#   ie. this is a known differnece in the toolkits. Is the inversion rate of sulfoxides similar to nitrogens? - no
#   Sulfoxides to no invert. (https://en.wikipedia.org/wiki/Esomeprazole). We do want stereo around Sulfoxide! 
#   Which toolkit is correct? 
# 

 



