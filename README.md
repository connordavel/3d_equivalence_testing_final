# 3d_equivalence_testing_final
Cheminformatics toolkits tend to give different results/chemical representations for equivalent inputs. This first project in the OpenForceField Initiative
hopes to find these differences, quantify approximately how often they occur, and offer possible solutions. 

## Testing Steps

The following tests were run to compare two parallel lists of `openff.toolkit.molecule.Molecule` objects.

1. comparing the output of `Molecule.to_smiles` 
    1. return 1 if smiles are identical, 0 if smiles are similar, -1 else
2. comparing the output of `Molecule.to_inchi`
    1. return 3 if inchi are identical, 2, 1, 0, or -1 else
3. comparing the results of `Molecule.are_isomorphic` with 3 different bond matching functions:
    1. match by aromaticity only `Molecule.are_isomorphic(mol1, mol2, bond_order_matching=False)`
    2. match by bond order only `Molecule.are_isomorphic(mol1, mol2, aromatic_matching=False)`
    3. match by both `Molecule.are_isomorphic(mol1, mol2)`
4. comparing the explicit hydrogen and heavy atoms as a double check on the previous tests. 

<img src="3d_equivalence_testing_final/3d_compare_flow_chart_new.png" alt="Drawing" style="width: 900px;"/>

## (Brief) Results

A full description of the results can be found in the 3d_equivalence_testing_final.ipynb notebook, but 
a short summary is given here. Errors can be placed into the following categories:

1. Different kekule structures around aromatic rings (Low Priority)
2. OpenEye including one or more extra protons on ambiguous molecules (High Priority)
3. Different stereochemistry in the sulfur of solfoxides (High Priority)
4. Difference in stereochemistry around pyramidal nitrogens (Low Priority - already accounted for) 
5. RDKit refusing to load in nitrogens with 4 bonds and ambiguous formal charge (Medium Priority)

## 1. Arom_kekule_diff
These errors are not very important and are handled by the internal logic of `Molecule.are_isomorphic()`. I have included some examples of small molecules below (larger molecules do are not scaled quite so well when graphing in networkX). On the left is the networkX graphs of each molecule. Red edges indicates double bonds. On the right is the direct pymol representation of the molecule.

<h2><center>Drug_Bank_3674</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/arom_kekule_diff/DrugBank_3674_DrugBank_3674_(117).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/arom_kekule_diff/Drug_Bank_3674.PNG" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>

<h2><center>Drug_Bank_4047</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/arom_kekule_diff/DrugBank_4047_DrugBank_4047_(162).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/arom_kekule_diff/Drug_Bank_4047.PNG" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>

<h2><center>Drug_Bank_4247</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/arom_kekule_diff/DrugBank_4247_DrugBank_4247_(197).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/arom_kekule_diff/Drug_Bank_4247.PNG" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>

<h2><center>Drug_Bank_6145</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/arom_kekule_diff/DrugBank_6145_DrugBank_6145_(122).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/arom_kekule_diff/Drug_Bank_6145.PNG" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>

In all of these cases, RDKit appears to adjust the kekule structures, while OpenEye says the most "faithful" to the original structure.
However, like I said before, this is a pretty minor issue and not worth spending any more time on. 

## proton_diff
These errors are mostly due to ambiguous sdf info blocks, specifically incomplete or nonsensical 
information about formal charge. RDKit does nothing to these molecules, while openeye adds a proton 
during the Molecule creation step. Here are a few indicative test cases. 
<h2><center>DrugBank_1564</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/proton_diff/DrugBank_1564_DrugBank_1564_(184).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/proton_diff/DrugBank_1564.png" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>
<h2><center>DrugBank_1700</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/proton_diff/DrugBank_1700_DrugBank_1700_(211).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/proton_diff/DrugBank_1700.png" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>
<h2><center>DrugBank_2148</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/proton_diff/DrugBank_2148_DrugBank_2148_(297).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/proton_diff/DrugBank_2148.png" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>
<h2><center>DrugBank_2210</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/proton_diff/DrugBank_2210_DrugBank_2210_(305).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/proton_diff/DrugBank_2210.png" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>
<h2><center>DrugBank_2642</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/proton_diff/DrugBank_2642_DrugBank_2642_(356).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/proton_diff/DrugBank_2642.png" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>
<h2><center>DrugBank_5418</center></h2>
<table><tr>
<td> <img src="3d_equivalence_testing_final/netx_figs/proton_diff/DrugBank_5418_DrugBank_5418_(11).png" alt="Drawing" style="width: 550px;"/> </td>
<td> <img src="3d_equivalence_testing_final/pymol_figs/proton_diff/DrugBank_5418.png" alt="Drawing" style="width: 550px;"/> </td>
</tr></table>

## Stereo Differences 
As seen above, printing the stereochemistry around sulfurs gives opposite results. As far as I know,
this cannot be blamed on ambiguity in the SDfile, since (I think) stereochemistry is interpretted by the 
toolkit upon gathering info form the file. To determine which toolkit appear to be "right" be the human eye, 
I included the following 4 cases with the lone pair oriented behind the 3 pyramidal bonds. 

<h2>Drug_Bank_1971 (oe->R, rd->S, human->R)</h2> 
<img src="3d_equivalence_testing_final/pymol_figs/stereo_diff/DrugBank_1971.png" alt="Drawing" style="width: 700px;"/>
<h2>Drug_Bank_2140 (oe->R, rd->S, human->R)</h2> 
<img src="3d_equivalence_testing_final/pymol_figs/stereo_diff/DrugBank_2140.png" alt="Drawing" style="width: 700px;"/>
<h2>Drug_Bank_2563 (oe->R, rd->S, human->R)</h2> 
<img src="3d_equivalence_testing_final/pymol_figs/stereo_diff/DrugBank_2563.png" alt="Drawing" style="width: 700px;"/>
<h2>Drug_Bank_4032 (oe->S, rd->R, human->S)</h2> 
<img src="3d_equivalence_testing_final/pymol_figs/stereo_diff/DrugBank_4032.png" alt="Drawing" style="width: 700px;"/>


In conclusion, openeye agrees with what a human would assign to the stereochemistry, assuming the lone pair is the least significant group.

# Brief Conclusion

Given the ambiguities in the SDFile format surrounding what constitutes an implicit hydrogen it is hard to prefer one toolkit
over the other without more information. So far, OpenEye appears to handle these issues in the most intuitive way given the input 
SDF block, and this applies to sulfoxide stereochemistry as well. Fixes for the future will make an attempt at remedying the differences 
between OpenEye with minimal deviation from the original behavior and providing low-cost error checking for the user. This could be implimented by 
streamlining some of the methods in this notebook. 
