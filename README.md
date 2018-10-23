# DR-SIP: Distance Restraints- and Cyclic Symmetry-Imposed Packing
DR-SIP contains tools and molecular docking protocols for predicting the quaternary structures of homo-oligomeric transmembrane proteins (HoTPs) by assuming the complex is cyclic (C<sub>n</sub>) symmetric and filtering docking poses with experimentally measured distance restraints between the monomers in the complex.

The DR-SIP package contains 4 Python modules:
1. drsip: Implements the molecular docking protocols. Current implementation accepts ZDOCK output files.
2. zdock-parser: Parses ZDOCK output files, generates and returns the coordinates of poses without first writing to PDB files.
3. docking-eval: Implements the CAPRI criteria to evaluate docking results with respect to a reference structure. 
4. drsip-common: Common functions used by the other modules.

# How to Install DR-SIP?
There are 2 ways to install the DR-SIP, using conda (recommended) or pip.

## Conda Installation
We recommend that you install the [Anaconda distribution](https://www.anaconda.com/download/) for Python 2.7.

To install DR-SIP:
```
conda install drsip
```

## Pip Installation
For pip:
```
pip install drsip
```

# How to Use?
```
./drsip_cli.py static_pdb_file mobile_pdb_file zdock_output_file -d distance_restraints_file --trans-helix "13-46, 69-98" -o DRSIP_results.xlsx -p top10/
```
The input for --trans-helix is a comma separated resids of each transmembrane helix. While the -o argument is for writing the results of top 10 poses to an Excel file and -p is the folder to write out the PDB files of the top 10 complexes generated from the top 10 poses.

Example:
```
./drsip_cli.py 2oar_static_marked.pdb 2oar_mobile_marked.pdb MscL_54000_ZDOCK.out -d MscL_FRET_Data.txt --trans-helix "13-46, 69-98" -o MscL/DRSIP_results.xlsx -p MscL/
```

For more details see the [Documentation](http://drsip.readthedocs.io/) or ./drsip_cli.py -h

# Documentation
See here for the [full documentation](http://drsip.readthedocs.io/)
