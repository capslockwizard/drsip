# DR-SIP: Distance Restraints- and Cyclic Symmetry-Imposed Packing
DR-SIP contains tools and molecular docking protocols for predicting the quaternary structures of homo-oligomeric transmembrane proteins (HoTPs) by assuming the complex is cyclic (C<sub>n</sub>) symmetric and filtering docking poses with experimentally measured distance restraints between the monomers in the complex.

The DR-SIP package contains 4 Python modules:
1. drsip: Implements the molecular docking protocols. Current implementation accepts ZDOCK output files.
2. zdock-parser: Parses ZDOCK output files, generates and returns the coordinates of poses without first writing to PDB files.
3. docking-eval: Implements the CAPRI criteria to evaluate docking results with respect to a reference structure. 
4. drsip-common: Common functions used by the other modules.

If you're using any of our packages please [cite us](#how-to-cite-us)

# How to Install DR-SIP?
There are 2 ways to install the DR-SIP, using conda (recommended) or pip.

## Conda Installation
We recommend that you install the [Anaconda distribution](https://www.anaconda.com/download/) for Python 2.7.

To install DR-SIP:
```
conda config --append channels conda-forge
conda config --append channels drsip
conda install drsip
```

## Pip Installation
For pip:
```
pip install drsip
```

# How to Use?
## Membrane Protein Docking Protocol
```
drsip membrane static-pdb-file mobile-pdb-file trans-helix zdock-output-file -d distance-restraints-file -o DRSIP-results.csv -p top20/
```
The trans-helix argument is a string containing comma separated resids of each transmembrane helix. Example: "17-46, 69-93", where the first transmembrane helix is from resid 17 to 46 and the second is from resid 69 to 93. Transmembrane helix assignments can be obtained from the [Orientations of Proteins in Membranes database](https://opm.phar.umich.edu/).

While the -o argument is for writing the results of top 20 poses to an CSV file and -p is the folder to write out the PDB files of the top 20 complexes generated from the top 20 poses.

The distance restraints file is optional, the filter will not be applied if there are no distance restraints. Each line in the distance restraints file contains a residue pair formatted as:
```
chainID1 resID1 chainID2 resID2 distance
```
NOTE: The columns are separated by tabs (tab-delimited).

Example:
```
drsip membrane 2oar-static_marked.pdb 2oar_mobile_marked.pdb "17-46, 69-93" MscL_54000_ZDOCK.out -d MscL_FRET_Data.txt -o MscL/DRSIP_results.csv -p MscL/
```

MscL_FRET_Data.txt contains:
```
A	40	B	40	5.033146272
B	25	A	25	4.528073406
...
```

For more details run "drsip membrane -h" or see the [documentation](http://drsip.readthedocs.io/).

## Soluble Protein Protocol
```
drsip soluble static-pdb-file mobile-pdb-file zdock-output-file distance-restraints-file -o DRSIP-results.csv -p top20/
```
Similar to running the membrane protein docking protocol except that the distance restraints file is required and there is no trans-helix argument.

Example:
```
drsip soluble 5ccg_SNARE_marked.pdb 2r83_aligned_domains_marked.pdb Syt1-SNARE_ZDOCK_54000.out Syt1-SNARE_FRET_Data.txt -o Syt1_SNARE/DRSIP_soluble_results.csv -p Syt1_SNARE/
```

For more details run "drsip soluble -h" or see the [documentation](http://drsip.readthedocs.io/).

# Documentation
Full documentation available [here](http://drsip.readthedocs.io/)

# How to Cite Us?
If you use any part of the DR-SIP package please cite us by:
```
Chan Justin, Chien Chi-Hong Chang, Zou Jinhao, Pan Rong-Long, Yang Lee-Wei.
(2019) DR-SIP: Protocols for Higher Order Structure Modeling with Distance
Restraints- and Cyclic Symmetry-Imposed Packing. Manuscript in preperation.
```

# References
The DR-SIP package uses the following packages:
1. [MDAnalysis](https://www.mdanalysis.org/pages/citations/)
2. [BioPython](https://biopython.org/wiki/Documentation#papers)
3. [NumPy](https://www.scipy.org/citing.html)
4. [SciPy](https://www.scipy.org/citing.html)
5. [Numba](https://numba.pydata.org/numba-doc/dev/user/faq.html#how-do-i-reference-cite-acknowledge-numba-in-other-work)
6. [Pandas](https://www.scipy.org/citing.html)
7. [u-msgpack-python](https://github.com/vsergeev/u-msgpack-python).