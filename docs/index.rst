.. DRSIP documentation master file, created by
   sphinx-quickstart on Thu May 31 15:20:54 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DRSIP: Distance Restraints- and Cyclic Symmetry-Imposed Packing
================================================================

DR-SIP contains tools and molecular docking protocols for predicting the quaternary structures of homo-oligomeric transmembrane proteins (HoTPs) by assuming the complex is cyclic (C\ :sub:`n`) symmetric and filtering docking poses with experimentally measured distance restraints between the monomers in the complex.

The DR-SIP package contains 4 Python modules:

1. drsip: Implements the molecular docking protocols. Current implementation accepts ZDOCK output files.
2. zdock-parser: Parses ZDOCK output files, generates and returns the coordinates of poses without first writing to PDB files.
3. docking-eval: Implements the CAPRI criteria to evaluate docking results with respect to a reference structure. 
4. drsip-common: Common functions used by the other modules.

If you're using any of our packages please :ref:`cite us <cite-us>`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   example
   cite
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
