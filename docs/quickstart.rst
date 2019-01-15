Install and Quick Start
================================================================

How to Install DR-SIP?
------------------------
There are 2 ways to install the DR-SIP, using conda (recommended) or pip.

Conda Installation
^^^^^^^^^^^^^^^^^^^^
We recommend that you install the `Anaconda distribution <https://www.anaconda.com/download/>`_ for Python 2.7.

To install DR-SIP::

    conda config --append channels conda-forge
    conda config --append channels drsip
    conda install drsip

Pip Installation
^^^^^^^^^^^^^^^^^^^^
For pip::

    pip install drsip

How to Use?
------------------------

Membrane Protein Docking Protocol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    drsip membrane static-pdb-file mobile-pdb-file trans-helix zdock-output-file -d distance-restraints-file -o DRSIP-results.csv -p top20/

The trans-helix argument is a string containing comma separated resids of each transmembrane helix. Example: "17-46, 69-93", where the first transmembrane helix is from resid 17 to 46 and the second is from resid 69 to 93. Transmembrane helix assignments can be obtained from the `Orientations of Proteins in Membranes database <https://opm.phar.umich.edu/>`_.
While the -o argument is for writing the results of top 20 poses to an CSV file and -p is the folder to write out the PDB files of the top 20 complexes generated from the top 20 poses.

The distance restraints file is optional, the filter will not be applied if there are no distance restraints. Each line in the distance restraints file contains a residue pair formatted as::

    chainID1 resID1 chainID2 resID2 distance

NOTE: The columns are separated by tabs (tab-delimited).

See an example :ref:`here <example-HoTP>`.

Soluble Protein Protocol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    drsip soluble static-pdb-file mobile-pdb-file zdock-output-file distance-restraints-file -o DRSIP-results.csv -p top20/

Similar to running the membrane protein docking protocol except that the distance restraints file is required and there is no trans-helix argument.

See an example :ref:`here <example-soluble>`.