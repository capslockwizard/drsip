Examples
========

Download the examples `here <https://github.com/capslockwizard/drsip/raw/master/examples.zip>`_ and extract the zip file.

First, change the current working directory to the examples folder::

    cd examples

.. _example-HoTP:

Then for the HoTP system, MscL::

    drsip membrane 2oar_static_marked.pdb 2oar_mobile_marked.pdb "17-46, 69-93" MscL_54000_ZDOCK.out -d MscL_FRET_Data.txt -o MscL/DRSIP_results.csv -p MscL/

where "17-46, 69-93" are the 2 transmembrane helices assigned by `OPM <https://opm.phar.umich.edu/proteins/35>`_. The results table and top 20 complexes will be written into the MscL folder.

While, MscL_FRET_Data.txt contains::

    A   40  B   40  5.033146272
    B   25  A   25  4.528073406
    ...

.. _example-soluble:

As for the soluble system, Syt1-SNARE::

    drsip soluble 5ccg_SNARE_marked.pdb 2r83_aligned_domains_marked.pdb Syt1-SNARE_ZDOCK_54000.out Syt1-SNARE_FRET_Data.txt -o Syt1_SNARE/DRSIP_results.csv -p Syt1_SNARE/

The results table and the top 20 poses will be written into the Syt1_SNARE folder.

See the MscL_ref and Syt1_SNARE_ref folders for the pre-computed results.