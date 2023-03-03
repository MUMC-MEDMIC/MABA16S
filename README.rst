=======
MABA16S
=======


.. image:: https://img.shields.io/pypi/v/maba16s.svg
        :target: https://pypi.python.org/pypi/maba16s

.. image:: https://img.shields.io/travis/Casperjamin/maba16s.svg
        :target: https://travis-ci.com/Casperjamin/maba16s

.. image:: https://readthedocs.org/projects/maba16s/badge/?version=latest
        :target: https://maba16s.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/Casperjamin/maba16s/shield.svg
     :target: https://pyup.io/repos/github/Casperjamin/maba16s/
     :alt: Updates



Oxford nanopore 16S for clinical samples


* Free software: MIT license
* Documentation: https://maba16s.readthedocs.io.

Quickstart
----------
As a quickstart to use this pipeline you need Python 3.6 or higher, conda environment manager  and snakemake.

Usage
^^^^^
.. code-block:: bash

    git clone https://github.com/MUMC-MEDMIC/MABA16S 
    cd MABA16S/maba16s
    python cli.py snakemake -i folders_containing_nanopore16s_reads -o my_output_directory --cores 1 

    # input are directories which hold your nanopore reads. Naming of the output will be done based on the names of these directories


How does it work?
-----------------
1. reads are classified on genus level using kraken2 and SILVA database
2. reads for each genus are extracted
3. each genus readset is mapped to the first species in the SILVA database of this genus
4. consensus sequence is extracted and BLASTed to the SILVA database to obtain a species ID
5. results are compiled and written to a spreadsheet

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
