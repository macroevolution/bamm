Quick-start guide to BAMM
=========================

To run ``bamm``, you must always specify a *contol* file. For example,
if your control file is named ``divcontrol.txt``, run the following::

    bamm -c divcontrol.txt

Data files and input
--------------------

BAMM comes with example control files (located in the directory examples).
Use these as templates to quickly get your data analyzed by BAMM.

Currently, BAMM can process two types of models:
speciation-extinction analysis or phenotypic evolution analysis.
The model type is specified by the parameter ``modeltype`` in the control file.

Specify the file path of the Newick-formatted tree with ``treefile``.
All file paths are relative to the directory in which bamm was executed.

Speciation-extinction analyses
******************************

(To be done.)

Incomplete taxon sampling
*************************

(To be done.)

Phenotypic evolution analyses
*****************************

Specify the file path of the trait data in ``traitfile``.

Setting up the control file
---------------------------

(To be done.)

Running the program
-------------------

If your control file is named ``divcontrol.txt``, run BAMM as follows::

    bamm -c divcontrol.txt

Understanding BAMM output
-------------------------

(To be done.)
