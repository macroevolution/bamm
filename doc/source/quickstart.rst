.. _quickstart:

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

The main results from a BAMM run are written to the file specified by the
``eventDataOutfile`` parameter in the *control* file.  If the ``outName``
parameter is also specified in the *control* it will be prefixed to the name
in ``eventDataOutfile``. So if ``outName = BAMM`` and ``eventDataOutfile = event_data.txt``
BAMM results are written to a file named *BAMM_event_data.txt*. 

Depending on the ``modeltype`` specification this is a comma-delimited text file with
six to eight columns.  Each row in the event data file contains information about an 
event. The first column will always specify the generation number of the MCMC chain at 
the point the data were written to the file. Each unique generation number is a sample
from the posterior distribution of the probability of the model given the likelihood of
the data and the prior probability of the model parameters. If duplicate generation 
numbers occur it means that the tree contains more than one event in that posterior
sample. Columns two and three are names of two taxa. The most recent common ancestor of 
these taxa locates the branch where the event occurred. If column three is NA it means the 
event occurred on the branch leading to the taxon in column two. The fourth column 
indicates the absolute time on the branch where the event occurred. Time zero occurs at 
the root of the tree. The remaining columns specify the parameters of the event. See
:ref:`bammtheory` to understand the meaning of these parameters. The companion R package
BAMMtools has a function that can extract all the relevant information from this file, so
you do not need to work with it directly. For more information about BAMMtools functionality
see :ref:`bammtools`.

Another important file is the ``mcmcOutfile``. This contains details on the state of
the Markov chain and will be important for :ref:`convergence`. For both diversification
and trait analysis the ``mcmcOutfile`` is a six column comma-delimited text file. Each row
represents a sample from the Markov chain. The first column is the generation in the chain at
the point the sample was taken. The second column is the number of shifts on the tree at that
point in the chain. A value of zero means there are no shifts so only a single event, the root
event, is present on the tree. Columns three and four contain the log prior probability of the
current model parameters and the log likelihood of the data given the current model parameters.
Column five contains the current parameter controlling the rate at which shifts occur on
the tree. BAMM models rate shifts using a compound Poisson process, so this parameter is
simply the rate parameter of a poisson distribution (:ref:`bammtheory`). The last column
is the acceptance rate of the chain, or how frequently proposals to change the state of the
chain are accepted. Ideally this should be around 0.234.

BAMM also generates a ``run_info`` file that contains the time, random number seed, and model
settings when the program was executed. 
