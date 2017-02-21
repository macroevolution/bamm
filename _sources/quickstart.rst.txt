.. _quickstart:

Quick-start guide to BAMM
=========================

This section assumes you have a compiled version of **BAMM** in a directory on your computer (see :ref:`Setting Up BAMM<bammsetup>` for installation). This is a **quick-start** guide: we will provide guidance for some general parameter settings that should get BAMM running on your dataset, but you should explore the `Configuration <configuration.html>`_ page in more detail to optimize performance.

BAMM can be used to model speciation-extinction rates and phenotypic evolutionary rates across phylogenetic trees. To run a speciation-extinciont analysis on your dataset, you need the following (easiest if all in the same directory): 

* A time-calibrated phylogenetic tree
* A *control* file
* The BAMM program

To use BAMM to analyze phenotypic evolution, you also need a file containing your phenotypic (trait) data. 

Your phylogenetic tree must be **ultrametric**, it must be **fully bifurcating** (no polytomies), and all branch lengths must be **greater than 0**. BAMM checks for these things, but it's also good to do a quick check in R using the *ape* package. You can do this as follows::

	library(ape)
	v <- read.tree("myPhylogeny.tre")
	is.ultrametric(v)
	is.binary.tree(v)
	# Now to check min branch length:
	min(v$edge.length)

If all of those checks are good, we can move on.

Control file
------------

To run ``bamm``, you must always specify a *control file*. The control file contains all of the settings necessary to run the program on your dataset, including the name(s) of the input files you seek to analyze. The easiest way to run BAMM is to place the control file, and all files to be analyzed (e.g., the phylogeny) in the same directory as the **BAMM** application. If your control file is named ``myControlFile.txt``, you would run BAMM as follows (on the OSX operating system)::

    ./bamm -c myControlFile.txt

On Windows machines, you would run::

    bamm -c myControlFile.txt

You may download example control files for :download:`diversification<filesForDownload/template_diversification.txt>`  or :download:`trait<filesForDownload/template_trait.txt>` analyses. Use these as a template for setting up your own analyses. The control file is simply a text file with a set of parameter names, followed by the "equals" sign (=), followed by the parameter value. Anything on a line of the control file to the right of a pound sign (#) will be ignored by the program (e.g., it is considered a *comment*). Part of your control file might look like this::

	# This line is a comment. It will not be read by BAMM
	modeltype = speciationextinction        
	treefile = myPhylogeny.tre                                 
	runInfoFilename = run_info.txt
	# This line is also a comment and not read
	runMCMC = 1                           
	useGlobalSamplingProbability = 1        
	globalSamplingFraction = 0.81            
	seed = 12345
	expectedNumberOfShifts = 1.0 # Another comment here

There are many possible settings that can be tweaked in BAMM. The next two sections give you a simple recipe for running the program on your data, depending on whether you are analyzing speciation-extinction rates or phenotypic evolutionary rates. **There is no guarantee that these settings will work for your dataset**.

.. _speciationextinctionquick:

Speciation-extinction analyses
------------------------------

You must have an ultrametric phylogenetic tree. For optimal performance with the *quick-start* settings, the tree should be calibrated in units of millions of years (e.g, a branch length of 1.0 implies 1.0 million years). As a template, use the example file :download:`template_diversification.txt<filesForDownload/template_diversification.txt>`. The default values in this file work reasonably well for most trees calibrated on million-year timescales but may not work for your data. Here's an example :download:`phylogenetic tree of whales<filesForDownload/whaletree.txt>` that is used elsewhere in this documentation.

If you open the template file, you'll see that there are extensive comments. For each parameter in the BAMM control file, we've included a brief description on the line following the parameter. For example: ::

	modeltype = speciationextinction        
	# Specify "speciationextinction" or "trait" analysis
   
This tells us what parameter ``modeltype`` does. It specifies the type of analysis (here, a speciation-extinction analysis). If we wanted to do a phenotypic evolutionary analysis, we would have set ``modeltype = trait``.

There are only a handful of parameters in the template file that you need to set explicitly in order to run BAMM. These are currently defined with the following symbol: ``%%%%``. For example, you see the following immediately after ``modeltype``::

    treefile = %%%%
    # File name of the phylogenetic tree to be analyzed
	
This is where you specify the name of your phylogenetic tree. For example, ``treefile = mytree.tre``. The other parameters we will force you to define explicitly have to do with output. This parameter block looks like this in the template::

	numberOfGenerations = %%%%

	mcmcWriteFreq = %%%%

	eventDataWriteFreq = %%%%
	
	printFreq = %%%%

	acceptanceResetFreq = %%%%

``numberOfGenerations`` is the number of simulation steps you want in your MCMC analysis. ``printFreq`` is the frequency that BAMM will write some simple information to the screen so you can track the progress of the run. ``mcmcWriteFreq`` and ``eventDataWriteFreq`` tells BAMM how often to write the two basic types of output to file.

BAMM generates two types of output. The first is a file containing basic attributes of the MCMC chain as you sample it at particular timepoints. This includes: the log-likelihood of the data under the current parameters, the number of diversification shifts in the current state, the log-prior density of the current parameters, and a moving-window average of the acceptance rate for the MCMC simulation. The second type of output is the *event data*. This is the real stuff of interest. It contains all parameters associated with the macroevolutionary rate regimes and is used for all the subsequent analyses of evolutionary rates. 

You can set these parameters to whatever you want. However, please remember that you will be working with the *event data file* in R, which is a bit limited on memory. As a rough guide, we suggest choosing a value for ``eventDataWriteFreq`` that gives at least 1000 samples from the posterior, but we also don't see much advantage to having more than 5000. 

``acceptanceResetFreq`` specifies the frequency (in generations) in which
to reset the acceptance and rejection counts for reporting.
Its value should be a divisor of ``mcmcWriteFreq``.

For starters, you should try a simple run with settings like this::

    numberOfGenerations = 5000
    mcmcWriteFreq = 1000
    eventDataWriteFreq = 1000
    printFreq = 100
    acceptanceResetFreq = 1000
	
You'll want to increase all of these once you are sure the program is correctly loading your data etc, but it's a good first check. 

One other block of parameters can be critical to BAMM performance: the priors that you place on your evolutionary rate parameters. The prior block in your control file looks similar to this (ignoring most comments in the template file)::

	expectedNumberOfShifts = 1.0
	lambdaInitPrior = 1.0
	lambdaShiftPrior = 0.05
	muInitPrior = 1.0

These priors may work for your dataset. They may also be extremely inadequate. To this end, we have included a function in the BAMMtools package to help you choose appropriate prior values. The function, ``setBAMMpriors``, will automatically generate a prior block as a text file that you can copy and paste over the prior block in the template file. To do this, you need to install BAMMtools (see `here <postprocess.html>`_), and you need your phylogenetic tree. Assuming you have a phylogenetic tree file ``my_tree.tre``, you can generate the prior block with::
	
	> library(BAMMtools) # Assuming you have installed BAMMtools!
	> setBAMMpriors(read.tree("my_tree.tre"))
	
and the relevant output file will be generated in your working directory. See the help file (``?setBAMMpriors``) for more information. To be clear: this does not optimize priors to your dataset. It simply chooses a set of priors that we have found to be reasonable for most datasets and scales the distributions based on the age (root depth) of your tree. A more complete explanation :ref:`can be found here<ratepriors>`.

Incomplete taxon sampling
*************************

For speciation-extinction analyses BAMM can analytically account for incomplete taxon sampling that might otherwise bias results. You can even correct for *non-random* taxon sampling. An explanation of how to account for both random and non-random taxon sampling is found :ref:`here<incompsampling>`.

.. _phenotypicquick:

Phenotypic evolution
--------------------

This section is redundant with the preceding section on **speciation-extinction**, with a few differences.

You must have an ultrametric phylogenetic tree. For optimal performance with the *quick-start* settings, the tree should be calibrated in units of millions of years (e.g, a branch length of 1.0 implies 1.0 million years). As a template, use the example file :download:`template_trait.txt<filesForDownload/template_trait.txt>`. The default values in this file work reasonably well for most trees calibrated on million-year timescales but may not work for your data.

If you open the template file, you'll see that there are extensive comments. For each parameter in the BAMM control file, we've included a brief description on the line following the parameter. For example: ::

	modeltype = trait        
   
This tells us what parameter `modeltype` does. It specifies the type of analysis (here, a phenotypic evolution analysis). If we wanted to do a speciation-extinction analysis, we would have set `modeltype = speciationextinction`.

There are only a handful of parameters in the template file that you need to set explicitly in order to run BAMM. These are currently defined with the following symbol: `%%%%`. For example, you see the following immediately after `modeltype`::

	treefile = %%%%
	
This is where you specify the name of your phylogenetic tree. For example, ``treefile = mytree.tre``. Since we are analyzing phenotypes, we also need to specify the location of the trait data, which we do here::

	traitfile = %%%%

The trait file should consist of a 2 column text file, with species name followed by a tab, followed by the relevant trait value. Here is an :download:`example file<filesForDownload/primates_logmass.txt>` of log-transformed primate body masses, and :download:`here<filesForDownload/primatetree.txt>` is the corresponding Newick format tree. You should be able to plug these into the control file and get BAMM to run.

The other parameters we will force you to define explicitly have to do with output. This parameter block looks like this in the template::

    numberGenerations = %%%%

    mcmcWriteFreq = %%%%

    eventDataWriteFreq = %%%%
	
    printFreq = %%%%

    acceptanceResetFreq = %%%%

``numberOfGenerations`` is the number of simulation steps you want in your MCMC analysis. ``printFreq`` is the frequency that BAMM will write some simple information to the screen so you can track the progress of the run. ``mcmcWriteFreq`` and ``eventDataWriteFreq`` tells BAMM how often to write the two basic types of output to file. BAMM generates two types of output. The first is a file containing basic attributes of the MCMC chain as you sample it at particular timepoints. This includes: the log-likelihood of the data under the current parameters, the number of diversification shifts in the current state, the log-prior density of the current parameters, and a moving-window average of the acceptance rate for the MCMC simulation. The second type of output is the *event data*. This is the real stuff of interest. It all parameters associated with the macroevolutionary rate regimes and is used for all the subsequent analyses of evolutionary rates.
``acceptanceResetFreq`` specifies the frequency (in generations) in which
to reset the acceptance and rejection counts for reporting.
Its value should be a divisor of ``mcmcWriteFreq``.

You can set these parameters to whatever you want. However, please remember that you will be working with the *event data file* in R, which is a bit limited on memory. As a rough guide, we suggest choosing a value for ``eventDataWriteFreq`` that gives at least 1000 samples from the posterior, but we also don't see much advantage to having more than 5000. 

For starters, you should try a simple run with settings like this::

    numberOfGenerations = 5000
    mcmcWriteFreq = 1000
    eventDataWriteFreq = 1000
    printFreq = 100
    acceptanceResetFreq = 1000
	
You'll want to increase most of these once you are sure the program is correctly loading your data etc, but it's a good first check. 

As for the speciation-extinction models, the priors you place on phenotypic evolutionary parameters can have a substantial impact on BAMM performance. The prior block in your (trait) template control file looks similar to this::

	expectedNumberOfShifts = 1
	betaInitPrior = 1
	betaShiftPrior = 0.05
	useObservedMinMaxAsTraitPriors = 1

These priors may work for your dataset, but they may also be very poor choices: it really depends on the scale of your tree (e.g., depth of the tree) and the variances in your trait values. The function ``setBAMMpriors`` (BAMMtools) will automatically generate a prior block as a text file that you can copy and paste over the prior block in the template file. This new set of priors is matched to the "scale" of your data. To do this, you need to install BAMMtools (see `here <postprocess.html>`_), and you need your phylogenetic tree. Assuming you have a phylogenetic tree file ``my_tree.tre`` and a trait dataset ``my_traitfile.txt``, you can generate the prior block with::
	
	> library(BAMMtools) # Assuming you have installed BAMMtools!
	> setBAMMpriors(phy = "my_tree.tre", traits = "my_traitfile.txt")
	
and the relevant output file will be generated in your working directory. See the help file (``?setBAMMpriors``) for more information. To be clear: this does not optimize priors to your dataset. It simply chooses a set of priors that we have found to be reasonable for most datasets and scales the distributions based on the age (root depth) of your tree and the variance of your trait data. A more complete explanation :ref:`can be found here<ratepriors>`.

As of BAMMtools v2.1, you can generate a customized controlfile from within R. Doing the following::

	> library(BAMMtools)
	> generateControlFile(file = 'divcontrol.txt', type = 'diversification')

will create a template controlfile similar to the one that is available for download from this webpage. Additionally, one can specify BAMM parameter values and have them directly supplied to the controlfile template. The parameters are supplied in the form of a list. Following the whales example used throughout the website, and using evolutionary rate parameter priors supplied by ``setBAMMpriors``, one can easily create a customized controlfile with the following::

	> library(BAMMtools)
	> generateControlFile('divcontrol.txt', type = 'diversification', params = list(
		treefile = 'whales.tre',
		globalSamplingFraction = '0.98',
		numberOfGenerations = '1000000',
		overwrite = '1',
		lambdaInitPrior = '1.889,
		lambdaShiftPrior = '0.032,
		muInitPrior = '1.889',
		expectedNumberOfShifts = '1'))


BAMM output: brief
------------------

BAMM generates multiple types of output files. These (usually) include:

* The ``run_info.txt`` file, containing a summary of your parameters/settings
* An ``mcmc_out.txt`` or equivalent file, containing raw MCMC information useful in diagnosing convergence
* An ``event_data.txt`` file or equivalent, containing all of evolutionary rate parameters and their topological mappings
* A ``prior.txt`` file or equivalent, giving the prior expectation on the number of shift events (this is optional and can be turned off).
* A ``chain_swap.txt`` file, containing data about each chain swap proposal
  (when a proposal occured, which chains might be swapped,
  and whether the swap was accepted).

In general, the post-BAMM workflow consists of:

#. Reading your MCMC file into R and testing whether your run appears to have converged. We advocate doing this using the ``coda`` package for R, which enables you to compute the *effective sample size* of your log-likelihoods and numbers of rate shifts sampled during the MCMC simulation.

#. Summarizing your posterior distribution on the number of rate shift events

#. Loading your ``event_data.txt`` file or equivalent into R using the **BAMMtools** package

#. Many potential downstream analyses, including summarizing mean evolutionary rates for clades, analyses of rate shift distributions, plotting model-averaged rate-through-time curves, and so on.

A more detailed description of BAMMtools workflows for postprocessing BAMM output can be found :ref:`here<bammtools>`.
