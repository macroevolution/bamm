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

The main results from a BAMM run are written to the file specified by ``eventDataOutfile`` in the control file.  If ``outName`` is also specified in the control file it will be prefixed to the name in ``eventDataOutfile``. So if ``outName = BAMM`` and ``eventDataOutfile = event_data.txt`` BAMM results are written to a file named *BAMM_event_data.txt*. 

Before diving too deeply into the results in this file you should first assess whether or not the Markov chain converged. The information needed to do that is written to the file specified by ``mcmcOutfile`` in the control file. There are several programs available for assessing convergence of Markov chains but an easy one to use is the coda package for R. Assuming you are running R you can install coda by simply typing::
	
	> install.packages("coda")

If that fails you can get the package directly from the website http://cran.r-project.org. Once installed you can use the functions in the coda package by loading it with::
	
	> library(coda) 

A quick way to assess convergence is to take a look at the effective sample size of the log likelihood or model parameters. Samples from a Markov chain are autocorrelated and this autocorrelation means samples are not independent. Effective sample size is the sample size after accounting for this non-independence. A low effective sample size means that samples are highly autocorrelated and the chain is not "mixing" well, which is an indication that it is has probably not converged to the posterior distribution you're hoping to sample from. We will calculate the effective sample size of the log likelihood and the number of regime shifts. First read in the MCMC outfile::
	
	> # this is a comment in R
	> # fn is a character string with the path to your mcmc_outfile text file, e.g. "mcmc_out.txt"
	> mcmc <- read.csv(fn)

To take a look at the head of this file do::
	
	> head(mcmc)

So the number of event shifts is the second column and the log likelihood is the fourth column. Before calculating the effective sample size we should discard the beginning of the chain to reduce autocorrelation from the starting configuration. It can be helpful to plot a trace of the log likelihood for deciding on a good cutoff value for this burn-in period::
	
	> plot(mcmc[,1], mcmc[,4])
	> # suppose we think ten percent is satisfactory
	>
	> nsamples <- nrow(mcmc)
	> postburn <- 0.1 * nsamples + 1
	> mcmc <- mcmc[postburn:nsamples, ]

We can now calculate effective sample size using the coda function::
	
	> effectiveSize(mcmc[,2]) 
	> effectiveSize(mcmc[,4])

In general these should be at least 200. Assessing convergence can be complicated and you are encouraged to research other methods.

Once you're satisfied about convergence you are ready to work with the event data file. To work with the data in this file use the utility functions in the BAMMtools for R. This can be downloaded just like in the example above.

Once you've loaded BAMMtools in your R session you can take a look at the main results::

	> # fn is character string specifying the path to your event data file, e.g. "event_data.txt"
	> # mytree is a phylogenetic tree is ape format. see ape documentation for the function read.tree
	>
	> edata <- getEventData(mytree, fn, burnin = 0.1, type = "diversification")
	>
	> # if you are working with BAMM trait data specify type = "trait"

Now that the event data is loaded we can take a look at what it contains. For an explanation of the R object that BAMMtools uses to work with the data simply type::

	> ?getEventData

To quickly summarize your data do::

	> summary(edata)

This will tell you how many posterior samples were analyzed as well as the number of shifts in the maximum shift credibility tree and the tipward node(s) (in ape format) of the branch(es) where those shifts occur. It will also print out the posterior distribution of the number of shifts so you can gauge the relative support for models with different numbers of events. Note that a value of zero means there are no shifts and the single root event describes the entire tree.

To visualize how speciation rates or rates of trait evolution vary through time and among lineages simply type::

	> plot(edata)

You can also plot a polarized version of the tree::

	> plot(edata, method = "polar")

This calculates the mean of the marginal posterior density of rates of speciation or trait evolution for many different points along each branch and maps those rates to colors such that cool colors represent slow rates and warm colors represent fast rates. If you want to take a look at just a single posterior sample rather than averaging over all posterior samples this is possible::

	> mysample <- 1
	> plot(edata, method = "polar", index = mysample)

If this posterior sample happens to contain shifts you can add these to the plotted tree::

	> addBAMMshifts(edata, method = "polar", index = mysample)

Many more types of analysis and visualization are available and you are encouraged to explore the documentation for BAMMtools.
