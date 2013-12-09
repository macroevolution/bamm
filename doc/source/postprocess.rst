Analyzing BAMM output with BAMMtools
===============================================

Brief
................
The R package **BAMMtools** contains almost everything you need to analyze and visualize evolutionary rate dynamics from BAMM output. You will need to install this package. To install BAMMtools, type the following in your R console::

	> install.packages("BAMMtools")
	
A comprehensive introduction to BAMMtools can be found **here** (**linked ref to BAMMtools latex**). The information below is a brief overview to get you started, but you should refer to the linked pdf documentation for details. 

BAMM output files
............................

BAMM generates three primary output files. The first is the *mcmc data file*, which contains several pieces of information about the MCMC simulation that may be useful in diagnosing convergence. The most important pieces of information from this file are the number of shift events, the log-likelihood of the data, and the log-prior probability of the data, for each sample from the posterior. 

The second is the *event data file*, which contains all of the actual model parameters. Each sample from the posterior can be described with complete knowledge of all of the shift events on the tree, including their location and evolutionary rate parameters. The *event data file* is just a long list of all the shift events and associated parameters that were sampled, as well as the MCMC generation in which they were sampled. **You will not directly do anything with this file**. BAMMtools has a function for extracting all the relevant information from the file and for mapping the rate shift configurations to phylogenetic trees. Finally, BAMM will (optionally) generate a second MCMC output file that contains the results of a prior-only simulation. This file can be used to reconstruct the prior distribution of the number of shift events and is important for the estimation of Bayes factors.


Diagnosing convergence
......................
The first question after running any MCMC simuluation should always be: *did my run converge?* While it may be difficult to prove convergence in an absolute sense, there are a few simple checks you can do. First, you can plot the log-likelihood trace of your MCMC output file::

	> mcmcout <- read.csv("mcmc_out.txt", header=T)
	> plot(mcmcout$logLik ~ mcmcout$generation)
	> 
	> # Also, remember that comments in R are lines that start with 
	> #     a pound sign. We'll occasionally use them here.
	> #     R will not execute these lines - they are only for reference!
	
This can give you a ballpark idea of whether your run has converged. The next step is to discard some as burnin. Here we'll discard the first 10% of samples as burnin::

	> burnstart <- floor(0.1 * nrow(mcmcout))
	> postburn <- mcmcout[burnstart:nrow(mcmcout), ]

And using these samples, it is good to check the *effective sample sizes* of the log-likelihood and the number of shift events present in each sample. We'll do this using the *coda* library for R::

	> library(coda)
	> effectiveSize(postburn$N_shifts)
	> effectiveSize(postburn$logLik)

In general, we want these to be at least 200 (and 200 is on the low side, but might be reasonable for very large datasets).
 
If you are having trouble with convergence, please see the section on troubleshooting convergence issues (**LINK** to FAQ). 

General BAMMtools workflow
..........................................

The general analyses described in this section apply both to **speciation-extinction** and **phenotypic evolution** studies. As such, they are not treated separately. The primary difference is the name of the parameters (:math:`\lambda` and :math:`\mu` for speciation and extinction, and :math:`\beta` for trait evolution).

The first step in the analysis of BAMM output is to ask some basic questions about the number of macroevolutionary rate regimes on our phylogenetic tree. We can do this directly using the post-burn *MCMC output file*. Here, we'll compute the posterior probabilities of models sampled using BAMM::

	> post_probs <- table(postburn$N_shifts) / nrow(postburn)

*post_probs* is now a vector of model posterior probabilities. We can look at these, we can plot them, we can compute posterior odds ratios, and so on. To see which models are part of the set that were sampled, you can just look at the names of this vector:: 

	> names(post_probs)
	
And to compute the posterior odds ratio for (say) two models 'X' and 'Y' (X and Y must be integers), you would do::

	> post_probs['X'] / post_probs['Y'] 

In general, any model that is not included in *names(post_probs)* was so lousy that it was never even sampled. Thus, if you fail to observe a model '0' in this set, this means that you have such overwhelming evidence for diversification rate heterogeneity in your data that this model probability is effectively 0 (bear in mind that a model with name '0' is model :math:`M_0`, or a model with no rate shifts). The probability of model '0' is the posterior probability of a model with just a single evolutionary rate dynamic (no rate shifts). We'll discuss the use of Bayes factors in gauging model support a little further down in this document (**INTERNAL LINK**), but posterior model probabilities are a valuable way of identifying the best model (or set of models). The BAMMtools vignette **#LINK#** shows other tricks, like identifying the 95% credibility set of models.

The next step in the analysis of BAMM output is to get the *event data file* into a more workable format. Using BAMMtools, we'll process this output into an R data object that will be much easier to work with. We also need our time-calibrated phylogenetic tree that we analyzed::
	
	> library(ape)  # Need the ape package!
	> mytree <- read.tree("my_example_tree.tre")
	> edata <- getEventData(mytree, eventfilename = "event_data.txt", burnin=0.1, verbose=T)

*edata* is now a "BAMM-data" object, which has all the attributes of a class "phylo" object, plus a few more. Please be patient with *getEventData* - this function can take some time to run for large datasets. 

Now we will focus on a few simple analyses that you can do with the 'BAMM-data' object. We won't go into details here about making particularly pretty plots (see the **LINK BAMM GRAPH GALLERY** for some ideas, as well as the **SWEAVE DOCUMENTATION**) but will use tools available in *BAMMtools* and the *ape* package.


Analyzing locations of rate shifts
----------------------------------
Once you have established there there is at least some evidence for heterogeneous evolutionary dynamics in your dataset, the obvious question is: where are these rate shifts? In the BAMM framework, this is a deceptively simple question, because BAMM does not generate a single *best* rate shift configuration. In the BAMM framework, many different shift configurations may be (more-or-less) equally plausible. BAMM samples shift configurations in proportion to their posterior probability. In principle, this means that each sample from your posterior contains a potentially unique configuration of regime shift events. 

A conceptual discussion of the meaning of rate shifts is included in this documentation and it is **strongly recommended** that you :ref:`read this section before continuing<bammtheory>`. Approaches that identify a single best shift configuration (e.g., stepwise AIC, or other approaches that simply maximize the likelihood) are inherently limited by their assumption that the model with the best information theoretic score (AIC etc) is *the* model, given the candidate set of models. However, for most real datasets, the best rate shift configuration is merely one of a large number of possible rate shift configurations that have similar probabilities. The BAMM philosophy is largely oriented around addressing this. 

The following instructions for visualizing rate shifts assume that you have read the relevant :ref:`documentation<bammtheory>` on the topic and that you understand the difference between the *marginal shift probabilities*, the *cumulative shift probability tree*, and the *maximum shift credibility configuration*. 
 
Starting at the beginning, let's re-load our event data file::

	> mytree <- read.tree("my_example_tree.tre")
	> edata <- getEventData(mytree, eventfilename = "event_data.txt", burnin=0.1)
	> # How many samples from the posterior are in edata?
	> # Here is a quick check:
	>
	> length(edata$eventData)
	
This let's us see how many post-burnin samples from the posterior are included in the *edata* object. Each of these samples is associated with a potentially unique shift configuration. Let's visualize just a single rate shift configuration from our *bamm-data* object::

	> mysample <- 1  # this is the sample we'll plot
	> nrow(edata$eventData[[ mysample ]]) 
	
Will give us the total number of rate rate regimes on our tree for the i'th sample from the posterior. If there is only 1, then you have no rate shifts: the single rate regime starts at the root and describes the entire tree. Assuming you have more than 1, we can get the node numbers (in *ape* format), as follows::

	> shiftnodes <- getShiftNodesFromIndex(edata, index = mysample)	
 
And we can plot these nodes on the tree like this::

	> plot.phylo(mytree)
	> nodelabels(node = shiftnodes, pch=21, col="red", cex=1.5)
	
This highlights the *downstream node* (e.g., "tipwards", as opposed to "rootwards") at the end of each branch on which a shift occurs in the specified sample. You should be able to repeat this exercise again with a different value for *mysample*, and sooner or later, you should be able to see that different shift configurations will "light up" on your tree. Note that if there are no shifts in a given sample, there are no nodes to plot, which will lead to an error message. Let's actually view the actual configuration of evolutionary rates across our tree, as well as the location where the shifts are inferred to occur. To do this, we will use BAMM's function *plot.dtrates*, which will plot rates at multiple points in time along every branch of the phylogeny, for the specified sample::

	> ##### THIS IS NOT DONE
	> mysample <- 1
	> subdata <- subsetEventData(edata, index = mysample)
	> plot.dtrates(mytree, subdata)

Which should generate a nice plot showing rate dynamics. And we can even visualize the actual location of the regime shift events themselves, like this::

	> ##### THIS IS NOT DONE
	> shiftnodes <- getShiftNodesFromIndex(subdata, index=1)
	> plotShiftsOnTree(node = shiftnodes, pch=21, col="white", cex=2.5)
 
We'll come back to the function *plot.dtrates* when we discuss :ref:`branch-specific rates<bammtoolsRTT>`.

One of the first things to look at is the marginal shift probabilities on individual branches. This is nothing more than the marginal probability that each branch contains a shift event (see :ref:`here<bammtheory>` for why these can be difficult to interpret). The next few lines of code will compute the marginal shift probabilities for each branch, then plot a new phylogenetic tree where the branch lengths are scaled by the probability that they contain a shift event::

	> marg_probs <- marginalShiftProbsTree(edata)
	> plot.phylo(marg_probs)
	
The variable *marg_probs* becomes a copy of our phylogenetic tree, but where each branch length has been transformed into the corresponding marginal shift probability. The marginal shift probabilities can be a little misleading, because we might have relatively low confidence in precisely which branch a shift occurred on, but nonetheless have extremely high confidence that a shift occurred *somewhere* in the vicinity. The *cumulative shift probability tree* shows the cumulative probability, on each branch, that a shift occurred somewhere between the focal branch and the root of the tree. The occurrence of such a shift implies that evolutionary dynamics on the focal branch are decoupled from the "background" diversification or trait evolutionary process at the root of the tree. We can compute and plot the cumulative shift probability tree as follows::

	> cst <- cumulativeShiftProbsTree(edata)
	> plot.phylo(cst)

Or, showing shift probs in color::

	> cst <- cumulativeShiftProbsTree(edata)
	> edgecols <- rep('black', length(mytree$edge.length))
	> is_highprobshift <- cst$edge.length >= 0.95
	> edgecols[ is_highprobshift ] <- "red"
	> plot.phylo(mytree, edge.color = edgecols)
	
And this should plot your tree (*mytree*) such that all branches with cumulative shift probabilities of 0.95 or higher are identified in red. See also the example in the :ref:`BAMM graph gallery<cst>`.  	

Another tool for visualizing shift configurations is to plot the *maximum credibility shift configuration*. This is the sample from the posterior that maximizes the marginal probability of the shift events (more info :ref:`here<bammtheory>`). We can find the index of a sample from the posterior with the highest marginal probability as follows::

	> msc <- maximumShiftCredibilityTree(edata)

The variable *msc* contains a bit of information, including the marginal shift probabilities of each sample in the posterior. There will often be ties. For example, if the best shift configuration is sampled multiple times, every sample with this configuration has exactly the same combined marginal probability. We thus select a single representative and plot the shift nodes::

	> samp <- msc$sampleindex
	> shiftnodes <- getShiftNodesFromIndex(edata, samp)
	> plot(mytree)
	> nodelabels(mytree, node = shiftnodes, cex=2, bg="red", pch=21)

And we can also view the rate-through-time dynamics implied by this sample::

	> subdata <- subsetEventData(edata, index = samp)
	> plot.dtrates(subdata.....)

  

Estimating clade-specific rates
-------------------------------

**Under construction**


Branch-specific evolutionary rates
----------------------------------

**Under construction**

.. _bammtoolsRTT:

Plotting rate-through-time curves
---------------------------------

Plotting rate-through-time curve (example **here (link)**) is trivial. BAMM's built-in function *plotRateThroughTime* makes it easy to generate plots of rates through time::

	> plotRateThroughTime(edata, ratetype="speciation")
	
should produce a plot with density shading on confidence regions for your speciation-through-time curve. See help on this function for more details about tweaking the plot. This function can take awhile to run, because it generates a rate-through-time matrix that includes all samples in the posterior distribution. 

You can also use *plotRateThroughTime* to plot speciation through time curves for just portion of your phylogeny. We can do this by feeding a node number in to *plotRateThroughTime*, and the function will just compute and plot the rates for this subtree. To find a particular node number for your tree, you can plot the tree (using ape), and then plot your node numbers directly on the tree, like this::

	> mytree <- read.tree("example_tree.tre")
	> plot.phylo(mytree)
	> nodelabels(mytree)
	
Another way of doing this is to extract the most recent common ancestor (MRCA) node for your clade, by specifying the names of 2 descendant species from the clade that span the focal clade::

	> species1 <- "Homo_sapiens"
	> species2 <- "Ctenotus_pantherinus"
	
Now to get the *tip node numbers* in ape format::	
	
	> tipnode1 <- which(mytree$tip.label == species1)
	> tipnode2 <- which(mytree$tip.label == species2)
	
And now the MRCA node::

	> mrca <- getMRCA(mytree, tip = c(tipnode1, tipnode2))
	
Now we feed this in to *plotRateThroughTime*::

	> plotRateThroughTime(edata, node = mrca, nodetype="include")
	
And we can also plot the entire rate through time curve after we **exclude** this clade (as in: just plot the background rates, without the focal clade)::

	> plotRateThroughTime(edata, node = mrca, nodetype = "exclude")
	
There are many other options available through *plotRateThroughTime*, so please see the R help on this function::

	> ?plotRateThroughTime
	
That's the quick and dirty way of plotting rates through time. Often, you will want more control over the plotting process. The core BAMM operation for plotting a rate-through-time curve involves the generation of a rate-through-time matrix, like this::

	> rtt <- getRateThroughTimeMatrix(edata)

which returns a list of rate-through-time matrices plus a vector of the time points at which the rates were computed. If your rate matrix was for trait evolution, you will have a component *rtt$beta* in your rtt object (components *rtt$lambda* and *rtt$mu* if you are modeling speciation-extinction). To get the mean rates at any point in time::

	> meanTraitRate <- colMeans(rtt$beta)
	
and to do a simple no-frills plot::

	> plot(meanTraitRate ~ rtt$times)
	
You can also include- and exclude nodes from the calculation of the rate-through-time matrix (assuming you know the node to exclude or include)::

	> rtt_subtree <- getRateThroughTimeMatrix(edata, node = mynode)
	
Please see code underlying some BAMM graph gallery plots for more on working with these objects. For example, the code linked **here (#LINK#)** demonstrates how you can directly work with the rate matrices for extremely flexible plotting options.


Computing Bayes factors
----------------------------------
BAMMtools makes it easy to compute Bayes factor evidence in favor of one model relative to another. The disadvantage of Bayes factors is that they provide a measure of pairwise model support and don't necessarily identify a single best model (this isn't necessarily bad: *is* there a single best model?). An advantage of Bayes factors as that they allow model comparisons to be made *independent of the prior on the model*. In BAMM, you specified a hyperprior distribution on the number of shift regimes, and this will have some effect on your posterior model probabilities, so it can be useful to look at the Bayes factor matrix for model comparisons.

This analysis assumes that you have generated an *MCMC output file* involving simulation from the **prior only**. BAMMtools will need to perform explicit comparisons of the prior and posterior model probabilities. Assuming you have files *prior_mcmc_out.txt* and *post_mcmc_out.txt* for your analysis, you can compute a pairwise Bayes factor matrix as::

	> postfile <- "post_mcmc_out.txt"
	> priorfile <- "prior_mcmc_out.txt"
	> computeBayesFactors(postfile, priorfile, burnin=0.1)
	
and this will return a pairwise matrix of Bayes factors. It is very important to recognize that model probabilities for rarely sampled models are likely to be inaccurate. Hence, BAMM will return a matrix with missing values (NA) if a given model was insufficiently sampled to estimate posterior or prior odds (see the *threshold* argument in ?computeBayesFactors). Also keep in mind that any model sampled too infrequently to estimate model odds is also a model that is highly improbable given the data, so the missing Bayes factors aren't really something to worry about. Please see the analysis detailed here (**#LINK# to bird analysis file in graph gallery**) for analysis and visualization of pairwise Bayes factors for a large set of candidate models.


 
BAMMtools workflow: speciation-extinction
..........................................


BAMMtools workflow: phenotypic evolution
..........................................
