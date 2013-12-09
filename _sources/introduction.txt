.. _bammfunction: 
Introduction
===============
  
BAMM (Bayesian analysis of macroevolutionary mixtures) is a program for modeling complex dynamics of speciation, extinction, and trait evolution on phylogenetic trees. The program is oriented entirely towards detecting and quantifying heterogeneity in evolutionary rates. BAMM uses reversible jump Markov chain Monte Carlo to automatically explore a vast universe of candidate models of lineage diversification and trait evolution. BAMM has been described and extended in several publications (*citations*). BAMM is a command line program written in C++. Post-run analysis and visualization is performed using the R package **BAMMtools**. In the `Brief Overview of BAMM Functionality`_, you will find a brief summary of key features of BAMM.


.. _whalefig: 
.. figure:: figs/xIntroFig_whalerates.png
   :width: 640
   :align: center

   Figure 1: Speciation rates in whales, reconstructed using BAMM. BAMM can estimate the posterior density of evolutionary rates (speciation, extinction, trait evolution) at any point in time along any branch in a phylogenetic tree. Major burst in rates at top of tree corresponds to the origin of the dolphin clade. See the :ref:`BAMM Graph Gallery<bammgraphs>` for more examples involving BAMM. Inset histogram shows posterior density of speciation rates and allows color-based interpretation of rates on the whale phylogeny. Figure was generated using the function ### from the *BAMMtools* package in R.   Dataset in this example is taken from M. Steeman *et al.*, *Syst. Biol.* xxx:xxxx, 2009). 

 
Brief Overview of BAMM Functionality
------------------

Allows modeling of rates *through time* and *among clades*
....................................

BAMM abandons the simple notion of *rate shifts* involving constant-rate diversification processes. A rate shift in the BAMM environment is a shift in the evolutionary dynamics of a system and can involve shifts to new time-varying or diversity-dependent evolutionary rate regimes.


Detect evolutionary rate heterogeneity
..............................

BAMM allows explicit tests of the hypothesis that a single process of diversification can account for a given phylogenetic pattern. The program automatically explores many candidate models of diversification and provides a variety of tools for interpretation. BAMM enables estimation of Bayesian credible intervals on the number of distinct macroevolutionary regimes, calculation of Bayes factors for pairwise evaluation of model fit, and calculation of model posterior probabilities. The user can easily determine the evidence favoring rate variation within phylogenies.



Detect key innovations and diversity-dependence
............................

BAMM is oriented towards the automatic detection of key innovations, rate shifts, and diversity-dependence on phylogenetic trees. BAMM finds locations for shifts in evolutionary dynamics that are maximally supported by the data, with no *a priori* specification as to where these shifts in dynamics might have occurred.


Analysis and visualization with BAMMtools
.............................
BAMMtools is a comprehensive R library for the analysis and visualization of macroevolutionary dynamics. The package includes a variety of functions for analyzing and visualizing BAMM output. Examples of BAMMtools functionality includes: visualization of distinct evolutionary regimes on phylogenetic trees, calculation of Bayes factors, plotting evolutionary rates through time, estimating clade-specific average rates, and visualizing rate variation along the branches of individual phylogenetic trees (as shown :ref:`here<whalefig>`). 


Account for *non-random* taxon sampling
............................
Speciation-extinction calculations in BAMM account for incomplete taxon sampling analytically. The program is designed to work with datasets that contain large numbers of missing species, and taxon sampling can be phylogenetically non-random. The only requirement is that the user is able to specify how taxon sampling varies across the tree. For example, you can allow individual clades (such as genera or families) in a large phylogenetic tree to have different sampling probabilities. See XXX crossreference XXX for more information.

Better than stepwise AIC
......................
Many methods for modeling evolutionary dynamics use stepwise AIC-based approaches that are limited to identifying a single best model. These approaches are inherently limited, because many distinct combinations of evolutionary shift regimes might have roughly equal probabilities. Rather than identifying a single best configuration of rate shifts, BAMM samples rate shift configurations in proportion to their posterior probability.



Fast C++ implementation
......................
BAMM's underlying C++ core allows rapid modeling of evolutionary dynamics in phylogenetic datasets that would be too large to handle with R alone.
 