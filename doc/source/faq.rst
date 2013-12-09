Frequently Asked Questions
==========================
 
General
.............

Citing BAMM
-----------

Currently, a note describing the BAMM software and BAMMtools R library is in review. For now, you can cite one of several papers that initially described the model and its implementation::

	Rabosky et al 
	
	Rabosky and Matute



BAMM extensions
-----------------------

We are currently testing extensions of BAMM that allow modeling evolutionary dynamics under a much greater range of phenotypic evolutionary scenarios, as well as the incorporation of paleontological data. Bleeding edge releases of BAMM can be obtained from here (**link**).


Problems running BAMM
......................


I can't get BAMM to run on my system
-----------------------------------------

**Under construction**

Brief list of possible things to try.


I get an error message when I try to run BAMM
----------------------------------------------

**Under construction**

This section needs to explain common error messages obtained with BAMM. For example,::
	
	ERROR: No control file specified

implies that BAMM is having a hard time reading your control file. Are you spelling the control file name correctly (including the extension), and/or are you sure it is in the correct directory?

Stepwise troubleshooting
-----------------------------------------

The simplest way to localize potential problems is to first try to load your data without doing anything else. You can do this by setting the following parameter in your controlfile::

	initializeModel = 0 

This tells BAMM to simply try reading the tree (and associated trait or taxon sampling datafiles). If this works fine, your problem is probably elsewhere. The next step is to try initializing your MCMC simulation, without actually running it. This step involves actually computing likelihoods and priors, so may reveal problems with those steps::

	initializeModel = 1 
	runMCMC = 0

If these issues fail to resolve the problem, and if you cannot 

 
I get a system error (segfault or equivalent) when running BAMM
---------------------------------------------------------------------

Other possible errors may arise, not discussed above. This may include errors that mention "Segmentation fault", or "Floating point error". Please contact Carlos Anderson (**link**) or Dan Rabosky (**link**) with information about system errors. Please send us as much information as possible, so that we can try to replicate the problem. You should include (if possible) the data files that led to the problem, your control file, and as much information as possible about your operating system and computer architecture. On unix/linux/OSX, you should be able to obtain most of this information using::
	
	uname -a

If you compiled BAMM on your system, please send details about the compiler used. You can determine the compiler you've used to build BAMM with:: 

	How to do this?

If BAMM actually begins to perform an analysis, then it will also generated a detailed run_info file that contains information about your analysis. Please include this file in your bug description. 


Troubleshooting BAMM runs
..........................

This section addresses a few of the most commonly-encountered problems with BAMM.

Resolving convergence problems
------------------------------

**Under construction**


The model and approximations
....................................

How can BAMM detect diversity-dependent changes in speciation rates?
------------------------------------------------------------------------

BAMM models the dynamics of speciation and extinction within rate regimes using an exponential change function. The speciation rate :math:`\lambda` at any point in time is modeled as

.. math::
	\lambda(t) = \lambda_{0}e^{k t}

where :math:`\lambda` is the initial speciation rate at the start of the rate regime, :math:`k` is a parameter that controls the dynamics of rate change through time, and :math:`t` is the elapsed time since the start of the rate regime. Theoretically, a linear *diversity-dependent* change in speciation rates through time leaves a signal in molecular phylogenies that is virtually indistinguishable from an exponential *time-dependent* change in rates. Our analyses of simulated datasets suggest that these two types of models are not distinguishable in practice. 

We have conducted extensive performance evaluations where we have simulated datasets under formal *diversity-dependent* scenarios, then used BAMM to reconstruct the number of macroevolutionary rate regimes as well as the dynamics of speciation and extinction through time. Our simulations indicate that BAMM can estimate both the number of distinct macroevolutionary regimes, as well as the underlying evolutionary rates, even though we are using the exponential approximation to the diversity-dependent process. We have published these results **here** (non-functional link).
 
It is (vastly) more efficient computationally to work with the exponential change model than the formal diversity-dependent model, and calculations of single likelihoods on phylogenies can be many orders of magnitude faster with the exponential approximation than with the formal diversity-dependent model. The multi-process explorations of macroevolutionary dynamics that are possible with BAMM wouldn't really be feasible without the ability to quickly compute likelihoods. 
 
As an aside, the user is encouraged to remember that all analytically tractable models of diversity-dependence (e.g., Rabosky & Lovette, *Proc. R. Soc. B.*, 2008; or Etienne *et. al.*, *Proc. R. Soc. B*, 2011) are models that we are imposing on the data: there is no reason why a true diversity-dependent process need follow a linear model.
 
 
Other questions
.................

What about the joint diversification-trait evolution process implemented in the 2013 *Nature Communications* paper?
---------------------------------------------------------------------------------------------------------------------------------------------------------
 
**Under construction** 
