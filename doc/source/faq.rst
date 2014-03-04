Frequently Asked Questions
==========================
 
General
.......

Citing BAMM
-----------

The primary methodological description of the BAMM model was published in *PLoS ONE*. A link to the paper is `available here <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0089543>`_ 

This paper contains a description of the model, the reversible jump MCMC implementation, a comprehensive performance evaluation, and an empirical application.


How much data does it need?
---------------------------
**Not much**. There is no general rule here, but if your dataset is large enough to consider doing any other sort of diversification analysis, then it is probably large enough for BAMM. The whale diversification analysis shown in the :ref:`graph gallery<bammgraph>` uses a time-calibrated tree with 89 tips. We've had good success using BAMM on trees that are considerably smaller than this. 


How often should I sample the chain?
------------------------------------
Your first consideration should be: for how long should I run the chain? The answer to this question has a lot to do with how quickly the chain converges and how well it samples the posterior. We suggest that there is little advantage in going to more than 5000 samples from the posterior, and potentially large costs: for large datasets, sampling too frequently with BAMM can literally generate gigabytes of highly autocorrelated and mostly-unusable output. Moreover, due to memory issues with R, BAMMtools cannot handle output files of that size. So, pick a chain length (e.g., :math:`10^7` generations), and then specify a sample frequency that should give you somewhere between 500 and 5000 samples from the posterior after burnin. 



BAMM extensions
---------------

We are currently testing extensions of BAMM that allow modeling evolutionary dynamics under a much greater range of phenotypic evolutionary scenarios, as well as the incorporation of paleontological data. Bleeding edge releases of BAMM can be obtained from here (**link**).


How can BAMM detect diversity-dependent changes in speciation rates?
--------------------------------------------------------------------

BAMM models the dynamics of speciation and extinction within rate regimes using an exponential change function. The speciation rate :math:`\lambda` at any point in time is modeled as

.. math::
	\lambda(t) = \lambda_{0}e^{k t}

where :math:`\lambda` is the initial speciation rate at the start of the rate regime, :math:`k` is a parameter that controls the dynamics of rate change through time, and :math:`t` is the elapsed time since the start of the rate regime. Theoretically, a linear *diversity-dependent* change in speciation rates through time leaves a signal in molecular phylogenies that is virtually indistinguishable from an exponential *time-dependent* change in rates. Our analyses of simulated datasets suggest that these two types of models are not distinguishable in practice. 

We have conducted extensive performance evaluations where we have simulated datasets under formal *diversity-dependent* scenarios, then used BAMM to reconstruct the number of macroevolutionary rate regimes as well as the dynamics of speciation and extinction through time. Our simulations indicate that BAMM can estimate both the number of distinct macroevolutionary regimes, as well as the underlying evolutionary rates, even though we are using the exponential approximation to the diversity-dependent process. We have published these results **here** (non-functional link).
 
It is (vastly) more efficient computationally to work with the exponential change model than the formal diversity-dependent model, and calculations of single likelihoods on phylogenies can be many orders of magnitude faster with the exponential approximation than with the formal diversity-dependent model. The multi-process explorations of macroevolutionary dynamics that are possible with BAMM wouldn't really be feasible without the ability to quickly compute likelihoods. 
 
As an aside, the user is encouraged to remember that all analytically tractable models of diversity-dependence (e.g., Rabosky & Lovette, *Proc. R. Soc. B.*, 2008; or Etienne *et. al.*, *Proc. R. Soc. B*, 2011) are models that we are imposing on the data: there is no reason why a true diversity-dependent process need follow a linear model.
