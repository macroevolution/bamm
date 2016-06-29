Frequently Asked Questions
==========================
 
General
.......

How do I cite BAMM?
-------------------------------

The primary methodological description of the full BAMM model was published in *PLoS ONE*. A link to the paper is `available here <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0089543>`_. 

This paper contains a description of the model, the reversible jump MCMC implementation, a comprehensive performance evaluation, and an empirical application. An earlier paper used a related implementation that assumed constant rates of diversification through time within evolutionary rate regimes (`Rabosky et al. 2013 <http://www.nature.com/ncomms/2013/130606/ncomms2958/full/ncomms2958.html>`_). 

Macroevolutionary cohort analyses were described and published in `Systematic Biology 2014 <http://sysbio.oxfordjournals.org/content/63/4/610>`_. The MC3 implementation of BAMM, and the usage of Bayes factors for assessing evidence in favor of rate shifts, were published in `Evolution 2015 <http://onlinelibrary.wiley.com/doi/10.1111/evo.12681/abstract>`_. 

If you use BAMMtools for analysis and visualization of BAMM results, you should cite the `BAMMtools package <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12199/abstract>`_.

How much data does it need?
---------------------------
**Not much**. There is no general rule here, but if your dataset is large enough to consider doing any other sort of diversification analysis, then it is probably large enough for BAMM. The whale diversification analysis shown in the :ref:`graph gallery<bammgraphs>` uses a time-calibrated tree with 89 tips. We've had good success using BAMM on trees that are considerably smaller than this. 



.. _howmanyevents:

How often should I sample the chain?
------------------------------------

Your first consideration should be: for how long should I run the chain? The answer to this question has a lot to do with how quickly the chain converges and how well it samples the posterior. We suggest that there is little advantage in going to more than 5000 samples from the posterior, and potentially large costs: for large datasets, sampling too frequently with BAMM can literally generate gigabytes of highly autocorrelated and mostly-unusable output. Moreover, due to memory issues with R, BAMMtools cannot handle output files of that size. So, pick a chain length (e.g., :math:`10^7` generations), and then specify a sample frequency that should give you somewhere between 500 and 5000 samples from the posterior after burnin. 


How many generations should I expect to run BAMM for?
------------------------------------------------------

This section is just to provide some ballpark numbers and general guidelines. For relatively small datasets (< 1000 taxa), you may observe good performance with less than 10 million generations. For datasets in the thousands of taxa, don't be surprised if you need 100 million or more generations to adequately sample the posterior. BAMM 2.0 includes Metropolis coupling (MC3), which we have found to considerably boost performance. See documentation on this feature to ensure that you are using it to its full potential.


What extensions are underway for BAMM?
--------------------------------------

We are currently testing extensions of BAMM that allow modeling evolutionary dynamics under a much greater range of phenotypic evolutionary scenarios, as well as the incorporation of paleontological data. Bleeding edge releases of BAMM can be obtained from our `GitHub page <https://github.com/macroevolution/bamm>`_.

What is the actual model implemented in BAMM for time-varying rates?
--------------------------------------------------------------------

The model originally described for time-varying macroevolutionary rate regimes (e.g., `Rabosky 2014
<http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0089543>`_) used a simple exponential change function for rates through time:

.. math::
	\lambda(t) = \lambda_{0}e^{k t}

We have since come to recognize that this model is not ideal, because it is highly asymmetric with respect to increasing (k > 0) and decreasing (k < 0) rates through time. Specifically, the prior we place on *k* is symmetric (normal), but the rates are not. The problem, put simply, is that negative values of *k* can only decrease the evolutionary rate to zero. However, positive values - even of equivalent absolute value - can lead to extremely large values. Just consider the values of the exponential function above evaluated at :math:`\lambda_{0} = 1`, :math:`t = 1`, :math:`k = -10` and :math:`k = 10`.  With *k = -10*, the value of the rate function approaches zero, but with *k = +10*, it is 22026.47. Yet both of these *k* values have equal prior densities. 

We thus modified the exponential change function to make it symmetric. For *k < 0*, we use the exponential change function defined above. For *k > 0*, we use the function

.. math::
	\lambda(t) = \lambda_{0} (2 - e^{-k t})
	
If *k > 0*, this function asymptotically approaches the limiting value :math:`2\lambda_{0}`. 

This has the nice property that the function :math:`\lambda(t)` is symmetric about the line :math:`x = \lambda_{0}`, for positive and negative values of k with the same absolute value.

Note that in BAMM v2.0, the model allows "flipping" between time-constant and time-varying diversification modes. The theory underlying this new proposal class is described :ref:`here <timefliptheory>`. Hence, the model does not assume that diversification shifts result in a new time-varying diversification process, but allows for such processes if they are supported by the data.

	

How can BAMM detect diversity-dependent changes in speciation rates?
--------------------------------------------------------------------

BAMM models the dynamics of speciation and extinction within rate regimes using an exponential change function described in the preceding section. Theoretically, a linear *diversity-dependent* change in speciation rates through time leaves a signal in molecular phylogenies that is virtually indistinguishable from an exponential *time-dependent* change in rates. Our analyses of simulated datasets suggest that these two types of models are not distinguishable in practice. 

We have conducted extensive performance evaluations where we have simulated datasets under formal *diversity-dependent* scenarios, then used BAMM to reconstruct the number of macroevolutionary rate regimes as well as the dynamics of speciation and extinction through time. Our simulations indicate that BAMM can estimate both the number of distinct macroevolutionary regimes, as well as the underlying evolutionary rates, even though we are using the exponential approximation to the diversity-dependent process. We have published these results `here <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0089543>`_.
 
It is (vastly) more efficient computationally to work with the exponential change model than the formal diversity-dependent model, and calculations of single likelihoods on phylogenies can be many orders of magnitude faster with the exponential approximation than with the formal diversity-dependent model. The multi-process explorations of macroevolutionary dynamics that are possible with BAMM wouldn't really be feasible without the ability to quickly compute likelihoods. 
 
As an aside, the user is encouraged to remember that all analytically tractable models of diversity-dependence (e.g., Rabosky & Lovette, *Proc. R. Soc. B.*, 2008; or Etienne *et. al.*, *Proc. R. Soc. B*, 2011) are models that we are imposing on the data: there is no reason why a true diversity-dependent process need follow a linear model.



Why does my phylorate plot not show variation when I know there are significant rate shifts in my dataset?
----------------------------------------------------------------------------------------------------------

It is possible that the distribution of rates across your phylogeny is highly skewed. See :ref:`this page <colorbreaks>` for more information. 


