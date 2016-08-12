
.. _replication: 

Replication tests of empirical results from Moore et al (2016)
================================================================================

This page is not our full or formal response to the Moore et al (2016) `critique <http://www.pnas.org/content/early/2016/08/09/1518659113.full>`_ of BAMM, but addresses replicability of empirical results from the paper.

Summary
----------
Empirical results from Moore et al (MEA) cannot be replicated with standard usage of BAMM:

* MEA added a **hidden and undocumented** setting to their BAMM control files and modified its default value. This was not acknowledged in their manuscript.
* This setting was a *developer-only* option that changes how the extinction probabilities are calculated in BAMM
* MEA's modification reactivated an implementation concern that they `pointed out <https://github.com/macroevolution/bamm/issues/137>`_ in 2015
* We had fixed this issue in the November 2015 release of BAMM (v2.5)
* However, we released v2.5 with a hidden shortcut for experimentation, which was purposefully activated by MEA in all analyses (it cannot be activated by accident, as the parameter is not included in any template files and the setting they used is not documented)
* When BAMM v2.5 is used with the correct (default) setting, empirical pathologies observed by MEA cannot be replicated.

**We strongly recommend that researchers do not alter the values of hidden/undocumented parameters in BAMM, most of which are intended for internal use by BAMM developers**. Altering these settings may induce pathological behaviors relative to the default/recommended value. 

Details
-----------

Moore et al (MEA) raise a number of important theoretical issues regarding macroevolutionary rate shift models in their recently published critique of BAMM. However, the empirical results from their paper **cannot be replicated without modifying hidden, developer-only settings**. Researchers using BAMM *as recommended by the developers* will be unable to replicate key empirical results from their paper, although we encourage users to try for themselves using the information below.

MEA made a critical modification to standard usage of BAMM in all analyses from their paper for which they have released their input files. In all control files provided, MEA activated a hidden developer setting, ``combineExtinctionAtNodes``, and specified a value that was not intended for routine data analysis. Specifically, the last line of every BAMM control file they provide on `Dryad <http://datadryad.org/resource/doi:10.5061/dryad.mb0sd>`_ to accompany their paper includes the following command::

	combineExtinctionAtNodes = random

This setting is a hidden option that we included to allow developers to easily explore likelihood calculations in previous implementations of BAMM. The default value for this parameter, which was *not visible on any template files*, was::

	combineExtinctionAtNodes = if_different

Our last major update to BAMM (v2.5; November 2015) included a number of major improvements and also addressed several implementation concerns. One significant implementation concern was raised by MEA coauthor Hoehna in his GitHub post `here <https://github.com/macroevolution/bamm/issues/137>`_ from October 2, 2015. Hoehna pointed out a potential concern with the flow of extinction calculations in BAMM during the traversal of the tree from the tips to the root. Importantly, these calculations in previous versions of BAMM could, in some instances, be affected by arbitrary rotation of left and right descendant branches from a given node. We performed extensive testing of several algorithms for likelihood calculations under rate-shift models and identified a `solution <http://bamm-project.org/likelihoodmodel.html#extinction-calculations-at-nodes>`_ that worked well in practice. The solution we have implemented is similar to that used by Etienne and Haegeman in their model for subclade shifts involving diversity-dependent `processes <http://www.journals.uchicago.edu/doi/10.1086/667574>`_.
 
The hidden (and undocumented) non-default settings for ``combineExtinctionAtNodes`` were intended to be developer-only options to allow us, internally, to replicate aspects of the computational algorithm used in earlier (pre-2.5) versions of BAMM. Hence, the control files provided with MEA's Dryad submission override the default setting and force BAMM to perform likelihood calculations in a manner that is related to the calculation flow in earlier implementations. By setting ``combineExtinctionAtNodes`` to (non-default) values of ``random``, ``left``, or ``right``, you are essentially **recreating the implementation concern identified by Hoehna in his GitHub post.** However, the ``random`` setting used by MEA actually leads to *worse performance* than you would have observed for previous versions of BAMM, as this option destabilizes BAMM's likelihood calculations (details :ref:`here<combineExtinctionRandom1>`). As such, we expect that the ``random`` option used by MEA leads to pathologies that would never have been observed for any previous (pre-v2.5) versions of BAMM.
 
To reiterate:
 
* The authors reported bugs and/or implementation concerns to us via GitHub in October 2015

* We made a number of major changes to BAMM and released BAMM v2.5, addressing (among other issues) the concern about extinction probability handling at internal nodes

* All results in MEA were obtained with BAMM v2.5, but required changing a hidden setting that not only recreates but **exacerbates** performance issues due to the implementation concern raised in Hoehna's GitHub post and that are otherwise solved with the default implementation of v2.5.

* This use of a non-standard and undocumented BAMM algorithm was not described in the MEA article.

* No studies of which we are aware, other than MEA, have changed the default setting for ``combineExtinctionAtNodes``
 
We have never documented the setting ``combineExtinctionAtNodes = random``. The only reference to this is buried deep within BAMM's C++ code as a comment to developers (see `these lines <https://github.com/macroevolution/bamm/blob/master/src/SpExModel.cpp#L463-L471>`_ from SpExModel.cpp). In :ref:`this section<empirical1>`, we demonstrate that the empirical analyses from Moore et al are greatly compromised by their use of ``combineExtinctionAtNodes = random``. 
 
.. _empirical1:

Several theoretical concerns raised by MEA cannot be evaluated  
--------------------------------------------------------------------------- 
MEA raise several important theoretical considerations concerning rate-shift models, which will be treated in a more comprehensive response. However, the MEA calculations that purport to demonstrate error in BAMM's extinction probabilities (MEA Fig. 2) and likelihoods (MEA Fig. 3) are not valid, because their calculations are compromised by their ``combineExtinctionAtNodes = random`` setting. If the BAMM likelihood is incorrect, Fig. 2 and Fig. 3 are incapable of illustrating this, because results in those figures are based on compromised analyses that activated a known and documented implementation issue.


Reanalysis of empirical data from Moore et al (2016)
--------------------------------------------------------------------------- 
 
We obtained all input files as `posted <http://datadryad.org/resource/doi:10.5061/dryad.mb0sd>`_ to Dryad by MEA. Every control file used for their BAMM analyses terminates with the line::

	combineExtinctionAtNodes = random

We used BAMM v2.5 to repeat all analyses with the control files from MEA exactly as published by the authors. We then performed a second set of analyses where we re-analyzed the same control files, but where we deleted the command ``combineExtinctionAtNodes = random``, thus ensuring that the program used the correct (default) value (``if_different``).  
 
One of the key results of MEA is their finding that the posterior is overly sensitive to the prior. For every dataset in their article, we find dramatic differences in the shape of the posterior distributions obtained with the incorrect ``random`` (MEA) and the correct ``if_different`` (BAMM default) value of ``combineExtinctionAtNodes``. For example, here are the results for the cetacean dataset that is distributed as an example file with ``BAMMtools``. Each plot shows the posterior distribution on the number of shifts for a given prior parameterization (:math:`\gamma`; this is the same as the ``expectedNumberOfShifts`` parameter in BAMM). Each analysis was performed using MEA's control files **exactly as given** in all respects but one: the top row is with the default value for ``combineExtinctionAtNodes``, and the bottom row shows results using the MEA modification (``random``). 

.. _cetaceans1:  
.. figure:: replication/rr_cetaceans.png
   :width: 1100
   :align: center
 
Compare these results to MEA, Fig. S21. For unmodified BAMM, the posterior is well-behaved with respect to the prior. However, with the MEA modification, the posterior is highly sensitive to the prior. This is especially apparent for large values of :math:`\gamma`. Here is another example dataset, *Adelpha*, corresponding to MEA Figure S19:

.. _adelpha1:  
.. figure:: replication/rr_adelpha.png
   :width: 1100
   :align: center
 
Again, we see the same pattern as with the cetaceans: the posterior is poorly behaved with the MEA modification to the BAMM algorithm. These differences between the BAMM default behavior and the MEA modification can be extreme. Here is *Senna* (MEA fig S29), where the posterior is quite poorly behaved with the MEA modification but where there is virtually no sensitivity to the prior under the default value of ``combineExtinctionAtNodes``:

.. _senna1:  
.. figure:: replication/rr_senna.png
   :width: 1100
   :align: center

And here is another clade, the *Terebinthaceae* (MEA Figure S30):

.. _tere1:  
.. figure:: replication/rr_tere.png
   :width: 1100
   :align: center

This gets a bit repetitive, so we present other empirical results in a separate :ref:`section<empirical1>` below. We turn now to MEA's analysis of constant-rate trees, which was presented in the main text of their article (Figure 4). 

.. _constant1:

Reanalysis of constant-rate trees from Moore et al (2016)
---------------------------------------------------------------------------

MEA simulated constant-rate phylogenies of the same size as the cetacean phylogeny (87 taxa) and found that the posterior was strikingly sensitive to the prior. Here is our reanalysis of their input files, which recreates the results they presented in Figure 4:


.. _cr_mea_mea1:  
.. figure:: replication/cr_mea_mea.png
   :width: 1100
   :align: center

You can see in the figure above that the posterior is virtually identical to the prior for :math:`\gamma` = 10. However, there are two important differences between the way MEA used BAMM and the way the program would typically be used. First, MEA set ``combineExtinctionAtNodes`` to ``random``, as discussed above. Second, for this particular analysis (but not their empirical analyses), MEA placed very strong priors on the rate parameters for speciation and extinction. The next figure considers the effects only of the MEA modification to ``combineExtinctionAtNodes``. When we restore this setting to its default values, using the MEA control files, we see that the posterior is markedly less sensitive to the prior than reported in MEA Figure 4:

.. _cr_2.5_mea1:  
.. figure:: replication/cr_2.5_mea.png
   :width: 1100
   :align: center

There is considerably *less* sensitivity, but there is still *some* sensitivity. This may or may not be of concern. However, we now consider the additional impact of the strong priors placed on the rate parameters by MEA. There are 3 rate priors in the BAMM model: the prior on the speciation rate (``lambdaInitPrior``), the prior on the extinction rate (``muInitPrior``), and the prior on the rate-shift parameter (``lambdaShiftPrior``). The first 2 priors are exponential, the rate shift prior is Gaussian. The function ``setBAMMpriors`` from ``BAMMtools`` matches the scale of these prior distributions to the scale of the tree, such that inferences on relative diversification rates across the tree are independent of tree scale (e.g., you can multiply the branch lengths by 0.001 or 1000 and get identical posterior distributions for *relative* rates, if you use ``setBAMMpriors``). Most researchers either use the default rate priors in BAMM or they use the BAMMtools recommendation (from ``setBAMMpriors``). However, the priors used by MEA on constant-rate trees are much more restrictive than the more liberal priors we recommend. Here is a set of pairwise plots for the recommended versus MEA priors for the 3 rate parameters in BAMM, for their set of 100 constant rate (87 taxon) phylogenies:
 
.. _ratepriors1:  
.. figure:: replication/rate_priors.png
   :width: 900
   :align: center

We acknowledge that there are outstanding issues to be addressed with respect to the sensitivity of BAMM inferences to the underlying priors on speciation and extinction rates. However, the figure above makes clear that all prior values used by MEA are significantly mismatched relative to the recommended values. For ``lambdaInitPrior`` and ``muInitPrior``, the MEA values impose much stronger constraints on the rates (e.g., greater rate parameter of the exponential distribution, as plotted above). The plot above does not adequately convey the range of values for ``muInitPrior``, some of which exceeded 200 (max = 1098) and were not plotted. Moreover, MEA did not attempt to scale the ``lambdaShiftPrior`` and simply used a fixed value of 0.05, despite the dramatic differences in values used for the other priors.

Here is a re-analysis of the constant-rate phylogenies using ``BAMMtools`` recommended values for each dataset (which is how most researchers would use BAMM):

.. _cr_2.5_bammtools:  
.. figure:: replication/cr_2.5_bammtools.png
   :width: 1100
   :align: center
 
We now see that there is even less sensitivity of the posterior to the prior, when using the default setting for ``combineExtinctionAtNodes`` and with the recommended priors from ``BAMMtools::setBAMMpriors``. 

As a final analysis, we considered the results that would be obtained if a researcher simply used off-the-shelf BAMM v2.5 with the default prior settings in the program (e.g., without using ``setBAMMpriors``):

.. _cr_2.5_defaults:  
.. figure:: replication/cr_2.5_defaults.png
   :width: 1100
   :align: center
 
Again, the posterior is not particularly sensitive to the prior when the defaults are used. Thus, claims that BAMM v2.5 is overly sensitive to the prior on the number of rate shifts appears untenable for constant-rate phylogenies.

 
.. _allempirical1:
  
All empirical datasets from Moore et al, reanalyzed
---------------------------------------------------------------------------
 
Reanalysis of the remaining empirical datasets from MEA, as above. The only modification made to MEA's control files was to delete their ``combineExtinctionAtNodes`` argument (top row); bottom row shows results from analyzing their control file exactly as published. Here is the *Byttneria* (MEA Figure S20):

.. _byt1:  
.. figure:: replication/rr_byt.png
   :width: 1100
   :align: center
 
And here is the *Ericaceae* (see MEA Fig S22):

.. _ericaceae1:  
.. figure:: replication/rr_ericaceae.png
   :width: 1100
   :align: center 
 
Here is the *Graphidaceae* (MEA Fig S23):

.. _graphid1:  
.. figure:: replication/rr_graphid.png
   :width: 1100
   :align: center  
 
Here is the *Paphiopedelum* (MEA Fig S25):

.. _paphio1:  
.. figure:: replication/rr_paphio.png
   :width: 1100
   :align: center  
  
Here is the *Parmeliaceae* (MEA Fig S26):

.. _parm1:  
.. figure:: replication/rr_parm.png
   :width: 1100
   :align: center  
   
Here is the *Pleopeltis* (MEA Fig S27):
   
.. _pleo1:  
.. figure:: replication/rr_pleo.png
   :width: 1100
   :align: center  
   
Here is the *Polygoneae* (MEA Fig S28):
   
.. _poly1:  
.. figure:: replication/rr_poly.png
   :width: 1100
   :align: center  
   
Here is the *Turnera* (MEA Fig S31):
   
.. _turn1:  
.. figure:: replication/rr_turn.png
   :width: 1100
   :align: center  
   
And the *Viburnum* (MEA Fig S32):
   
.. _viburnum1:  
.. figure:: replication/rr_viburnum.png
   :width: 1100
   :align: center  
   
 
.. _combineExtinctionRandom1: 

Appendix: The "random" option used by Moore et al is not valid
---------------------------------------------------------------------------

Moore et al (2016) changed a hidden developer setting in BAMM v2.5 for all analyses performed in their paper. Specifically, they changed the default value of the parameter ``combineExtinctionAtNodes`` to the value ``random``. This option is related to the implementation concern raised by Hoehna's GitHub post, but does not exactly re-create it. To recreate the precise flow of extinction calculations as used in earlier (pre v2.5) versions of BAMM, you would set::

	combineExtinctionAtNodes = left

The reasons for the current default value of ``if_different`` are explained in detail in this `subsection <http://bamm-project.org/likelihoodmodel.html#extinction-calculations-at-nodes>`_, although there is no documentation of *random* in any previous BAMM documentation. The ``random`` option is related to ``left``, but is even more arbitrary. Unfortunately, the option has the effect of destabilizing likelihood calculations with BAMM, because it stochastically assigns inheritance of extinction probabilities at internal nodes to favor the right or left descendant node. Hence, the likelihood calculations with ``random`` depend not only on arbitrary node rotations, but will be assigned randomly each time a tree is loaded into BAMM. 

This means that likelihoods computed for the same dataset, with the same parameters, **are not guaranteed to be identical**. In fact, if there are any rate shifts on the tree, it is likely that *independent calculations of the likelihood with the same parameters will yield different values*. As a simple example, we will use the ``random`` option to compute the likelihood a set of parameters from the posterior distribution of rate regimes simulated for the whale dataset distributed with ``BAMMtools``. 

``BAMMtools`` includes a simple likelihood calculator that allows us to compute tree likelihoods under the ``random`` option, although this option *is not documented anywhere*; the default, as for BAMM, is ``if_different``. The argument for this option in ``BAMMtools::BAMMlikelihood`` is ``e_prob_condition`` (this is identical to ``combineExtinctionAtNodes`` in BAMM). Here, we will pull out the last 50 generations from the posterior sampled with BAMM; for each generation, we will compute the likelihood twice. Now, the parameters sampled for a particular generation are constant: hence, **all calculations of the likelihood with those values must be identical**, or the method is not mathematically coherent. ::

	library(BAMMtools)
	data(whales, events.whales)
	ux <- unique(events.whales$generation)[1951:2000]
	ll <- matrix(NA, nrow=50, ncol=2)
	for (i in 1:50){
		cat(i, "\n")
		ll[i,1] <- BAMMlikelihood(whales, events.whales, gen= ux[i], e_prob_condition = "random")
		ll[i,2] <- BAMMlikelihood(whales, events.whales, gen= ux[i], e_prob_condition = "random")
 
	}

We now have a 2 column matrix, ``ll``, where the right and left columns correspond to independent calculations of the likelihood with the same parameter sets. We will now plot ``ll[ ,1]`` against ``ll[ ,2]``:

.. _randomCompare:  
.. figure:: replication/randomCompare.png
   :width: 400
   :align: center
 
The likelihoods from successive computations are not identical. They are correlated, but they *should be identical*, as the parameters and tree are identical. Hence, ``random`` is perhaps the worst of several suboptimal settings that should not be used.
 
Just to demonstrate that ``if_different`` does not show this pathology, here will we repeat this exercise but using ``combineExtinctionAtNodes = if_different`` (recall that this parameter in BAMMtools has the label ``e_prob_condition``) ::

	library(BAMMtools)
	data(whales, events.whales)
	ux <- unique(events.whales$generation)[1951:2000]
	ll <- matrix(NA, nrow=50, ncol=2)
	for (i in 1:50){
		cat(i, "\n")
		ll[i,1] <- BAMMlikelihood(whales, events.whales, gen= ux[i], e_prob_condition = "if_different")
		ll[i,2] <- BAMMlikelihood(whales, events.whales, gen= ux[i], e_prob_condition = "if_different")
 
	}

And, as *should* be the case for any mathematically consistent method, the likelihoods are identical:

.. _ifdiff:  
.. figure:: replication/if_differentCompare.png
   :width: 400
   :align: center
 

