:orphan:

Welcome
=======

BAMM (Bayesian Analysis of Macroevolutionary Mixtures) is a program for
modeling complex dynamics of speciation, extinction, and trait evolution on
phylogenetic trees.

BAMM is oriented entirely towards detecting and quantifying heterogeneity in
evolutionary rates. It uses reversible jump Markov chain Monte Carlo to
automatically explore a vast universe of candidate models of lineage
diversification and trait evolution. BAMM and associated methods have been described
and extended in several publications (`PLoS ONE 2014 <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0089543>`_ ,  `Nature Communications 2013 <http://www.nature.com/ncomms/2013/130606/ncomms2958/full/ncomms2958.html>`_ , `Systematic Biology 2014 <http://sysbio.oxfordjournals.org/content/63/4/610>`_, and `Evolution 2015 <http://onlinelibrary.wiley.com/doi/10.1111/evo.12681/abstract>`_). BAMM is a command line program written in C++. Post-run analysis and visualization is performed using
the R package `BAMMtools <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12199/abstract>`_.

- `Download BAMM and BAMMtools <download.html>`_ or go to our
  `GitHub page <https://github.com/macroevolution/bamm>`_
  to get the development source code.

- Explore the `Graph Gallery <bammgraph.html>`_ for a sample of analyses
  produced using BAMM and BAMMtools.

- Quickly start using and analyzing data with BAMM by reading the
  `Quick-start Guide to BAMM <quickstart.html>`_.

- Go to our `Frequently Asked Questions <faq.html>`_ page to see common
  questions and answers.

Recommended reading
===========================

* We strongly recommend that you read our documentation on how to (and how not to) interpret `rate shifts <rateshifts.html>`_ on phylogenies. This section addresses some of the most common pitfalls with interpreting BAMM analyses.  

* **What are the assumptions of the BAMM likelihood?** Find out `here <likelihoodmodel.html>`_

* **Is the BAMM likelihood computed correctly**, given the assumptions of the model? Learn more about testing it :ref:`here <testlikelihood>`



Recent changes
=================

We've made a number of changes - some major - to both BAMM and BAMMtools
that significantly improve performance and reliability.
If you have used a previous version of BAMM or BAMMtools,
we recommend installing the latest version.

* BAMMtools now computes the prior analytically for the compound Poisson process model; there is thus no need to simulate the prior distribution of the number of shifts. More on this :ref:`here<analyticalprior>`.

* BAMMtools 2.1 uses branch-specific marginal odds ratios to identify
  credible sets of :ref:`shift configurations <rateshifts>`. More about why we made this change :ref:`here<marginalodds>`.

* Users can now use BAMMtools to test whether BAMM is correctly computing likelihoods (see :ref:`here<testlikelihood>`). 

* Some important bug fixes are documented :ref:`here<changes>`.  

* Comprehensive overhaul of BAMM's C++ core for transparency
  and extensibility
  
* Metropolis coupled MCMC implemented by default to facilitate convergence.
  The MC3 is described :ref:`here <mc3>`.

`Take a look <colorbreaks.html>`_ at a new webpage that explains some of the intricacies of phylorate plot interpretation.

Please see the `Changes <changes.html>`_ page for more information.

Support
=======

The development of BAMM is funded by the National Science Foundation.

.. figure:: nsf-logo.gif
   :width: 58
