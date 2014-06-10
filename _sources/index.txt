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
and extended in several publications (`PLoS ONE 2014 <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0089543>`_ ,  `Nature Communications 2013 <http://www.nature.com/ncomms/2013/130606/ncomms2958/full/ncomms2958.html>`_ , and `Systematic Biology 2014 <http://sysbio.oxfordjournals.org/content/early/2014/04/01/sysbio.syu025>`_). BAMM is a command line program written in C++. Post-run analysis and visualization is performed using
the R package `BAMMtools <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12199/abstract>`_

- `Download BAMM and BAMMtools <download.html>`_ or go to our
  `GitHub page <https://github.com/macroevolution/bamm>`_
  to get the development source code.

- Explore the `Graph Gallery <bammgraph.html>`_ for a sample of analyses
  produced using BAMM and BAMMtools.

- Quickly start using and analyzing data with BAMM by reading the
  `Quick-start Guide to BAMM <quickstart.html>`_.

- Go to our `Frequently Asked Questions <faq.html>`_ page to see common
  questions and answers.


BAMM and BAMMtools have been updated!
=====================================

We've made a number of changes - some major - to both BAMM and BAMMtools that significantly improve performance and reliability. If you have used a previous version of BAMM or BAMMtools, we recommend installing the latest version. The documentation on this site is oriented towards BAMM v2.0 and BAMMtools v2.0.
  
Major changes 
========================== 

- Implementation of a flip model (see :ref:`timefliptheory`) for rate shifts. Rate shifts are not assumed to lead to either time-constant processes (as in the previous BAMM version) or time-constant processes (as in MEDUSA); rather, a new class of reversible jump moves allows diversification dynamics to toggle between time-constant and time-variable modes.
  **Note: The time-flip proposal is still being tested for reliability,
  so it is currently turned off in the template and example control files.**

- Metropolis coupling (MC3) implemented by default to facilitate convergence (automatically run in parallel if you have OpenMP). The MC3 is described :ref:`here <mc3>`.

- Comprehensive overhaul of BAMM's C++ core for transparency and extensibility

- BAMMtools 2.0 uses branch-specific Bayes factors to identify credible sets of :ref:`shift configurations <rateshifts>`


Support
=======

The development of BAMM is funded by the National Science Foundation.

.. figure:: nsf-logo.gif
   :width: 58
