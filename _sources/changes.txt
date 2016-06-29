:orphan:

.. _changes:

BAMM Changes
============

2.5.0
-----
*Nov 4, 2015*

New features and enhancements
.............................

* The prior distribution on the number of shifts is now specified with parameter ``expectedNumberOfShifts``. This is simply the inverse of the previously-used parameter ``poissonRatePrior``. BAMM control files will accept either parameter name but we find that ``expectedNumberOfShifts`` is more intuitive.

* No more prior simulation in BAMM. The relevant BAMMtools functions now accept an argument for ``expectedNumberOfShifts`` and compute the full prior distribution from that information.

* We made a major change to the algorithm for handling extinction probabilities :math:`E(t)` at internal nodes when the two branches descended from a node differ in their shift histories. We don't consider the previous approach a bug *per se*, but have convinced ourselves that the new way has a better theoretical justification. We explain why we made this change :ref:`here <extinctionNodes>`. 


BAMMtools enhancements
............................. 

* We have added a new function to BAMMtools, ``BAMMlikelihood``, which is an R-based tool for replicating the calculations implemented in BAMM. This is described :ref:`here <testlikelihood>`.

* Changes to BAMMtools rate plotting functions, including improved histograms, color break options, and more.

* New functions, including ``generateControlFile`` (to generate BAMM control files from within R), ``plotPrior`` (to visualize the prior and posterior distributions on the number of rate shifts), and ``ratesHistogram`` (to visualize the distribution of rates across a phylogeny).

* No more "branch specific Bayes factors" for distinguishing between core and non-core shifts. We have replaced this terminology with a related concept - the *marginal odds ratio*. Note that neither "branch specific Bayes factors" nor "marginal odds ratios" can be used for formal model selection; we have explained this in detail on our conceptual page about the interpretation of `rate shifts <rateshifts.html>`_. 

Bug fixes
............................. 

* the option ``validateEventConfiguration`` is an internal debugging option that should have had a default value of 0, but was instead set to 1 (this led to a bug in the Hastings ratio that would have affected a very small fraction of MCMC moves).
 
* Fixed bug introduced during programming of fossil BAMM that sometimes recomputed extinction probabilities :math:`E(t)`. An explanation for why extinction probabilities must account for downstream shift histories is found :ref:`here <whatprocess>` (this bug was not incorporated into the compiled version of BAMM that was distributed on the website but could have affected some analyses where users compiled BAMM code themselves).   


2.4.0
-----

*June 13, 2015*

* Caught a bug introduced during programming of fossil BAMM that allowed numerical overflow issues to occur to the extinction probability when E(t) approached 1.0 (this bug was not incorporated into the compiled version of BAMM that was distributed on the website but could have affected some analyses where users compiled BAMM code themselves). 

*May 26, 2015*

* Fixed bug in acceptance probability for MCMC moves that update the Poisson rate parameter. This is the most serious bug we have encountered in BAMM and **would have magnified the effects of the prior on the posterior distribution of the number of shifts**. If you previously used the default value of ``poissonRatePrior = 1.0``, you will probably notice little effect on inference, but *it may have had an impact in some cases*. We are grateful to Cécile Ané and Bret Larget for helping us solve this.  
 

2.2.0
-----

*September 5, 2014*

* Parallelism is now implemented using C++11 threads

  * No dependency on OpenMP libraries
  * Requires C++11 compiler (if building from source)

* Provide a tool that helps choosing a ``deltaT`` value
  (available from GitHub), and it is documented in :ref:`mc3`.

* Add ``--version`` option in ``bamm`` and update "Usage" message

* Fix bug in initialization of trait values

2.1.0
-----

*July 21, 2014*

* Add proposal to change the number of events
  by picking a random branch rather than a random location

  * To use this proposal, add the setting ``updateRateEventNumberForBranch``
    and give it a value greater than 0 (e.g., 0.1)
  * It is recommended that the setting ``updateRateEventNumber``
    be set to 0 when this new proposal is used

* Fix log-likelihood calculation when loading event data
  or adding initial events

2.0.0
-----

*June 11, 2014*

New features and enhancements
.............................

* Metropolis-coupled MCMC

  * Implementation from Altekar et al 2004
    (*Parallel Metropolis coupled Markov chain Monte Carlo
    for Bayesian phylogenetic inference*)
  * In the control file, specify the number of chains with ``numberOfChains``,
    the temperature increment parameter (delta T) with ``deltaT``,
    and the number of generations to wait for each chain swap proposal
    with ``swapPeriod``
  * Supports running chains in parallel with OpenMP
  * Output acceptance information for chain swaps to the file
    specified by ``chainSwapFileName``

* Time-flip proposal

  * Allows shift events to change from time-variable to time-constant
    and vice versa
  * Enable this proposal with ``updateRateLambdaTimeMode`` for the
    speciation/extinction model type and with ``updateRateBetaTimeMode``
    for the trait model type

* Command-line options

  * Change any setting at the command-line; for example, to change
    the number of generations to 1000, add ``--numberOfGenerations 1000``
    as a command-line option
  * Any setting in the control file can be changed, and the new value
    specified at the command-line has precedence over that in the control file

* Use different rate functions for the rate parameter *k* > 0 and for *k* < 0;
  when *k* > 0, use the symmetrical function of *k* < 0
* Accept tree files with internal node names
* Display better error messages related to reading a tree

Bug fixes
.........

* Properly prefix file names if it includes a directory
* Fix detection of negative branches
* Fix simulating prior distributions for small number of generations

Other changes
.............
* ``numberGenerations`` is now ``numberOfGenerations``
* ``updateBetaScale`` is now ``updateBetaInitScale``
* Update example control files for this version of BAMM

1.0.0
-----

*March 5, 2014*

* Initial release
