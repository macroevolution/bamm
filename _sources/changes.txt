:orphan:

BAMM Changes
============

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
