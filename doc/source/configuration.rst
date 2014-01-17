.. highlight:: none

Configuration
=============


Control File
------------

Configuration options and parameters in BAMM are specified in a *control file*,
a plain text file in which each line contains the name of the option or
parameter, an equal sign, and the value of the option or parameter::

    treefile = whaletree.tre
    runMCMC = 1

The path of the control file, relative to the directory in which ``bamm``
is called, is specified with the flag ``-c``::

    ./bamm -c divcontrol.txt


General Options and Parameters
------------------------------

The following describes the configuration options and parameters
that are required regardless of the specific model used.
For true or false values, 1 is used for true and 0 is used for false.
File paths are relative to the directory in which ``bamm`` is called.

General
.......

modeltype
    The type of model BAMM uses for analysis.
    If ``speciationextinction``, run with a speciation/extinction model.
    If ``trait,`` run with a phenotypic evolution model.

treefile
    The file path of the input tree (in Newick format).
    For speciation/extinction analyses, the tree must be ultrametric and fully bifurcating.

runInfoFilename
    The path of the file to output general information about the current run.

sampleFromPriorOnly
    If ``1``, run by sampling from the prior only
    (ignoring likelihood contribution to posterior).
    If ``0``, run the full analysis.

autotune
    Experimental option for tuning MCMC operators.

runMCMC
    If ``1``, run the MCMC sampler.
    If ``0``, just check to see if the data can be loaded correctly.

simulatePriorShifts
    If ``1``, simulate prior distribution of the number of shift events,
    given the hyperprior on the Poisson rate parameter.
    This is necessary to compute the Bayes factor.

loadEventData
    If ``1``, load the event data from a file.

eventDataInfile
    The file path of the event data (used only if ``loadEventData = 1``).

initializeModel
    If ``1``, initializes MCMC.
    If ``0``, just check parameter file and ensure that data can be read.

seed
    Seed for the random number generator.
    If ``-1``, the seed is obtained from the clock time.
  
overwrite
    If ``1``, overwrite output files if they already exist with identical
    filenames.
    If ``0``, do not run if output files already exist with identical filenames.

MCMC Simulation
...............

numberGenerations
    Number of MCMC generations to run.

mcmcWriteFreq
    Frequency (in generations) at which to print MCMC details
    to ``mcmcOutfile``.

eventDataWriteFreq
    Frequency (in generations) at which to print event details
    to ``eventDataOutfile``.

acceptWriteFreq
    Frequency (in generations) at which to print accept rate to the screen.

printFreq
    Frequency (in generations) at which to print output to the screen.
  
outName
    If present (may be commented out), prefixes output files with the given
    text.

mcmcOutfile
    The path of the file to which to write the MCMC output.

eventDataOutfile
    The path of the file to which to write the raw event data.
    All of the results are contained in this file, and all branch-specific
    speciation rates, shift positions, marginal distributions, etc.,
    can be reconstructed from this output. See :ref:`bammtools`
    for more information on working with this output format.
  
updateEventRateScale
  *Description not yet available.*

localGlobalMoveRatio
  *Description not yet available.*

acceptrateOutfile
  *Description not yet available.*

Parameter Update Rates
......................

updateRateEventNumber
  Frequency of updating the number of events (shifts) on the tree.

updateRateEventPosition
  Frequency of moving the position of a shift point.

updateRateEventRate
  Frequency of updating the rate at which events occur.

initialNumberEvents
  *Description not yet available.*

Speciation/Extinction Model
---------------------------

The following describes configuration options and parameters
specifically for speciation/extinction analyses.

General
.......

useGlobalSamplingProbability
  If true (1), will look for a global correction for incomplete sampling (globalSamplingProbability)
  If false (0), will look for a file that specifies clade-specific corrections for incomplete sampling (sampleProbsFilename).
  
globalSamplingProbability
  Percentage of total number of species sampled in given phylogeny (between 0 and 1).

sampleProbsFilename
  Specifies a file with clade-specific corrections for incomplete sampling.

Priors
......

lambdaInitPrior
  Mean of the exponential distribution prior on speciation.

lambdaShiftPrior
  Prior on the speciation rate change parameter.

muInitPrior
  Exponential prior on extinction.

muShiftPrior
  Prior on the extinction rate change parameter.

segLength
  *Description not yet available.*

General MCMC Simulation Settings & Output Options
..................................................

lambdaOutfile
  Branch-specific speciation rates will be written to this file as Newick-formatted trees.

muOutfile
  Branch-specific extinction rates will be written to this file as Newick-formatted trees.

lambdaNodeOutfile
  *Description not yet available.*

updateLambdaInitScale
  *Description not yet available.*

updateMuInitScale
  *Description not yet available.*

updateLambdaShiftScale
  *Description not yet available.*

updateMuShiftScale
  *Description not yet available.*

minCladeSizeForShift
  *Description not yet available.*

Starting Parameters
...................

lambdaInit0
  Starting initial speciation rate.

lambdaShift0
  Starting initial rate change parameter for speciation
  (if 0, speciation rates will not change through time).
  A negative value implies decreasing rates through time.

muInit0
  Starting Initial extinction rate.

muShift0
  Starting initial rate change parameter for extinction. Currently not implemented.

Parameter Update Rates
......................

updateRateLambda0
  Frequency in which to update the initial speciation rate for an event.

updateRateLambdaShift
  Frequency in which to update how speciation rates change through time.

updateRateMu0
  Frequency in which to update the initial extinction rate.


Phenotypic Evolution Model
--------------------------

The following describes the configuration options and parameters
specifically for the phenotypic evolution model in BAMM.
The parameter "beta" represents the rate of phenotypic evolution
at any point in time.

General
.......

traitfile
  File that names the trait data. Traits must be continuous characters.
  Each line must have a species name and the corresponding trait,
  separated by a tab.
  No header row is permitted.
  All species in the trait data file must be in the tree and vice versa.

MCMC Tuning
...........

updateBetaScale
  Controls the amount by which to change the value of beta
  at any step in the MCMC sampling.

updateNodeStateScale
  *Description not yet available.*

updateBetaShiftScale
  *Description not yet available.*

Starting Parameters
...................

betaInit
  Starting initial rate.

betaShiftInit
  Starting initial rate change parameter for phenotypic evolution.
  (if 0, then constant-rate).
  A negative value implies decreasing rates through time.

Priors
......

betaInitPrior
  *Description not yet available.*

betaShiftPrior
  *Description not yet available.*

useObservedMinMaxAsTraitPriors
  *Description not yet available.*

traitPriorMin
  *Description not yet available.*

traitPriorMax
  *Description not yet available.*

Output
......

betaOutfile
  The file name in which to write the phenotypic rates as newick-formatted trees where the branches are scaled to the rate of phenotypic evolution.

nodeStateOutfile
  *Description not yet available.*

Parameter Update Rates
......................

updateRateBeta0
  *Description not yet available.*

updateRateBetaShift
  *Description not yet available.*

updateRateNodeState
  Relative rate at which to update individual node state values.
  This value should, in general, be substantially higher
  than the other parameter values (recommended 25:1 or 50:1)
  because there are so many internal nodes states that need to be updated.
