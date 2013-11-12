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

The path of the control file (relative to the current directory) is specified
with the flag ``-control`` when running bamm::

    ./bamm -control divcontrol.txt


Global Options and Parameters
-----------------------------

The following describes the configuration options and parameters
that are required regardless of the specific model used.
For true or false values, 1 is used for true and 0 is used for false.
All file paths are relative to the current directory.

General
.......

treefile
  The file name of the input tree (in Newick format).

sampleFromPriorOnly
  If true (1), run BAMM by sampling from the prior only
  (ignoring likelihood contribution to posterior).
  If false (0), run the full analysis.

runMCMC
  If true (1), run the MCMC sampler.
  If false (0), just check to see if the data can be loaded correctly.

loadEventData
  *Description not yet available.*

eventDataInfile
  *Description not yet available.*

initializeModel
  Should always be true (1).

NumberGenerations
  Number of MCMC generations to run.

MCMC Tuning
...........

MeanSpeciationLengthFraction
  *Description not yet available.*

updateEventRateScale
  *Description not yet available.*

localGlobalMoveRatio
  *Description not yet available.*

Priors
......

targetNumber
  Expeced number of "events" or rate shifts on the tree,
  if there is no signal in the data.

Output
......

mcmcOutfile
  The file name in which to write the MCMC parameter output.

acceptrateOutfile
  *Description not yet available.*

eventDataOutfile
  *Description not yet available.*

treeWriteFreq
  Frequency (in generations) in which to write speciation/extinction rates.
  To avoid very large files, use a frequency of at least 10000.

mcmcWriteFreq
  *Description not yet available.*

eventDataWriteFreq
  *Description not yet available.*

acceptWriteFreq
  *Description not yet available.*

printFreq
  Frequency (in generations) in which to write to the screen.

Parameter Update Rates
......................

updateRateEventNumber
  Frequency in which to update the number of events (shifts) on the tree.

updateRateEventPosition
  Frequency in which to move the position of a shift point.

updateRateEventRate
  Frequency in which to update the rate at which events occur.

initialNumberEvents
  *Description not yet available.*


Speciation/Extinction Model
---------------------------

The following describes the configuration options and parameters
for the speciation/extinction model in BAMM.
To use this model, run the executable called ``bamm``.

General
.......

useGlobalSamplingProbability
  If true (1), use global correction for incomplete taxon sampling
  in the likelihood calculation.

globalSamplingFraction
  Fraction (0.0 - 1.0) of the total number of species in the clade
  that are in the tree being analyzed.

MCMC Tuning
...........

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
  Starting speciation rate.

lambdaShift0
  Starting rate change parameter for speciation
  (if 0, speciation rates will not change through time).

muInit0
  Initial extinction rate.

muShift0
  *Description not yet available.*

Priors
......

lambdaInitPrior
  Mean of the exponential distribution prior on speciation.

lambdaShiftPrior
  Prior on the speciation rate change parameter.

muInitPrior
  Exponential prior on extinction.

muShiftPrior
  *Description not yet available.*

segLength
  *Description not yet available.*

Output
......

lambdaOutfile
  The file name in which to write branch-specific speciation rates.

muOutfile
  The file name in which to write branch-specific extinction rates.

lambdaNodeOutfile
  *Description not yet available.*

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
for the phenotypic evolution model in BAMM.
The parameter "beta" represents the rate of phenotypic evolution
at any point in time.
To use this model, run the executable called ``bammtrait``.

General
.......

traitfile
  The file name the trait data, relative to the current directory.
  The traits must be continuous characters.
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
  Initial rate.

betaShiftInit
  Initial time-dependent shift.
  If negative, it implies a decrease in the rate of phenotypic evolution
  through time.

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
  The file name in which to write the phenotypic rates.

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
