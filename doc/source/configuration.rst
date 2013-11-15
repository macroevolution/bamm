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
with the flag ``-control`` when running bamm. For example::

    ./bamm speciationextinction -control divcontrol.txt


Global Options and Parameters
-----------------------------

The following describes the configuration options and parameters
that are required regardless of the specific model used.
For true or false values, 1 is used for true and 0 is used for false.
Files will be relative to the current directory unless otherwise
specified by the user.

General
.......

modeltype
   Should be "speciationextinction" or "trait".
   Defines which application of BAMM will un.

treefile
  The file name of the input tree (in Newick format). For diversification analyses, this should be ultrametric.

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
  Expected number of "events" or rate shifts on the tree,
  if there is no signal in the data.

Output
......

mcmcOutfile
  MCMC parameter output will be written to this file.

acceptrateOutfile
  *Description not yet available.*

eventDataOutfile
  Event details will be written to this file. 

treeWriteFreq
  Frequency (in generations) of writing newick-formatted trees where branch lengths are scaled to rate.
  Use caution when setting a high write frequency, as the file size can grow quickly for large trees.

mcmcWriteFreq
  Frequency (in generations) at which MCMC details will be written to the mcmcOutfile.

eventDataWriteFreq
  Frequency (in generations) at which event details are written to the eventDataOutfile. 

acceptWriteFreq
  *Description not yet available.*

printFreq
  Frequency (in generations) of printing output to the screen.

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
  Percentage of total number of species sampled in your phylogeny (between 0 and 1).

sampleProbsFilename
  Specifies a file with clade-specific corrections for incomplete sampling.

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
  Starting initial speciation rate.

lambdaShift0
  Starting initial rate change parameter for speciation
  (if 0, speciation rates will not change through time).
  A negative value implies decreasing rates through time.

muInit0
  Starting Initial extinction rate.

muShift0
  Starting initial rate change parameter for extinction. Currently not implemented.

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
  Branch-specific speciation rates will be written to this file as newick-formatted trees.

muOutfile
  Branch-specific extinction rates will be written to this file as newick-formatted trees.

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
