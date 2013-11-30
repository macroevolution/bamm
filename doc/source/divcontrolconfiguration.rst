.. highlight:: none

Parameters & Settings: speciation-extinction BAMM
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

    ./bamm -control divcontrol.txt

This document refers to options and parameters for the speciationextinction analysis in BAMM.

General Setup & Data Input
-----------------------------

The following describes the configuration options and parameters
that are required regardless of the specific model used.
For true or false values, 1 is used for true and 0 is used for false.
Files will be relative to the current directory unless otherwise
specified by the user.

General
.......

modeltype
   Should be "speciationextinction" for diversity analyses.
   Defines which application of BAMM will run.

treefile
  The file name of the input tree (in Newick format). For diversification analyses, this should be ultrametric.

sampleFromPriorOnly
  If true (1), run BAMM by sampling from the prior only
  (ignoring likelihood contribution to posterior).
  If false (0), run the full analysis.
  
autotune
  Experimental option for tuning MCMC operators.
  
runMCMC
  If true (1), run the MCMC sampler.
  If false (0), just check to see if the data can be loaded correctly.

loadEventData
  If true (1), load event configuration (see below, eventDataInfile).

eventDataInfile
  *Description not yet available.*

initializeModel
  If true (1), initializes MCMC. If false (0), will just check parameter file and ensure that data can be read.

useGlobalSamplingProbability
  If true (1), will look for a global correction for incomplete sampling (globalSamplingProbability)
  If false (0), will look for a file that specifies clade-specific corrections for incomplete sampling (sampleProbsFilename).
  
globalSamplingProbability
  Percentage of total number of species sampled in your phylogeny (between 0 and 1). If useGlobalSamplingFraction = 0, this will be ignored.

sampleProbsFilename
  Specifies a file with clade-specific corrections for incomplete sampling.

seed
  Seed for random number generator. If unspecified or -1, seed is obtained from clock time.
  
overwrite
  If true (1), will overwrite analysis result files if output files already exist with identical filenames. If false (0), analyses will not run if output files already exist with identical filenames.

Priors
......

rootPrior
  *Description not yet available.*
  
poissonRatePrior
  *Description not yet available.*
  
lambdaInitPrior
  Mean of the exponential distribution prior on speciation.

lambdaShiftPrior
  Prior on the speciation rate change parameter.

muInitPrior
  Exponential prior on extinction.

muShiftPrior
  *Description not yet available.*

MCMC Simulation Settings & Output Options
...........

numberGenerations
  Number of MCMC generations to run.

branchRatesWriteFreq
  Frequency (in generations) at which branch-specific mean rates are written to file. Can be computed directly from eventdata using BAMMtools functions.

mcmcWriteFreq
  Frequency (in generations) at which MCMC details will be written to the mcmcOutfile.

eventDataWriteFreq
  Frequency (in generations) at which event details are written to the eventDataOutfile. 

acceptWriteFreq
  Prints accept frequency.

printFreq
  Frequency (in generations) of printing output to the screen.
  
outName
  Prefix for output files.

mcmcOutfile
  MCMC parameter output will be written to this file.

eventDataOutfile
  Event details will be written to this file. Raw event data containing all of the results. See BAMMtools for working with this output file.

lambdaOutfile
  Branch-specific speciation rates will be written to this file as newick-formatted trees.

lambdaNodeOutfile
  *Description not yet available.*

muOutfile
  Branch-specific extinction rates will be written to this file as newick-formatted trees.

acceptrateOutfile
  *Description not yet available.*

MCMC Scaling Operators
......................

updateLambdaInitScale
  Scale parameter for updating the initial speciation rate for each process.
  
updateLambdaShiftScale
  Scale parameter for the exponential change parameter for speciation.
  
updateMuInitScale
  Scale parameter for updating the initial extinction rate for each process.
  
updateMuShiftScale
  Scale parameter for the exponential change parameter for extinction.
  
updateRateMuShift
  *Description not yet available.*

updateEventLocationScale
  Scale parameter for updating local moves of events on tree. This defines the width of the sliding window proposal.

updateEventRateScale
  Scale parameter for updating the rate parameter of the Poisson process.
  
MCMC Move Frequencies
......................

updateRateEventNumber
  Relative frequency of MCMC moves that change the number of events.
  
updateRateEventPosition
  Relative frequency of MCMC moves that change the location of an event on the tree.
  
updateRateEventRate
  Relative frequency of MCMC moves that change the rate at which events occur.
  
updateRateLambda0
  Relative frequency of MCMC moves that change the initial speciation rate associated with an event.
  
updateRateLambdaShift
  Relative frequency of MCMC moves that change the exponential shift parameter of a speciation rate associated with an event.
  
updateRateMu0
  Relative frequency of MCMC moves that change the extinction rate for a given event.
  
localGlobalMoveRatio
  Ratio of local to global moves of events.

Initial Parameter Values
...................

lambdaInit0
  Initial speciation rate at the root of the tree.

lambdaShift0
  Initial rate change parameter for speciation
  (if 0, speciation rates will not change through time).
  A negative value implies decreasing rates through time.

muInit0
  Initial extinction rate at the root of the tree.

muShift0
  Initial rate change parameter for extinction. Currently not implemented (vaue is 0).
  
initialNumberEvents
  Initial number of non-root processes.

Numerical & Other Parameters
......................

minCladeSizeForShift
  Set the minimum number of descendant tips a branch must have to be the location of a possible rate-change event. A value of 1 allows shifts to occur on any branch.

segLength
  *More information to come...*
