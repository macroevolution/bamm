.. highlight:: none

Configurations & Settings: BAMM for Trait Evolution
===================================================


Control File
------------

Configuration options and parameters in BAMM are specified in a *control file*,
a plain text file in which each line contains the name of the option or
parameter, an equal sign, and the value of the option or parameter::

    treefile = whaletree.tre
    runMCMC = 1

The path of the control file (relative to the current directory) is specified
with the flag ``-control`` when running bamm. For example::

    ./bamm -control traitcontrol.txt

This document refers to options and parameters for the phenotypic analysis in BAMM.

Global Options and Parameters
-----------------------------

The following describes the configuration options and parameters
that are required regardless of the specific model used.
For true or false values, 1 is used for true and 0 is used for false.
Files will be relative to the current directory unless otherwise
specified by the user.

General
.......

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
  
betaInitPrior
  *Description not yet available.*
  
betaShiftPrior
  *Description not yet available.*
  
MCMC Simulation Settings & Output Options
............................................

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

betaOutfile
  *Description not yet available.*
  
acceptrateOutFile
  *Description not yet available.*

eventDataOutfile
  Event details will be written to this file. Raw event data containing all of the results. See BAMMtools for working with this output file.

nodeStateOutfile
  *Description not yet available.*

MCMC Scaling Operators
......................

updateBetaScale
  *Description not yet available.*

updateBetaShiftScale
  *Description not yet available.*

updateNodeStateScale
  *Description not yet available.*

updateEventLocationScale
  Scale parameter for updating local moves of events on the tree
  This defines the width of the sliding window proposal.

updateEventRateScale
  Scale parameter (proportional shrinking/expanding) for updating
  the rate parameter of the Poisson process

MCMC Move Frequencies
......................

updateRateEventNumber
  Relative frequency of MCMC moves that change the number of events.
  
updateRateEventPosition
  Relative frequency of MCMC moves that change the location of an event on the tree.
  
updateRateEventRate
  Relative frequency of MCMC moves that change the rate at which events occur.

updateRateBeta0
  *Description not yet available.*
  
updateRateBetaShift
  *Description not yet available.*
  
updateRateNodeState
  *Description not yet available.*
  
localGlobalMoveRatio
  Ratio of local to global moves of events.

Initial Parameter Values
......................................

betaInit
  Initial Brownian motion rate parameter at the base of the tree.
  
betaShiftInit
  Initial rate change parameter for Brownian motion.
  
initialNumberEvents
  Initial number of non-root processes.
  
Numerical & Other Parameters
............................................

useObservedMinMaxAsTraitPriors
  *Description not yet available.*
  
traitPriorMin
  *Description not yet available.*
  
traitPriorMax
  *Description not yet available.*

