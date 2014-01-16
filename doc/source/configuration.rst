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

``modeltype``
    The type of model BAMM uses for analysis.
    If ``speciationextinction``, run with a speciation/extinction model.
    If ``trait,`` run with a phenotypic evolution model.

``treefile``
    The file path of the input tree (in Newick format).
    For speciation/extinction analyses, the tree must be ultrametric.

``runInfoFilename``
    The path of the file to output general information about the current run.

``sampleFromPriorOnly``
    If ``1``, run by sampling from the prior only
    (ignoring likelihood contribution to posterior).
    If ``0``, run the full analysis.

``autotune``
    Experimental option for tuning MCMC operators.

``runMCMC``
    If ``1``, run the MCMC sampler.
    If ``0``, just check to see if the data can be loaded correctly.

``simulatePriorShifts``
    If ``1``, simulate prior distribution of the number of shift events,
    given the hyperprior on the Poisson rate parameter.
    This is necessary to compute the Bayes factor.

``loadEventData``
    If ``1``, load the event data from a file.

``eventDataInfile``
    The file path of the event data (used only if ``loadEventData = 1``).

``initializeModel``
    If ``1``, initializes MCMC.
    If ``0``, just check parameter file and ensure that data can be read.

``seed``
    Seed for the random number generator.
    If ``-1``, the seed is obtained from the clock time.

``overwrite``
    If ``1``, overwrite output files if they already exist with identical
    filenames.
    If ``0``, do not run if output files already exist with identical filenames.

MCMC Simulation
...............

``numberGenerations``
    Number of MCMC generations to run.

``mcmcWriteFreq``
    Frequency (in generations) at which to print MCMC details
    to ``mcmcOutfile``.

``eventDataWriteFreq``
    Frequency (in generations) at which to print event details
    to ``eventDataOutfile``.

``acceptWriteFreq``
    Frequency (in generations) at which to print accept rate to the screen.

``printFreq``
    Frequency (in generations) at which to print output to the screen.

``outName``
    If present (may be commented out), prefixes output files with the given
    text.

``mcmcOutfile``
    The path of the file to which to write the MCMC output.

``eventDataOutfile``
    The path of the file to which to write the raw event data.
    All of the results are contained in this file, and all branch-specific
    speciation rates, shift positions, marginal distributions, etc.,
    can be reconstructed from this output. See :ref:`bammtools`
    for more information on working with this output format.

``updateEventRateScale``
    Scale parameter (proportional shrinking/expanding) for updating
    the rate parameter of the Poisson process.

``localGlobalMoveRatio``
    Ratio of local to global moves of events.

Parameter Update Rates
......................

``updateRateEventNumber``
    Relative frequency of MCMC moves that change the number of events.

``updateRateEventPosition``
    Relative frequency of MCMC moves that change the location of an event
    on the tree.

``updateRateEventRate``
    Relative frequency of MCMC moves that change the rate at which events occur.

``initialNumberEvents``
    Initial number of non-root processes.


Speciation/Extinction Model
---------------------------

The following describes the configuration options and parameters
that are specific to speciation/extinction analyses in BAMM.

General
.......

``useGlobalSamplingProbability``
    If ``1``, look for a global correction for incomplete sampling
    (globalSamplingProbability).
    If ``0``, look for a file that specifies clade-specific corrections
    for incomplete sampling (``sampleProbsFilename``).

``globalSamplingProbability``
    Percentage of total number of species sampled in the phylogeny
    (between 0 and 1).

``sampleProbsFilename``
    The path of a file containing clade-specific corrections for
    incomplete sampling.

Priors
......

``lambdaInitPrior``
    Prior on the inital lambda (rate parameter of the exponential distribution)
    for the speciation rate.

``lambdaShiftPrior``
    Prior on the the lambda shift parameter (standard deviation of the normal
    distribution) for the speciation rate. The mean of the distribution
    is fixed at zero, which is equal to a constant rate diversification process.

``muInitPrior``
    Prior on the extinction rate (rate paramater of the exponential
    distribution).

``segLength``
    The "grain" of the likelihood calculations. It approximates the 
    continuous-time change in diversification rates by breaking each branch
    into a constant-rate diversification segments, with each segment equal
    to ``segLength``. So, a branch of length 10 will have the exponential
    speciation-rate change approximated by 10 segments if ``segLength = 1.0``.
    If the value is greater than the branch length (e.g., ``segLength = 100``
    in this case) BAMM will not break the branch into segments but use the mean 
    rate across the entire branch.

MCMC Simulation
...............

``updateLambdaInitScale``
    Scale parameter for updating the initial speciation rate for each process.

``updateLambdaShiftScale``
    Scale parameter for the exponential change parameter for speciation.

``updateMuInitScale``
    Scale parameter for updating initial extinction rate for each process.

``minCladeSizeForShift``
    Allows you to constrain the location of possible rate-change events
    to occur only on branches with at least this many descendant tips.
    A value of ``1`` allows shifts to occur on all branches.

Starting Parameters
...................

``lambdaInit0``
    Initial speciation rate (at the root of the tree).

``lambdaShift0``
    Initial rate change parameter for speciation at the root.
    If ``0``, speciation rates will not change through time.
    A negative value implies decreasing rates through time.

``muInit0``
    Initial extinction rate at the root.

Parameter Update Rates
......................

``updateRateLambda0``
    Relative frequency of MCMC moves that change the initial speciation rate
    associated with an event.

``updateRateLambdaShift``
    Relative frequency of MCMC moves that change the exponential shift parameter
    of a speciation rate associated with an event.

``updateRateMu0``
    Relative frequency of MCMC moves that change the extinction rate for a given
    event.


Phenotypic Evolution Model
--------------------------

The following describes the configuration options and parameters
specific to the phenotypic evolution model in BAMM.
The parameter "beta" represents the rate of phenotypic evolution
at any point in time.

General
.......

``traitfile``
    The path to a file that contains the phenotypic trait data.
    Traits must be continuous characters.
    Each line must have a species name and the corresponding trait value,
    separated by a tab.
    A header row is **not** permitted.
    All species in the trait data file must be in the tree and vice versa.

MCMC Tuning
...........

``updateBetaScale``
    Scale operator for proportional shrinking/expanding move to update
    the initial phenotypic rate for rate regimes.

``updateNodeStateScale``
    Scale operator for sliding window move to update ancestral states
    at internal nodes.

``updateBetaShiftScale``
    Scale operator for sliding window move to update initial phenotypic rate.

Starting Parameters
...................

``betaInit``
    Initial value of the phenotypic evolutionary process at the root
    of the tree.

``betaShiftInit``
    Initial value of the exponential change parameter for the phenotypic
    evolutionary process (at the root of the tree).
    If ``0``, then the process has a constant rate.
    If negative, it implies decreasing rates through time.

Priors
......

``betaInitPrior``
    Parameter (rate) of the prior (exponential) on the inital phenotypic
    evolutionary rate associated with regimes.

``betaShiftPrior``
    Parameter (stdandard deviation) of the prior (normal) on the rate-change
    parameter.

``useObservedMinMaxAsTraitPriors``
    If ``1``, puts a uniform prior density on the distribution of ancestral
    character states, with upper and lower bonds determined by the min and max
    of the observed data.

``traitPriorMin``
    User-defined minimum value for the uniform density on the distribution
    of ancestral charater states. Only used if
    ``useObservedMinMaxAsTraitPriors = 0``.

``traitPriorMax``
    User-defined maximum value for the uniform density on the distribution
    of ancestral charater states. Only used if
    ``useObservedMinMaxAsTraitPriors = 0``.
    
Parameter Update Rates
......................

``updateRateBeta0``
    Relative freuency of moves that change the initial phenotypic rate
    associated with an event.

``updateRateBetaShift``
    Relative frequency of moves that change the exponential shift parameter
    of a phenotypic rate associated with an event.

``updateRateNodeState``
    Relative frequency of moves update the value of ancestral character stats.
    You have as many ancestral states as you have internal nodes in your tree,
    so there are a lot of parameters: this value should, in general,
    be substantially higher than the other parameter values
    (recommended 25:1 or 50:1) because there are so many internal nodes states
    that need to be updated.
