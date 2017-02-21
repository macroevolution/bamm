
# This script illustrates how you can simulate a tree using MEA's tree 
# simulator from a given starting seed.
# You should be able to repeat this analysis exactly, e.g., generate
#  the exact same tree that we do, given a particular random seed.
# The script also illustrates how to output the tree
#  and set up a BAMM run so you can test BAMM's performance 
#
 
################################################################################
library(ape)
 
library(geiger)
library(BAMMtools)
 
# load the tree simulation functions from MEA
source("R_test_BAMM_functions.R")

# -------------------------------------
# Parameter initialization block

# Here we will set up the input parameterizations to exactly
#   equal those used by MEA in their "variable rates" tree dataset

# define speciation & extinction functions for the tree simulation
#  these are exponential distributions with means of 0.15 and 0.05 
lamfx <- function() return(rexp(1, 1/0.15))
mufx <- function() return(rexp(1, 1/0.05))

# rate at which events occur along the phylogeny
trate <- 0.006

# -------------------------------------
# Tree simulation block

# choose an integer-valued seed 
seed <- 8
seed_stub <- paste("s", seed, sep="")
treename <- paste("tree_", seed_stub, ".tre", sep="")

# this is the important bit, which sets the seed!
set.seed(seed)

# Using MEA tree simulator to generate the tree:
tree <- SimulateCBDPTree(35.5, trate, lamfx, mufx, verbose=F, NMAX = 5e3, MAX_FAILS = 1)

# Map shifts onto tree
simTree <- CPBDPStochasticMap(tree)
 
# Prune extinct taxa 
prunedTree <- pruneCPBDPTree(simTree)

# Convert the regime mappings from the simulation to BAMM-style eventData format
edata <- mea_to_edata(prunedTree)
 
# Write the event data to file:
write.table(edata[,1:8], file = paste(seed_stub, "_true_eventdata.txt", sep=""), sep=",", quote=F, row.names=F)

# ------------------------------------
# Generate control file for BAMM analysis

# setBAMMpriors on the tree
priors <- setBAMMpriors(prunedTree, outfile = NULL)
 
write.tree(prunedTree, treename)

bammcontrolfile <- paste("control_", seed_stub, ".txt", sep="")

# This will generate a BAMM control file that can be run for
#   5 million generations
#   We will assume time-constant BAMM:

generateControlFile(file = bammcontrolfile, 
					type = "diversification", 
					params = 
					list(
					  treefile = treename,
					  numberOfGenerations = '5000000',
					  overwrite = '1',
					  expectedNumberOfShifts = '100',
					  lambdaInitPrior = as.numeric(priors['lambdaInitPrior']),
					  lambdaShiftPrior = '0',
					  muInitPrior = as.numeric(priors['muInitPrior']),
					  lambdaIsTimeVariablePrior = '0',
					  updateRateLambdaShift = '0',
					  updateRateEventPosition = '0.5',
					  numberOfChains = '1', 
					  outName = seed_stub
					 )
					)




# -------------------------------------------
# You are now ready to run BAMM on this tree
# 
# On OSX, simply move your control file that you generated 
#  above to the directory where BAMM is located 
#  and you can invoke the program with:
#
#  ./bamm -c controlfile.txt
#
# where bamm is the name of your bamm executable, and 
#  controfile.txt is the name of your controlfile 
#  (note that your controlfile will actually be control_s8.txt or similar)
# 














