# This script illustrates ascertainment bias in MEA's simulations
#  Script loops over a set of random starting seeds,
#  generates trees under MEAs simulator
#  
#  Note that MEA rejected all trees with fewer than 50 or more than 150 tips
#   As you can see from executing this script, trees larger than 150 tips 
#   are very common given their parameterization
#
#  We have made a simple modification to MEA's tree simulator such that it can easily be 
#    used to generate trees outside of their target size range,
#    such that we can see what types of trees they rejected 


library(geiger)
library(ape)

# load tree simulation functions
source("R_test_bamm_functions.R")

# -------------------------------------
# Parameter initialization block
#   these are MEA's parameters from their "variable rates" dataset

lamfx <- function() return(rexp(1, 1/0.15))
mufx <- function() return(rexp(1, 1/0.05))

# rate at which events occur along the phylogeny
trate <- 0.006


# -------------------------------------
# 

# Make vector of integer-valued seeds so this can be repeated:
# Here we will just look at the first 50 integers from 1 to 50
#  but you can modify this as desired

seedvec <- 1:20

treeset <- list()                  # this stores the raw trees
treeset_pruned <- list()           # This stores the trees without extinct tips
ntips <- rep(NA, length(seedvec))  # vector of tip counts for trees


for (i in 1:length(seedvec)){
	cat(i, "\n")
	
	set.seed(seedvec[i])
	
	tree <- SimulateCBDPTree(35.5, trate, lamfx, mufx, verbose=F, NMAX = 2000, MAX_FAILS = 1)
	
	treeset[[i]] <- tree
	
	if (length(tree) > 1){
		
		# Map shifts onto tree
		simTree <- CPBDPStochasticMap(tree)
		
		# Prune out extinct lineages 
		treeset_pruned[[i]] <- pruneCPBDPTree(simTree)
		ntips[ i ] <- length(treeset_pruned[[i]]$tip.label)
	}else{
		treeset_pruned[[i]] <- NA
	}
	
	
}

# You now have a list of trees that can be explored.
#
# For now, we will simply look at the most basic property of these trees:
#   their size.
# 
# Create a dataframe that holds taxon counts
#

dff <- data.frame(seed = seedvec, ntips = ntips)
write.table(dff, file = "MEA_treedata_unbiased.csv", sep="", quote=F, row.names=F)
 
# IMPORTANT: this output file stores as "NA" every tree simulation that failed.
#             tree simulations fail whenever the total number of tips exceeds the 
#              maximum that we specified. We set this max to 2000 for these simulations.
#              MEA set this to 150. You can see that many, many trees are generated that 
#              exceed 150 tips. 
#          










	

