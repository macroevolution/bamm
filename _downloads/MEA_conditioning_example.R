
# This script requires files from the MEA Dryad submission 
# found in directory
# supplementary_data/code/monte_carlo_simulation

library(ape)
library(Rcpp)
library(BAMMtools)

############# 
# This 

whale_tree  <- read.tree("cetaceans.tre")

# You must set stringsAsFactors = F to replicate the code below, as 
#  factor encoding causes problems when we make a smaller dummy dataset from this tree

event_data  <- read.table("cetaceans_event_data.txt",header=TRUE,sep=",", stringsAsFactors = F)


mcmc_output <- read.table("cetaceans_mcmc_out.txt",header=TRUE,sep=",")
source("likelihoodFunctions.R") #from Moore et al (2016)
sourceCpp("MonteCarloSimulator.cpp")


# This function includes a modified version of the MEA likelihood calculator
# that simply returns the parameters of the process at the root along with 
# corresponding extinction probability computed by the MEA functions

source("MonteCarloLikelihood_MOD.R")   
 

 
generation       <- 300
these_parameters <- event_data[event_data$generation == generation,]
transition_rate  <- mcmc_output$eventRate[mcmc_output$generation == generation] / sum(whale_tree$edge.length)

# these_parameters is a BAMM event data object that has a root process (row 1)
# and 2 rate shifts.

# You can verify everything below using the analysis script exactly as presented in the 
# MEA html code using_the_likelihood_function.html
#
# However, it runs slowly and we can speed things up dramatically by making a much smaller
# input phylogeny (the goal here is not to actually correctly compute the likelihood of the cetacean
# phylogeny, merely to demonstrate the form of conditioning used by MEA)

# Here, we subset the whale tree down to just a 5 taxon tree that includes all branches with rate shifts:

tmpwhales <- whale_tree
taxset <- setdiff(tmpwhales$tip.label, as.vector(c(these_parameters$leftchild, these_parameters$rightchild)))
tmpwhales <- drop.tip(tmpwhales, taxset)


# We will also make these_parameters favor a really high extinction rate 
# at the root of the tree, to make this incomplete conditioning obvious,
# and we will set extinction on branches after the rate shifts to zero.
# This increases the effects of the theoretical error in MEA's calculations and,
# in the extreme, causes likelihoods to become infinite.

these_parameters$muinit <- c(0.25, 0, 0)
these_parameters$abstime <- c(0, 1, 1)

# look at our modifications to the data frame:
these_parameters

res  <- MonteCarloLikelihood_MOD(tmpwhales,
                                       these_parameters,
                                       transition_rate = transition_rate,
                                       lambda_prior = 1,
                                       mu_prior = 1)

# look at results:

res

# The root extinction probability is 0.98 or so, and the only parameters used for this are the 
# $node_process argument, which is just the paraemeter set at the root of the tree.

# We can now use MEA's MonteCarlo extinction simulator to verify that their calculations 
#    were performed with the parameters only at the root. 
#  

# We will just plug the root parameters in to their function and compute the relevant extinction probability.
# The true extinction probability should be low, because the data augmentation 
# involves rate shifts that occur very early in the clade history that lead to zero extinction.

root_age <- 35.85784
sim_E_root <- replicate(50000,SimulateCPBDP(time = root_age,init_lambda=these_parameters$lambdainit[1],init_mu=these_parameters$muinit[1],
                                                             transition_rate=transition_rate,lambda_prior=1,
                                                             mu_prior=1))

mean(sim_E_root < root_age)
# gives approximately 0.985

# this gives almost identical results to the MEA likelihood calculator:
res$root_e_prob 
# approx 0.985



