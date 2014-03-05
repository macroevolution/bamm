
library(BAMMtools)

data(whales, events.whales)

# make the bammdata object:
edata <- getEventData(whales, events.whales, burnin=0.1)

shiftprobs <- summary(edata)

# View the dataframe with shift posterior probabilities:
shiftprobs

#########################################
#########################################
########   Mean phylorate plot

plot.bammdata(edata, lwd=2)
# see helpfile for other color palette examples

plot.bammdata(edata, lwd=2, legend=T)


## Plotting a random sample from posterior:
index <- 25
e2 <- subsetEventData(edata, index = index)
plot.bammdata(e2, lwd=2)
addBAMMshifts(e2, cex=2)

#########################################
#########################################
########   Bayesian credible set of shifts

data(prior.whales)
priorshifts <- getBranchShiftPriors(whales, prior.whales)
dss <- distinctShiftConfigurations(edata, priorshifts)

css <- credibleShiftSet(edata, threshold = priorshifts)

# Number of shift configurations in the credible set:
css$number.distinct

summary(css)

# Plot the 9 most-frequently sampled shift configurations
# 	as phylorate plots
plot.credibleshiftset(css, plotmax=9, legend=F)


#########################################
#########################################
########   Best shift configuration

best <- getBestShiftConfiguration(edata, threshold = priorshifts);
plot.bammdata(best, lwd=2)
addBAMMshifts(best, cex=2.5)


first <- subsetEventData(css, index=1)
second <- subsetEventData(css, index = 2)

plot.bammdata(second);
addBAMMshifts(second, cex=2)


#########################################
#########################################
########   maximum shift credibility configuration

msc.set <- maximumShiftCredibility(edata, maximize='product');
msc.config <- subsetEventData(edata, index = msc.set$sampleindex);
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)

#########################################
#########################################
########   some random shift configurations

dsc <- distinctShiftConfigurations(edata, threshold=0.01)
# Here is one random sample with the BEST shift configuration
plot.bammshifts(dsc, edata, rank=1, legend=F)
# Here is another (read the label text):
plot.bammshifts(dsc, edata, rank=1, legend=F)
 
plot.bammshifts(dsc, edata, rank=2, legend=F)
plot.bammshifts(dsc, edata, rank=2, legend=F)
  


#########################################
#########################################
########   shift plotting using plot.phylo from ape:

mysample <- 25  # this is the sample we'll plot
	
nrow(edata$eventData[[ mysample ]]) 
 
shiftnodes <- getShiftNodesFromIndex(edata, index = mysample)	

plot.phylo(whales)
nodelabels(node = shiftnodes, pch=21, bg="red", cex=1.5)


#########################################
########   marginal shift probabilities:

marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, show.tip.label=F)

# cumulative shift probs tree
cst <- cumulativeShiftProbsTree(edata)
plot.phylo(cst)


cst <- cumulativeShiftProbsTree(edata)
edgecols <- rep('black', length(whales$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"
plot.phylo(whales, edge.color = edgecols)


#########################################
########   clade-specific rates:


library(BAMMtools)
data(whales)
data(events.whales)
edata <- getEventData(whales, events.whales, burnin=0.1)
#and here we get the rates
allrates <- getCladeRates(edata)

# Mean and 95% HPD on whale speciation rates:
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))
	
# Rates for dolphins only:
dolphinrates <- getCladeRates(edata, node=140)	
mean(dolphinrates$lambda)	
	
# Rates for all NON-DOLPHINS
nondolphinrate <- getCladeRates(edata, node = 140, nodetype = "exclude")
mean(nondolphinrate$lambda)
quantile(nondolphinrate$lambda, c(0.05, 0.95))	

#########################################
########  branch and tip specific rates

edata$meanTip


#########################################
########  rate through time analysis

# Finding node labels by inspecting the tree:
plot.phylo(whales)
nodelabels()
	
# Extracting by finding MRCA
species1 <- "Tursiops_truncatus"
species2 <- "Orcinus_orca"
	
#Now to get the *tip node numbers* in ape format::	
	
tipnode1 <- which(whales$tip.label == species1)
tipnode2 <- which(whales$tip.label == species2)
	
#And now the MRCA node::

mrca <- getMRCA(whales, tip = c(tipnode1, tipnode2))
	
# Now we feed this in to plotRateThroughTime::
plotRateThroughTime(edata, node = mrca, nodetype="include", ylim=c(0,0.7))
	
# Plotting for background after excluding dolphins:
plotRateThroughTime(edata, node = mrca, nodetype = "exclude", ylim=c(0, 0.7))

rtt_subtree <- getRateThroughTimeMatrix(edata, node = mrca)

###########





