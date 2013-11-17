

library(ape);

# This script assumes you have done some analysis of the whale tree
# and that the event data from that run is stored in 
#	file beventdata.txt


v <- read.tree('whaletree.tre');
ed <- getEventData(v, 'beventdata.txt', burnin=0.2, verbose=T, header=F);



###########################################
### The marginal shift probability tree:

mst <-marginalShiftProbsTree(ed);

quartz.options(height=10, width=8);
plot.phylo(ladderize(mst), show.tip.label=T, cex=0.4);
add.scale.bar(length = 0.5); #shows cumulative shift prob of 0.5 


###########################################
# The cumulative shift probability tree.
cum_shifttree <- cumulativeShiftProbsTree(ed);

# Plotting it.:
cum_shifttree <- ladderize(cum_shifttree);

quartz.options(height=10, width=8);
plot.phylo(cum_shifttree, cex=0.65);
add.scale.bar(length = 1); #shows cumulative shift prob of 1.0 

plot(cum_shifttree, type = 'f', show.tip.label=F);

################### 
# The maximum shift credibility tree.

bestsample <- maximumShiftCredibilityTree(ed);

quartz.options(height=10, width=8);
plot.phylo(v, cex=0.65);

# Get shift nodes (excluding the root)
shiftnodes <- bestsample$bestShiftConfig$node;
# exclude the root
shiftnodes <- shiftnodes[shiftnodes > (length(v$tip.label)+1)];

nodelabels(node = shiftnodes, pch=21, cex=2.2, bg ='red', col='black')

#################




