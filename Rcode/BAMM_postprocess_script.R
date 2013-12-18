
library(ape);

#########################################################
#	Testing convergence of the MCMC
#		simplest way is to look at effective sample sizes
#		 of log-likelihoods and numbers of processes

library(coda);
mcmc_res <- read.csv('bmcmcout.txt', header=F);

# BAMM generates this file by default -
# 1st column: generation
# 2nd column: log-likelihood
# 3rd column: log-prior
# 4th column: Number of non-root processes on tree

plot(mcmc_res[,4] ~ mcmc_res[,1]);

burn_fraction <- 0.1; # fraction of sample to discard as burnin
burn_start <- floor(burn_fraction * nrow(mcmc_res));
keepset <- mcmc_res[burn_start:nrow(mcmc_res), ];

# compute effective sample sizes:
effectiveSize(mcmc_res[,2]); # number of processes
effectiveSize(mcmc_res[,4]); # log-likelihood


#########################################################
#	Summarize the BAMM analyses
#	First, we will generate an object of class 'bammdata'
#	using the event data output file:

source('/Users/danrabosky/DanWork/bamm/devel/testrun/process_BAMM_fxns.R');

# Here we read in the model tree - the ultrametric tree that we analyzed 
#	with BAMM:

modeltree <- read.tree('oz_spheno_mcc216.tre');

bdata <- getEventData(modeltree, 'beventdata.txt', burnin=0.25, nsamples=30, verbose=T);

# Once this is processed, we can do all sorts of things with this data object.
# For example, suppose we want to get the AVERAGE rate of speciation 
#	across the entire tree

mean_rates <- getCladeRates(bdata);

# Now, what if we want it for just a particular clade? 
# Let us try this for lizards in the genus Ctenotus. First,
#	we find which node corresponds to the Ctenotus MRCA:
ctenotus_species <- modeltree$tip.label[grep('ctenotus', modeltree$tip.label)];


span_pair <- getSpanningTaxonPair(modeltree, ctenotus_species);

#Get the MRCA node:

ctenotus_node <- getMRCA(modeltree, which(modeltree$tip.label %in% span_pair));

# Now we can get the mean rate for just the Ctenotus subtree:
source('/Users/danrabosky/DanWork/bamm/devel/testrun/process_BAMM_fxns.R');

ctenotus_rate <- getCladeRates(bdata, node=ctenotus_node, nodetype = 'include');

# Now we'll get the background rate excluding Ctenotus and its 
#	sister clade, Lerista:

parentnode <- modeltree$edge[,1][modeltree$edge[,2] == ctenotus_node];

background_rate <- getCladeRates(bdata, node=parentnode, nodetype = 'exclude');

#######################################

# Now we will generate a rate-through-time curve

rtt.matrix <- getRateThroughTimeMatrix(bdata);

# look at what is in the returned object:
attributes(rtt.matrix);

# Now we sweep out the mean rates and make a simple plot:
meanLambda <- colMeans(rtt.matrix$lambda);
meanMu <- colMeans(rtt.matrix$mu);

# now we will make a nice plots of speciation and extinction through time:

quartz.options(height=6, width=6, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(6,6,1,1));
plot.window(xlim=c(0, max(rtt.matrix$time)), ylim=c(0 , 0.8));
lines(x=rtt.matrix$time, y=meanLambda, lwd=3, col='red');
lines(x=rtt.matrix$time, y=meanMu, lwd=3, col='blue');
axis(at=seq(-5, 25, by=5), cex.axis=1.2, side=1);
axis(at=seq(-0.2, 1.2, by=0.2), las=1, cex.axis=1.2, side=2);
mtext(side=1, text='Time from start of radiation', line =3, cex=1.4);
mtext(side=2, text='Rate', line =3, cex=1.4);

#######################

# We could also sweep out the 95% confidence intervals on speciation 
#	and extinction from the rate-through-time matrix:

mm <- apply(rtt.matrix$lambda, MARGIN = 2, quantile, c(0.025, 0.5, 0.975));

# Here we will plot with a confidence polygon for speciation:

quartz.options(height=6, width=6, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(6,6,1,1));
plot.window(xlim=c(0, max(rtt.matrix$time)), ylim=c(0 , 0.8));

xvals <- c(rtt.matrix$time, rev(rtt.matrix$time));
yvals <- c(mm[1,], rev(mm[3,]));
polygon(x=xvals, y=yvals, col='gray70', border=F);
lines(x=rtt.matrix$time, y=meanLambda, lwd=3, col='red');
axis(at=seq(-5, 25, by=5), cex.axis=1.2, side=1);
axis(at=seq(-0.2, 1.2, by=0.2), las=1, cex.axis=1.2, side=2);
mtext(side=1, text='Time from start of radiation', line =3, cex=1.4);
mtext(side=2, text='Speciation', line =3, cex=1.4);



######################################
# We can also pull out the rates for any clade like this (here we 
# are pulling out rates for ctenotus and lerista)

rtt_ctenotus_lerista <- getRateThroughTimeMatrix(bdata, node=parentnode, nodetype='include');

# and here is the "background" rate: non-ctenotus, non-lerista
rtt_background <- getRateThroughTimeMatrix(bdata, node=parentnode, nodetype = 'exclude');

# we will make these rates into "net diversification rates" 
#	by subtracting extinction (mu) from speciation (lambda)

ndr_CL <- rtt_ctenotus_lerista$lambda - rtt_ctenotus_lerista$mu;
ndr_back <- rtt_background$lambda - rtt_background$mu;

# These are matrices, so we will sweep out the means:
mean_CL <- colMeans(ndr_CL);
mean_back <- colMeans(ndr_back);

# now we plot all together:

quartz.options(height=6, width=6, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(6,6,1,1));
plot.window(xlim=c(0, max(rtt_background$time)), ylim=c(0 , 0.8));

lines(x=rtt_background$time, y=mean_back, lwd=3, col='gray50');
lines(x=rtt_ctenotus_lerista$time, y=mean_CL, lwd=3, col='red');
axis(at=seq(-5, 25, by=5), cex.axis=1.2, side=1);
axis(at=seq(-0.2, 1.2, by=0.2), las=1, cex.axis=1.2, side=2);
mtext(side=1, text='Time from start of radiation', line =3, cex=1.4);
mtext(side=2, text='Rate', line =3, cex=1.4);

#####################################################
# Now we move on to another issue: plotting trees and 
#	extracting means of marginal rate distributions for branches

# Two ways of doing this with BAMM output.
# The slower way is to pull these out of the 'bammdata' object

ratemat <- getMarginalBranchRateMatrixDiversification(bdata);

# this gives speciation and extinction matrices: each
# column corresponds to branch rates from a single sample
# from the posterior

# we can sweep out the means like this:

mean_branch_rates <- rowMeans(ratemat$lambda_branch_matrix);

# and this matches the 'edge.length' component
#	of the phylogenetic tree. This lets us do a couple of things. We can
#	make a new tree, where branches are scaled proportional to means of 
#	marginal distributions on individual branches:

# Reading in the model tree:
modeltree <- read.tree('oz_spheno_mcc216.tre');

lambda_tree <- modeltree;
lambda_tree$edge.length <- mean_branch_rates;
lambda_tree <- ladderize(lambda_tree); #make pretty

# Plotting:
quartz.options(height=10, width=10, dpi=72);
plot.phylo(lambda_tree, cex=0.3);

# What if we want to plot the original tree, but color the branches 
# based on the reconstructed rates?

# Lots of possible colors schemes - my favorite is rich.colors
#	from gplots

mbr <- rowMeans(ratemat$lambda_branch_matrix);

library(gplots);
NCOLORS <- 100;
mypalette <- rich.colors(NCOLORS);
colset <- numeric(nrow(modeltree$edge));

# get color breakpoints:
bks <- quantile(mean_branch_rates, seq(0, 1, length.out=(NCOLORS+1)));

for (i in 2:length(bks)){
	if (i == 2){
		colset[mbr < bks[2]] <- mypalette[1];
	}else if (i == length(bks)){
		colset[mbr >= bks[length(bks)-1]] <- mypalette[length(bks) - 1];
 	}else{
 		colset[mbr >= bks[i-1] & mbr < bks[i]] <- mypalette[i-1];
 	}
}

quartz.options(height=12, width=8, dpi=72);
plot.phylo(modeltree, edge.color=colset, edge.width=1, cex=0.3);

## 


































