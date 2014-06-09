
data(whales, events.whales, prior.whales)
edata <- getEventData(whales, events.whales, burnin=0.1)

set.seed(5)
tipset <- sample(whales$tip.label, size=50)
edata <- subtreeBAMM(edata, tips = tipset)

set.seed(1);
sset <- sample(1:length(edata$eventData), size=20)

dset <- subsetEventData(edata, index=sset)


png(height=700, width=600, file = 'distinctshiftconfigs1.png')
par(mfrow=c(5,4))
for (i in 1:20){
	par(mar=c(1,1,1,1))
	tmp <- subsetEventData(edata, index=sset[i])
	plot(as.phylo.bammdata(tmp), no.margin=F, show.tip.label=F)
	addBAMMshifts(tmp, bg='red', cex=2, par.reset=F)
}
dev.off();



########### plot whale phylo with 
data(whales, events.whales, prior.whales)
edata <- getEventData(whales, events.whales, burnin=0.1)
mprobs <- marginalShiftProbsTree(edata)
nodes <- mprobs$edge[,2];
probs <- mprobs$edge.length

nodes <- nodes[mprobs$edge.length > 0]
probs <- probs[mprobs$edge.length > 0]

minx <- 1;
maxx <- 3;
minp <- min(probs)
maxp <- max(probs)

mm <- (maxx - minx) / (maxp - minp)

cexvec <- (probs - minp) * mm + minx;

png(height = 700, width=600, file = 'distinctshiftconfigs2.png')
plot.phylo(whales, show.tip.label=F, edge.width=1.3)

nodelabels(node=nodes, cex=1.3, pch=21, bg='blue') 

dev.off()


### 
data(prior.whales)
prior <- getBranchShiftPriors(whales, prior.whales)
bf <- bayesFactorBranches(edata, prior)

nodeset <- bf$edge[,2][bf$edge.length >= 3]

png(height = 700, width=600, file = 'distinctshiftconfigs3.png')
plot.phylo(whales, show.tip.label=F, edge.width=1.3)

nodelabels(node=nodeset, cex=1.6, pch=21, bg='red') 

dev.off()

############


data(prior.whales)
prior <- getBranchShiftPriors(whales, prior.whales)
css <- credibleShiftSet(edata, prior, BFcriterion=3) 
plot(css)
 
png(height = 600, width=600, file = 'distinctshiftconfigs4.png')
plot(css)
dev.off()














