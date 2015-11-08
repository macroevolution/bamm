setwd('~/bamm/doc/source/figures')

library(BAMMtools)
library(TeachingDemos)


##################################
## RATESHIFTS PAGE
##################################

# -------------------------------
# x_whale_marginals.png

library(BAMMtools)
data(whales, mcmc.whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.1)

mst <- marginalShiftProbsTree(ed)

n1 <- getMRCA(whales, c("Delphinus_delphis", "Orcaella_brevirostris"))
n2 <- getMRCA(whales, c("Delphinus_delphis", "Orcinus_orca"))
n3 <- getMRCA(whales, c("Delphinus_delphis", "Phocoena_sinus"))

edges <- round(mst$edge.length[mst$edge[,2] %in% c(n1,n2,n3)], 2)

#par(mar=c(1,1,1,7))
par(oma=c(1,1,1,1))
par(mar=c(1,1,1,10))

png(height=400, width=700, file = "x_whale_marginals.png")

#quartz.options(height = 6, width = 10)
plot.new()

plot.bammdata(ed, lwd=2, breaksmethod = "quantile", tau = 0.002)

lines(x=c(25.8, 23.4), y=c(61.2, 77.9), lwd=1.5, col = "gray40")

text(x = 21, y = 82, label = "0.14", cex = 1.8, pos=4, font = 2)
text(x = 12.8, y = 54.82, label = "0.29", cex = 1.8, pos=4, font=2)
text(x = 19.5, y = 58.9, label = "0.46", cex = 1.8, pos=4, font=2)

text(x = 37.1, y = 85, label = "Dolphins", cex = 1.8, pos=2, font = 2,  srt=90)
 
dev.off()


# -------------------------------
# xx_whale_postshiftsamples.png

png(height=700, width=700, file = "xx_whale_postshiftsamples.png")
par(mfrow=c(2,2))
par(oma=c(1,1,1,1))
par(mar=c(1,1,1,1))
lx <- 1.7
ii <- 1
etmp <- subsetEventData(ed, ii)
z1 <- plot.bammdata(etmp, breaksmethod="quantile", tau = 0.004, lwd=lx)
addBAMMshifts(etmp, cex = 2, par.reset=F)
text(x=10, y = 70, label = "gen 1", cex = 1.5, font = 2)
ii <- 11
etmp <- subsetEventData(ed, ii)
plot.bammdata(etmp, breaksmethod="quantile", tau = 0.004, lwd=lx)
addBAMMshifts(etmp, cex = 2, par.reset=F)
text(x=10, y = 70, label = "gen 11", cex = 1.5, font = 2)

ii <- 21
etmp <- subsetEventData(ed, ii)
plot.bammdata(etmp, breaksmethod="quantile", tau = 0.004, lwd=lx)
addBAMMshifts(etmp, cex = 2, par.reset=F)
text(x=10, y = 70, label = "gen 21", cex = 1.5, font = 2)

ii <- 31
etmp <- subsetEventData(ed, ii)
plot.bammdata(etmp, breaksmethod="quantile", tau = 0.004, lwd=lx)
addBAMMshifts(etmp, cex = 2, par.reset=F)
text(x=10, y = 70, label = "gen 31", cex = 1.5, font = 2)


dev.off()


# ---------------------------------------
# node 140:
# node 141

# library(BAMMtools)
# data(whales, mcmc.whales, events.whales)
# ed <- getEventData(whales, events.whales, burnin=0.1)

# nodemat <- matrix(0, nrow=length(ed$eventData), ncol=max(ed$edge))

 
# for (i in 1:length(ed$eventData)){
	# tmp <- ed$eventData[[i]]
	# nodemat[i,tmp$node] <- 1
# }


# plot(nodemat[,146])

# ss <- spectrum0.ar(nodemat[,140])


# plot.new()
# par(mfrow=c(3,1))
# plot(nodemat[,132], type = "l")
# plot(nodemat[,140], type = "l")
# plot(nodemat[,141], type = "l")




# ---------------------------------------
# x_whale_priors.png

data(whales, events.whales)
edata <- getEventData(whales, events.whales, burnin=0.1)
margprobs <- marginalShiftProbsTree(edata)
branch_priors <- getBranchShiftPriors(whales, expectedNumberOfShifts = 1)


par(oma=c(1,1,1,1))
par(mar=c(1,1,1,10))

png(height=400, width=700, file = "x_whale_priors.png")

#quartz.options(height = 6, width = 10)
plot.new()


plot.bammdata(ed, lwd=2, breaksmethod = "quantile", tau = 1, mask = ed$edge[,2], mask.color="black")

lines(x=c(25.8, 23.4), y=c(61.2, 77.9), lwd=1.5, col = "gray40")


text(x = 21, y = 82, label = "0.002", cex = 1.8, pos=4, font = 2)
text(x = 12.8, y = 54.82, label = "0.01", cex = 1.8, pos=4, font=2)
text(x = 19.5, y = 58.9, label = "0.009", cex = 1.8, pos=4, font=2)


#text(x = 37.1, y = 85, label = "Dolphins", cex = 1.8, pos=2, font = 2,  srt=90)
 
dev.off()

# ----------------------------------------
# x_whale_marginalodds.png

data(whales, events.whales)
edata <- getEventData(whales, events.whales, burnin=0.1)
branch_priors <- getBranchShiftPriors(whales, expectedNumberOfShifts = 1)
mo <- marginalOddsRatioBranches(edata, expectedNumberOfShifts = 1)

par(oma=c(1,1,1,1))
par(mar=c(1,1,1,10))

png(height=400, width=700, file = "x_whale_marginalodds.png")

#quartz.options(height = 6, width = 10)
plot.new()


plot.bammdata(ed, lwd=2, breaksmethod = "quantile", tau = 1, mask = ed$edge[,2], mask.color="black")

lines(x=c(25.8, 23.4), y=c(61.2, 77.9), lwd=1.5, col = "gray40")


text(x = 21, y = 82, label = "91.9", cex = 1.8, pos=4, font = 2)
text(x = 12.8, y = 54.82, label = "53.1", cex = 1.8, pos=4, font=2)
text(x = 19.5, y = 58.9, label = "29.7", cex = 1.8, pos=4, font=2)


#text(x = 37.1, y = 85, label = "Dolphins", cex = 1.8, pos=2, font = 2,  srt=90)
 
dev.off()


# ------------------------------------
# x_whale_credibleshiftset.png

png(height=700, width=700, file = "x_whale_credibleshiftset.png")

data(whales, events.whales)
edata <- getEventData(whales, events.whales, burnin=0.1)
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=3) 
plot.credibleshiftset(css, tau = 0.002, breaksmethod="quantile", lwd=1.5, cex=2)
dev.off()





# -------------------------------
## xprimates_shiftconfigs.png

library(BAMMtools)
library(TeachingDemos)
data(primates, events.primates)
edata_primates <- getEventData(primates, events.primates, burnin=0.1, nsamples=1000, type='trait')

png(height=700, width=700, file = "xprimates_shiftconfigs.png")
plot.new()
par(mfrow=c(2,2), mar=c(0,0.5,0,0.5), bg='white', family='Courier');
 
plot(primates, show.tip.label = FALSE)
addBAMMshifts(edata_primates, index=1, col='white', bg='red', cex=1.5, par.reset=FALSE);
text('Sample 1', x=15, y=210, cex=2)

plot(primates, show.tip.label = FALSE)
addBAMMshifts(edata_primates, index=50, col='white', bg='red', cex=1.5, par.reset=FALSE);
text('Sample 50', x=15, y=210, cex=2)

plot(primates, show.tip.label = FALSE)
addBAMMshifts(edata_primates, index=100, col='white', bg='red', cex=1.5, par.reset=FALSE);
text('Sample 100', x=15, y=210, cex=2)

q <- plot.bammdata(edata_primates, tau=0.001, lwd=1.2, breaksmethod='jenks')

rect(-5, 130, 35, 235, col='white', border='black', xpd=NA, lwd=0.5)
subplot(ratesHistogram(q, plotBrks=FALSE, xlab='Phenotypic rate', ylab='Density', cex=0.8), x=17, y=195, size=c(2,1.2))

dev.off()


# ---------------------------------------
# rateshifts_prior.png

library(BAMMtools)

prob.k <- function(k, poissonRatePrior=1)	{
	Denom <- (poissonRatePrior + 1)^(k+1)
	Prob <- poissonRatePrior / Denom
	return(Prob)
}

png(file = 'rateshifts_prior.png', width=450, height=450);

obsK <- seq(from=0, to=19, by=1)
expectedNumberofShifts <- 1

priorD <- sapply(obsK, prob.k, poissonRatePrior=1/expectedNumberofShifts)

par(mar=c(5,5,1,1))
plot(1, 1, xlim=c(0,max(obsK)), ylim=c(0,0.5), type='n', xlab="Number of shifts", ylab="Probability")
points(obsK, priorD, pch=21, bg='red', type='b')

dev.off();


# ---------------------------------------
# rateshifts_prior_0.1.png

png(file = 'rateshifts_prior_0.1.png', width=450, height=450);


prob.k <- function(k, poissonRatePrior=1)	{
	Denom <- (poissonRatePrior + 1)^(k+1)
	Prob <- poissonRatePrior / Denom
	return(Prob)
}

obsK <- seq(from=0, to=19, by=1)
expectedNumberofShifts <- 10

priorD <- sapply(obsK, prob.k, poissonRatePrior=1/expectedNumberofShifts)

par(mar=c(5,5,1,1))
plot(1, 1, xlim=c(0,max(obsK)), ylim=c(0,0.5), type='n', xlab="Number of shifts", ylab="Probability")
points(obsK, priorD, pch=21, bg='red', type='b')

dev.off();

# -----------------------------------------
# distinctshiftconfigs1.png

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



# ---------------------------------------
# distinctshiftconfigs2.png

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



# ---------------------------------------
# cohort_whales_illustrated.png

data(whales, events.whales)
edata_whales <- getEventData(whales, events.whales, burnin=0.1)

cmat <- getCohortMatrix(edata_whales)

png(height=1000, width=1000, file = "cohort_whales_illustrated.png")

cohorts(cmat, edata_whales, lwd=3, use.plot.bammdata = TRUE, pal='temperature')

text(0.355, 0.67, 'A', cex=7)
text(0.75, 0.67, 'B', cex=7)
text(0.355, 0.22, 'C', cex=7)
text(0.75, 0.22, 'D', cex=7)
text(0.568, 0.465, 'E', cex=3.5)
text(0.355, 0.465, 'F', cex=3.5)
text(0.568, 0.22, 'F', cex=3.5)
text(0.75, 0.465, 'G', cex=3.5)
text(0.568, 0.67, 'G', cex=3.5)

dev.off()










