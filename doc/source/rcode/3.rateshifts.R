
### R code for generating figures on rateshifts: theoretical background page
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

######
cst <- cumulativeShiftProbsTree(ed)



### sample shift configurations


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

##########
# node 140:
# node 141

library(BAMMtools)
data(whales, mcmc.whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.1)

nodemat <- matrix(0, nrow=length(ed$eventData), ncol=max(ed$edge))

 
for (i in 1:length(nn)){
	tmp <- ed$eventData[[i]]
	nodemat[i,tmp$node] <- 1
}



plot(nodemat[,146])

ss <- spectrum0.ar(nodemat[,140])



plot.new()
par(mfrow=c(3,1))
plot(nodemat[,132], type = "l")
plot(nodemat[,140], type = "l")
plot(nodemat[,141], type = "l")




#############
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

########
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


########

png(height=700, width=700, file = "x_whale_credibleshiftset.png")

	data(whales, events.whales)
	edata <- getEventData(whales, events.whales, burnin=0.1)
	css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=3) 
	plot.credibleshiftset(css, tau = 0.002, breaksmethod="quantile", lwd=1.5, cex=2)
dev.off()






