
data(primates, events.primates, prior.primates)
prior <- getBranchShiftPriors(primates, prior.primates)
ed <- getEventData(primates, events.primates, burnin=0.1, type ='trait')
best <- getBestShiftConfiguration(ed, prior=prior, BFcriterion=25)
plot.bammdata(best, lwd=1.25)
addBAMMshifts(best, cex=2)

data(whales, events.whales, prior.whales)
prior <- getBranchShiftPriors(whales, prior.whales)
ed <- getEventData(whales, events.whales, burnin=0.1)
best <- getBestShiftConfiguration(ed, prior=prior, BFcriterion=3)
plot.bammdata(best, lwd=1.25)
addBAMMshifts(best, cex=2)



png(height=1000, width=1000, file = 'overallbestshiftconfig1.png')
plot.bammdata(best, lwd=1.75)
addBAMMshifts(best, cex=2)
dev.off()

