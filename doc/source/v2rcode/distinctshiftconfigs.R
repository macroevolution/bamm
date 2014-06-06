
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


 




