library(BAMMtools)
data(whales, events.whales, prior.whales)
edata_whales <- getEventData(whales, events.whales, burnin=0.1)
plot.bammdata(edata_whales, lwd=2, method="polar", pal="temperature")


png(height=500, width=500, file = "bammgraphgallery1.png");
par(mar=c(1,1,1,1))
plot.bammdata(edata_whales, lwd=3, method="polar", pal="temperature")
dev.off();


##### #########################################
# bammgraphgallery2
##### #########################################

 
ixx <- rep(c(10, 30, 40), 2);

png(height=700, width=1000, file = "bammgraphgallery2.png")
plot.new()
par(mfrow=c(2,3));
 
# Here we set up the color schemes 
colschemes <- list();
colschemes[1:3] <- 'temperature';
colschemes[4:6] <- list(c('blue', 'gray', 'red'));

for (i in 1:length(ixx)) {
	par(mar=c(0,0,0,0))
	index <- ixx[i];
	eventsub <- subsetEventData(edata_whales, index=index);
	plot.bammdata(eventsub, method='polar', pal= colschemes[[i]], par.reset=FALSE, lwd=3);
	addBAMMshifts(eventsub, method='polar', index=1, col='white', bg='black', cex=5, par.reset=FALSE);
}
dev.off();


##### #########################################
# bammgraphgallery3
##### #########################################

data(prior.whales)
pset <- getBranchShiftPriors(whales, prior.whales)
cset <- credibleShiftSet(edata_whales, pset, BFcriterion=3)
png(height=700, width=700, file = "bammgraphgallery3.png")
plot.credibleshiftset(cset, lwd=2.5)
dev.off()


##### #########################################
# bammgraphgallery4
##### #########################################


png(height=1000, width=1000, file = "bammgraphgallery4.png")

cmat <- getCohortMatrix(edata_whales)
cohorts(cmat, edata_whales, lwd=3, pal="temperature")

dev.off()


##### #########################################
# bammgraphgallery5
##### #########################################

png(height=700, width=700, file = "bammgraphgallery5.png")


data(primates, events.primates, prior.primates)
ed_prim <- getEventData(primates, events.primates, burnin=0.1, type = "trait")
pprior <- getBranchShiftPriors(primates, prior.primates)
css_prim <- credibleShiftSet(ed_prim, pprior)
plot.credibleshiftset(css_prim, lwd=1.7, plotmax=4)

dev.off()


 
##### #########################################
# bammgraphgallery6
##### #########################################

png(width=1000, height=350, file = "bammgraphgallery6.png")

par(mfrow=c(1,3))
st <- max(branching.times(whales))
plotRateThroughTime(edata_whales, intervalCol="red", avgCol="red", start.time=st, ylim=c(0,1), cex.axis=2)
text(x=30, y= 0.8, label="All whales", font=4, cex=2.0, pos=4)
plotRateThroughTime(edata_whales, intervalCol="blue", avgCol="blue", start.time=st, node=141, ylim=c(0,1),cex.axis=1.5)
text(x=30, y= 0.8, label="Dolphins only", font=4, cex=2.0, pos=4)

plotRateThroughTime(edata_whales, intervalCol="darkgreen", avgCol="darkgreen", start.time=st, node=141, nodetype = "exclude", ylim=c(0,1), cex.axis=1.5)
text(x=30, y= 0.8, label="Non-dolphins", font=4, cex=2.0, pos=4)



dev.off()



##### #########################################
# bammgraphgallery7
##### #########################################

png(width=1000, height=350, file = "bammgraphgallery7.png")
 
par(mfrow=c(1,3))
st <- max(branching.times(whales))
plotRateThroughTime(edata_whales, avgCol="black", start.time=st, ylim=c(0,1), cex.axis=2, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="All whales", font=4, cex=2.0, pos=4)
plotRateThroughTime(edata_whales, avgCol="black", start.time=st, node=141, ylim=c(0,1),cex.axis=1.5,intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="Dolphins only", font=4, cex=2.0, pos=4)

plotRateThroughTime(edata_whales, avgCol="black", start.time=st, node=141, nodetype = "exclude", ylim=c(0,1), cex.axis=1.5,intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="Non-dolphins", font=4, cex=2.0, pos=4)



dev.off()


##### #########################################
# bammgraphgallery8
##### #########################################













