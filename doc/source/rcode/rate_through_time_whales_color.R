

library(BAMMtools);
data(events.whales);
data(whales);

## Process the event data
ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=200, verbose=F);

## Generate the rate through time matrices
##		Do for ALL lineages, and for dolphins only, 
##		and for BACKGROUND lineages

bmat <- getRateThroughTimeMatrix(ed);
bmat.dolphin <- getRateThroughTimeMatrix(ed, node=140, nodetype='include');
bmat.background <- getRateThroughTimeMatrix(ed, node=140, nodetype='exclude');



quartz.options(height=4, width=10);
par(oma=c(1,1,2,1));

par(mfrow=c(1,3));

plotRateThroughTime(bmat, intervalCol='red', avgCol='red', ylim=c(0,1));
mtext(side=3, line=1, "All Cetaceans", cex=1.2, font=4);

plotRateThroughTime(bmat.dolphin, intervalCol='blue', avgCol='blue', ylim=c(0, 1));
mtext(side=3, line=1, "Dolphins", cex=1.2, font=4);

plotRateThroughTime(bmat.background, intervalCol='darkgreen', avgCol='darkgreen', ylim=c(0, 1));
mtext(side=3, line=1, "Background (non-dolphin)", cex=1.2, font=4);

