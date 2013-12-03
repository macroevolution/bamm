

eventfname <- '/Users/danrabosky/DanWork/bamm/devel/bamm/doc/sweave/data/event_data.txt';
whaletree <- read.tree('/Users/danrabosky/DanWork/bamm/devel/bamm/doc/sweave/data/whaletree.tre')
source('/Users/danrabosky/DanWork/bamm/devel/bamm/Rcode/Load_BAMMtools.R')


## Process the event data
ed <- getEventData(whaletree, eventfname, burnin=0.25, nsamples=200, verbose=T);

## Generate the rate through time matrices
##		Do for ALL lineages, and for dolphins only, 
##		and for BACKGROUND lineages

bmat <- getRateThroughTimeMatrix(ed);
bmat.dolphin <- getRateThroughTimeMatrix(ed, node=140, nodetype='include');
bmat.background <- getRateThroughTimeMatrix(ed, node=140, nodetype='exclude');



quartz.options(height=10, width=10);
par(mfrow=c(2,2));

plotRateThroughTime(bmat, ratetype='speciation', intervalCol='red', avgCol='red', start.time = 0, end.time=35);


plotRateThroughTime(bmat, ratetype='speciation', intervalCol='red', avgCol='red', start.time = NULL, end.time=NULL);


