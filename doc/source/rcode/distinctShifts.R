
library(BAMMtools);
data(whales);
data(events.whales);
data(prior.whales)
ed <- getEventData(whales, events.whales, burnin=0.1);
priordist <- getBranchShiftPriors(whales, prior.whales)
cset <- credibleShiftSet(ed, threshold = priordist)

quartz.options(height=5, width=8);
plot.credibleshiftset(cset, plotmax=6, lwd=1.1);



