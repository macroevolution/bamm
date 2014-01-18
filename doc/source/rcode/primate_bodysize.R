
library(BAMMtools);
data(primates);
data(events.primates);

ed <- getEventData(primates, events.primates, burnin=0.2, type='trait');

quartz.options(height=8, width=8, dpi=72); # on mac osx 

plot.bammdata(ed, lwd=1.35, pal = 'temperature');

# Compute maximum shift credibility tree
msctree <- maximumShiftCredibilityTree(ed);

addBAMMshifts(ed, method='phylogram', index=2, cex=2);



 