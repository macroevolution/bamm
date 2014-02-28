
library(BAMMtools);
data(whales);
data(events.whales);
ed <- getEventData(whales, events.whales, burnin=0.1);

dsc <- distinctShiftConfigurations(ed, threshold=0.04);

quartz.options(height=5, width=8);
plot.bammshifts(dsc, ed, plotmax=6, lwd=1.1);



