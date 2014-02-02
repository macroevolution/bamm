

# Make sure BAMMtools package is installed.

library(BAMMtools);

data(whales);
data(events.whales);

ed <- getEventData(phy=whales, eventdata = events.whales, burnin=0.1);

# This plots the 
plot.bammdata(ed, method='polar', pal='temperature', lwd=2);



