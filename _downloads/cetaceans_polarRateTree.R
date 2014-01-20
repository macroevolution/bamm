

# Make sure BAMMtools package is installed.

library(BAMMtools);

data(cetaceans);
data(events.cetaceans);

ed <- getEventData(phy=cetaceans, eventdata = events.cetaceans, burnin=0.1);

# This plots the 
plot.bammdata(ed, method='polar', pal='temperature', lwd=2);



