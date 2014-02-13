




library(BAMMtools);
data(whales);
data(events.whales);
ed <- getEventData(whales, events.whales, burnin=0.2);


configs <- distinctShiftConfigurations(ed, threshold=0.05);
 
 
plot.bammshifts(configs, ed);
 
 
 