
library(BAMMtools)
library(RColorBrewer)

data(whales, events.whales)

ed <- getEventData(whales, events.whales);

plot.bammdata(ed, method='phylogram', pal= c('blue', 'gray', 'red'))




