setwd('~/bamm/doc/source/figures')
library(BAMMtools)
library(TeachingDemos)

#################################
## INTRODUCTION PAGE
#################################

## xIntroFig_whalerates.png
library(BAMMtools)
data(whales, events.whales)
edata_whales <- getEventData(whales, events.whales, burnin=0.1)

png(height=500, width=500, file = "xIntroFig_whalerates.png");
par(mar=c(1,1,1,1))
q <- plot.bammdata(edata_whales, lwd=3, pal="temperature")

subplot(ratesHistogram(q, plotBrks=FALSE), x=7, y=75, size=c(2, 2))

dev.off();
