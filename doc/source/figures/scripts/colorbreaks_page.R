setwd('~/bamm/doc/source/figures')

library(BAMMtools)
library(TeachingDemos)


#####################################
## COLOR BREAKS PAGE
#####################################

# ------------------------------------------------------------------
# rateHistograms.png: rates histograms figure for colorbreaks page

require(BAMMtools)
data(primates, events.primates)
ed <- getEventData(primates, events.primates, burnin=0.25, type = 'trait')

png('rateHistograms.png', height=12, width=8, units='in', res=300)
par(mfrow=c(5,1))

q <- plot.bammdata(ed, breaksmethod='linear', show=FALSE)
ratesHistogram(q, xlab='', ylab='')    
title(main='linear', cex=2)

q <- plot.bammdata(ed, breaksmethod='linear', show=FALSE, logcolor=TRUE)
ratesHistogram(q, xlab='', ylab='')   
title(main='linear - log', cex=2)

q <- plot.bammdata(ed, breaksmethod='linear', show=FALSE, color.interval=c(NA, 0.12))
ratesHistogram(q, xlab='', ylab='') 
title(main='linear - color.interval', cex=2)

q <- plot.bammdata(ed, breaksmethod='quantile', show=FALSE)
ratesHistogram(q, xlab='', ylab='') 
title(main='quantile', cex=2)

q <- plot.bammdata(ed, breaksmethod='jenks', show=FALSE)
ratesHistogram(q, xlab='', ylab='') 
title(main='jenks', cex=2)

dev.off()


#----------------------------------------
# breaksmethods example for colorbreaks page

png('breaksmethodPhylorates.png', width=12, height=20, units='in', res=300)

par(mfrow=c(3,2), xpd=T)
q <- plot.bammdata(ed, tau=0.001, breaksmethod='linear', lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=2)
title(main='linear',cex.main=2)
addBAMMlegend(q, location=c(0,1,140, 250))
q <- plot.bammdata(ed, tau=0.001, breaksmethod='linear', logcolor=T, lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=2)
title(main='linear - log',cex.main=2)
addBAMMlegend(q, location=c(0,1,140, 250))
q <- plot.bammdata(ed, tau=0.001, breaksmethod='linear', color.interval=c(NA,0.12), lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=2)
title(main='linear - color.interval',cex.main=2)
addBAMMlegend(q, location=c(0,1,140, 250))
q <- plot.bammdata(ed,tau=0.001, breaksmethod='quantile', lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=2)
title(main='quantile',cex.main=2)
addBAMMlegend(q, location=c(0,1,140, 250))
q <- plot.bammdata(ed,tau=0.001, breaksmethod='jenks', lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=2)
title(main='jenks',cex.main=2)
addBAMMlegend(q, location=c(0,1,140, 250))

dev.off()
