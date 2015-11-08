setwd('~/bamm/doc/source/figures')

library(BAMMtools)
library(TeachingDemos)

################################
## GRAPH GALLERY
################################

#--------------------------------------------------
# whales_polar.png: polar phyloplot of whales

library(BAMMtools)
data(whales, events.whales)
edata_whales <- getEventData(whales, events.whales, burnin=0.1)

png(height=500, width=500, file = "whales_polar.png");
par(mar=c(1,1,1,1))
plot.bammdata(edata_whales, lwd=3, method="polar", pal="temperature")
dev.off();

# ----------------------------------------------------------------
# bammgraphgallery2.png: primates example of colorbreaks methods

data(primates, events.primates)
ed <- getEventData(primates, events.primates, burnin=0.25, nsamples=1000, type = 'trait')

png("bammgraphgallery2.png", height=5, width=7, units='in', res=300)
par(mfrow=c(1,3), mar=c(1, 0.5, 0.5, 0.5), xpd=TRUE)

q <- plot.bammdata(ed, tau=0.001, breaksmethod='linear', lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=1)
title(sub='linear',cex.sub=2, line=-1)
addBAMMlegend(q, location=c(0, 1, 140, 220))

q <- plot.bammdata(ed, tau=0.001, breaksmethod='linear', color.interval=c(NA,0.12), lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=1)
title(sub='linear - color.interval',cex.sub=2, line=-1)
addBAMMlegend(q, location=c(0, 1, 140, 220))

q <- plot.bammdata(ed, tau=0.001, breaksmethod='jenks', lwd=2)
addBAMMshifts(ed, par.reset=FALSE, cex=1)
title(sub='jenks',cex.sub=2, line=-1)
addBAMMlegend(q, location=c(0, 1, 140, 220))

dev.off()


# -------------------------------------------------------------
# whales_sepRateShiftConfigs.png: separate rate configurations

data(whales, events.whales)
edata_whales <- getEventData(whales, events.whales, burnin=0.1)

ixx <- rep(c(10, 30, 40), 2);

png(height=700, width=1000, file = "whales_sepRateShiftConfigs.png")
plot.new()
par(mfrow=c(2,3));
 
# Here we set up the color schemes 
colschemes <- list();
colschemes[1:3] <- 'temperature';
colschemes[4:6] <- list(c('blue', 'gray', 'red'));

for (i in 1:length(ixx)) {
	par(mar=c(0,0,0,0))
	index <- ixx[i];
	eventsub <- subsetEventData(edata_whales, index=index);
	plot.bammdata(eventsub, method='polar', pal= colschemes[[i]], par.reset=FALSE, lwd=3);
	addBAMMshifts(eventsub, method='polar', index=1, col='white', bg='black', cex=5, par.reset=FALSE);
}
dev.off();


# -------------------------------------------------------------
# whales_distinctShiftConfigs.png: distinct shift configurations

pset <- getBranchShiftPriors(whales, expectedNumberOfShifts = 1)
cset <- credibleShiftSet(edata_whales, expectedNumberOfShifts = 1, threshold=3)
png(height=700, width=700, file = "whales_distinctShiftConfigs.png")
plot.credibleshiftset(cset, lwd=2.5)
dev.off()



# -------------------------------------------------------------
# whales_cohort.png: cohort matrix for whales

png(height=1000, width=1000, file = "whales_cohort.png")

data(whales, events.whales)
edata_whales <- getEventData(whales, events.whales, burnin=0.1)

cmat <- getCohortMatrix(edata_whales)
cohorts(cmat, edata_whales, lwd=3, pal="temperature", use.plot.bammdata = TRUE)

dev.off()


# ------------------------------------------------------
# primates_credShiftSet.png: trait credible shift sets

png(height=700, width=700, file = "primates_credShiftSet.png")

data(primates, events.primates)
ed_prim <- getEventData(primates, events.primates, burnin=0.1, type = "trait", nsamples=1000)
css_prim <- credibleShiftSet(ed_prim, expectedNumberOfShifts = 1)
plot.credibleshiftset(css_prim, lwd=1.7, plotmax=4)

dev.off()



# --------------------------------------------------------------------
# whales_RatesThroughTime.png: separated rate through time plots for whales

data(whales, events.whales)
edata_whales <- getEventData(whales, events.whales, burnin=0.1)

png(width=1000, height=350, file = "whales_RatesThroughTime.png")

par(mfrow=c(1,3))
st <- max(branching.times(whales))
plotRateThroughTime(edata_whales, intervalCol="red", avgCol="red", start.time=st, ylim=c(0,1), cex.axis=2)
text(x=30, y= 0.8, label="All whales", font=4, cex=2.0, pos=4)
plotRateThroughTime(edata_whales, intervalCol="blue", avgCol="blue", start.time=st, node=140, ylim=c(0,1),cex.axis=1.5)
text(x=30, y= 0.8, label="Dolphins only", font=4, cex=2.0, pos=4)

plotRateThroughTime(edata_whales, intervalCol="darkgreen", avgCol="darkgreen", start.time=st, node=140, nodetype = "exclude", ylim=c(0,1), cex.axis=1.5)
text(x=30, y= 0.8, label="Non-dolphins", font=4, cex=2.0, pos=4)

dev.off()


# --------------------------------------------------------------------------------
# whales_RatesThroughTimeBW.png: separated rate through time plots for whales, grayscale

png(width=1000, height=350, file = "whales_RatesThroughTimeBW.png")
 
par(mfrow=c(1,3))
st <- max(branching.times(whales))
plotRateThroughTime(edata_whales, avgCol="black", start.time=st, ylim=c(0,1), cex.axis=2, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="All whales", font=4, cex=2.0, pos=4)
plotRateThroughTime(edata_whales, avgCol="black", start.time=st, node=140, ylim=c(0,1),cex.axis=1.5,intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="Dolphins only", font=4, cex=2.0, pos=4)

plotRateThroughTime(edata_whales, avgCol="black", start.time=st, node=140, nodetype = "exclude", ylim=c(0,1), cex.axis=1.5,intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="Non-dolphins", font=4, cex=2.0, pos=4)

dev.off()

