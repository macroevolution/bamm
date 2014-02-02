# This file illustrates how you can plot a rate-through-time curve
#	with confidence intervals as shaded polygons.
#	
#  We avoid using the function plot(...) directly here, as it 
#	is generally easier to make high-quality figures by building a plot piecewise
#	and specifying all the components individually 
#	(e.g., axis labels, margins, legends, plot window, etc)

library(BAMMtools);
data(whales);
data(events.whales);


## Process the event data
ed <- getEventData(whales, events.whales, burnin=0.25, nsamples=200);

## Generate the rate through time matrices
##		Do for ALL lineages, and for dolphins only, 
##		and for BACKGROUND lineages

bmat <- getRateThroughTimeMatrix(ed);
bmat.dolphin <- getRateThroughTimeMatrix(ed, node=140, nodetype='include');
bmat.background <- getRateThroughTimeMatrix(ed, node=140, nodetype='exclude');

### Now to just extract speciation rate quantiles from these distributions:

lambda_all <- apply(bmat$lambda, 2, quantile, seq(0.1, 0.9, by=0.1));
lambda_dolphin <- apply(bmat.dolphin$lambda, 2, quantile, seq(0.1, 0.9, by=0.1));
lambda_background <- apply(bmat.background$lambda, 2, quantile, seq(0.1, 0.9, by=0.1));

###############################
########## Now the plotting

cols <- c('gray90', 'gray75', 'gray60', 'gray45');


quartz.options(height=4.5, width=12);
par(mfrow=c(1,3));

########### First, all whales

tvec <- max(bmat$times) - bmat$times;

plot.new();
par(mar=c(6,6,1,1));
plot.window(xlim=c(35, 0), ylim=c(0, 0.75));

for (i in 1:4){
	xcoords <- c(tvec, rev(tvec));
	ycoords <- c(lambda_all[i,], rev(lambda_all[10-i, ]));
	polygon(x = xcoords, y = ycoords, col=cols[i], border=F)
}

lines(x=tvec, y=colMeans(lambda_all), lwd=3, col='black');
axis(1, at=seq(40, -5, by=-5), cex.axis=1.3);
axis(2, at=seq(-0.2, 0.8, by=0.2), cex.axis=1.3, las=1);
mtext(side = 1, text="Time before present (mya)", line=3.5, cex=1.2)
mtext(side = 2, text="Speciation", line=3.5, cex=1.2)
text(x=35, y=0.6, label="All whales", font=4, cex=2, pos=4)

################## Now the dolphins

tvec <- max(bmat.dolphin$times) - bmat.dolphin$times;
plot.new();
par(mar=c(6,6,1,1));
plot.window(xlim=c(35, 0), ylim=c(0, 0.75));

for (i in 1:4){
	xcoords <- c(tvec, rev(tvec));
	ycoords <- c(lambda_dolphin[i,], rev(lambda_dolphin[10-i, ]));
	polygon(x = xcoords, y = ycoords, col=cols[i], border=F)
}

lines(x=tvec, y=colMeans(lambda_dolphin), lwd=3, col='black');
axis(1, at=seq(40, -5, by=-5), cex.axis=1.3);
axis(2, at=seq(-0.2, 0.8, by=0.2), cex.axis=1.3, las=1);
mtext(side = 1, text="Time before present (mya)", line=3.5, cex=1.2)
mtext(side = 2, text="Speciation", line=3.5, cex=1.2)
text(x=35, y=0.6, label="Dolphins only", font=4, cex=2, pos=4)


######## Now the "non-dolphin" background rate:

tvec <- max(bmat.background$times) - bmat.background$times;
plot.new();
par(mar=c(6,6,1,1));
plot.window(xlim=c(35, 0), ylim=c(0, 0.75));

for (i in 1:4){
	xcoords <- c(tvec, rev(tvec));
	ycoords <- c(lambda_background[i,], rev(lambda_background[10-i, ]));
	polygon(x = xcoords, y = ycoords, col=cols[i], border=F)
}

lines(x=tvec, y=colMeans(lambda_background), lwd=3, col='black');
axis(1, at=seq(40, -5, by=-5), cex.axis=1.3);
axis(2, at=seq(-0.2, 0.8, by=0.2), cex.axis=1.3, las=1);
mtext(side = 1, text="Time before present (mya)", line=3.5, cex=1.2)
mtext(side = 2, text="Speciation", line=3.5, cex=1.2)
text(x=35, y=0.6, label="Background: no dolphins", font=4, cex=2, pos=4)








