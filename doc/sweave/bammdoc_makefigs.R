
source("BAMMtools.R");
require("ape");
require("diversitree");

whaletree <- read.tree("data/whaletree.tre")
eventfile <- "data/event_data.txt"
bammdata <- getEventData(whaletree, eventfile, burnin=0.1, verbose=T)




### My general purpose histogram plotting function.
 

histFx <- function(x, breaks='sturges', col='black', scale=0.9, border=FALSE)
{
	z <- hist(x, breaks=breaks, plot=F);
	if (!is.na(scale)){
		sf <- scale/max(z$counts);
		z$counts <- z$counts*sf;		
		for (i in 1:(length(z$breaks)-1))	
			polygon(x=c(rep(z$breaks[i],2), rep(z$breaks[i+1],2)), y=c(0, z$counts[i], z$counts[i], 0), col=col, border=border);	
	
	}else{
		# plot density:
		for (i in 1:(length(z$breaks)-1))	
			polygon(x=c(rep(z$breaks[i],2), rep(z$breaks[i+1],2)), y=c(0, z$density[i], z$density[i], 0), col=col, border=border);			
		
	}


}



################################################

pdf('mst.pdf', height=4, width=9)
mst <- marginalShiftProbsTree(bammdata);
quartz.options(height=6, width=8, dpi=72);
plot(mst, show.tip.label=F, no.margin=T);
add.scale.bar(length=0.5, x= 0.25, y=20, lcol="gray40", lwd=2, cex=1.5);
n1 <- 140;
n2 <- 141;
nodelabels(node=141, bg="red", cex=3, pch=21);
nodelabels(node=140, bg="blue", cex=3, pch=21);

dev.off()

################################################
################################################

pdf('cst.pdf', height=4, width=9)
cst <- cumulativeShiftProbsTree(bammdata);

quartz.options(height=6, width=8, dpi=72);
plot(cst, show.tip.label=F, no.margin=T);
nodelabels(node=141, bg="red", cex=3, pch=21);
nodelabels(node=140, bg="blue", cex=3, pch=21);
add.scale.bar(length=1.0, x= 5, y=45, lcol="darkgreen", lwd=3, cex=1.5);

 
dev.off()

################################################
################################################

tax1 <- "Orcinus_orca"
tax2 <- "Delphinus_delphis"
tipnode1 <- which(whaletree$tip.label == tax1)
tipnode2 <- which(whaletree$tip.label == tax2)
dolphin_node <- getMRCA(whaletree, c(tipnode1, tipnode2))

rates_all <- getCladeRates(bammdata);


rates_dolphin<- getCladeRates(bammdata, node=dolphin_node);

rates_background <- getCladeRates(bammdata, node=dolphin_node, nodetype = "exclude");


pdf('dolphins_ratedistributions.pdf', height=5, width=10)

plot.new();
par(mar=c(6,6,1,1));
plot.window(xlim=c(0, 0.6), ylim=c(0, 25));
histFx(rates_background$lambda, col='white', border=T, breaks=10, scale=NA);
histFx(rates_dolphin$lambda, col='coral', border=T, breaks=40, scale=NA);
axis(1, at=seq(-0.1, 0.7, by=0.1), cex.axis=1.1)
axis(2, at=seq(-5, 25, by=5), labels=F, cex.axis=1.1)
mtext(side=1, text="Speciation rate", line=3, cex=1.2)
mtext(side=2, text="Density", line=1.5, cex=1.2)

dev.off()

################################################
################################################

pdf('timepoints.pdf', height=6, width=12)

plot(whaletree, show.tip.label=F, no.margin=F);
axisPhylo()
mx <- max(branching.times(whaletree));
mn <- 0;
sq <- seq(mn, mx, length.out=30);
for (i in sq){
	lines(c(i,i), c(0, 90), lwd=1, col='gray60')
}
mtext(side=1, 'Time before present', cex=1.5, line=3)

dev.off()

################################################
################################################
rtt.all <- getRateThroughTimeMatrix(bammdata)
rtt.dolphin <- getRateThroughTimeMatrix(bammdata, node=dolphin_node)
rtt.background <- getRateThroughTimeMatrix(bammdata, node=dolphin_node, nodetype="exclude")

pdf('cetaceanRTT.pdf', height=4, width=7)

lam_rates.all <- colMeans(rtt.all$lambda)
plot.new()
par(mar=c(6, 6, 1,1))
plot.window(xlim=c(0, 36), ylim=c(0, 0.6))
lines(x=rtt.all$times, y=lam_rates.all, lwd=3, col="black")
axis(1, at=seq(-5, 35, by=5));
axis(2, at=seq(-0.1, 0.7, by=0.1), las=1)
mtext(side=1, text="Time from start of radiation", line =3.5, cex=1.2)
mtext(side=2, text="Speciation rate", line =3.5, cex=1.2) 
dev.off()

################################################
################################################
pdf('dolphinssplitRTT.pdf', height=4, width=7)

lam_rates.dolphin <- colMeans(rtt.dolphin$lambda)
lam_rates.background <- colMeans(rtt.background$lambda)
plot.new()
par(mar=c(6, 6, 1,1))
plot.window(xlim=c(0, 36), ylim=c(0, 0.6))

lines(x=rtt.dolphin$times, y=lam_rates.dolphin, lwd=3, col='red')
lines(x=rtt.background$times, y=lam_rates.background, lwd=3, col='blue')
axis(1, at=seq(-5, 35, by=5));
axis(2, at=seq(-0.1, 0.7, by=0.1), las=1)
mtext(side=1, text="Time from start of radiation", line =3.5, cex=1.2)
mtext(side=2, text="Speciation rate", line =3.5, cex=1.2)
 
dev.off()

################################################
################################################




