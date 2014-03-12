

library(BAMMtools);
data(whales);
v <- whales;

data(events.whales);

ed <- getEventData(v, events.whales, burnin=0.1);
 
msc <- maximumShiftCredibility(ed);

quartz.options(height=5, width=10, dpi=72);
par(mfrow=c(1,2));

plot(v, show.tip.label=F, edge.width=1.1);
addBAMMshifts(ed, method='phylogram', index=msc$sampleindex, pch=21, bg='red', cex=1.7, par.reset=FALSE)
 
ee <- rep('black', nrow(v$edge));
ee[v$edge[,2] %in% c(132, 140, 129, 141)] <- 'red';

plot(v, edge.color=ee, show.tip.label=F, edge.width=1.1);

##########

mm <- marginalShiftProbsTree(ed)

data(primates);
data(events.primates);
ep <- getEventData(primates, events.primates,burnin=0.1, type ='trait');

mm <- marginalShiftProbsTree(ep)


######
v2 <- drop.tip(v, tip = v$tip.label[8:55]);

### 
library(RColorBrewer)
quartz.options(height=6, width=8)
plot.bammdata(ed, pal=list(c('blue', 'gray', 'red')));


plot(v, show.tip.label=F)


####
quartz.options(height=6, width=8)
library(BAMMtools)
library(RColorBrewer)

data(whales, events.whales)

ed <- getEventData(whales, events.whales);


plot.bammdata(ed, method='phylogram', pal= c('blue', 'gray', 'red'), lwd=1.5, tau=0.001)

#####
mm <- marginalShiftProbsTree(ed)
cst <- cumulativeShiftProbsTree(ed)

######
data(primates)
data(events.primates)
ed <- getEventData(primates, events.primates, burnin=0.1, type = 'trait')
mm <- marginalShiftProbsTree(ed)

xx <- computeJointShiftCorrelations(ed, threshold=0.05)

#236,281, 287 vs 235

mm$edge.length[mm$edge[,2] == 235]
mm$edge.length[mm$edge[,2] == 236]

plot.bammdata(ed, pal='temperature');

########
quartz.options(height=8, width=10, dpi=72);
par(oma=c(0,0,0,0));
par(mar=c(0,0,0,0));
par(mfrow=c(2,2));
plot.phylo(primates, show.tip.label=F, edge.color='gray30');
addBAMMshifts(ed, method='phylogram', index=1, cex=1.75, par.reset=FALSE);

plot.phylo(primates, show.tip.label=F, edge.color='gray30');
addBAMMshifts(ed, method='phylogram', index=50, cex=1.75, par.reset=FALSE);

plot.phylo(primates, show.tip.label=F, edge.color='gray30');
addBAMMshifts(ed, method='phylogram', index=100, cex=1.75, par.reset=FALSE);

plot.bammdata(ed, pal=c("darkgreen","yellow2","red"), lwd=1, tau=0.001, legend=F)

plot.bammdata(ed, pal=c("darkgreen","yellow2","red"), lwd=1, tau=0.001, legend=T)


cst <- cumulativeShiftProbsTree(ed);

######## Plot color matrix of pairwise shifts:

xx <- computeJointShiftCorrelations(ed, threshold=0.10)
dx <- nrow(xx$p.value);
colmat <- matrix('white', nrow=dx, ncol=dx);

quartz.options(height=8, width=8, dpi=72);
plot.new();
plot.window(xlim=c(0, dx+1), ylim=c(0, dx+1));

for (i in 1:(dx-1)){
	for (j in (i+1):dx){
		xv <- c(i, i, i+1, i+1);
		yv <- c(j, j+1, j+1, j);
		if (xx$p.value[i,j] < 0.05 & xx$phi[i,j] < 0){
			col <- 'blue';
		}else if (xx$p.value[i,j] < 0.05 & xx$phi[i,j] > 0){
			col <- 'red';
		}else{
			col <- 'white';
		}
		polygon(x=xv, y=yv, border='black', col=col);
	}
}


#########

data(events.primates, primates);
ed <- getEventData(primates, events.primates, type = 'trait', burnin=0.1)
msc <- maximumShiftCredibility(ed)


quartz.options(height=5, width=8)
plot(primates, show.tip.label=F, edge.width=0.7)
addBAMMshifts(ed, method='phylogram', index=msc$sampleindex, cex=2)

 










