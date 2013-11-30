
# library('BAMMtools');
# library(diversitree);

# ed <- getEventData(....);
# birdtree <- read.tree(....);

# Compute the marginal shift probabilities on branches
marg_tree <- marginalShiftProbsTree(ed);

# Compute the cumulative shift probabilities for branches:
cst <- cumulativeShiftProbsTree(ed);

# Compute the maximum shift credibility configuration:
# This is the joint distribution of shift events that maximizes the 
#	marginal probability of the data.
msc <- maximumShiftCredibilityTree(ed);

# get the relevant rate shift nodes from the 
#	maximum shift credibility configuration
nodes <- getShiftNodesFromIndex(ed, index = msc$sampleindex);



######### get scaling and colors for shift nodes based on
#			their importance

library(gplots);
pal <- rich.colors(10);

 
cexmin <- 2;
cexmax <- 6;

size <- cexmin + probvec*(cexmax - cexmin);

 
# Get marginal probabilities associated with each node in the the 
#		maximum shift credibility configuration
probvec <- numeric(length(nodes));

for (i in 1:length(nodes)){
	probvec[i] <- marg_tree$edge.length[marg_tree$edge[,2] == nodes[i] ];
}


## Get colors for nodes based on marginal shift probs
colvec <- rep(pal[length(pal)], length(nodes));

cuts <- seq(0, 1, length.out=length(pal)+1); 

for (i in 1:(length(cuts) - 1)){
	for (k in 1:length(probvec)){
		if (probvec[k] > cuts[i] & probvec[k] <= cuts[i+1]){
			colvec[k] <- pal[i];
		}
	}	
}


lmat <- matrix(c(1,1,1,1,1,1,1, 2), nrow=1);


edgewid <- rep(0.5, length(birdtree$edge.length));
edgecol <- rep('gray60', length(birdtree$edge.length));

edgewid[cst$edge.length > 0.9] <- 0.7;
edgecol[cst$edge.length > 0.9] <- pal[10];
 

## First, the marginal shift probability tree:
quartz.options(height=5, width=17);
plot.new();
layout(lmat);

plot.phylo(ed, show.tip.label=F, edge.width=0.8, edge.color='gray40', no.margin=T, type = 'p', direction='upwards');
nodelabels(node=nodes, cex=size, pch=21, bg=colvec);

plot.new();
par(mar=c(3,0,3,6));
plot.window(xlim=c(0,4), ylim=c(-0.3, 1.3));
for (i in 1:(length(cuts)-1)){
	xco <- c(0,2,2, 0);
	yco <- c(cuts[i], cuts[i], cuts[i+1], cuts[i+1]);
	polygon(x=xco, y=yco, col=pal[i]);
}
axis(side=4,at=seq(0,1, by=0.2), cex.axis=2, font=2, las=1, pos=2.2)
mtext(text='Shift probability', side=4, line=1, cex=1.8) 


##############################################
###### THe cumulative shift probability tree:
#	
#  Highlight all branches (red) with a probability > 0.95 of 
#	having rate dynamics different from the root of the tree.   

edgecol <- rep('gray40', length(birdtree$edge.length));
edgecol[cst$edge.length > 0.95] <- 'red';

quartz.options(height=5, width=17);
plot.new();
layout(lmat);
 
plot.phylo(ed, show.tip.label=F, edge.width=0.8, edge.color=edgecol, no.margin=T, type = 'p', direction='upwards');
nodelabels(node=nodes, cex=size, pch=21, bg=colvec);
 
plot.new();
par(mar=c(3,0,3,6));
plot.window(xlim=c(0,4), ylim=c(-0.3, 1.3));
for (i in 1:(length(cuts)-1)){
	xco <- c(0,2,2, 0);
	yco <- c(cuts[i], cuts[i], cuts[i+1], cuts[i+1]);
	polygon(x=xco, y=yco, col=pal[i]);
}
axis(side=4,at=seq(0,1, by=0.2), cex.axis=2, font=2, las=1, pos=2.2)
mtext(text='Shift probability', side=4, line=1, cex=1.8) 









