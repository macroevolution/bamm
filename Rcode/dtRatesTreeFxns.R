############################################
#	dtRate(ephy,tau)
#
#	A function to calculate approximations of
#	mean instantaneous speciation rates or
#	phenotypic rates along each branch. Still has 
#	some bugs, could be sped up, too.
#
#	Arguments: ephy = a bammdata object
#	           tau = fraction of tree height for approximation (e.g. 0.01)
#
#	Returns: an ephy object with a list appended containing a vector of branch
#			 rates and the step size used for calculation.
dtRates = function(ephy,tau)
{
	if(!'bamm-data' %in% class(ephy))
	{
		stop('Function requires a bammdata object');
	}
	phy = as.phylo.bammdata(ephy);
	phy = getStartStopTimes(phy);
	if(attributes(phy)$order != 'cladewise')
	{
		phy = reorder(phy,'cladewise');
	}
	tH = max(branching.times(phy));
	
	segs = segMap(phy$edge[,2],phy$begin/tH,phy$end/tH,tau);
	rates = numeric(nrow(segs));
	
	nsamples = length(ephy$eventBranchSegs);
	for(k in 1:nsamples)
	{
		eventSegs = ephy$eventBranchSegs[[k]];
		eventData = ephy$eventData[[k]];
		for(j in 1:nrow(eventSegs))
		{
			node = eventSegs[j,1];
			ev = eventSegs[j,4];
			Start = eventData[eventData$index == ev,]$time;
			lam1 = eventData[eventData$index == ev,]$lam1;			
			lam2 = eventData[eventData$index == ev,]$lam2;
			
			isGoodSeg = segs[,1] == node;
			isGoodStart = segs[,2]*tH >= eventSegs[j,2];
			isGoodEnd = segs[,3]*tH <= eventSegs[j,3];
			
			relStart = segs[isGoodSeg & isGoodStart & isGoodEnd, 2]*tH - Start;
			relEnd = segs[isGoodSeg & isGoodStart & isGoodEnd, 3]*tH - Start;
			if(length(relStart))
			{
				rates[isGoodSeg & isGoodStart & isGoodEnd] = 
						rates[isGoodSeg & isGoodStart & isGoodEnd] + 
							branchMeanRateExponential(relStart,relEnd,lam1,lam2)/nsamples;		
			}
			else		#straddler, fudge for now
			{
				isGoodStart = segs[,2]*tH <= eventSegs[j,3];
				isGoodEnd = segs[,3]*tH >= eventSegs[j,3];
				relStart = segs[isGoodSeg & isGoodStart & isGoodEnd, 2]*tH - Start;
				relEnd = segs[isGoodSeg & isGoodStart & isGoodEnd, 3]*tH - Start;
				rates[isGoodSeg & isGoodStart & isGoodEnd] = 
						rates[isGoodSeg & isGoodStart & isGoodEnd] + 
							branchMeanRateExponential(relStart,relEnd,lam1,lam2)/nsamples;
			}
		}
	}
	ephy$dtrates = list(tau=tau,rates=rates);
	return(ephy);	
}

######################################
#	Internal function called by dtRates
#
#
segMap = function(nodes,begin,end,tau)
{
	foo = function(x,tau)
	{
		ret = seq(x[2],x[3],by=tau);
		if(length(ret) == 1) return(matrix(x,nrow=1));
		ret = seq(x[2],x[3],length.out=length(ret));
		ret = rep(ret,each=2); ret=ret[2:(length(ret)-1)];
		ret = matrix(ret,ncol=2,byrow=TRUE);
		return(cbind(matrix(rep(as.integer(x[1]),nrow(ret)),ncol=1), ret));
	}
	times = cbind(nodes,begin,end);
	ret = apply(times,1,foo,tau);
	return(do.call(rbind,ret));	
}


###########################################
#	polartree(...)
#
#	A function to plot polarized trees
#	May be used to plot trees with branches
#	colored by a continuously changing value.
#	Lots of improvements could be made to how
#	it does this.
#
#	Arguments: 
#	ephy = bammdata object
#	open_angle = angle in degrees to prevent overplotting first and last tips 	
#	rbf = fraction of tree height that the length of root branch should be
#		  Use >0 if you want an arc connecting the first two branches
#	lwd = line width
#	edge.color = edge colors
#	xlim,ylim = plotting window
#	labels = should the tip labels be plotted
#	show = should the tree be plotted
#	colorize = should the branches be colored by a continuously changing value
#
#	Returns invisibly: coordinates of the beginning and end of each branch and their theta value
#
polartree = function(ephy,open_angle=10,rbf=0.001,lwd=1,edge.color=1,xlim=c(-1,1),ylim=c(-1,1),labels=FALSE,show=TRUE,colorize=FALSE)
{
	
	phy = as.phylo.bammdata(ephy);
	if(attributes(phy)$order != 'cladewise')
	{
		phy = reorder.phylo(phy,order='cladewise');
	}
	if(colorize)
	{
		if(!'dtrates' %in% names(ephy))
		{
			colorize = FALSE;
			warning('No values to color with. Plotting tree normally.');
		}
		else
		{
			require(gplots);
			cols = rich.colors(32);
			tau = ephy$dtrates$tau;
			edge.color = cols[1 + 31*ephy$dtrates$rates/max(ephy$dtrates$rates)];
		}
	}
	
	phy = getStartStopTimes(phy);
	tH = max(branching.times(phy));
	open_angle = open_angle*(pi/180);
	theta_step = (2*pi-open_angle)/phy$Nnode;
	
	theta = matrix(0,nrow(phy$edge),3);
	for(node in tree.traverse(phy, phy$Nnode+2,'postorder'))
	{
		if(node <= phy$Nnode+1)
		{
			theta[phy$edge[,2]==node,] = (node-1)*theta_step;
		}
		else
		{
			isChild = phy$edge[,2] %in% phy$edge[phy$edge[,1] == node,2];
			dth = sum(theta[isChild,1])/2;
			if(node == phy$Nnode+2)
			{
				root = c(dth,theta[isChild,1])
				next;
			}
			theta[phy$edge[,2] == node,] = c(dth,theta[isChild,1]);
		}
	}
	rb = tH*rbf;
	theta = rbind(root,theta);
	x0 = c(rb,rb+(phy$begin/tH))*cos(theta[,1]);
	y0 = c(rb,rb+(phy$begin/tH))*sin(theta[,1]);
	x1 = c(rb,rb+(phy$end/tH))*cos(theta[,1]);
	y1 = c(rb,rb+(phy$end/tH))*sin(theta[,1]);
	ret = cbind(x0[-1],y0[-1],x1[-1],y1[-1],theta[-1,1]);
	if(colorize)
	{
		p = cbind(x0[-1],y0[-1],x1[-1],y1[-1], phy$edge[,2]);
		p = apply(p,1,matrify,tau);
		p = do.call(rbind, p);
		x0 = c(x0[1],p[,1]);x1=c(x1[1],p[,2]);y0=c(y0[1],p[,3]);y1=c(y1[1],p[,4]);
	}
	if(show)
	{
		plot.new();
		plot.window(xlim=xlim+c(-rb,rb),ylim=ylim+c(-rb,rb),asp=1);
		if(length(edge.color)==1)
		{
			edge.color = rep(edge.color,length(x0)-1);
		}
		offset = table(p[,5])[as.character(unique(p[,5]))];
		arc.color = c(edge.color[1],edge.color[match(unique(p[,5]),p[,5])+offset-1]);
		
		arc(0,0,theta[,2],theta[,3],c(rb,rb+phy$end/tH),border=arc.color,lwd=lwd);
		segments(x0,y0,x1,y1,col=c(edge.color[1],edge.color),lwd=lwd,lend=2);
		
		if(labels)
		{
			for(k in 1:length(phy$tip.label))
			{
				text(ret[phy$edge[,2]==k,3],ret[phy$edge[,2]==k,4],phy$tip.label[k],cex=.5,
				srt = (180/pi)*theta[-1,][phy$edge[,2]==k,1],adj=c(0,NA));	
			}
		}
	}
	rownames(ret) = phy$edge[,2];
	colnames(ret) = c('x0','y0','x1','y1','theta');
	invisible(ret);
}



##################################
#	Functions called by polartree
#
#

##################################
#	x,y = coordinates of center of curvature of arc, e.g. (0,0)
#	theta1 = initial theta of arc (radians)
#	theta2 = ending theta of arc (radians)
#	rad = radius of arc
arc = function(x,y,theta1,theta2,rad,border,...)
{
	step = (theta2-theta1)/100;
	xv = x+rad*cos(seq(theta1,theta2,step));
	yv = y+rad*sin(seq(theta1,theta2,step));
	polygon(c(xv,rev(xv)),c(yv,rev(yv)),border=border,...);
}

arc = Vectorize(arc);

##################################
#	D.Rabosky
#	get vector of beginning and ending times
#	for each branch. time 0 at root
#
#
getStartStopTimes = function(phy)
{
	bt = branching.times(phy);
	bt = max(bt) - bt;
	begin = bt[as.character(phy$edge[,1])];
	end = begin + phy$edge.length;
	phy$begin = as.numeric(begin);
	phy$end = as.numeric(end);
	return(phy);
}

##################################
#	traverse a tree in 'phylo' format from specified node. 
#	if tips.only = TRUE will return tips only
#	if internal.only = TRUE will return internal nodes only
tree.traverse = function(phy,node,order='preorder',tips.only=FALSE, internal.only=FALSE,path=NULL)
{
	if(order == 'preorder')
	{
		path = c(path,node);
	}
	if(length(which(phy$edge[,1] == node)))
	{
		for(child in phy$edge[phy$edge[,1] == node,2])
		{
			path = tree.traverse(phy,child,path=path,order=order);	
		}
	}
	if(order == 'postorder')
	{
		path = c(path,node);
	}
	if(tips.only)
	{
		return(path[which(path <= length(phy$tip.label))]);	
	}
	else if(internal.only)
	{
		return(path[which(path > length(phy$tip.label))]);	
	}
	else
	{
		return(path);
	}
}

##################################
#	Internal function called by polartree
#	to plot approximation of continuous rate
#
matrify = function(x,tau)
{
	bn = sqrt((x[3]-x[1])^2 + (x[4]-x[2])^2);
	len = bn/tau;
	
	j = seq(x[1],x[3],length.out=len);
	if(length(j) == 1) return(matrix(x[c(1,3,2,4,5)],nrow=1));
		
	k = seq(x[2],x[4],length.out = len);
	
	j = rep(j,each=2); j = j[2:(length(j)-1)];
	j = matrix(j,ncol=2,byrow=TRUE);	
	k = rep(k,each=2); k = k[2:(length(k)-1)];
	k = matrix(k,ncol=2,byrow=TRUE);	
	l = matrix(rep(x[5],nrow(j)),ncol=1);
	return(cbind(j,k,l));	
}

#####################################
#	D.Rabosky
#
#
as.phylo.bammdata <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}		
	
	newphylo <- list();
	newphylo$edge <- ephy$edge;
	newphylo$Nnode <- ephy$Nnode;
	newphylo$tip.label <- ephy$tip.label;
	newphylo$edge.length <- ephy$edge.length;
	class(newphylo) <- 'phylo';
	attributes(newphylo)$order <- 'cladewise';
	return(newphylo);
}
