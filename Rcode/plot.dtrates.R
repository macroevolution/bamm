##################################
#	plot.dtrates(...)
#
#	A function to plot dynamic rates through time
#	onto a phylogeny
#
#	Arguments: ephy = a bammdata object.
#	           method = method used to plot the tree.
#	                    May be 'polar' or 'phylogram'.
#	           tau = fraction of tree height for approximation (e.g. 0.01).
#	                 This is the step size used for calculating rate changes 
#	                 along branches, so 0.01 is a step size equal to 1% of tree height.
#	           index = index of posterior sample(s). Currently may be NULL or 
#	                   a vector of integer values.  if NULL the function will use all 
#	                   posterior samples, otherwise it will use only
#	                   the samples specified in index.
#	           show = TRUE or FALSE. If TRUE the tree will plot.
#	           labels = TRUE or FALSE. If TRUE the tip labels will plot.
#	           hrates = TRUE or FALSE. If TRUE a histogram is plotted in the same
#	                    device for interpreting the meaning of plotted colors.
#	                    You will be asked to supply an anchor point for plotting.
#	           lwd = The line width used for plotting the tree.
#	           cex = Character expansion for plotting tip labels.
#	           ncolors = The number of color bins for mapping rates to colors.
#	           pal = A string or vector of strings to specify colors for mapping to rates. 
#	                 Currently this may be one of the named diverging palette options
#	                 in the RColorBrewer package documented in description of brewer.pal,
#	                 a vector of 3 valid named colors, or the string 'temperature'. 
#	                 The first and last options require RColorBrewer and gplots packages,
#	                 respectively
#	           ... = arguments passed to function used for polar tree plotting.
#	                 These should be 'vtheta', which specifies the angle in degrees
#	                 separating the first and last tips, and 'rbf', which specifies
#	                 the length of the root branch as a fraction of the total tree
#	                 height.  Values rbf > 0 will allow an arc between the first
#	                 two descendant branches from the root. rbf = 0.001 seems to be 
#	                 a good first choice

plot.dtrates = function(ephy, method='phylogram', tau=0.01, index=NULL, show=TRUE, labels=FALSE, hrates=TRUE, lwd=3, cex=1, ncolors=64, pal='temperature', ...)
{
	if ('bamm-data' %in% class(ephy)) phy = as.phylo.bammdata(ephy) else stop('Trying to work with a non-bammdata object');
	if (!('dtrates' %in% names(ephy))) ephy = dtRates(ephy, tau, index);
	
	if (method == 'polar')
	{
		if(!hasArg(vtheta) & !hasArg(rbf)) stop('vtheta and rbf are needed to plot in polar format');
		phy = getStartStopTimes(phy);
		vtheta = list(...)$vtheta; rbf = list(...)$rbf;
		ret = setPolarTreeCoords(phy,vtheta,rbf);
		tH = max(branching.times(phy));
		rb = tH*rbf;
	}	
	else if (method == 'phylogram')
	{
		ret = setPhyloTreeCoords(phy);
	}
	else
	{
		stop('Unimplemented method');
	}
	x0 = ret$segs[,1];y0=ret$segs[,2];x1=ret$segs[,3];y1=ret$segs[,4];
	
	tau = ephy$dtrates$tau;
	edge.color = colorMap(ephy$dtrates$rates,pal,ncolors);
	p = cbind(x0[-1],y0[-1],x1[-1],y1[-1],phy$edge[,2]);
	p = apply(p,1,matrify,tau);
	p = do.call(rbind, p);
	x0 = c(x0[1],p[,1]);x1=c(x1[1],p[,2]);y0=c(y0[1],p[,3]);y1=c(y1[1],p[,4]);
	offset = table(p[,5])[as.character(unique(p[,5]))];
	arc.color = c(edge.color[1],edge.color[match(unique(p[,5]),p[,5])+offset]);
	edge.color = c(edge.color[1],edge.color);
	
	if (show)
	{
		plot.new(); ofs = 0;
		if (labels)
		{
			ofs = max(nchar(phy$tip.label) * 0.03 * cex);
		}
		
		if (method == 'polar') 
		{
			plot.window(xlim=c(-1,1)+c(-rb,rb)+c(-ofs,ofs),ylim=c(-1,1)+c(-rb,rb)+c(-ofs,ofs),asp=1);
			segments(x0,y0,x1,y1,col=edge.color,lwd=lwd,lend=2);	
			arc(0,0,ret$arcs[,1],ret$arcs[,2],c(rb,rb+phy$end/tH),border=arc.color,lwd=lwd);
			if(labels)
			{
				for(k in 1:length(phy$tip.label))
				{
					text(ret$segs[-1,][phy$edge[,2]==k,3],ret$segs[-1,][phy$edge[,2]==k,4],phy$tip.label[k],cex=cex, srt = (180/pi)*ret$arcs[-1,][phy$edge[,2]==k,1],adj=c(0,NA));	
				}
			}
		}
		if (method == 'phylogram')
		{
			plot.window(xlim=c(0,1+ofs),ylim=c(0,phy$Nnode*1/(phy$Nnode+1)),asp=1);
			segments(x0,y0,x1,y1,col=edge.color,lwd=lwd,lend=2);
			isTip = phy$edge[,2] <= phy$Nnode+1; isTip = c(FALSE,isTip);
			segments(ret$arcs[!isTip,1],ret$arcs[!isTip,2],ret$arcs[!isTip,3],ret$arcs[!isTip,4],col=arc.color[!isTip],lwd=lwd);
			if(labels)
			{
				text(ret$segs[-1,][phy$edge[,2] <= phy$Nnode+1,3],ret$segs[-1,][phy$edge[,2] <= phy$Nnode+1,4], phy$tip.label, cex=cex, pos=4, offset = 0.25);
			}
		}
		if(hrates)
		{
			histRates(ephy$dtrates$rates,pal,ncolors);
		}
	}
	index = order(as.numeric(rownames(ret$segs)));
	if (method == 'phylogram')
	{
		assign("last_plot.phylo", list(Ntip = phy$Nnode+1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,3], yy = ret$segs[index,4]), envir = .PlotPhyloEnv);
	}
	else if (method == 'polar')
	{
		assign("last_plot.phylo", list(Ntip = phy$Nnode+1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,3], yy = ret$segs[index,4], theta = ret$segs[index,5], rb = rb), envir = .PlotPhyloEnv);
	}
	invisible(ret$segs[-1,]);
}

setPolarTreeCoords = function(phy,vtheta,rbf)
{
	phy = getStartStopTimes(phy);
	tH = max(branching.times(phy));
	vtheta = vtheta*(pi/180);
	theta_step = (2*pi-vtheta)/phy$Nnode;
	
	theta = matrix(0,nrow(phy$edge),3);
	for (node in tree.traverse(phy, phy$Nnode+2,'postorder'))
	{
		if (node <= phy$Nnode+1)
		{
			theta[phy$edge[,2]==node,] = (node-1)*theta_step;
		}
		else
		{
			isChild = phy$edge[,2] %in% phy$edge[phy$edge[,1] == node,2];
			dth = sum(theta[isChild,1])/2;
			if (node == phy$Nnode+2)
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
	ret = cbind(x0,y0,x1,y1,theta[,1]);
	rownames(ret) = c(phy$edge[1,1],phy$edge[,2]); colnames(ret) = c('x0','y0','x1','y1','theta');
	return(list(segs = ret, arcs = theta[,2:3]) );	
}

setPhyloTreeCoords = function(phy)
{
	ntips = length(phy$tip.label);
	tH = max(branching.times(phy));
	
	dy = 1/ntips;
	xy = matrix(0,nrow=nrow(phy$edge),ncol=4);
	bar = matrix(0,nrow=nrow(phy$edge),ncol=4);
	for (node in tree.traverse(phy,phy$Nnode+2,'postorder'))
	{
		bl = phy$edge.length[phy$edge[,2] == node]/tH;
		if (node <= phy$Nnode + 1)
		{
			xy[phy$edge[,2]==node,3:4] = c(1,(node-1)*dy);
			xy[phy$edge[,2]==node,1:2] = c(1-bl,(node-1)*dy);
			bar[phy$edge[,2]==node,] = 0;
		}
		else
		{
			isChild = phy$edge[,2] %in% phy$edge[phy$edge[,1] == node,2];
			if (node == phy$Nnode+2)
			{
				root = c(xy[isChild,1][1],xy[isChild,4][1],xy[isChild,1][1],xy[isChild,4][2]);
				next;
			}
			xy[phy$edge[,2] == node,4] = sum(xy[isChild,4])/2;
			xy[phy$edge[,2] == node,2] = sum(xy[isChild,4])/2;	
			xy[phy$edge[,2] == node,3] = xy[isChild,1][1];
			xy[phy$edge[,2] == node,1] = xy[isChild,1][1] - bl;
			
			bar[phy$edge[,2] == node,c(2,4)] = xy[isChild,4];
			bar[phy$edge[,2] == node,c(1,3)] = xy[isChild,1];
		}
	}
	xy = rbind(c(xy[1,1],sum(xy[1:2,4])/2,xy[1,1],sum(xy[1:2,4])/2),xy);
	bar = rbind(root,bar);
	rownames(xy) = c(phy$edge[1,1],phy$edge[,2]); colnames(xy) = c('x0','y0','x1','y1');
	return(list (segs = xy, arcs = bar ) );	
}

histRates = function(rates,pal,NCOLORS)
{
	opar = par(no.readonly = TRUE);
	fx = density(rates);
	dpal = c('BrBG','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral');
	if(length(pal) == 3)
	{
		rate.colors = colorRampPalette(pal,space='Lab')(NCOLORS);	
	}
	else if(pal %in% dpal)
	{
		rate.colors = colorRampPalette(rev(brewer.pal(3,pal)),space='Lab')(NCOLORS);
	}
	else if(pal == 'temperature')
	{
		rate.colors = rich.colors(NCOLORS);	
	}
	qx = quantile(rates,seq(0,1,length.out = NCOLORS+1));
	cat("Click once where you want the lower left corner of the figure\n");
	cxy = locator(n=1);
	xc = (grconvertX(cxy$x,to='ndc'));
	yc = (grconvertY(cxy$y,to='ndc'));
	ofs = min(1 - xc, 1 - yc);
	fig = c(xc,xc+ofs,yc,yc+ofs);
	par(fig = fig, new=TRUE, xpd=TRUE, mar=c(1.5,1.5,0.25,0.25));
	plot.new();
	plot.window(xlim=c(0,max(fx$x)),ylim=c(0,max(fx$y)));
	for(i in 1:length(fx$x))
	{
		index = which(qx > fx$x[i])[1];
		if(is.na(index)) break;
		if(index > 1) index = index - 1;
		bcol = rate.colors[index];
		segments(fx$x[i],fx$y[i],fx$x[i],0,lend=2,col=bcol);
	}
	axis(1,round(seq(0,max(fx$x),length.out=5),1),pos=0,cex.axis=0.75,tcl=NA,mgp=c(0,0.25,0));
	axis(2,round(seq(0,max(fx$y),length.out=3),0),las=1,pos=0,cex.axis=0.75,tcl=NA,mgp=c(0,0.25,0));
	mtext('Evolutionary Rate',1,line=1,cex=0.75);
	mtext('Density',2,line=1,cex=0.75);
	par(opar);
}


############################################
#	dtRates(ephy,tau)
#
#	A function to calculate approximations of
#	mean instantaneous speciation rates or
#	phenotypic rates along each branch. Still has 
#	some bugs, could be sped up, too.
#
#	Arguments: ephy = a bammdata object
#	           tau = fraction of tree height for approximation (e.g. 0.01)
#	           ism = index of posterior sample(s). Currently may be NULL or 
#	                 a vector of integer values.  if NULL the function will use all 
#	                 posterior samples, otherwise it will use only
#	                 the samples corresponding to the indices in ism,
#	                 e.g. 50, e.g. 50:100.
#
#	Returns: an ephy object with a list appended containing a vector of branch
#			 rates and the step size used for calculation.
dtRates = function(ephy, tau, ism = NULL)
{
	if(!'bamm-data' %in% class(ephy))
	{
		stop('Function requires a bammdata object');
	}
	
	ephy$eventBranchSegs = lapply(ephy$eventBranchSegs, function(x) x[order(x[,1]), ]); 
	
	phy = as.phylo.bammdata(ephy);
	phy = getStartStopTimes(phy);
	if(attributes(phy)$order != 'cladewise')
	{
		phy = reorder(phy,'cladewise');
	}
	tH = max(branching.times(phy));
	
	segmat = segMap(phy$edge[,2],phy$begin/tH,phy$end/tH,tau);
	segmat[,2] = segmat[,2] * tH;
	segmat[,3] = segmat[,3] * tH;
	
	tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
	
	if (storage.mode(segmat) != "double") stop('Exiting');
	if (storage.mode(tol) != "double") stop('Exiting');
	if (storage.mode(ephy) != "list") stop('Exiting');
	
	if (is.null(ism)) ism = as.integer(1:length(ephy$eventBranchSegs)) else ism = as.integer(ism);
	if (ism[length(ism)] > length(ephy$eventBranchSegs) )
	{
		warning("Sample index out of range");
		ism = as.integer(1:length(ephy$eventBranchSegs));
	}
	
	index = 1:nrow(segmat)
	rownames(segmat) = index;
	segmat = segmat[order(segmat[,1]),];
	dtrates = .Call("dtrates", ephy, segmat, tol, ism);
	names(dtrates) = rownames(segmat);
	dtrates = dtrates[as.character(index)];
	names(dtrates) = NULL;

	ephy$dtrates = list(tau = tau, rates = dtrates);
	return(ephy);
}

# dtRates = function(ephy,tau)
# {
	# if(!'bamm-data' %in% class(ephy))
	# {
		# stop('Function requires a bammdata object');
	# }
	# phy = as.phylo.bammdata(ephy);
	# phy = getStartStopTimes(phy);
	# if(attributes(phy)$order != 'cladewise')
	# {
		# phy = reorder(phy,'cladewise');
	# }
	# tH = max(branching.times(phy));
	
	# segs = segMap(phy$edge[,2],phy$begin/tH,phy$end/tH,tau);
	# rates = numeric(nrow(segs));
	
	# nsamples = length(ephy$eventBranchSegs);
	# for(k in 1:nsamples)
	# {
		# eventSegs = ephy$eventBranchSegs[[k]];
		# eventData = ephy$eventData[[k]];
		# for(j in 1:nrow(eventSegs))
		# {
			# node = eventSegs[j,1];
			# if(j < nrow(eventSegs)) nnode = eventSegs[j+1,1];
			# ev = eventSegs[j,4];
			# Start = eventData[eventData$index == ev,]$time;
			# lam1 = eventData[eventData$index == ev,]$lam1;			
			# lam2 = eventData[eventData$index == ev,]$lam2;
			
			# isGoodSeg = segs[,1] == node;
			# #isGoodStart = segs[,2]*tH >= eventSegs[j,2];
			# #isGoodEnd = segs[,3]*tH <= eventSegs[j,3];
			# isGoodStart = safeCompare(segs[,2]*tH,eventSegs[j,2],">=",1*10^-decimals(eventSegs[j,2]));
			# isGoodEnd = safeCompare(segs[,3]*tH,eventSegs[j,3],"<=",1*10^-decimals(eventSegs[j,3]));
			
			# if(sum(isGoodSeg & isGoodStart & isGoodEnd))
			# {
				# relStart = segs[isGoodSeg & isGoodStart & isGoodEnd, 2]*tH - Start;
				# relEnd = segs[isGoodSeg & isGoodStart & isGoodEnd, 3]*tH - Start;
				# rates[isGoodSeg & isGoodStart & isGoodEnd] = 
						# rates[isGoodSeg & isGoodStart & isGoodEnd] + 
							# branchMeanRateExponential(relStart,relEnd,lam1,lam2)/nsamples;		
			# }
			# if(node == nnode)
			# {
				# isGoodStart = segs[,2]*tH < eventSegs[j,3];
				# isGoodEnd = segs[,3]*tH > eventSegs[j,3];
				# if(sum(isGoodSeg & isGoodStart & isGoodEnd) == 1)
				# {		
					# relStart = segs[isGoodSeg & isGoodStart & isGoodEnd, 2]*tH - Start;
					# relEnd = eventSegs[j,3] - Start;
					# leftshift = timeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
					
					# relStart = 0;
					# relEnd = segs[isGoodSeg & isGoodStart & isGoodEnd,3]*tH - eventSegs[j,3];
					# lam1 = eventData[eventData$index==eventSegs[j+1,4],]$lam1;
					# lam2 = eventData[eventData$index==eventSegs[j+1,4],]$lam2;
					# rightshift = timeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
					
					# shift = (leftshift+rightshift)/(segs[isGoodSeg & isGoodStart & isGoodEnd,3]*tH - segs[isGoodSeg & isGoodStart & isGoodEnd,2]*tH);
					
					# rates[isGoodSeg & isGoodStart & isGoodEnd] = rates[isGoodSeg & isGoodStart & isGoodEnd] + shift/nsamples;	
				# }
			# }	
		# }
	# }
	# ephy$dtrates = list(tau=tau,rates=rates);
	# return(ephy);	
# }

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

safeCompare = function(vec,val,FUN,tol=1e-4)
{
	if(FUN == ">=")
	{
		ret = rep(FALSE,length(vec));
		ret[(vec-val) >= 0] = TRUE;
		ret[abs(vec-val) <= tol] = TRUE;
		return(ret);
	}
	if(FUN == "<=")
	{
		ret = rep(FALSE,length(vec));
		ret[(vec-val) <= 0] = TRUE;
		ret[(vec-val) <= tol] = TRUE;
		return(ret);
	}
}

decimals = function(x)
{
	if(x%%1 != 0)
	{
		return(nchar(strsplit(as.character(x),".",fixed=TRUE)[[1]][[2]]));
	}	
	else
	{
		return(10);
	}
}


##################################
#	Internal functions called by plot.dtrates
#
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
	if(step) polygon(c(xv,rev(xv)),c(yv,rev(yv)),border=border,...);
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

matrify = function(x,tau)
{
	bn = sqrt((x[3]-x[1])^2 + (x[4]-x[2])^2);
	len = bn/tau; if (len %% 1 == 0) len = len + 1;
	
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

colorMap = function(x, pal, NCOLORS)
{
	dpal = c('BrBG','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral');
	colset = numeric(length(x));
	if(length(pal) == 3)
	{
		colpalette = colorRampPalette(pal,space='Lab')(NCOLORS);	
	}
	else if(pal %in% dpal)
	{
		colpalette = colorRampPalette(rev(brewer.pal(3,pal)),space='Lab')(NCOLORS);
	}
	else if(pal == 'temperature')
	{
		colpalette = rich.colors(NCOLORS);	
	}
	
	bks = quantile(x, seq(0,1,length.out=(NCOLORS+1)));
	for(i in 2:length(bks))
	{
		if(i == 2)
		{
			colset[x < bks[2]] = colpalette[1];	
		}
		else if(i == length(bks))
		{
			colset[x >= bks[length(bks)-1]] = colpalette[length(bks)-1];
		}
		else
		{
			colset[x >= bks[i-1] & x < bks[i]] = colpalette[i-1];
		}
	}
	return(colset);
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
	attributes(newphylo)$order <- attributes(ephy)$order;
	return(newphylo);
}
