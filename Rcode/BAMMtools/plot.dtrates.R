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
#	           vtheta = specifies the angle in degrees separating the first and last
#	                    tips to prevent over plotting. Ignored if method = 'phylogram'.
#	           rbf = specifies the length of the root branch as a fraction of the 
#	                 total tree height. rbf > 0 will cause an arc to connect the immediate
#	                 descendants of the root branch. Ignored if method = 'phylogram'. 
#	           show = TRUE or FALSE. If TRUE the tree will plot.
#	           labels = TRUE or FALSE. If TRUE the tip labels will plot.
#	           multi = TRUE or FALSE for multipanel plotting.
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
#	           ... = further arguments passed to par to control plotting, e.g. mar.
#
plot.dtrates = function(ephy, method='phylogram', tau=0.01, index=NULL, vtheta=5, rbf=0.001, show=TRUE, labels=FALSE, multi = FALSE, hrates=FALSE, lwd=3, cex=1, ncolors=64, pal='Spectral', ...)
{
	if ('bamm-data' %in% class(ephy)) phy = as.phylo.bammdata(ephy) else stop('Trying to work with a non-bammdata object');
	if (!is.binary.tree(phy)) stop('Function requires fully bifurcating tree.');
	if (any(phy$edge.length == 0)) warning('Tree contains zero length branches. Rates for these will be NA and coerced to zero');
	if (!('dtrates' %in% names(ephy))) ephy = dtRates(ephy, tau, index);
	if(sum(is.na(ephy$dtrates$rates))) 
	{
		warning(sprintf("Found %s NA values. Coercing to zero.",sum(is.na(ephy$dtrates$rates))));
		ephy$dtrates$rates[is.na(ephy$dtrates$rates)] = 0;
	}
	tH = max(branching.times(phy));
	phy$begin = ephy$begin; phy$end = ephy$end;
	if (method == 'polar')
	{
		ret = setPolarTreeCoords(phy,vtheta,rbf);
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
	#p = apply(p,1,mkdtsegs,tau);
	p = apply(p,1,mkdtsegs,tau,phy,tH);
	p = do.call(rbind, p);
	x0 = c(x0[1],p[,1]);x1=c(x1[1],p[,2]);y0=c(y0[1],p[,3]);y1=c(y1[1],p[,4]);
	offset = table(p[,5])[as.character(unique(p[,5]))];
	arc.color = c(edge.color[1],edge.color[match(unique(p[,5]),p[,5])+offset]);
	edge.color = c(edge.color[1],edge.color);
	
	if (show)
	{
		if(length(list(...)))
		{
			op = par(no.readonly=TRUE);
			par(...);
		}
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
			segments(ret$arcs[!isTip,1],ret$arcs[!isTip,2],ret$arcs[!isTip,3],ret$arcs[!isTip,4],col=arc.color[!isTip],lwd=lwd,lend=2);
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
		assign("last_plot.phylo", list(Ntip = phy$Nnode+1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,3], yy = ret$segs[index,4], pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
	}
	else if (method == 'polar')
	{
		assign("last_plot.phylo", list(Ntip = phy$Nnode+1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,3], yy = ret$segs[index,4], theta = ret$segs[index,5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
	}
	if(!multi)
		if(length(list(...))) par(op);
	invisible(ret$segs[-1,]);
}
