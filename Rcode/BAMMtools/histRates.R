############################################
#	Internal function called by plot.dtrates(...)
#
#
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