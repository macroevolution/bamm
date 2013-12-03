############################################
#	Internal function called by plot.dtrates(...)
#
#
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