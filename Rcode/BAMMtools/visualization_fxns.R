





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







