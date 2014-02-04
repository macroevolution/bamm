## show.all.nodes: should plot core and "floater" nodes with different colors 
#	that can be specified by the user.
#    argument to specify layout , e.g., c(2,2) would specify a 2x2 panel?
#	 maybe print the number of shift configurations above the threshold that were 
#	 not plotted along with their frequency?
#	
#	add.freq.text = argument to add text to each plot indicating the frequency of the 
#		shift configuration

plotDistinctShiftConfigurations <- function(sc, ephy, plotmax=9, method='phylogram', pal = 'temperature', cex=2){
	if (class(sc) != 'bammshifts'){
		stop('arg sc must be of class "bammshifts"');
	}
	if (class(ephy) != 'bammdata'){
		stop('arg ephy must be of class "bammdata"');
	}
	
	plot.new();
	if (length(sc$frequency) <= 2){
		par(mfrow=c(1,2));
	}else if (length(sc$frequency) <= 4){
		par(mfrow = c(2,2));
	}else if (length(sc$frequency) <= 6){
		par(mfrow=c(2,3));
	}else{
		par(mfrow=c(3,3));
	}
	
	mm <- min(c(length(sc$frequency), plotmax));
	
	for (i in 1:mm){
		tmp <- subsetEventData(ephy, index=sc$sampleset[[i]][1]);
		plot.bammdata(tmp, method, pal=pal);
		addBAMMshifts(tmp, method, index=1, cex=cex, multi=T);	
	}
	
	
}

