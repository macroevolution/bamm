

credibleShiftSet <- function(ephy, set.limit=0.95, threshold=0.01){
	
	dsc <- distinctShiftConfigurations(ephy, threshold = threshold);
	cfreq <- cumsum(dsc$frequency);
	cut <- min(which(cfreq >= set.limit));
	nodeset <- dsc$marg.probs[dsc$marg.probs >= threshold];
 	
 	
	shiftnodes <- dsc$shifts[1:cut];
	indices <- dsc$samplesets[1:cut];
	frequency <- dsc$frequency[1:cut];
	cumulative <- cumsum(dsc$frequency)[1:cut];

	obj <- list(shiftnodes=shiftnodes, indices=indices, frequency=frequency, cumulative=cumulative, number.distinct=cut, set.limit=set.limit);
	class(obj) <- 'credibleshiftset'
	return(obj);
}






