getMeanTipRateStates <-
function(ephy, use.names = FALSE) {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}

	if (!use.names) {
		return(ephy$meanTipLambda - ephy$meanTipMu);
	} else {
		ret <- ephy$meanTipLambda - ephy$meanTipMu;
		names(ret) <- ephy$tip.label;
		return(ret);
	}
}
