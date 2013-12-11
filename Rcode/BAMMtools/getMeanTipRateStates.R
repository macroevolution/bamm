#############################################################
#
#	getMeanTipRateStates(....)
#
# returns vector of mean tip states (either net diversification rate or brownian motion
#	rate parameter depending on the 'type' of ephy).
#	use.names=TRUE returns a named vector

getMeanTipRateStates <- function(ephy, use.names = FALSE){
	
	if ('bammdata' != class(ephy)){
		stop("Object ephy must be of class bammdata\n");
	}

	if (use.names == FALSE){
		return(ephy$meanTipLambda - ephy$meanTipMu);
	}
	else{
		ret <- ephy$meanTipLambda - ephy$meanTipMu;
		names(ret) <- ephy$tip.label;
		return(ret);
	}
}
