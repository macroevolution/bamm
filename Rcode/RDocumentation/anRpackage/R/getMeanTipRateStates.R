getMeanTipRateStates <-
function(ephy, use.names = FALSE){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
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
