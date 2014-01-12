#############################################################
#
#	getMeanTipRateStates(....)
#
# returns vector of mean tip states (either speciation, extinction, net diversification or trait rate
#	use.names=TRUE returns a named vector (sp names)

getMeanTipRateStates <- function(ephy, ratetype, use.names = TRUE) {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}
	
	if (!ratetype %in% c('speciation','extinction','netdiv','trait')) {
		stop("Ratetype must be either speciation, extinction, netdiv or trait.");
	}

	if ((ratetype == 'trait' & ephy$type == 'diversification') | (ratetype != 'trait' & ephy$type == 'trait')) {
		stop("Ratetype does not correspond to bamm analysis type.");
	}
	
	if (ratetype == 'speciation' | ratetype == 'trait') {
		res <- ephy$meanTipLambda;
	}
	
	if (ratetype == 'extinction') {
		res <- ephy$meanTipMu;
	}
	
	if (ratetype == 'netdiv') {
		res <- ephy$meanTipLambda - ephy$meanTipMu;
	}
	
	if (!use.names) {
		return(res);
	} else {
		names(res) <- ephy$tip.label;
		return(res);
	}
}
