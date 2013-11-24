#############################################################
#
#	getBySpeciesRateMatrix <- function(...)
#
#	Start time: how many time units before present to include
#	nbreaks: how many time points to compute the rate
#	ndr: boolean, should net diversification rate be returned, ignored if ephy is of type "traits".
#	node: if NULL, all species returned, if defined, species descending from that node are returned. 
#
#	This version computes for all species
#	


getBySpeciesRateMatrix <- function(ephy, start.time, nbreaks, ndr=TRUE, node=NULL){
	
	if (!is.null(node)){
		if (node >= Ntip(ephy)){
			spset <- tips(ephy,node)
		}
	} else{
		spset <- ephy$tip.label;
	}
	
	if (ephy$type == 'traits'){
		ndr <- FALSE
	}
	
	mm <- t(sapply(spset,function(x) getSpeciesRateThroughTime(ephy,start.time, nbreaks, ndr, species=x)))
	
	return(mm);
}
