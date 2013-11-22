#############################################################
#
#	getBySpeciesRateMatrix <- function(...)
#
#	Start time: how many time units before present to include
#	nbreaks: how many time points to compute the rate
#
#	This version computes for all species
#	
#	node argument will just compute for the subtree descended from "node"

getBySpeciesRateMatrix <- function(ephy, start.time, nbreaks, ndr=TRUE, node){
	
	
	spset <- ephy$tip.label;
	
	mm <- matrix(NA, nrow=length(spset), ncol=nbreaks);
	for (i in 1:nrow(mm)){
		#cat(spset[i], '\n')
		mm[i,] <- getSpeciesRateThroughTime(ephy, start.time, nbreaks, ndr, species=spset[i]);
	}
	
	rownames(mm) <- spset;
	return(mm);
}
