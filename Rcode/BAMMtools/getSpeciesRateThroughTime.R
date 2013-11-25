##################################################################
#
#	getSpeciesRateThroughTime <- function(...)
#
#	Start time: how many time units before present to include
#	nbreaks: how many time points to compute the rate
#	ndr: boolean, should net diversification rate be returned (automatically adjusted if ephy$type = traits)
#	species: name of a single species
#	returnAll: boolean, whether or not to return all values from posterior, rather than averaging.
#
#
getSpeciesRateThroughTime <- function(ephy, start.time, nbreaks, ndr, species,returnAll = F){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	if (!species %in% ephy$tip.label | length(species) > 1){
		stop("Species must be a single species whose name matches a tip in the phylogeny.")
	}
	
	if (ephy$type == 'traits'){
		ndr <- FALSE;
	}
		
	tend <- max(branching.times(ephy))*0.999;
	tstart <- max(branching.times(ephy)) - tend;
	tseq <- seq(tstart, tend, length.out=nbreaks);
 
	res <- numeric(nbreaks);
	
	index <- which(ephy$tip.label == species);
	
	path <- getPathToRoot(ephy, node=index);
	
	resMat <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=nbreaks);
	
	for (k in 1:length(ephy$eventBranchSegs)){

		ed <- ephy$eventData[[k]];
		
		resVec <- vector(mode='numeric', length=nbreaks);
		for (z in 1:length(tseq)){
			isGoodNode <- ephy$eventBranchSegs[[k]][,1] %in% path;
			isGoodStart <- ephy$eventBranchSegs[[k]][,2] <= tseq[z];
			isGoodEnd <- ephy$eventBranchSegs[[k]][,3] >= tseq[z];
				
			ev <- ephy$eventBranchSegs[[k]][,4][isGoodNode & isGoodStart & isGoodEnd];
				
			if (ndr){
				lam <- exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
				mu <- ed$mu1[ev];
				resVec[z] <- lam - mu;	
			}else{
				resVec[z] <- exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
			}	
		}		
	resMat[k,] <- resVec;
	}
	if (returnAll == TRUE){
		return(resMat);
	}
	if (returnAll == FALSE){
		return(colMeans(resMat));
	}
}
