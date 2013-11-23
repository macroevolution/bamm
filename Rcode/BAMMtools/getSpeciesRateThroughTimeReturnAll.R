##################################################################
#
#	getSpeciesRateThroughTimeReturnAll <- function(...)
#
#	Primarily intended as an internal function to support plotSpeciesRatesThroughTime()
#	Identical to getSpeciesRateThroughTime() except that this function returns all values from posterior, rather than averaging.
#
#
getSpeciesRateThroughTimeReturnAll <- function(ephy, start.time, nbreaks, ndr, species){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
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
	return(resMat);
}
