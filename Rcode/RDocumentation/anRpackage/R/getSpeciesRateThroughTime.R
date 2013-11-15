getSpeciesRateThroughTime <-
function(ephy, start.time, nbreaks=10, ndr=TRUE, species){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
		
	tend <- max(branching.times(ephy))*0.999;
	tstart <- start.time #tend - start.time;
	tseq <- seq(tstart, tend, length.out=nbreaks);
 
	res <- numeric(nbreaks);
	
	index <- which(ephy$tip.label == species);
	
 
	path <- getPathToRoot(ephy, node=index);
		
	for (k in 1:length(ephy$eventBranchSegs)){
			
		ed <- ephy$eventData[[k]];
			
		for (z in 1:length(tseq)){
			isGoodNode <- ephy$eventBranchSegs[[k]][,1] %in% path;
			isGoodStart <- ephy$eventBranchSegs[[k]][,2] <= tseq[z];
			isGoodEnd <- ephy$eventBranchSegs[[k]][,3] >= tseq[z];
				
			ev <- ephy$eventBranchSegs[[k]][,4][isGoodNode & isGoodStart & isGoodEnd];
				
			if (ndr){
				lam <- exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
				mu <- ed$mu1[ev];
				res[z] <- res[z] + (lam - mu);	
			}else{
				res[z] <- res[z] + exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
			}
				
				
		}
			
	}
 	res <- res / length(ephy$eventBranchSegs);
 
 	return(res);
	
}
