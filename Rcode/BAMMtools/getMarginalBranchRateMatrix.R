#############################################################
#
#	getMarginalBranchRateMatrix(....)
#
#	get matrix of marginal rates on each branch for each sample from posterior
# 	
#	This function can handle either a bamm-data object or a multiphylo object (i.e., list of trees)

getMarginalBranchRateMatrix <- function(obj, verbose=FALSE){
	
	if (!'bamm-data' %in% class(obj) & class(obj[[1]]) != 'phylo'){
		stop("Object must be of class bamm-data or a list of trees.\n");
	}
	
	if ('bamm-data' %in% class(obj)){
		lammat <- matrix(0, ncol=length(obj$eventBranchSegs), nrow=nrow(obj$edge));
		mumat <- matrix(0, ncol=length(obj$eventBranchSegs), nrow=nrow(obj$edge));
		
		for (i in 1:length(obj$eventBranchSegs)){
			if (verbose){
				cat('Processing sample ', i, '\n');
			}
			esegs <- obj$eventBranchSegs[[i]];
			events <- obj$eventData[[i]];
			events <- events[order(events$index), ];			
			
			# relative start time of each seg, to event:
			relsegmentstart <- esegs[,2] - obj$eventData[[i]]$time[esegs[,4]];
			relsegmentend <- esegs[,3] - obj$eventData[[i]]$time[esegs[,4]];
			lam1 <- obj$eventData[[i]]$lam1[esegs[,4]];
			lam2 <-  obj$eventData[[i]]$lam2[esegs[,4]];
			mu1 <-  obj$eventData[[i]]$mu1[esegs[,4]];
			mu2 <-  obj$eventData[[i]]$mu2[esegs[,4]];
	 
			lamint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, lam1, lam2);
			muint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, mu1, mu2);
			seglengths <- esegs[,3] - esegs[,2];	
					
			for (k in 1:nrow(obj$edge)){
				isRightBranch <- esegs[,1] == obj$edge[k,2];
				lammat[k, i] <- sum(lamint[isRightBranch]) / sum(seglengths[isRightBranch]);
				mumat[k, i] <- sum(muint[isRightBranch]) / sum(seglengths[isRightBranch]);
				
			}
		
		}
		
		if (obj$type == 'diversification'){
			return(list(lambda_branch_matrix = lammat, mu_branch_matrix = mumat));
		}
		if (obj$type == 'traits'){
			return(list(beta_branch_matrix = lammat));
		}
	}
	else if (class(obj[[1]]) == 'phylo'){
		return(sapply(obj,with,edge.length));
	}
}








