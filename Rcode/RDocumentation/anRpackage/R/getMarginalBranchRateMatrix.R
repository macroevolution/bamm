getMarginalBranchRateMatrix <-
function(ephy, verbose=FALSE){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	lammat <- matrix(0, ncol=length(ephy$eventBranchSegs), nrow=nrow(ephy$edge));
	mumat <- matrix(0, ncol=length(ephy$eventBranchSegs), nrow=nrow(ephy$edge));
	
	for (i in 1:length(ephy$eventBranchSegs)){
		if (verbose){
			cat('Processing sample ', i, '\n');
		}
		esegs <- ephy$eventBranchSegs[[i]];
		events <- ephy$eventData[[i]];
		events <- events[order(events$index), ];			
		
		# relative start time of each seg, to event:
		relsegmentstart <- esegs[,2] - ephy$eventData[[i]]$time[esegs[,4]];
		relsegmentend <- esegs[,3] - ephy$eventData[[i]]$time[esegs[,4]];
		lam1 <- ephy$eventData[[i]]$lam1[esegs[,4]];
		lam2 <-  ephy$eventData[[i]]$lam2[esegs[,4]];
		mu1 <-  ephy$eventData[[i]]$mu1[esegs[,4]];
		mu2 <-  ephy$eventData[[i]]$mu2[esegs[,4]];
 
		lamint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, lam1, lam2);
		muint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, mu1, mu2);
		seglengths <- esegs[,3] - esegs[,2];	
				
		for (k in 1:nrow(ephy$edge)){
			isRightBranch <- esegs[,1] == ephy$edge[k,2];
			lammat[k, i] <- sum(lamint[isRightBranch]) / sum(seglengths[isRightBranch]);
			mumat[k, i] <- sum(muint[isRightBranch]) / sum(seglengths[isRightBranch]);
			
		}
	
	}
	
	return(list(lambda_branch_matrix = lammat, mu_branch_matrix = mumat));
	
}
