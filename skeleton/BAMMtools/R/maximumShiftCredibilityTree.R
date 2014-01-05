maximumShiftCredibilityTree <-
function(ephy, maximize = 'product') {

	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}			
	
	probvec <- numeric(length(ephy$eventData));
	
	mtree <- marginalShiftProbsTree(ephy);
	px <- mtree$edge.length;
	
	for (i in 1:length(ephy$eventData)) {
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		branchprobs <- (hasShift)*px  + (!hasShift)*(1 - px) ;
		if (maximize == 'product') {
			probvec[i] <- sum(log(branchprobs));
		} else if (maximize == 'sum') {
			probvec[i] <- sum(branchprobs);
		} else {
			stop("Unsupported optimize criterion in maximumShiftCredibilityTree");
		}
	}
	
	best <- which(probvec == max(probvec));
	
	# Now test for multiple trees with same log-prob:
	bestconfigs <- list();
		
	index <- 0;	
	while (length(best) > 0) {
		index <- index + 1;	
		lv <- logical(length = length(best));
		for (i in 1:length(best)) {
			lv[i] <- areEventConfigurationsIdentical(ephy, best[1], best[i]);
		}
		bestconfigs[[index]] <- best[lv];
		best <- best[!lv];
	}
	
	sampleindex <- numeric(length(bestconfigs));
	for (i in 1:length(bestconfigs)) {
		sampleindex[i] <- bestconfigs[[i]][1];
	}
	
	obj <- list();
	obj$bestconfigs <- bestconfigs;
	obj$scores <- probvec;
	obj$optimalityType = maximize;
 	obj$sampleindex <- sampleindex;
	return(obj);
}
