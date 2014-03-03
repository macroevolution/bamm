


# Drop all nodes from event data with marginal probs < 0.05
# test is unique
#  if so, add to list
# Returns:
#	$marg.probs = marginal probs for nodes
#	$shifts = unique shift sets
#	$samplesets = list of sample indices that reduce to each of the unique shift sets
#	$frequency = vector of frequencies of each shift configuration
#	$threshold = marginal prob threshold for nodes (ignores all nodes less than this value)
#	
#	Results are sorted by frequency. 
#	$frequency[1] gives the most common shift config sampled
#	$shifts[[1]] gives the corresponding node indices for that configuration
#	$samplesets[[1]] gives the indices of samples with this configuration



distinctShiftConfigurations <- function(ephy, threshold) {
	mm <- marginalShiftProbsTree(ephy);
	
	if (class(threshold) == 'branchprior'){
		goodnodes <- mm$edge[,2][mm$edge.length >= threshold$edge.length];
	}else if (class(threshold) == 'numeric'){
		goodnodes <- mm$edge[,2][mm$edge.length >= threshold];		
	}else{
		stop('Threshold is of wrong class\n');
	}
	
	

	xlist <- list();
	for (i in 1:length(ephy$eventData)) {
		xlist[[i]] <- intersect(goodnodes, ephy$eventData[[i]]$node);
	}

	ulist <- list();
	treesets <- list();
	
	ulist[[1]] <- xlist[[1]];
	treesets[[1]] <- 1;
	
	for (i in 2:length(xlist)) {
		lx <- length(ulist);
		#cat(lx, '\n')
		for (k in 1:lx) {
			if (areShiftSetsEqual(ulist[[k]], xlist[[i]])){
				treesets[[k]] <- c(treesets[[k]], i);
				break;	
			} else {
				if (k == length(ulist)){
					xlen <- length(ulist);
					ulist[[xlen + 1]] <- xlist[[i]];
					treesets[[xlen + 1]] <- i;
				}
			}
		}
	}
	
	freqs <- unlist(lapply(treesets, length));
	freqs <- freqs / sum(freqs);
	
	ord <- order(freqs, decreasing=TRUE);
	
	obj <- list();
	obj$marg.probs <- mm$edge.length;  
	names(obj$marg.probs) <- mm$edge[,2]; 
	obj$shifts <- ulist[ord]; 
	obj$samplesets <- treesets[ord];
	obj$frequency <- freqs[ord];
	obj$threshold <- threshold;

	class(obj) <- 'bammshifts';
	
	return(obj);
}






