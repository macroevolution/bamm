
#############################################################
#
#	marginalShiftProbsTree(....)
#
#	Args: ephy	=	object of class 'bamm-data'
#	
#	Returns:		a phylogenetic tree, but where each 
#	             	branch length (edge length) is equal to the
#					marginal probability of shift occuring 
#					on that particular branch. 
#							

marginalShiftProbsTree <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	shiftvec <- numeric(length(ephy$edge.length));
 	rootnode <- length(ephy$tip.label) + 1;
 
 
	for (i in 1:length(ephy$eventData)){
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		shiftvec[hasShift] <- shiftvec[hasShift] + rep(1, sum(hasShift));
	}
	
	shiftvec <- shiftvec / length(ephy$eventData);	
	
	newphy <- as.phylo.bammdata(ephy);
	newphy$edge.length <- shiftvec;
	return(newphy);
}


#############################################################
#
#	cumulativeShiftProbsTree(....)
#
#	Args: ephy	=	object of class 'bamm-data'
#	
#	Returns:		a phylogenetic tree, but where each 
#	             	branch length (edge length) is equal to the
#					cumulative probability of a shift somewhere 
#					between the focal branch and the root of the 
#					tree. The branch length itself does not tell 
#					you where the shifts occur, but they tell 
#					you which clades/lineages have diversification
#					dynamics that are decoupled from the root of the tree 
#							

cumulativeShiftProbsTree <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}	
			
	shiftvec <- numeric(length(ephy$edge.length));
 	rootnode <- length(ephy$tip.label) + 1;
 
 
	for (i in 1:length(ephy$eventData)){
		
		snodes <- unique(ephy$eventBranchSegs[[i]][,1][ephy$eventBranchSegs[[i]][,4] != 1]);
		hasShift <- ephy$edge[,2] %in% snodes;
		shiftvec[hasShift] <- shiftvec[hasShift] + rep(1, sum(hasShift));
	}	
	
	shiftvec <- shiftvec / length(ephy$eventData);		
	newphy <- as.phylo.bammdata(ephy);
	newphy$edge.length <- shiftvec;
	return(newphy);	
	
	
}

#############################################################
#
#	maximumShiftCredibilityTree(....)
#
#	Args: ephy	=	object of class 'bamm-data'
#	
#
#	maximize		=	'sum', = sum of branch probabilities for each tree
#						'product', = product of branch probabilities for each tree
#	
#	Returns: 		- bestconfigs: a list of length equal the number of 
#						unique shift configurations in the maximum shift
#						credibility set. Each element is a vector of sample
#						indices from the BAMM-data object with identical
#						shift configurations.
#					  
#					- A vector of optimality scores for all other samples 
#						in posterior from the bamm-data object.
#
#					- sampleindex: a representative index for samples from each 
#						set of unique shift configurations. The length of this vector
#						is equal to the length of the bestconfigs list. If this vector was
#						sampleindex = c(2, 20, 50), this would mean that there are 3 distinct
#						sets of shift configurations with equal credibility under the optimality 
#						criterion. More commonly, a single shift configuration will be dominant, and 
#						although the length of bestconfigs[[1]] may be greater than 1, the sampleindex
#						vector will contain a single representative event from that set.
#
#	See example file.
#   This is analogous to the maximum clade credibility tree from a 
#		Bayesian phylogenetic analysis.

maximumShiftCredibilityTree <- function(ephy, maximize = 'product'){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}			
	
	probvec <- numeric(length(ephy$eventData));
	
	mtree <- marginalShiftProbsTree(ephy);
	px <- mtree$edge.length;
	

	for (i in 1:length(ephy$eventData)){
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		branchprobs <- (hasShift)*px  + (!hasShift)*(1 - px) ;
		if (maximize == 'product'){
			probvec[i] <- sum(log(branchprobs));
		}else if (maximize == 'sum'){
			probvec[i] <- sum(branchprobs);
		}else{
			stop("Unsupported optimize criterion in maximumShiftCredibilityTree");
		}
	}
	
	best <- which(probvec == max(probvec));
	
	# Now test for multiple trees with same log-prob:
	bestconfigs <- list();
		
	index <- 0;	
	while (length(best) > 0){
		index <- index + 1;	
		lv <- logical(length = length(best));
		for (i in 1:length(best)){
			lv[i] <- areEventConfigurationsIdentical(ephy, best[1], best[i]);
		}
		bestconfigs[[index]] <- best[lv];
		best <- best[!lv];
	}
	
	
	sampleindex <- numeric(length(bestconfigs));
	for (i in 1:length(bestconfigs)){
		sampleindex[i] <- bestconfigs[[i]][1];
	}
	
	obj <- list();
	obj$bestconfigs <- bestconfigs;
	obj$scores <- probvec;
	obj$optimalityType = maximize;
 	obj$sampleindex <- sampleindex;
	return(obj);
}

areEventConfigurationsIdentical <- function(ephy, index1, index2){
	
	nodeset <- c(ephy$eventData[[index1]]$node, ephy$eventData[[index2]]$node);
	diffs <- sum(table(nodeset) == 1);
	return(diffs == 0);	
}

#############################################################
#
#	getShiftNodesFromIndex (....)
#
#	Args: ephy	=	object of class 'bamm-data'
#		  index =   the index of the sample you wish to view, e.g.,
#					 if index = 5, this will give you the nodes subtending
#                    all branches with rate shifts for the 5th sample
#					 from the posterior in your BAMM-data object.								 
#
#	
#	Returns: 		- a vector of the nodes where shifts occurred, excluding the root.
#					Note: if NO shifts occured, this will return a 
#							numeric object of length zero
# 

getShiftNodesFromIndex <- function(ephy, index){
	
	if (index > length(ephy$eventData)){
		cat("Error. Attempting to access non-existent element from BAMM-data object\n");
		cat("You have << ", length(ephy$eventData), " >>> samples in your BAMM-data object\n");
		cat("Value of index must be no greater than the number of samples\n\n");
		stop();
	}
	root <- length(ephy$tip.label) + 1;
	nodes <- ephy$eventData[[index]]$node;
	nodes <- nodes[nodes != root];

	return(nodes);
}


as.phylo.bammdata <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}		
	
	newphylo <- list();
	newphylo$edge <- ephy$edge;
	newphylo$Nnode <- ephy$Nnode;
	newphylo$tip.label <- ephy$tip.label;
	newphylo$edge.length <- ephy$edge.length;
	class(newphylo) <- 'phylo';
	return(newphylo);
}











