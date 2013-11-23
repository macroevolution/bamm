
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
#	Returns: 		- a dataframe with the event data for the "best" shift tree
#					- A vector of optimality scores for all other samples 
#						in posterior from the bamm-data object.
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
	
	obj <- list();
	obj$bestShiftConfig <- ephy$eventData[[best]];
	obj$scores <- probvec;
	obj$optimalityType = maximize;
	obj$best_index <- best;
	return(obj);
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
	attributes(newphylo)$order <- 'cladewise';
	return(newphylo);
}











