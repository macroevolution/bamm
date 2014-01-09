#############################################################
#
#	getRecursiveSequence(....)
#
#	Private function, called by getEventDataDiversification

getRecursiveSequence <- function(phy){
		
	root.node <- length(phy$tip.label) + 1;
	
	phy$downseq <- root.node
	phy$lastvisit <- numeric(length(unique(phy$edge[,2])));
 
 	phy <- getSequenceForwardTraversal(phy, root.node);

	return(phy);
}
