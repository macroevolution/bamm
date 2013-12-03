getPathToRoot <-
function(phy, node){
	
	root <- length(phy$tip.label) + 1;
	nset <- node;
	while (node != root){
		node <- phy$edge[,1][phy$edge[,2] == node];
		nset <- c(nset, node);
	}
	return(nset);
}
