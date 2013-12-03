recursivelySetNodeStates <-
function(phySV, node, state) {
	
	phySV$statevec[phySV$edge[,2] == node] <- state;
	if (sum(phySV$edge[,1] == node) > 0){
		# node is internal
		dset <- phySV$edge[,2][phySV$edge[,1] == node];
		phySV <- recursivelySetNodeStates(phySV, dset[1], state);
		phySV <- recursivelySetNodeStates(phySV, dset[2], state);
	}
	
	return(phySV);
}
