getShiftNodesFromIndex <-
function(ephy, index) {
	
	if (index > length(ephy$eventData)) {
		cat("Error. Attempting to access non-existent element from 'bammdata' object\n");
		cat("You have << ", length(ephy$eventData), " >>> samples in your 'bammdata' object\n");
		cat("Value of index must be no greater than the number of samples\n\n");
		stop();
	}
	root <- length(ephy$tip.label) + 1;
	nodes <- ephy$eventData[[index]]$node;
	nodes <- nodes[nodes != root];

	return(nodes);
}
