as.phylo.bammdata <- function(ephy) {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}		
	
	newphylo <- list();
	newphylo$edge <- ephy$edge;
	newphylo$Nnode <- ephy$Nnode;
	newphylo$tip.label <- ephy$tip.label;
	newphylo$edge.length <- ephy$edge.length;
	class(newphylo) <- 'phylo';
	attributes(newphylo)$order = attributes(ephy)$order;
	if (attributes(newphylo)$order != "cladewise") {
		newphylo = reorder(newphylo);
	}
	return(newphylo);
}
