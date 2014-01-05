dtRates <-
function(ephy, tau, ism = NULL) {
	if (!'bammdata' %in% class(ephy)) {
		stop('Object ephy must be of class bammdata');
	}
	
	ephy$eventBranchSegs = lapply(ephy$eventBranchSegs, function(x) x[order(x[,1]), ]); 
	
	phy = as.phylo.bammdata(ephy);
	phy = getStartStopTimes(phy);
	if (attributes(phy)$order != 'cladewise') {
		phy = reorder(phy,'cladewise');
	}
	tH = max(branching.times(phy));
	
	segmat = segMap(phy$edge[,2],phy$begin/tH,phy$end/tH,tau);
	segmat[,2] = segmat[,2] * tH;
	segmat[,3] = segmat[,3] * tH;
	
	tol = max(1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]),1*10^-decimals(ephy$eventBranchSegs[[1]][1,3]));
	
	if (storage.mode(segmat) != "double") stop('Exiting');
	if (storage.mode(tol) != "double") stop('Exiting');
	if (storage.mode(ephy) != "list") stop('Exiting');
	
	if (is.null(ism)) ism = as.integer(1:length(ephy$eventBranchSegs)) else ism = as.integer(ism);
	if (ism[length(ism)] > length(ephy$eventBranchSegs)) {
		warning("Sample index out of range");
		ism = as.integer(1:length(ephy$eventBranchSegs));
	}
	
	index = 1:nrow(segmat)
	rownames(segmat) = index;
	segmat = segmat[order(segmat[,1]),];
	dtrates = .Call("dtrates", ephy, segmat, tol, ism);
	names(dtrates) = rownames(segmat);
	dtrates = dtrates[as.character(index)];
	names(dtrates) = NULL;

	ephy$dtrates = list(tau = tau, rates = dtrates);
	return(ephy);
}
