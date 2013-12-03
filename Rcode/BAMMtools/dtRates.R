############################################
#	dtRates(ephy,tau)
#
#	A function to calculate approximations of
#	mean instantaneous speciation rates or
#	phenotypic rates along each branch.
#
#	Arguments: ephy = a bammdata object
#	           tau = fraction of tree height for approximation (e.g. 0.01).
#	                 This is the step size over which rates are calculated along
#	                 a branch, 0.01 corresponds to a step size of 1% of tree height.
#	           ism = index of posterior sample(s). Currently may be NULL or 
#	                 a vector of integer values.  if NULL the function will use all 
#	                 posterior samples, otherwise it will use only
#	                 the samples corresponding to the indices in ism,
#	                 e.g. 50, e.g. 50:100.
#
#	Returns: an ephy object with a list appended containing a vector of branch
#			 rates and the step size used for calculation.
dtRates = function(ephy, tau, ism = NULL)
{
	if(!'bamm-data' %in% class(ephy))
	{
		stop('Function requires a bammdata object');
	}
	
	ephy$eventBranchSegs = lapply(ephy$eventBranchSegs, function(x) x[order(x[,1]), ]); 
	
	phy = as.phylo.bammdata(ephy);
	phy = getStartStopTimes(phy);
	if(attributes(phy)$order != 'cladewise')
	{
		phy = reorder(phy,'cladewise');
	}
	tH = max(branching.times(phy));
	
	segmat = segMap(phy$edge[,2],phy$begin/tH,phy$end/tH,tau);
	segmat[,2] = segmat[,2] * tH;
	segmat[,3] = segmat[,3] * tH;
	
	tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
	
	if (storage.mode(segmat) != "double") stop('Exiting');
	if (storage.mode(tol) != "double") stop('Exiting');
	if (storage.mode(ephy) != "list") stop('Exiting');
	
	if (is.null(ism)) ism = as.integer(1:length(ephy$eventBranchSegs)) else ism = as.integer(ism);
	if (ism[length(ism)] > length(ephy$eventBranchSegs) )
	{
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
