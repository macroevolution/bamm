#############################################################
#
#	getTipRateStates(....)
#
# Returns a list with:
# lambda = matrix of tip rates where rows are species and columns are posterior samples, 
# mu if ephy$type == 'diversification',
# beta if ephy$type='trait',
# lambda.avg, mu.avg, beta.avg: named vector of average tip rates. 

getTipRateStates <- function(ephy, statistic='mean') {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}
	if (!statistic %in% c('mean','median')) {
		stop("statistic must be either 'mean' or 'median'.");
	}

	obj <- list();
	if (ephy$type == 'diversification') {
		obj$lambda <- do.call(cbind, ephy$tipLambda);
		rownames(obj$lambda) <- as.phylo.bammdata(ephy)$tip.label;
		
		obj$mu <- do.call(cbind, ephy$tipMu);
		rownames(obj$mu) <- as.phylo.bammdata(ephy)$tip.label;
		
		if (statistic == 'mean') {
			obj$lambda.avg <- rowMeans(obj$lambda);
			obj$mu.avg <- rowMeans(obj$mu);
		}
		if (statistic == 'median') {
			obj$lambda.avg <- apply(obj$lambda, 1, median);
			obj$mu.avg <- apply(obj$mu, 1, median);
		}
	}
	if (ephy$type == 'trait') {
		obj$beta <- do.call(cbind, ephy$tipLambda);
		rownames(obj$beta) <- as.phylo.bammdata(ephy)$tip.label;
		
		if (statistic == 'mean') {
			obj$beta.avg <- rowMeans(obj$lambda);
		}
		if (statistic == 'median') {
			obj$betabeta.avg <- apply(obj$beta, 1, median);
		}
	}
	return(obj);
}









