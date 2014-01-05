getRateThroughTimeMatrix <-
function(ephy, start.time=NULL, end.time=NULL, nslices=100, node=NULL, nodetype = 'include') {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class 'bammdata'\n");
	}
	
	if (is.null(node)) {
		nodeset <- ephy$edge[,2];
	} else if (!is.null(node) & nodetype == 'include') {
		#nodeset <- getDesc(ephy, node)$desc_set;
		nodeset <- unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set))
	} else if (!is.null(node) & nodetype == 'exclude') {
		nodeset <- setdiff( ephy$edge[,2], unlist(sapply(node, function(x) getDesc(ephy, x)$desc_set)));
	} else {
		stop('error in getRateThroughTimeMatrix\n');
	}
	
	bt <- branching.times(as.phylo.bammdata(ephy));
	maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[,1]))]);

	if (is.null(start.time)) {
		start.time <- max(bt) - maxpossible;
	}
	if (is.null(end.time)) {
		end.time <- max(bt);
	}
 
	tvec <- seq(start.time, end.time, length.out= nslices);
	tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
	mm <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));
	mumat <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));
	
	# this is much faster than 'safeCompare'.
	goodTime <- function (vec, val, tol) {
		(vec[,2] - val <= tol) & (val - vec[,3] <= tol);
	}
	
	for (i in 1:nrow(mm)) {
		es <- ephy$eventBranchSegs[[i]];
		events <- ephy$eventData[[i]];
		
		isGoodNode <- rep(TRUE, nrow(es));
		
		for (k in 1:length(tvec)) {
			isGoodTime <- goodTime(es, tvec[k], tol=tol);
			if (!(is.null(node))) { # only enter this if not the root. otherwise, only have to set once per i.
				isGoodNode <- es[,1] %in% nodeset;	
			}
			
			estemp <- es[isGoodTime & isGoodNode, ];
			tvv <- tvec[k] - events$time[estemp[,4]];
			rates <- exponentialRate(tvv, events$lam1[estemp[,4]], events$lam2[estemp[,4]]);
			mm[i, k] <- mean(rates);
			mumat[i,k] <- mean(events$mu1[estemp[,4]]);
		}
		
	}
	
	obj <- list();
	if (ephy$type == 'diversification') {
		obj$lambda <- mm;
		obj$mu <- mumat;
	}
	if (ephy$type == 'trait') {
		obj$beta <- mm;
	}
	obj$times <- tvec;
	
	class(obj) <- 'bamm-ratematrix';
	if (ephy$type=='diversification') {
		obj$type = 'diversification';
	} else {
		obj$type = 'trait';	
	}
	return(obj);
}
