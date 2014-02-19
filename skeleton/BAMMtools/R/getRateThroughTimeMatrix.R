#############################################################
#
#	getRateThroughTimeMatrix(....)
#	
#	include or exclude options:
#	start.time		=	 start time (in units before present)
#						 if NULL, starts at root
#	end.time		=	 end time 
#						 if NULL, ends at present
#	nslices			=	 number of time cuts
#	Return
#	list with components:
#					$times	= the time vector
#					$lambda = speciation rate matrix
#					$mu 	= extinction rate matrix
# 					$type   = diversification or trait (needs to be extended to trait data)		
# returns object of class bamm-ratematrix
#	 

getRateThroughTimeMatrix <- function(ephy, start.time=NULL, end.time=NULL, nslices=100, node=NULL, nodetype = 'include') {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class 'bammdata'\n");
	}
	
	if (is.null(node)) {
		nodeset <- c(length(ephy$tip.label)+1, ephy$edge[,2]);
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
	#tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
	tol <- 0.00001
	
	goodTime <- function (vec, val, tol) {
		(vec[,2] - val <= tol) & (val - vec[,3] <= tol);
	}
	
	getRates <- function(time, es, events, isGoodNode) {
		isGoodTime <- goodTime(es, time, tol=tol);
		if (!(is.null(node))) { # only enter this if not the root. otherwise, only have to set once per i.
			isGoodNode <- es[,1] %in% nodeset;	
		}
		estemp <- es[isGoodTime & isGoodNode, ];
		if (is.vector(estemp)) {
			index <- estemp[4];
		} else {
			index <- estemp[,4];
		}
		tvv <- time - events$time[index];
		rates <- exponentialRate(tvv, events$lam1[index], events$lam2[index]);
		return(list(rates,index));
	}

	bySample <- function(counter, ephy) {
		es <- ephy$eventBranchSegs[[counter]];
		events <- ephy$eventData[[counter]];
		isGoodNode <- rep(TRUE, nrow(es));
		ret <- lapply(tvec, function(x) getRates(time = x, es, events, isGoodNode));
		mmRow <- unlist(lapply(ret, function(x) mean(x[[1]])));
		mumatRow <- unlist(lapply(ret, function(x) mean(events$mu1[x[[2]]])));
		return(list(mmRow,mumatRow));
	}

	
	ret <- lapply(1:length(ephy$eventBranchSegs), function(y) bySample(y, ephy));
	mm <- lapply(ret, function(x) x[[1]]);
	mm <- do.call(rbind, mm);
	mumat <- lapply(ret, function(x) x[[2]]);
	mumat <- do.call(rbind, mumat);
		
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
