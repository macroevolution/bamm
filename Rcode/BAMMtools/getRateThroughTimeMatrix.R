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

getRateThroughTimeMatrix <- function(ephy, start.time=NULL, end.time=NULL, nslices=100, node=NULL, nodetype = 'include'){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	if (is.null(node)){
		nodeset <- ephy$edge[,2];
	}else if (!is.null(node) & nodetype == 'include'){
		nodeset <- getDesc(ephy, node)$desc_set;
	}else if (!is.null(node) & nodetype == 'exclude'){
		nodeset <- setdiff( ephy$edge[,2],  getDesc(ephy, node)$desc_set);
	}else{
		stop('error in getRateThroughTimeMatrix\n');
	}
	
	bt <- branching.times(ephy);
		
	maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[,1]))]);

	if (is.null(start.time)){
		start.time <- max(bt) - maxpossible;
	}
	
	if (is.null(end.time)){
		end.time <- max(bt);
	}
 
	tvec <- seq(start.time, end.time, length.out= nslices);
	tol = 1*10^-decimals(ephy$eventBranchSegs[[1]][1,2]);
		
	mm <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));
	mumat <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));
		
	for (i in 1:nrow(mm)){
 	
		es <- ephy$eventBranchSegs[[i]];
		events <- ephy$eventData[[i]];
				
		for (k in 1:length(tvec)){
			isGoodTime <- safeCompare(es[,2],tvec[k],'<=',tol=tol) & safeCompare(es[,3],tvec[k],'>=',tol=tol)
			
			if (is.null(node)){ 
				isGoodNode <- rep(TRUE, nrow(es));
			}else{
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
	if (ephy$type == 'diversification'){
		obj$lambda <- mm;
		obj$mu <- mumat;
	}
	if (ephy$type == 'traits'){
		obj$beta <- mm;
	}
	obj$times <- tvec;
	
	class(obj) <- 'bamm-ratematrix';
	if(ephy$type=='diversification'){
		obj$type = 'diversification';
	}
	else{
		obj$type = 'traits';	
	}
	return(obj);
}
