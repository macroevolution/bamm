getMeanBranchLengthTree <-
function(obj,ndr=TRUE) {
	
	if ('bammdata' == class(obj)) {
		#v <- list(edge = obj$edge, Nnode = obj$Nnode,tip.label = obj$tip.label, edge.length = obj$edge.length);
		#attributes(v) <- list(names = c("edge","Nnode","tip.label","edge.length"),class="phylo",order="cladewise");
		v <- as.phylo.bammdata(obj);
		
		if (ndr) {
			if (obj$type == 'diversification') {
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$lambda_branch_matrix) - rowMeans(obj$mu_branch_matrix);
			} else if (obj$type == 'trait') {
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$beta_branch_matrix);
			}
		} else {
			if (obj$type == 'diversification') {
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$lambda_branch_matrix);
			} else if (obj$type == 'trait') {
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$beta_branch_matrix);
			}	
		}
		v$edge.length <- el;
		obj <- list();
		obj$phy <- v;
		obj$median <- median(el);
		obj$mean <- mean(el);			
	} else if (class(obj[[1]]) == 'phylo') {
		el = rowMeans(sapply(obj,with,edge.length));

		#v <- obj[[1]];
		#meanvec <- numeric(length(obj));
 
		#el <- v$edge.length;
		#meanvec[1] <- mean(v$edge.length);
		#for (i in 2:length(obj)) {
		#	el <- el + obj[[i]]$edge.length;
		#	meanvec[i] <- mean(obj[[i]]$edge.length);
		#}
		#el <- el / length(obj);
		#v$edge.length <- el;
		
		obj <- list();
		obj$phy <- obj[[1]];
		obj$median <- median(el);
		obj$mean <- mean(el);
	} else {
		stop("Method not implemented for supplied object class");
	}
	return(obj);		
}
